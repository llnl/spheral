from math import *
from Spheral import *
import mpi
import distributeNodesGeneric
from NodeGeneratorBase import ConstantRho

from spheralDimensions import spheralDimensions
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# A factory string pattern for making distributors based on different algorithms
#-------------------------------------------------------------------------------
def makeDistributeNodesMethod(distribute, SpheralModule):
    assert type(distribute) == str

    def f(*listOfNodeTuples):
        "Domain decompose using the %s method)" % distribute
        import mpi
        if mpi.procs > 1:
            if not hasattr(SpheralModule, distribute):
                raise RuntimeError("ERROR: %s is not available in the Spheral module" % distribute)
            distributeNodesGeneric.distributeNodesGeneric(listOfNodeTuples,
                                                          SpheralModule.DataBase,
                                                          SpheralModule.globalNodeIDs,
                                                          getattr(SpheralModule, distribute))
        else:
            distributeNodesGeneric.distributeNodesGeneric(listOfNodeTuples,
                                                          SpheralModule.DataBase,
                                                          SpheralModule.globalNodeIDs,
                                                          None)

    return f

# For backwards compatibility we provide a default distributor here
from PeanoHilbertDistributeNodes import *

#-------------------------------------------------------------------------------
# Old method for handling 1-D distributions.  Deprecated now.
# Distribute sets of nodes evenly in the given ranges (1d)
# Input:  a list of tuples with (nodeList, (global) numNodes, (global) xRange)
#-------------------------------------------------------------------------------
def distributeNodesInRange1d(listOfNodeTuples,
                             nPerh,
                             reverse = False):

    # Count how many nodes total we need to distribute
    numNodesPerNodeList = []
    for nodeTuple in listOfNodeTuples:
        assert len(nodeTuple) == 2 or len(nodeTuple) == 4
        if len(nodeTuple) == 2:
            numNodesPerNodeList.append(sum([n for n, rho, (x1, x2) in nodeTuple[1]]))
        else:
            numNodesPerNodeList.append(nodeTuple[1])
    totalNumNodes = sum(numNodesPerNodeList)

    # Determine the global node ID range for this domain.
    nNodesPerDomain = totalNumNodes // mpi.procs
    minNodeID = mpi.rank*nNodesPerDomain
    if mpi.rank == mpi.procs - 1:
        maxNodeID = totalNumNodes
    else:
        maxNodeID = minNodeID + nNodesPerDomain

    # Iterate over each NodeList and assign its nodes.
    nCumulative = 0
    for nodeTuple in listOfNodeTuples:
        if len(nodeTuple) == 2:
            nodes, initialConditions = nodeTuple
            pos = nodes.positions()
            mass = nodes.mass()
            H = nodes.Hfield()
            rho = nodes.massDensity()
            for n, rho0, (x0, x1) in initialConditions:
                if type(rho0) == float:
                    rhof = ConstantRho(rho0)
                else:
                    rhof = rho0
                if minNodeID <= (nCumulative + n) and maxNodeID >= nCumulative:
                    iglobal0 = max(minNodeID, nCumulative)
                    iglobal1 = min(maxNodeID, nCumulative + n)
                    numNewNodes = iglobal1 - iglobal0
                    nodeOffset = iglobal0 - nCumulative
                    indexOffset = nodes.numInternalNodes
                    nodes.numInternalNodes += numNewNodes
                    dx = (x1 - x0)/max(1, n)
                    Hi = SymTensor1d(1.0/max(nodes.hmin, min(nodes.hmax, nPerh*dx)))
                    for i in range(numNewNodes):
                        if reverse:
                            xi = x1 - (nodeOffset + i + 0.5)*dx
                        else:
                            xi = x0 + (nodeOffset + i + 0.5)*dx
                        pos[indexOffset + i].x = xi
                        mass[indexOffset + i] = dx*rhof(xi)
                        rho[indexOffset + i] = rhof(xi)
                        H[indexOffset + i] = Hi
                nCumulative += n
        else:
            nodes, n, rho0, (x0, x1) = nodeTuple
            if type(rho0) == float:
                rhof = ConstantRho(rho0)
            else:
                rhof = rho0
            pos = nodes.positions()
            mass = nodes.mass()
            H = nodes.Hfield()
            rho = nodes.massDensity()
            if minNodeID <= (nCumulative + n) and maxNodeID >= nCumulative:
                iglobal0 = max(minNodeID, nCumulative)
                iglobal1 = min(maxNodeID, nCumulative + n)
                numNewNodes = iglobal1 - iglobal0
                assert numNewNodes > 0
                nodeOffset = iglobal0 - nCumulative
                indexOffset = nodes.numInternalNodes
                nodes.numInternalNodes += numNewNodes
                dx = (x1 - x0)/max(1, n)
                Hi = SymTensor1d(1.0/max(nodes.hmin, min(nodes.hmax, nPerh*dx)))
                for i in range(numNewNodes):
                    if reverse:
                        xi = x1 - (nodeOffset + i + 0.5)*dx
                    else:
                        xi = x0 + (nodeOffset + i + 0.5)*dx
                    pos[indexOffset + i].x = xi
                    mass[indexOffset + i] = dx*rhof(xi)
                    rho[indexOffset + i] = rhof(xi)
                    H[indexOffset + i] = Hi
            nCumulative += n

    return

#-------------------------------------------------------------------------------
# Handle psuedo-1D spherical node distribution similarly to above.
# Input:  a list of tuples with (nodeList, (global) numNodes, (global) xRange)
#-------------------------------------------------------------------------------
def distributeNodesInSphericalRange3d(listOfNodeTuples,
                                      nPerh = 2.01,
                                      reverse = False):
    
    # Count how many nodes total we need to distribute.
    globalNumNodes = 0
    for (nodeList, n, rho, range) in listOfNodeTuples:
        assert n >= 0
        globalNumNodes += n

    # Determine the global node ID range this process will cover.
    nNodesPerDomain = globalNumNodes // mpi.procs
    minNodeID = mpi.rank*nNodesPerDomain
    if mpi.rank == mpi.procs - 1:
        maxNodeID = globalNumNodes
    else:
        maxNodeID = minNodeID + nNodesPerDomain
    
    # Set the size of each node list, and assign node positions, masses, and
    # H's.
    cumulativeN = 0
    for (nodeList, n, rho, (r0, r1)) in listOfNodeTuples:
        dr = (r1 - r0)/max(1, n)
        minGlobalNodeListID = max(0, min(n, minNodeID - cumulativeN))
        maxGlobalNodeListID = max(0, min(n, maxNodeID - cumulativeN))
        nNodesThisDomain = maxGlobalNodeListID - minGlobalNodeListID
        assert minGlobalNodeListID >= minGlobalNodeListID
        nodeList.numInternalNodes = nNodesThisDomain
        cumulativeN += n
        Hxx = 1.0/max(nodeList.hmin, min(nodeList.hmax, nPerh*dr))
        Hi = SymTensor3d(Hxx, 0.0, 0.0,
                         0.0, Hxx, 0.0,
                         0.0, 0.0, Hxx)

        vol = 4.0/3.0*pi*(r1**3 - r0**3)
        neff = vol/dr**3
        mi = vol*rho/neff

        for globalNodeID in range(minGlobalNodeListID, maxGlobalNodeListID):
            localNodeID = globalNodeID - minGlobalNodeListID
            if reverse:
                nodeList.positions()[localNodeID].x = r1 - (globalNodeID + 0.5)*dr
            else:
                nodeList.positions()[localNodeID].x = r0 + (globalNodeID + 0.5)*dr
            ri = nodeList.positions()[localNodeID].x

            nodeList.mass()[localNodeID] = mi
            nodeList.massDensity()[localNodeID] = rho
            nodeList.Hfield()[localNodeID] = Hi

    return

