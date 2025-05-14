# SpheralTestUtilities -- small helper functions used in Spheral unit tests.

import sys
from math import *
from collections.abc import Iterable

#-------------------------------------------------------------------------------
# Helper method, sort a set of lists by the first one.
#-------------------------------------------------------------------------------
def multiSort(*args):
    # All the lists have to be the same length.
    for l in args:
        assert len(l) == len(args[0])

    # This is the obscure zip trick!
    result = list(zip(*sorted(zip(*args))))

    # Copy the sorted stuff back to the input arguments.
    for ilist in range(len(args)):
        for i in range(len(args[ilist])):
            args[ilist][i] = result[ilist][i]

    return

#-------------------------------------------------------------------------------
# Get the frame of the caller of a routine
# Based on Pat Millers getCurrentFrame()
#-------------------------------------------------------------------------------
def previousFrame():
    try:
        raise None
    except:
        type, value, traceback = sys.exc_info()
    currentFrame = traceback.tb_frame
    previousFrame = currentFrame.f_back
    return previousFrame

#-------------------------------------------------------------------------------
# Get the outermost (presumably global) frame.
# Based on Pat Millers getOuterFrame()
#-------------------------------------------------------------------------------
def globalFrame():
    currentFrame = previousFrame()
    globalFrame = currentFrame
    while currentFrame:
        globalFrame = currentFrame
        currentFrame = currentFrame.f_back
    return globalFrame

#-------------------------------------------------------------------------------
# Provide a standard title text
#-------------------------------------------------------------------------------
def title(titleText, lineLength=80):
    fillerText = "-"*((lineLength - len(titleText))//2)
    print(fillerText, titleText, fillerText)

#-------------------------------------------------------------------------------
# Provide a standard title text
#-------------------------------------------------------------------------------
def output(cmd, dict=None):
    if not dict:
        frame = globalFrame()
        dict = frame.f_globals
    print(cmd, " : ", eval(cmd, dict))

#-------------------------------------------------------------------------------
# Functions to help testing neighbor selection.
#-------------------------------------------------------------------------------
def findNeighborNodes(pos1, H1, radius, nodes,
                      potentials = None):
    pos = nodes.positions()
    H = nodes.Hfield()
    if potentials is None:
        from Spheral import vector_of_int
        master, coarse, refine = vector_of_int(), vector_of_int(), vector_of_int()
        nodes.neighbor().setMasterList(pos1, H1, master, coarse)
        nodes.neighbor().setRefineNeighborList(pos1, H1, coarse, refine)
        potentials = list(refine)
    return [i for i in potentials
            if (min((H1*(pos[i] - pos1)).magnitude(),
                    (H[i]*(pos[i] - pos1)).magnitude()) <= radius)]

def findGatherNeighborNodes(pos1, H1, radius, nodes,
                            potentials = None):
    pos = nodes.positions()
    H = nodes.Hfield()
    if potentials is None:
        from Spheral import vector_of_int
        master, coarse, refine = vector_of_int(), vector_of_int(), vector_of_int()
        nodes.neighbor().setMasterList(pos1, H1, master, coarse)
        nodes.neighbor().setRefineNeighborList(pos1, H1, coarse, refine)
        potentials = list(refine)
    return [i for i in potentials
            if (H1*(pos1 - pos[i])).magnitude() <= radius]

def findScatterNeighborNodes(pos1, radius, nodes,
                             potentials = None):
    pos = nodes.positions()
    H = nodes.Hfield()
    if potentials is None:
        from Spheral import vector_of_int
        master, coarse, refine = vector_of_int(), vector_of_int(), vector_of_int()
        nodes.neighbor().setMasterList(pos1, master, coarse)
        nodes.neighbor().setRefineNeighborList(pos1, coarse, refine)
        potentials = list(refine)
    return [i for i in potentials
            if (H[i]*(pos1 - pos[i])).magnitude() <= radius]

def findOverlapNeighbors(pos1, H1, radius, db):
    pos = db.globalPosition
    H = db.globalHfield
    result = [findNeighborNodes(pos1, H1, radius, nodes) for nodes in db.nodeLists]
    for iNL, inodes in enumerate(db.nodeLists):
        mygather = findGatherNeighborNodes(pos1, H1, radius, inodes)
        for i in mygather:
            for jNL, jnodes in enumerate(db.nodeLists):
                result[jNL] += findScatterNeighborNodes(pos(iNL, i), radius, jnodes)
    result = [list(set(x)) for x in result]
    return result

def findOverlapRegion(pos1, H1, pos2, H2, radius, nodes):
    pos = nodes.positions()
    result = [i for i in range(nodes.numInternalNodes)
              if ((H1*(pos[i] - pos1)).magnitude() <= radius and
                  (H2*(pos[i] - pos2)).magnitude() <= radius)]
    return result

def checkNeighbors(neighborList, answer):
    return len([x for x in answer if x not in neighborList]) == 0

#-------------------------------------------------------------------------------
# Print statistic about the H tensors for a set of NodeLists.
#-------------------------------------------------------------------------------
def hstats(nodeSet):
    import mpi
    for nodes in nodeSet:
        hmin, hmax, havg = 1e50, -1e50, 0.0
        hratmin, hratmax, hratavg = 1e50, -1e50, 0.0
        for H in nodes.Hfield().internalValues():
            Heigen = H.eigenValues()
            hmin = min(hmin, 1.0/(Heigen.maxElement()))
            hmax = max(hmax, 1.0/(Heigen.minElement()))
            havg += sum([1.0/x for x in Heigen])/len(Heigen)
            hrati = Heigen.minElement()/Heigen.maxElement()
            hratmin = min(hratmin, hrati)
            hratmax = max(hratmax, hrati)
            hratavg += hrati
        n = max(1, mpi.allreduce(nodes.numInternalNodes, mpi.SUM))
        hmin = mpi.allreduce(hmin, mpi.MIN)
        hmax = mpi.allreduce(hmax, mpi.MAX)
        havg = mpi.allreduce(havg, mpi.SUM)/n
        hratmin = mpi.allreduce(hratmin, mpi.MIN)
        hratmax = mpi.allreduce(hratmax, mpi.MAX)
        hratavg = mpi.allreduce(hratavg, mpi.SUM)/n
        print("%20s : %g %g %g %g %g %g" % (nodes.name, hmin, hmax, havg, hratmin, hratmax, hratavg))

#-------------------------------------------------------------------------------
# Print statistic about the numbers of neighbors per node.
#-------------------------------------------------------------------------------
def neighborStats(connectivityMap):
    import mpi
    for inodes in range(connectivityMap.numNodeLists()):
        nodes = connectivityMap.nodeList(inodes)
        n = [connectivityMap.numNeighborsForNode(nodes, i) for i in range(nodes.numInternalNodes)]
        print("%20s : %g %g %g" % (nodes.name,
                                   mpi.allreduce(min(n + [1e10]), mpi.MIN),
                                   mpi.allreduce(max(n + [0]), mpi.MAX),
                                   mpi.allreduce(sum(n), mpi.SUM)/max(1, mpi.allreduce(len(n), mpi.SUM))))

#-------------------------------------------------------------------------------
# Grab the node set values (internal, ghost, or all) for FieldLists.
#-------------------------------------------------------------------------------
def internalValues(fieldList):
    result = []
    for f in fieldList:
        vals = f.internalValues()
        result.extend(vals)
    return result

def ghostValues(fieldList):
    result = []
    for f in fieldList:
        vals = f.ghostValues()
        result.extend(vals)
    return result

def allValues(fieldList):
    result = []
    for f in fieldList:
        vals = f.allValues()
        result.extend(vals)
    return result

#-------------------------------------------------------------------------------
# Fuzzy comparisons.
#-------------------------------------------------------------------------------
def fuzzyEqual(lhs, rhs,
               fuzz = 1.0e-5):
    if isinstance(lhs, Iterable):
        assert isinstance(rhs, Iterable) and len(lhs) == len(rhs)
        return min([fuzzyEqual(x, y, fuzz) for (x, y) in zip(lhs, rhs)])
    else:
        return abs(lhs - rhs)/max(1.0, abs(lhs) + abs(rhs)) < fuzz;

def fuzzyLessThanOrEqual(lhs, rhs,
                         fuzz = 1.0e-5):
    return lhs < rhs or fuzzyEqual(lhs, rhs, fuzz);

def fuzzyGreaterThanOrEqual(lhs, rhs,
                            fuzz = 1.0e-5):
    return lhs > rhs or fuzzyEqual(lhs, rhs, fuzz)

def distinctlyLessThan(lhs, rhs,
                       fuzz = 1.0e-5):
    return lhs < rhs and not fuzzyEqual(lhs, rhs, fuzz)

def distinctlyGreaterThan(lhs, rhs,
                          fuzz = 1.0e-5):
    return lhs > rhs and not fuzzyEqual(lhs, rhs, fuzz)

#-------------------------------------------------------------------------------
# sign functions.
#-------------------------------------------------------------------------------
def sgn(x):
    if x > 0.0:
        return 1.0
    else:
        return -1.0

def sgn0(x):
    if x == 0.0:
        return 0.0
    elif x > 0.0:
        return 1.0
    else:
        return -1.0

#-------------------------------------------------------------------------------
# Return a sequence with just the unique elements (also sorted).
# I found this at http://www.peterbe.com/plog/uniqifiers-benchmark
#-------------------------------------------------------------------------------
def uniqueSequence(seq, idfun=None): 
    # order preserving
    if idfun is None:
        def idfun(x): return x
    seen = {}
    result = []
    for item in seq:
        marker = idfun(item)
        # in old Python versions:
        # if seen.has_key(marker)
        # but in new ones:
        if marker in seen: continue
        seen[marker] = 1
        result.append(item)
    return result

#-------------------------------------------------------------------------------
# Test the communicated info for consistency between domains in parallel.
#-------------------------------------------------------------------------------
def testParallelConsistency(mesh, xmin, xmax):
    from Spheral import hashPosition, quantizedPosition
    import mpi
    rank = mpi.rank
    numDomains = mpi.procs

    # Convert the parallel info to a convenient form.
    neighborDomains = [int(x) for x in mesh.neighborDomains]
    sharedNodes, sharedFaces = [], []
    for ll in mesh.sharedNodes:
        sharedNodes.append([int(x) for x in ll])
    for ll in mesh.sharedFaces:
        sharedFaces.append([int(x) for x in ll])
    assert len(neighborDomains) == len(mesh.sharedNodes)
    assert len(neighborDomains) == len(mesh.sharedFaces)

    # Bounding box scale for hashing positions.
    boxInv = xmax - xmin
    for i in range(len(boxInv)):
        boxInv[i] = 1.0/boxInv[i]

    # Local method to help accumulating the error messages in parallel.
    def allReduceMsg(msg):
        if msg == "ok":
            badProc = numDomains
        else:
            badProc = rank
        badProc = mpi.allreduce(badProc, mpi.MIN)
        if badProc != numDomains:
            msg = mpi.bcast(msg, root=badProc)
        return msg

    # Generic method to check communicated indices between domains.
    def checkConsistentCommInfo(myHashes, sharedIDs):
        msg = "ok"
        positions = [quantizedPosition(hashi, xmin, xmax) for hashi in myHashes]
        for sendProc in range(numDomains):
            numChecks = mpi.bcast(len(neighborDomains), root=sendProc)
            assert mpi.allreduce(numChecks, mpi.MIN) == mpi.allreduce(numChecks, mpi.MAX)
            for k in range(numChecks):
                if rank == sendProc:
                    ksafe = k
                else:
                    ksafe = 0
                recvProc = mpi.bcast(neighborDomains[ksafe], root=sendProc)
                recvHashes = mpi.bcast([myHashes[i] for i in sharedIDs[ksafe]], root=sendProc)
                recvPos = mpi.bcast([positions[i] for i in sharedIDs[ksafe]], root=sendProc)
                if rank == recvProc:
                    assert sendProc in neighborDomains
                    kk = neighborDomains.index(sendProc)
                    assert kk < len(sharedIDs)
                    if not ([myHashes[i] for i in sharedIDs[kk]] == recvHashes):
                        msg = ("Shared indices don't match %i %i\n   %s != %s\n    %s\n    %s" %
                               (rank, sendProc, 
                                str([myHashes[i] for i in sharedIDs[kk]]),
                                recvHashes,
                                [str(positions[i]) for i in sharedIDs[kk]],
                                [str(xi) for xi in recvPos]))
                msg = allReduceMsg(msg)
                if msg != "ok":
                    return msg
        return msg

    # Common code for checking hashed positions so we can do both nodes and faces.
    def checkHashConsistency(myHashes, sharedIDs, allowEmptyComm):
        msg = "ok"
        myHashSet = set(myHashes)
        for sendProc in range(numDomains):
            theirHashSet = mpi.bcast(myHashSet, root=sendProc)
            if sendProc != mpi.rank:
                commonHashes = myHashSet.intersection(theirHashSet)
                if len(commonHashes) > 0 and (not sendProc in neighborDomains):
                    msg = "Missed a neighbor domain : %i %i : %i" % (mpi.rank, sendProc, len(commonHashes))
    
                elif len(commonHashes) == 0 and (sendProc in neighborDomains) and not allowEmptyComm:
                    msg = "Erroneously communicating between domains : %i %i" % (mpi.rank, sendProc)
                
                elif len(commonHashes) > 0:
                    k = neighborDomains.index(sendProc)
                    if len(commonHashes) != len(sharedIDs[k]):
                        msg = "Size of shared elements does not match: %i %i : %i %i" % (mpi.rank, sendProc,
                                                                                         len(commonHashes),
                                                                                         len(sharedIDs[k]))
                    else:
                        sharedHashes = set([myHashes[i] for i in sharedIDs[k]])
                        if sharedHashes != commonHashes:
                            msg = ("Set of common hashes does not match: " + 
                                   str(sharedHashes) + " != " +
                                   str(commonHashes))
            msg = allReduceMsg(msg)
            if msg != "ok":
                return msg
        return msg

    # Check that the communicated mesh nodes.
    nodeHashes = [hashPosition(mesh.node(i).position, xmin, xmax, boxInv) for i in range(mesh.numNodes)]
    msg = checkConsistentCommInfo(nodeHashes, sharedNodes)
    if msg != "ok":
        return "Node failure : " + msg
    msg = checkHashConsistency(nodeHashes, sharedNodes, False)
    if msg != "ok":
        return "Node failure : " + msg

    # Check that the communicated mesh faces.
    faceHashes = [hashPosition(mesh.face(i).position, xmin, xmax, boxInv) for i in range(mesh.numFaces)]
    msg = checkConsistentCommInfo(faceHashes, sharedFaces)
    if msg != "ok":
        return "Face failure : " + msg
    msg = checkHashConsistency(faceHashes, sharedFaces, True)
    if msg != "ok":
        return "Face failure : " + msg

    return msg

#-------------------------------------------------------------------------------
# Estimate number of radial and axial SPH nodes.
#-------------------------------------------------------------------------------
def num_3d_cyl_nodes(rmin, rmax, zmin, zmax, thetamin, thetamax, nr0, nz0, Ntot):
    """
    This routine estimates values for nr and nz for a 3D cylindrical distribution
    to closely match a total desired number of SPH nodes

    Parameters:
    rmin, rmax: Inner/outer radial location
    zmin, zmax: Lower/upper axial locations
    thetamin, thetamax: Lower/upper polar angles
    nr0, nz0: Initial guesses for number of radial and axial SPH nodes
    Ntot: Desired total number of nodes
    Returns:
    nr, nz: Optimal number of radial and axial SPH nodes
    """
    import numpy as np
    from scipy.optimize import basinhopping

    class takestep:
        def __init__(self, x0):
            self.x0 = int(x0)
        def __call__(self, x):
            x += self.x0
            return x
    x0 = 1 # Optimal step size
    ig = [nr0, nz0]
    zlen = zmax - zmin
    if nr0*nz0 > Ntot:
        raise Exception("Initial guesses are too high")
    rlen = rmax - rmin
    Dtheta = thetamax - thetamin
    mintheta = int(Dtheta/(0.5*np.pi) + 0.5)
    def calc_ntot(x):
        nr = int(x[0])
        nz = int(x[1])
        dz = zlen / nz
        dr = rlen / nr
        maxd = max(dr, dz)
        new_tot = 0
        for i in range(nr):
            ri = rmin + (i + 0.5)*dr
            n_th = max(mintheta, int(ri*Dtheta/maxd))
            new_tot += n_th
        return abs(new_tot*nz - Ntot)
    result = basinhopping(calc_ntot, ig, take_step=takestep(x0))
    nr, nz = result.x[0], result.x[1]
    return int(nr), int(nz)
