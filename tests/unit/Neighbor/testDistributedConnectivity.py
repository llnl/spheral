#ATS:t0 = test(      SELF,       label="test distributed connectivity (1 proc)")
#ATS:t2 = testif(t0, SELF, np=2, label="test distributed connectivity (2 proc)")
#ATS:t4 = testif(t0, SELF, np=4, label="test distributed connectivity (4 proc)")

from Spheral import *
from SpheralTestUtilities import *
import os, shutil, time, sys
import mpi

import random
random.seed(4599281940)

title("distributed connectivity")

commandLine(
    # Discretization options
    dimension = 1,
    nx = 32,
    x0 = 0.0,
    x1 = 1.0,
    nPerh = 2.01,

    # Optionally randomize the node positions
    randomizeNodes = False,
    ranfrac = 0.2,
    
    # Testing options
    testGhosts = True,         # Check all the nodes, not just internal
    testOverlap = False,       # Check the overlap connectivity
    printErrors = True,        # Print the differences in the connectivity
    standardGlobalIDs = False, # This should work in 1D, but probably not otherwise
    testSuperset = False,      # Test if parallel connectivity is superset of serial; otherwise they must be equal

    # Base directory for storing connectivity
    baseTestDir = "data-distributed")

if dimension == 1:
    from Spheral1d import *
elif dimension == 2:
    from Spheral2d import *
else:
    from Spheral3d import *

#-------------------------------------------------------------------------------
# This test only has any meaning if we're testing an MPI enabled build
#-------------------------------------------------------------------------------
if mpi.is_fake_mpi():
    sys.exit(0)

#-------------------------------------------------------------------------------
# Set up the output directories
#-------------------------------------------------------------------------------
dataDir = os.path.join(baseTestDir,
                       "dim={}".format(dimension),
                       "nx={}".format(nx),
                       "x0,x1={},{}".format(x0, x1),
                       "nPerh={}".format(nPerh),
                       "overlap={}".format(testOverlap),
                       "glIDs={}".format(standardGlobalIDs),
                       "rand={}".format(randomizeNodes if not randomizeNodes else "{},{}".format(randomizeNodes, ranfrac)))
if mpi.rank == 0:
    if not os.path.exists(dataDir):
        os.makedirs(dataDir)
        
#-------------------------------------------------------------------------------
# Interpolation kernels
#-------------------------------------------------------------------------------
WT = TableKernel(WendlandC4Kernel(), 1000)
kernelExtent = WT.kernelExtent
output("WT")
output("kernelExtent")

delta = (x1 - x0) / nx
hmid = delta * nPerh * kernelExtent
hmin = hmid * 1.e-3
hmax = hmid * 1.e3

#-------------------------------------------------------------------------------
# Material properties
#-------------------------------------------------------------------------------
units = MKS()
gamma = 5./3.
mu = 1.0
eos = GammaLawGas(gamma, mu, units)

#-------------------------------------------------------------------------------
# Make the NodeList
#-------------------------------------------------------------------------------
nodes = makeFluidNodeList("nodes", eos, 
                          hmin = hmin,
                          hmax = hmax,
                          nPerh = nPerh,
                          kernelExtent = kernelExtent)
output("nodes")
output("nodes.hmin")
output("nodes.hmax")
output("nodes.nodesPerSmoothingScale")
    
#-------------------------------------------------------------------------------
# Seed the nodes
#-------------------------------------------------------------------------------
rho0 = 1.0
if dimension == 1:
    from DistributeNodes import distributeNodesInRange1d
    distributeNodesInRange1d([(nodes, nx, rho0, (x0, x1))],
                             nPerh = nPerh)
elif dimension == 2:
    from GenerateNodeDistribution2d import *
    generator = GenerateNodeDistribution2d(distributionType="lattice",
                                           nRadial = nx, nTheta = nx,
                                           xmin = (x0, x0),
                                           xmax = (x1, x1),
                                           rho = rho0,
                                           nNodePerh = nPerh)
    if mpi.procs > 1:
        from VoronoiDistributeNodes import distributeNodes2d
    else:
        from DistributeNodes import distributeNodes2d
    distributeNodes2d((nodes, generator))
else:
    from GenerateNodeDistribution3d import *
    generator = GenerateNodeDistribution3d(distributionType="lattice",
                                           n1 = nx, n2 = nx, n3 = nx,
                                           xmin = (x0, x0, x0),
                                           xmax = (x1, x1, x1),
                                           rho=rho0,
                                           nNodePerh = nPerh)
    if mpi.procs > 1:
        from VoronoiDistributeNodes import distributeNodes3d
    else:
        from DistributeNodes import distributeNodes3d
    distributeNodes3d((nodes, generator))

#-------------------------------------------------------------------------------
# Randomize the node positions
#-------------------------------------------------------------------------------
if randomizeNodes:
    delta = (x1 - x0) / nx
    pos = nodes.positions()
    for i in range(nodes.numInternalNodes):
        if dimension == 1:
            pos[i].x += ranfrac * delta * random.uniform(-1.0, 1.0)
        elif dimension == 2:
            pos[i].x += ranfrac * delta * random.uniform(-1.0, 1.0)
            pos[i].y += ranfrac * delta * random.uniform(-1.0, 1.0)
        elif dimension == 3:
            pos[i].x += ranfrac * delta * random.uniform(-1.0, 1.0)
            pos[i].y += ranfrac * delta * random.uniform(-1.0, 1.0)
            pos[i].z += ranfrac * delta * random.uniform(-1.0, 1.0)

#-------------------------------------------------------------------------------
# Make the DataBase
#-------------------------------------------------------------------------------
dataBase = DataBase()
dataBase.appendNodeList(nodes)
output("dataBase")
output("dataBase.numNodeLists")
output("dataBase.numFluidNodeLists")
output("dataBase.numNodes")
output("dataBase.numInternalNodes")

#-------------------------------------------------------------------------------
# Get some fields from DataBase
#-------------------------------------------------------------------------------
mass = dataBase.globalMass
position = dataBase.globalPosition
h = dataBase.globalHfield
fieldLists = [mass, position, h]

#-------------------------------------------------------------------------------
# Iterate h
#-------------------------------------------------------------------------------
domainbc = TreeDistributedBoundary.instance()
bounds = vector_of_Boundary([domainbc])
method = SPHSmoothingScale(IdealH, WT)
iterateIdealH(dataBase,
              vector_of_Physics([method]),
              bounds,
              100, # max h iterations
              1.e-8) # h tolerance
dataBase.updateConnectivityMap(testGhosts, testOverlap)

#-------------------------------------------------------------------------------
# Create the distributed boundary and compute the connectivity
#-------------------------------------------------------------------------------
for NN in dataBase.nodeLists:
    NN.numGhostNodes = 0
    NN.neighbor().updateNodes()
domainbc.setAllGhostNodes(dataBase)
domainbc.finalizeGhostBoundary()
for fl in fieldLists:
    domainbc.applyFieldListGhostBoundary(fl)
domainbc.finalizeGhostBoundary()
for NN in dataBase.nodeLists:
    NN.neighbor().updateNodes()
dataBase.updateConnectivityMap(testGhosts, testOverlap)
connectivity = dataBase.connectivityMap(testGhosts, testOverlap)
output("connectivity")
output("dataBase.numNodes")
output("dataBase.numInternalNodes")

#-------------------------------------------------------------------------------
# Compute the global node IDs
#-------------------------------------------------------------------------------
if standardGlobalIDs:
    if dimension > 1:
        print("warning: test will likely fail, as global IDs depends on mpi.procs")
    globalIndicesFL = globalNodeIDsAll(dataBase)
else:
    globalIndicesFL = mortonOrderIndices(dataBase)
    # globalIndicesFL = peanoHilbertOrderIndices(dataBase) # not working
domainbc.applyFieldListGhostBoundary(globalIndicesFL)
domainbc.finalizeGhostBoundary()
output("globalIndicesFL")

#-------------------------------------------------------------------------------
# Get connectivity in terms of global indices
#-------------------------------------------------------------------------------
globalIndices = []
globalConnectivity = []
for ni, nodeList in enumerate(dataBase.nodeLists):
    for i in range(nodeList.numNodes if testGhosts else nodeList.numInternalNodes):
        globali = globalIndicesFL(ni, i)
        globalConnectivityi = []
        connectivityi = connectivity.overlapConnectivityForNode(ni, i) if testOverlap else connectivity.connectivityForNode(ni, i)
        for nj in range(dataBase.numNodeLists):
            for j in connectivityi[nj]:
                globalj = globalIndicesFL(nj, j)
                globalConnectivityi.append(globalj)
        globalIndices.append(globali)
        globalConnectivity.append(set(globalConnectivityi))
# if mpi.procs > 1:
#     globalIndices = mpi.allreduce(globalIndices)
#     globalConnectivity = mpi.allreduce(globalConnectivity)
globalData = dict(list(zip(globalIndices, globalConnectivity)))

for key, vals in list(globalData.items()):
    print(key, " : ", vals)

#-------------------------------------------------------------------------------
# Save connectivity if serial, load if parallel
#-------------------------------------------------------------------------------
import pickle
dataFilename = os.path.join(dataDir, "connectivity")
if mpi.procs == 1:
    with open(dataFilename, 'wb') as f:
        pickle.dump(globalData, file = f)
else:
    if os.path.exists(dataFilename):
        with open(dataFilename, 'rb') as f:
            serialData = pickle.load(f)
    else:
        raise IOError("need to run in serial first")

#-------------------------------------------------------------------------------
# Check connectivity
#-------------------------------------------------------------------------------e
if mpi.procs > 1:
    # Go through serial connectivity and make sure parallel connectivity is the same
    # Note that we do this on every processor so that we are checking the ghost nodes
    # for connectivity as well
    numFailures = 0
    localKeys = set(globalData.keys())
    for i, cp in list(globalData.items()):
        if i in serialData:
            cs = serialData[i]

            # We only expect to find the nodes that are on this process, so reduce
            # to that set
            answer = cs.intersection(localKeys)

            # The parallel connectivity contains more points than the serial
            # (is this expected?), so we check whether the parallel connectivity
            # contains the complete serial connectivity
            if testSuperset:
                if not cp.issuperset(cs):
                    if printErrors:
                        diff = cs.difference(cp)
                        print("missing for point {}: {}".format(i, diff))
                    numFailures += 1
            else:
                diff = answer.symmetric_difference(cp)
                if diff:
                    if printErrors:
                        print("diff for point {}: {}".format(i, diff))
                    numFailures += 1
        else:
            numFailures += 1
            print("point {} not found in serial data".format(i))
    if numFailures > 0:
        raise ValueError("num points with connectivity differences: {}".format(numFailures))
    else:
        print("test passed")
        
    

    
