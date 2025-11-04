#ATS:test(SELF, "--dimension 1", np=1, label="Flat connectivity, 1d 1p")
#ATS:test(SELF, "--dimension 2 --reflecting True", np=4, label="Flat connectivity, 2d")
#ATS:test(SELF, "--dimension 3 --periodic True", np=8, label="Flat connectivity, 3d")

#-------------------------------------------------------------------------------
# Set up a system of equations and solve them with the linear solver operator
#-------------------------------------------------------------------------------
from Spheral import *
from SpheralTestUtilities import *
import numpy as np
import time, mpi
title("Flat connectivity test")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(
    # Spatial stuff
    dimension = 1,
    nx = 32,
    x0 = 0.0,
    x1 = 1.0,

    # Interpolation kernel
    nPerh = 4.01,

    # Include boundary conditions
    reflecting = False,
    periodic = False,

    # Test options
    verboseErr = False,
    ignoreErr = False,
)    
exec("from Spheral%id import *" % dimension, globals())

#-------------------------------------------------------------------------------
# Create nodes
#------------------------------------------------------------------------------
units = MKS()
output("units")

eos = GammaLawGas(5.0/3.0, 1.0, units)
output("eos")

WT = TableKernel(WendlandC4Kernel(), 1000)
kernelExtent = WT.kernelExtent
output("WT")

delta = (x1 - x0) / nx
hmid = delta * nPerh * kernelExtent
hmin = hmid * 1.e-3
hmax = hmid * 1.e3
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
output("nodes.numNodes")

#-------------------------------------------------------------------------------
# Create DataBase
#-------------------------------------------------------------------------------
dataBase = DataBase()
dataBase.appendNodeList(nodes)
output("dataBase")

#-------------------------------------------------------------------------------
# Construct boundary conditions.
#-------------------------------------------------------------------------------
limits = [x0, x1]
bounds = []
assert(not (reflecting and periodic))
if reflecting or periodic:
    for d in range(dimension):
        planes = []
        for n in range(2):
            point = Vector.zero
            point[d] = limits[n]
            normal = Vector.zero
            normal[d] = 1.0 if n == 0 else -1.0
            planes.append(Plane(point, normal))
            if reflecting:
                bounds.append(ReflectingBoundary(planes[n]))
        if periodic:
            bounds.append(PeriodicBoundary(planes[0], planes[1]))
if mpi.procs > 1:
    bounds.append(TreeDistributedBoundary.instance())
output("bounds")

#-------------------------------------------------------------------------------
# Iterate h
#-------------------------------------------------------------------------------
method = SPHSmoothingScale(IdealH, WT)
iterateIdealH(dataBase,
              [method],
              bounds,
              100, # max h iterations
              1.e-4) # h tolerance

#-------------------------------------------------------------------------------
# Finish setting up connectivity
#-------------------------------------------------------------------------------
nodes.numGhostNodes = 0
nodes.neighbor().updateNodes()
for bc in bounds:
    bc.setAllGhostNodes(dataBase)
    bc.finalizeGhostBoundary()
    nodes.neighbor().updateNodes()
dataBase.updateConnectivityMap()
conn = dataBase.connectivityMap()
output("conn")

#-------------------------------------------------------------------------------
# Get global indices directly
#------------------------------------------------------------------------------
globalInd = globalNodeIDsAll(dataBase)

for bc in bounds:
    bc.applyFieldListGhostBoundary(globalInd)
    bc.finalizeGhostBoundary()
output("globalInd")

#-------------------------------------------------------------------------------
# Create map using generated points
#------------------------------------------------------------------------------
fc = FlatConnectivity()
fc.computeIndices(dataBase)
fc.computeGlobalIndices(dataBase, bounds)
fc.computeBoundaryInformation(dataBase, bounds)
output("fc")

#-------------------------------------------------------------------------------
# Test function
#------------------------------------------------------------------------------
def check(s1, s2, label):
    if s1 != s2:
        message = "fail: {}\n\tcalculated={}\n\t  expected={}".format(label, s1, s2)
        if ignoreErr:
            print(message)
        else:
            raise ValueError(message)
    elif verboseErr:
        print("pass: {}".format(label))
        
#-------------------------------------------------------------------------------
# Check that things are initialized
#------------------------------------------------------------------------------
check(fc.indexingInitialized(), True, "indexing initialized")
check(fc.overlapIndexingInitialized(), False, "overlap indexing not initialized")
check(fc.globalIndexingInitialized(), True, "global indexing initialized")
check(fc.surfaceIndexingInitialized(), False, "surface indexing not initialized")
check(fc.boundaryInformationInitialized(), True, "boundary information initialized")

#-------------------------------------------------------------------------------
# Calculate information directly in Python
#------------------------------------------------------------------------------
class PythonFlatConnectivity:
    def __init__(self):
        self.numNodes = dataBase.numNodes
        self.numInternalNodes = dataBase.numInternalNodes
        self.numGlobalNodes = mpi.allreduce(self.numInternalNodes, mpi.SUM)
        numNodesPerProc = mpi.allreduce([self.numInternalNodes], mpi.SUM)
        self.firstGlobalIndex = sum(numNodesPerProc[:mpi.rank])
        self.lastGlobalIndex = self.firstGlobalIndex + self.numInternalNodes - 1
        self.numBoundaryNodes = 0
        self.neighbors = [conn.connectivityForNode(0, i) for i in range(self.numInternalNodes)]
        return
    
    def nodeToLocal(self, nodeListi, nodei):
        return nodei
    def localToNode(self, locali):
        return [0, locali]
    def localToGlobal(self, locali):
        return self.firstGlobalIndex + locali
    def numNeighbors(self, locali):
        return len(self.neighborIndices(locali))
    def numConstNeighbors(self, locali):
        return 0
    def numNonConstNeighbors(self, locali):
        return self.numNeighbors(locali)
    def isConstantBoundaryNode(self, locali):
        return False
    def neighborIndices(self, locali):
        return [locali] + list(self.neighbors[locali][0])
    def constNeighborIndices(self, locali):
        return []
    def nonConstNeighborIndices(self, locali):
        return self.neighborIndices(locali)
    def globalNeighborIndices(self, locali):
        return [globalInd[0][i] for i in self.neighborIndices(locali)]
    def uniqueNeighborIndices(self, locali):
        v1 = self.nonConstNeighborIndices(locali)
        v2 = self.globalNeighborIndices(locali)
        v3 = []
        globalMap = {}
        i = 0
        for j in range(len(v2)):
            if v2[j] in globalMap:
                continue
            globalMap[v2[j]] = i
            i += 1
        v4 = [globalMap[j] for j in v2]
        return [len(globalMap),
                v1, v2,
                v3, v4]
    
    # # Cop out on these for now
    # def localToFlat(self, locali, localj):
    #     return fc.localToFlat(locali, localj)
    # def flatToLocal(self, locali, flatj):
    #     return fc.flatToLocal(locali, flatj)

fcp = PythonFlatConnectivity()
output("fcp")

#-------------------------------------------------------------------------------
# Grab the function for both classes and call them using the provided args
#------------------------------------------------------------------------------
def compare(name, *args, **kwargs):
    """Call the same-named method on fc and fcp with given arguments."""
    f1 = getattr(fc, name)
    f2 = getattr(fcp, name)

    r1 = f1(*args, **kwargs)
    r2 = f2(*args, **kwargs)
    try:
        r1 = list(r1)
        r2 = list(r2)
    except:
        pass
    # print(r1, r2)
    argstr = ", ".join([repr(a) for a in args] + ["{}={}".format(k, v) for k, v in kwargs.items()])
    return check(r1, r2, name + " " + argstr)

#-------------------------------------------------------------------------------
# Compare the Python implementation to the C++ one
#------------------------------------------------------------------------------
numInternalNodes = nodes.numInternalNodes

for name in ["nodeToLocal"]:
    for i in range(numInternalNodes):
        compare(name, 0, i)
    if not verboseErr:
        print("{} passed".format(name))

for name in ["localToNode", "localToGlobal", "numNeighbors",
             "numConstNeighbors", "numNonConstNeighbors", "isConstantBoundaryNode",
             "neighborIndices", "constNeighborIndices", "nonConstNeighborIndices",
             "globalNeighborIndices", "uniqueNeighborIndices"]:
    for i in range(numInternalNodes):
        compare(name, i)
    if not verboseErr:
        print("{} passed".format(name))

# for name in ["localToFlat", "flatToLocal"]:
#     for i in range(numInternalNodes):
#         numNeighbors = conn.numNeighborsForNode(0, i)
#         for j in range(numNeighbors):
#             compare(name, i, j)
#     if not verboseErr:
#         print("{} passed".format(name))
