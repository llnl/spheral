#ATS:test(SELF, "--dimension 1 --nx 64 --preconditioner 'none' --numpyTest True", np=1, label="Linear solver, 1d GMRES")
#ATS:test(SELF, "--dimension 2 --nx 16 --preconditioner 'ilu' --reflecting True --numpyTest True", np=4, label="Linear solver, 2d ILU")
#ATS:test(SELF, "--dimension 3 --nx 8 --preconditioner 'amg' --periodic True --numpyTest True --nPerh 4.01", np=8, label="Linear solver, 3d AMG")
# #ATS:test(SELF, "--dimension 3 --nx 50 --preconditioner 'amg' --reflecting True --nPerh 4.01 --numTests 10 --numSolves 5 --caliperConfig 'spot,runtime-report'", np=32, label="Linear solver, timing check")

#-------------------------------------------------------------------------------
# Set up a system of equations and solve them with the linear solver operator
#-------------------------------------------------------------------------------
from Spheral import *
from SpheralTestUtilities import *
import numpy as np
import time, mpi
title("Linear solver test")

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
    nPerh = 2.01,

    # Solver options
    preconditioner = "amg",

    # Test reinitialization
    numTests = 4, # Multiple tests, same graph but different data
    numSolves = 2, # Solve the same thing multiple times for timing

    # Include boundary conditions
    reflecting = False,
    periodic = False,

    # Test options
    tolerance = 1e-4,
    numpyTest = False,
    seed = None,
    saveSystem = False, # Save hypre matrices
    exactGuess = False,
)    
exec("from Spheral%id import *" % dimension, globals())

if seed is None:
    seed = np.random.randint(0, 2**8)
output("seed")
np.random.seed(seed)

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
numLocal = nodes.numInternalNodes
output("numLocal")

#-------------------------------------------------------------------------------
# Create DataBase
#-------------------------------------------------------------------------------
dataBase = DataBase()
dataBase.appendNodeList(nodes)
output("dataBase")

#-------------------------------------------------------------------------------
# Set initial conditions
#------------------------------------------------------------------------------
inp = dataBase.newFluidScalarFieldList(0.0, "input")
out = dataBase.newFluidScalarFieldList(0.0, "output")
inp0 = dataBase.newFluidScalarFieldList(0.0, "initial input")
out0 = dataBase.newFluidScalarFieldList(0.0, "initial output")
for i in range(numLocal):
    inp[0][i] = np.random.uniform(0.0, 2.0)
    inp0[0][i] = inp[0][i]
    out[0][i] = np.random.uniform(-2.0, 2.0)
    out0[0][i] = out[0][i]
pos = dataBase.fluidPosition

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
for bc in bounds:
    for fl in [inp, out, inp0, out0]:
        bc.applyFieldListGhostBoundary(fl)
for bc in bounds:
    bc.finalizeGhostBoundary()

#-------------------------------------------------------------------------------
# Create map using generated points
#------------------------------------------------------------------------------
flatConnectivity = FlatConnectivity()
flatConnectivity.computeIndices(dataBase)
flatConnectivity.computeGlobalIndices(dataBase, bounds)
flatConnectivity.computeBoundaryInformation(dataBase, bounds)
output("flatConnectivity")
numGlobal = flatConnectivity.numGlobalNodes()
output("numGlobal")

#-------------------------------------------------------------------------------
# Test functions
#------------------------------------------------------------------------------
maxErr = 0.0
def norm(s1, s2, eps = 1.0e-16):
    v1 = np.array(s1)
    v2 = np.array(s2)
    return np.amax(2.0 * np.abs(v1 - v2) / (np.abs(v1) + np.abs(v2) + eps))
def check(s1, s2, label, tol = tolerance, eps = 1.0e-16):
    err = norm(s1, s2, eps)
    global maxErr
    if err > maxErr:
        maxErr = err
    if err > tol:
        message = "{}\n\tcalculated={} does not agree with\n\texpected={}\n\terr={}".format(label, s1, s2, err)
        if ignoreErr:
            print(message)
        else:
            raise ValueError(message)
        
#-------------------------------------------------------------------------------
# Generate matrix
#------------------------------------------------------------------------------
globRowInd = []
numColsPerRow = []
globColInd = []
colVal = [[] for t in range(numTests)]
multDirect = [np.zeros(numLocal) for t in range(numTests)]
if numpyTest:
    spRow = []
    spCol = []
    spDat = [[] for t in range(numTests)]
for i in range(numLocal):
    # Also check unique neighbors as a bonus since this function is designed for matrices
    globali = flatConnectivity.localToGlobal(i)
    numNeighbors = flatConnectivity.numNonConstNeighbors(i)
    numConst = flatConnectivity.numConstNeighbors(i)
    numUnique, localCols, globalCols, constCols, indexMap = flatConnectivity.uniqueNeighborIndices(i)
    check(localCols, flatConnectivity.nonConstNeighborIndices(i), "local cols")
    check(globalCols, flatConnectivity.globalNeighborIndices(i), "global cols")
    if numConst > 0:
        check(constCols, flatConnectivity.constNeighborIndices(i), "const cols")

    # Set row information
    numColsPerRow.append(numUnique)
    globRowInd.append(globali)

    # Initialize the column data
    firstInd = len(globColInd)
    globColInd.extend([-1 for j in range(numUnique)])
    for t in range(numTests):
        colVal[t].extend([0.0 for j in range(numUnique)])
    
    for j in range(numNeighbors):
        ind = localCols[j]
        globInd = globalCols[j]
        uniqueInd = indexMap[j]
        flatInd = firstInd + uniqueInd # Index for globColInd and colVal
        globColInd[flatInd] = globInd
        if numpyTest:
            spRow.append(globali)
            spCol.append(globInd)
        for t in range(numTests):
            val = -1.0 * np.random.uniform(1.5, 2.5) if ind == i else np.random.uniform(0.5, 1.5) / numNeighbors
            colVal[t][flatInd] += val
            multDirect[t][i] += val * inp0[0][ind]
            if numpyTest:
                spDat[t].append(val)
if numpyTest:
    import scipy.sparse as sps
    globDat = [mpi.allreduce(spDat[t]) for t in range(numTests)]
    globRow = mpi.allreduce(spRow)
    globCol = mpi.allreduce(spCol)
    npMat = [sps.coo_matrix((globDat[t], (globRow, globCol)), shape=(numGlobal, numGlobal)).tocsr() for t in range(numTests)]
    lhsVec = [inp0[0][i] for i in range(numLocal)]
    lhsVec = np.array(mpi.allreduce(lhsVec))
    multNumpy = [npMat[t] @ lhsVec for t in range(numTests)]

#-------------------------------------------------------------------------------
# Initialize solver
#------------------------------------------------------------------------------
options = HypreOptions()
options.maxNumberOfIterations = 1000
options.saveLinearSystem = saveSystem
if preconditioner == "amg":
    options.preconditionerType = HypreOptions.AMGPreconditioner
elif preconditioner == "ilu":
    options.preconditionerType = HypreOptions.ILUPreconditioner
    options.factorLevelILU = 5
else:
    options.preconditionerType = HypreOptions.NoPreconditioner
solver = HypreLinearSolver(options);
output("solver")

start = time.perf_counter()
firstGlobInd = flatConnectivity.firstGlobalIndex()
mpi.barrier()
solver.initialize(numLocal, firstGlobInd, vector_of_unsigned(numColsPerRow))
hypre_init_time = mpi.allreduce(time.perf_counter() - start, mpi.MAX)
output("hypre_init_time")

#-------------------------------------------------------------------------------
# Run tests
#------------------------------------------------------------------------------
for t in range(numTests):
    mpi.barrier()
    start = time.perf_counter()
    solver.beginFill()
    solver.setMatRows(numLocal,
                      vector_of_unsigned(numColsPerRow),
                      vector_of_unsigned(globRowInd),
                      vector_of_unsigned(globColInd),
                      vector_of_double(colVal[t]))
    hypre_fill_time = mpi.allreduce(time.perf_counter() - start, mpi.MAX)
    output("hypre_fill_time")
    
    mpi.barrier()
    start = time.perf_counter()
    solver.assemble()
    hypre_assemble_time = mpi.allreduce(time.perf_counter() - start, mpi.MAX)
    output("hypre_assemble_time")

    mpi.barrier()
    start = time.perf_counter()
    solver.finalize()
    hypre_finalize_time = mpi.allreduce(time.perf_counter() - start, mpi.MAX)
    output("hypre_finalize_time")

    #-------------------------------------------------------------------------------
    # Multiply through to get RHS
    #------------------------------------------------------------------------------
    # output("inp0[0][0]")
    
    lhsVec = vector_of_double([inp0[0][i] for i in range(numLocal)])
    mpi.barrier()
    solver.set(LinearSolver.LHS, numLocal, firstGlobInd, lhsVec)
    mpi.barrier()
    solver.multiply()
    
    rhsVec = solver.get(LinearSolver.RHS, numLocal, firstGlobInd)
    for i in range(numLocal):
        out[0][i] = rhsVec[i]

    for i in range(numLocal):
        if abs(out[0][i] - out0[0][i]) < 1e-16:
            raise ValueError("output has not been changed")
    
    for bc in bounds:
        for fl in [inp, out]:
            bc.applyFieldListGhostBoundary(fl)
    for bc in bounds:
        bc.finalizeGhostBoundary()
    
    #-------------------------------------------------------------------------------
    # Check that the multiplication result is as expected
    #------------------------------------------------------------------------------
    if numpyTest:
        firstGlobal = flatConnectivity.localToGlobal(0)
        for i in range(numLocal):
            globali = flatConnectivity.localToGlobal(i)
            if abs(multNumpy[t][globali] - multDirect[t][i]) > tolerance:
                print("warning: numpy and direct mult differ: np={}, dir={}".format(multNumpy[t][globali], multDirect[t][globali]))
            if abs(multNumpy[t][globali] - out[0][i]) > tolerance:
                raise ValueError("numpy comparison failed for test {}, point {}: exp={}, calc={}".format(t, i, multDirect[t][i], out[0][i]))
        print("numpy comparison successful: {} -> {}, {}".format(lhsVec[0], multDirect[t][0], multNumpy[t][firstGlobal]))
        
    for i in range(numLocal):
        if abs(multDirect[t][i] - out[0][i]) > tolerance:
            raise ValueError("multiplication incorrect for test {}, point {}: exp={}, calc={}".format(t, i, multDirect[t][i], out[0][i]))
    
    #-------------------------------------------------------------------------------
    # Solve the linear system to get back the LHS
    #------------------------------------------------------------------------------
    if exactGuess:
        for i in range(numLocal):
            inp[0][i] = inp0[0][i]
    else:
        for i in range(numLocal):
            inp[0][i] = np.random.uniform(0.0, 1.0)

    for i in range(numSolves):
        lhsVec = vector_of_double([inp[0][i] for i in range(numLocal)])
        rhsVec = vector_of_double([out[0][i] for i in range(numLocal)])
        
        start = time.perf_counter()
        mpi.barrier()
        solver.set(LinearSolver.LHS, numLocal, firstGlobInd, lhsVec)
        solver.set(LinearSolver.RHS, numLocal, firstGlobInd, rhsVec)
        mpi.barrier()
        solver.solve()
        hypre_solve_time = mpi.allreduce(time.perf_counter() - start, mpi.MAX)
        output("hypre_solve_time")
    
    mpi.barrier()
    lhsVec = solver.get(LinearSolver.LHS, numLocal, firstGlobInd)
    for i in range(numLocal):
        inp[0][i] = lhsVec[i]

    # output("inp[0][0]")

    #-------------------------------------------------------------------------------
    # Check the results
    #------------------------------------------------------------------------------
    for i in range(numLocal):
        if abs(inp[0][i] - inp0[0][i]) > tolerance:
            raise ValueError("error too high for test {} point {}: begin={}, end={}".format(t + 1, i, inp0[0][i], inp[0][i]))
    print("test {} passed".format(t+1))

