#ATS:dim = 3
#ATS:dimstr = f"{dim}d"
#ATS:ntotals = [50, 60, 70, 80]
#ATS:for nxv in ntotals:
#ATS:    nx = int(nxv**dim)
#ATS:    test_name = f"EVALDERIV_{nxv}"
#ATS:    cali_name = f"{test_name}.cali"
#ATS:    inputs = f"--raja True --ntotal {nx} --testDim {dimstr} --adiakData 'test_name: {test_name}' --caliperFilename {cali_name}"
#ATS:    test(SELF, label=test_name, clas=inputs, ngpu=1, np=1, nt=1, caliper_filename=cali_name, independent=False)

#-------------------------------------------------------------------------------
# Isolated evaluateDerivatives for performance testing.
#-------------------------------------------------------------------------------
from Spheral import *
from SpheralTestUtilities import *

title("Evaluate derivatives performance test")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(
    # Parameters for seeding nodes.
    ntotal = 192000,
    L = 1.0,
    rho1 = 1.0,
    eps = 0.0,
    nPerh = 2.01,
    hmin = 0.0001, 
    hmax = 10.0,

    # What hydro operator should we test?
    HydroChoice = "SPH",
    solid = True,
    raja = False,

    # Should we randomly perturb the positions?
    ranfrac = 0.2,
    seed = 14892042,

    # What test problem are we doing?
    testDim = "3d",
    testCase = "linear",

    # The fields we're going to interpolate.
    # Linear coefficients: y = y0 + m0*x
    y0 = 1.0,
    m0 = 1.0,

    # Quadratic coefficients: y = y2 + m2*x^2
    y2 = 1.0,
    m2 = 0.5,

    gamma = 5.0/3.0,
    mu = 1.0,

    # Test parameters
    velCorrection = False,
    numIters = 10,
    doJitter = False,
    initVel = False,

    # Parameters for iterating H.
    iterateH = False,
    maxHIterations = 200,
    Htolerance = 1.0e-4
)

testDim = testDim.lower()
assert testCase in ("linear", "quadratic", "step")
assert testDim in ("1d", "2d", "3d", "spherical")
x0 = 0.
x1 = L
x2 = L
x3 = L

#-------------------------------------------------------------------------------
# Appropriately set generic object names based on the test dimensionality.
#-------------------------------------------------------------------------------
if testDim == "spherical":
    from SphericalSpheral import *
    nx1 = ntotal
else:
    exec("from Spheral%s import *" % testDim, globals())
    dval = int(testDim[0])
    nx1 = int((ntotal+1)**(1/dval))
nx2 = nx1
nx3 = nx1
print(f"N: ({nx1}, {nx2}, {nx3})")
HydroConstructor = eval(HydroChoice)

#-------------------------------------------------------------------------------
# Create a random number generator.
#-------------------------------------------------------------------------------
import random
random.seed(seed)

#-------------------------------------------------------------------------------
# Material properties.
#-------------------------------------------------------------------------------
eos = GammaLawGasMKS(gamma, mu)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
if testDim == "spherical":
    WT = TableKernel3d(BSplineKernel3d(), 1000)
else:
    WT = TableKernel(BSplineKernel(), 1000)
output("WT")
kernelExtent = WT.kernelExtent

#-------------------------------------------------------------------------------
# Make the NodeList.
#-------------------------------------------------------------------------------
if solid:
    nodes1 = makeSolidNodeList("nodes1", eos,
                               hmin = hmin,
                               hmax = hmax,
                               nPerh = nPerh)
else:
    nodes1 = makeFluidNodeList("nodes1", eos,
                               hmin = hmin,
                               hmax = hmax,
                               nPerh = nPerh)
output("nodes1")
output("nodes1.hmin")
output("nodes1.hmax")
output("nodes1.nodesPerSmoothingScale")

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
if testDim == "1d":
    from DistributeNodes import distributeNodesInRange1d
    distributeNodesInRange1d([(nodes1, [(nx1, rho1, (x0, x1))])], nPerh = nPerh)
elif testDim == "2d":
    from DistributeNodes import distributeNodes2d
    from GenerateNodeDistribution2d import GenerateNodeDistribution2d
    from CompositeNodeDistribution import CompositeNodeDistribution
    gen = GenerateNodeDistribution2d(nx1, nx2, rho1,
                                     distributionType = "lattice",
                                     xmin = (x0, x0),
                                     xmax = (x1, x2),
                                     nNodePerh = nPerh,
                                     SPH = True)
    distributeNodes2d((nodes1, gen))

elif testDim == "3d":
    from DistributeNodes import distributeNodes3d
    from GenerateNodeDistribution3d import GenerateNodeDistribution3d
    from CompositeNodeDistribution import CompositeNodeDistribution
    gen = GenerateNodeDistribution3d(nx1, nx2, nx3, rho1,
                                     distributionType = "lattice",
                                     xmin = (x0, x0, x0),
                                     xmax = (x1, x2, x3),
                                     nNodePerh = nPerh,
                                     SPH = True)
    distributeNodes3d((nodes1, gen))

elif testDim == "spherical":
    from PeanoHilbertDistributeNodes import distributeNodes1d
    from GenerateSphericalNodeDistribution1d import GenerateSphericalNodeDistribution1d
    from CompositeNodeDistribution import CompositeNodeDistribution
    gen = GenerateSphericalNodeDistribution1d(nr = nx1,
                                              rho = rho1,
                                              rmin = x0,
                                              rmax = x1,
                                              nNodePerh = nPerh)
    distributeNodes1d((nodes1, gen))

else:
    raise ValueError("Only tests cases for 1d, 2d, 3d, and Spherical.") 

output("nodes1.numNodes")

#-------------------------------------------------------------------------------
# Optionally randomly jitter the node positions.
#-------------------------------------------------------------------------------
if doJitter:
    # Set node properties.
    epsv = nodes1.specificThermalEnergy()    
    dx = (x1 - x0)/nx1
    dy = (x2 - x0)/nx2
    dz = (x3 - x0)/nx3
    pos = nodes1.positions()
    for i in range(nodes1.numInternalNodes):
        epsv[i] = eps
        if testDim in ("1d", "spherical"):
            pos[i].x += ranfrac * dx * random.uniform(-1.0, 1.0)
        elif testDim == "2d":
            pos[i].x += ranfrac * dx * random.uniform(-1.0, 1.0)
            pos[i].y += ranfrac * dy * random.uniform(-1.0, 1.0)
        elif testDim == "3d":
            pos[i].x += ranfrac * dx * random.uniform(-1.0, 1.0)
            pos[i].y += ranfrac * dy * random.uniform(-1.0, 1.0)
            pos[i].z += ranfrac * dz * random.uniform(-1.0, 1.0)

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase()
output("db")
output("db.appendNodeList(nodes1)")
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Initialize the velocity.
#-------------------------------------------------------------------------------
if initVel:
    f = nodes1.velocity()
    pos = nodes1.positions()
    for i in range(nodes1.numInternalNodes):
        for j in range(db.nDim):
            x = pos[i][j]
            if testCase == "linear":
                f[i][j] = (y0 + m0*x)
            elif testCase == "quadratic":
                f[i][j] = (y2 + m2*x*x)
            elif testCase == "step":
                if x < x1:
                    f[i][j] = y0
                else:
                    f[i][j] = 2*y0

#-------------------------------------------------------------------------------
# Build a hydro
#-------------------------------------------------------------------------------
hydro_args = dict(dataBase = db, W = WT, correctVelocityGradient = velCorrection)
if raja:
    hydro_args.update(dict(RAJA = True))

if HydroChoice not in ("PSHP", "PASPH"):
    hydro_args.update(dict(gradhCorrection = True))

hydro = HydroConstructor(**hydro_args)

#-------------------------------------------------------------------------------
# Build a controller
#-------------------------------------------------------------------------------
# The controller has non-trivial logic to get our packages initialized and
# organized correctly, so we just use it here
integrator = SynchronousRK2Integrator(db, [hydro])
control = SpheralController(integrator,
                            kernel = WT,
                            iterateInitialH = iterateH)

state = State(db, integrator.physicsPackages())
derivs = StateDerivatives(db, integrator.physicsPackages())
for i in range(numIters):
    print(f"Iter {i}")
    integrator.preStepInitialize(state, derivs)
    integrator.initializeDerivatives(0.0, 1.0, state, derivs)
    derivs.Zero()
    integrator.evaluateDerivatives(0.0, 1.0, db, state, derivs)
    integrator.finalizeDerivatives(0.0, 1.0, db, state, derivs)
