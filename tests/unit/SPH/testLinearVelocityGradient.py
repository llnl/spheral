#ATS:for testDim in ("1d", "2d", "3d"):
#ATS:    for HydroChoice in ("SPH", "ASPH", "PSPH", "PASPH"):
#ATS:        for solid in (False, True):
#ATS:            test(SELF, "--graphics False --nx1 10 --nx2 10 --testCase linear --testDim %s --HydroChoice %s" % (testDim, HydroChoice), 
#ATS:                 label="%s linear gradient correction test -- %s (serial)" % (HydroChoice, testDim))
#-------------------------------------------------------------------------------
# Unit test of the linear velocity gradient correction for SPH.
#-------------------------------------------------------------------------------
from Spheral import *
from SpheralTestUtilities import *

title("SPH velocity gradient linear correction test")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(
    # Parameters for seeding nodes.
    nx1 = 50,     
    nx2 = 50,
    rho1 = 1.0,
    rho2 = 1.0,
    eps1 = 0.0,
    eps2 = 0.0,
    x0 = 0.0,
    x1 = 0.5,
    x2 = 1.0,
    nPerh = 2.01,
    hmin = 0.0001, 
    hmax = 10.0,

    # What hydro operator should we test?
    HydroChoice = "SPH",
    solid = True,

    # Should we randomly perturb the positions?
    ranfrac = 0.2,
    seed = 14892042,

    # What test problem are we doing?
    testDim = "1d",
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

    # Parameters for iterating H.
    iterateH = True,
    maxHIterations = 200,
    Htolerance = 1.0e-4,

    # Parameters for passing the test
    derivativeTolerance = 5.0e-5,

    graphics = True,
    plotKernels = False,
    plotSPH = True,
)

testDim = testDim.lower()
assert testCase in ("linear", "quadratic", "step")
assert testDim in ("1d", "2d", "3d", "spherical")

#-------------------------------------------------------------------------------
# Appropriately set generic object names based on the test dimensionality.
#-------------------------------------------------------------------------------
if testDim == "spherical":
    from SphericalSpheral import *
else:
    exec("from Spheral%s import *" % testDim, globals())
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
    distributeNodesInRange1d([(nodes1, [(nx1, rho1, (x0, x1)),
                                        (nx2, rho2, (x1, x2))])], nPerh = nPerh)
elif testDim == "2d":
    from DistributeNodes import distributeNodes2d
    from GenerateNodeDistribution2d import GenerateNodeDistribution2d
    from CompositeNodeDistribution import CompositeNodeDistribution
    gen1 = GenerateNodeDistribution2d(nx1, nx1 + nx2, rho1,
                                      distributionType = "lattice",
                                      xmin = (x0, x0),
                                      xmax = (x1, x2),
                                      nNodePerh = nPerh,
                                      SPH = True)
    gen2 = GenerateNodeDistribution2d(nx2, nx1 + nx2, rho2,
                                      distributionType = "lattice",
                                      xmin = (x1, x0),
                                      xmax = (x2, x2),
                                      nNodePerh = nPerh,
                                      SPH = True)
    gen = CompositeNodeDistribution(gen1, gen2)
    distributeNodes2d((nodes1, gen))

elif testDim == "3d":
    from DistributeNodes import distributeNodes3d
    from GenerateNodeDistribution3d import GenerateNodeDistribution3d
    from CompositeNodeDistribution import CompositeNodeDistribution
    gen1 = GenerateNodeDistribution3d(nx1, nx1 + nx2, nx1 + nx2, rho1,
                                      distributionType = "lattice",
                                      xmin = (x0, x0, x0),
                                      xmax = (x1, x2, x2),
                                      nNodePerh = nPerh,
                                      SPH = True)
    gen2 = GenerateNodeDistribution3d(nx2, nx1 + nx2, nx1 + nx2, rho2,
                                      distributionType = "lattice",
                                      xmin = (x1, x0, x0),
                                      xmax = (x2, x2, x2),
                                      nNodePerh = nPerh,
                                      SPH = True)
    gen = CompositeNodeDistribution(gen1, gen2)
    distributeNodes3d((nodes1, gen))

elif testDim == "spherical":
    from PeanoHilbertDistributeNodes import distributeNodes1d
    from GenerateSphericalNodeDistribution1d import GenerateSphericalNodeDistribution1d
    from CompositeNodeDistribution import CompositeNodeDistribution
    gen1 = GenerateSphericalNodeDistribution1d(nr = nx1,
                                               rho = rho1,
                                               rmin = x0,
                                               rmax = x1,
                                               nNodePerh = nPerh)
    gen2 = GenerateSphericalNodeDistribution1d(nr = nx2,
                                               rho = rho2,
                                               rmin = x1,
                                               rmax = x2,
                                               nNodePerh = nPerh)
    gen = CompositeNodeDistribution(gen1, gen2)
    distributeNodes1d((nodes1, gen))

else:
    raise ValueError("Only tests cases for 1d, 2d, 3d, and Spherical.") 

output("nodes1.numNodes")

# Set node properties.
eps = nodes1.specificThermalEnergy()
for i in range(nx1):
    eps[i] = eps1
for i in range(nx2):
    eps[i + nx1] = eps2

#-------------------------------------------------------------------------------
# Optionally randomly jitter the node positions.
#-------------------------------------------------------------------------------
dx1 = (x1 - x0)/nx1
dx2 = (x2 - x1)/nx2
dy = (x2 - x0)/(nx1 + nx2)
dz = (x2 - x0)/(nx1 + nx2)
pos = nodes1.positions()
for i in range(nodes1.numInternalNodes):
    if pos[i].x < x1:
        dx = dx1
    else:
        dx = dx2
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
if HydroChoice in ("PSPH", "PASPH"):
    hydro = HydroConstructor(dataBase = db,
                             W = WT,
                             correctVelocityGradient = False)
else:
    hydro = HydroConstructor(dataBase = db,
                             W = WT,
                             gradhCorrection = True,
                             correctVelocityGradient = False)

#-------------------------------------------------------------------------------
# Build a controller
#-------------------------------------------------------------------------------
# The controller has non-trivial logic to get our packages initialized and
# organized correctly, so we just use it here
integrator = SynchronousRK2Integrator(db, [hydro])
control = SpheralController(integrator,
                            kernel = WT,
                            iterateInitialH = iterateH)

#-------------------------------------------------------------------------------
# Invoke the SPH evaluateDerivatives, which will put velocity gradients in the 
# derivatives state object.
#-------------------------------------------------------------------------------
# Compute the uncorrected derivatives.
hydro.correctVelocityGradient = False
state = State(db, integrator.physicsPackages())
derivs = StateDerivatives(db, integrator.physicsPackages())
integrator.preStepInitialize(state, derivs)
integrator.initializeDerivatives(0.0, 1.0, state, derivs)
derivs.Zero()
integrator.evaluateDerivatives(0.0, 1.0, db, state, derivs)
integrator.finalizeDerivatives(0.0, 1.0, db, state, derivs)
dfSPH0_fl = derivs.tensorFields(HydroFieldNames.velocityGradient)
dfSPH0_fl.copyFields()
dfSPH0 = dfSPH0_fl[0]

# Compute the corrected derivatives.
hydro.correctVelocityGradient = True
derivs.Zero()
integrator.preStepInitialize(state, derivs)
integrator.initializeDerivatives(0.0, 1.0, state, derivs)
integrator.evaluateDerivatives(0.0, 1.0, db, state, derivs)
integrator.finalizeDerivatives(0.0, 1.0, db, state, derivs)
dfSPH1_fl = derivs.tensorFields(HydroFieldNames.velocityGradient)
dfSPH1_fl.copyFields()
dfSPH1 = dfSPH1_fl[0]

#-------------------------------------------------------------------------------
# Prepare variables to accumulate the test values.
#-------------------------------------------------------------------------------
position_fl = db.fluidPosition
weight_fl = db.fluidMass
H_fl = db.fluidHfield

# Extract the field state for the following calculations.
positions = position_fl[0]
H = H_fl[0]

#-------------------------------------------------------------------------------
# Prepare the answer to check against.
#-------------------------------------------------------------------------------
xans = [positions[i].x for i in range(nodes1.numInternalNodes)]
dyans = TensorField("derivative answer", nodes1)
for i in range(nodes1.numInternalNodes):
    if testCase == "linear":
        dyans[i] = m0 * Tensor.one
    elif testCase == "quadratic":
        dyans[i] = 2*m2*xans[i] * Tensor.one
    elif testCase == "step":
        dyans[i] = Tensor.zero

#-------------------------------------------------------------------------------
# Check our answers accuracy.
#-------------------------------------------------------------------------------
errdySPH0 =  ScalarField("uncorrected SPH derivative error", nodes1)
errdySPH1 =  ScalarField("corrected SPH derivative error", nodes1)
for i in range(nodes1.numInternalNodes):
    errdySPH0[i] =  (dfSPH0[i] - dyans[i]).selfDoubledot()
    errdySPH1[i] =  (dfSPH1[i] - dyans[i]).selfDoubledot()

maxdySPHerror0 = max([abs(x) for x in errdySPH0])
maxdySPHerror1 = max([abs(x) for x in errdySPH1])

print("Maximum error in uncorrected SPH: %g" % maxdySPHerror0)
print("Maximum error in   corrected SPH: %g" % maxdySPHerror1)

#-------------------------------------------------------------------------------
# Plot the things.
#-------------------------------------------------------------------------------
if graphics:
    from SpheralMatplotlib import *

    p3 = plotField(dyans,
                   yFunction = "%s.xx",
                   plotStyle = "k-",
                   lineTitle = "Answer")
    plotField(dfSPH0,
              yFunction = "%s.xx",
              plotStyle = "r*",
              lineTitle = HydroChoice + " (uncorrected)",
              plot = p3)
    plotField(dfSPH1,
              yFunction = "%s.xx",
              plotStyle = "ko",
              lineTitle = HydroChoice + " (corrected)",
              kwords = {"fillstyle" : "none"},
              plot = p3)
    p3.set_title("Derivative values")
    p3.legend(loc = "best")
              
    p4 = plotField(errdySPH0,
                   plotStyle = "r*",
                   lineTitle = HydroChoice + " (uncorrected)")
    plotField(errdySPH1,
              plotStyle = "ko",
              lineTitle = HydroChoice + " (corrected)",
              kwords = {"fillstyle" : "none"},
              plot = p4)
    p4.set_title("Error in derivatives")
    p4.legend(loc = "best")

    # # If we're in 2D dump a silo file too.
    # if testDim == "2d":
    #     from siloPointmeshDump import siloPointmeshDump
    #     siloPointmeshDump("testLinearVelocityGradient_%s_2d" % testCase,
    #                       fields = [dfSPH0, dfSPH1,
    #                                 dyans,
    #                                 errdySPH0, errdySPH1])

#-------------------------------------------------------------------------------
# Check the maximum corrected SPH error and fail the test if it's out of bounds.
#-------------------------------------------------------------------------------
if maxdySPHerror1 > derivativeTolerance:
    raise ValueError("corrected SPH derivative error out of bounds: %g > %g" % (maxdySPHerror1, derivativeTolerance))
