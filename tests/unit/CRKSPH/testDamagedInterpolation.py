#ATS:test(SELF, "--graphics False --nx1 10 --nx2 10 --testDim 1d", label="CRKSPH interpolation test -- 1D (serial)")
#ATS:test(SELF, "--graphics False --nx1 10 --nx2 10 --testDim 2d", label="CRKSPH interpolation test -- 2D (serial)")
#ATS:test(SELF, "--graphics False --nx1 10 --nx2 10 --testDim 3d", label="CRKSPH interpolation test -- 3D (serial)")
#-------------------------------------------------------------------------------
# A set of tests to compare how different meshless methods interpolate fields.
#-------------------------------------------------------------------------------
from Spheral import *
from SpheralTestUtilities import *

title("Interpolation tests")

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
    nPerh = 1.25,
    hmin = 0.0001, 
    hmax = 10.0,

    # Radius to damage nodes in.
    breakRadius = 0.04,

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

    numGridLevels = 20,
    topGridCellSize = 20.0,

    # Parameters for iterating H.
    iterateH = True,
    maxHIterations = 200,
    Htolerance = 1.0e-4,

    # Parameters for passing the test
    interpolationTolerance = 5.0e-7,
    derivativeTolerance = 5.0e-5,

    graphics = True,
    plotKernels = False,
)

assert testCase in ("linear", "quadratic", "step")
assert testDim in ("1d", "2d", "3d")

FacetedVolume = {"1d" : Box1d,
                 "2d" : Polygon,
                 "3d" : Polyhedron}[testDim]

#-------------------------------------------------------------------------------
# Appropriately set generic object names based on the test dimensionality.
#-------------------------------------------------------------------------------
exec("from Spheral%s import *" % testDim, globals())

## import Spheral
## for name in [x for x in Spheral.__dict__ if testDim in x]:
##     exec("%s = Spheral.__dict__['%s']" % (name.replace(testDim, ""), name))

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
WT = TableKernel(BSplineKernel(), 1000)
output("WT")
kernelExtent = WT.kernelExtent

#-------------------------------------------------------------------------------
# Make the NodeList.
#-------------------------------------------------------------------------------
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
    gen1 = GenerateNodeDistribution2d(nx1, nx1, rho1,
                                      distributionType = "lattice",
                                      xmin = (x0, x0),
                                      xmax = (x1, x2),
                                      nNodePerh = nPerh,
                                      SPH = True)
    gen2 = GenerateNodeDistribution2d(nx2, nx2, rho2,
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
    gen1 = GenerateNodeDistribution3d(nx1, nx1, nx1, rho1,
                                      distributionType = "lattice",
                                      xmin = (x0, x0, x0),
                                      xmax = (x1, x1, x2),
                                      nNodePerh = nPerh,
                                      SPH = True)
    gen2 = GenerateNodeDistribution3d(nx2, nx2, nx2, rho2,
                                      distributionType = "lattice",
                                      xmin = (x1, x0, x0),
                                      xmax = (x2, x1, x2),
                                      nNodePerh = nPerh,
                                      SPH = True)
    gen = CompositeNodeDistribution(gen1, gen2)
    distributeNodes3d((nodes1, gen))

else:
    raise ValueError("Only tests cases for 1d,2d and 3d.") 

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
for i in range(nodes1.numInternalNodes):
    if i < nx1:
        dx = dx1
    else:
        dx = dx2
    nodes1.positions()[i].x += ranfrac * dx * random.uniform(-1.0, 1.0)

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase()
output("db")
output("db.appendNodeList(nodes1)")
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Iterate the h to convergence if requested.
#-------------------------------------------------------------------------------
if iterateH:
    bounds = vector_of_Boundary()
    method = SPHSmoothingScale()
    iterateIdealH(db,
                  bounds,
                  WT,
                  method,
                  maxHIterations,
                  Htolerance)

#-------------------------------------------------------------------------------
# Create the damage field.
#-------------------------------------------------------------------------------
D_fl = db.newFluidSymTensorFieldList(SymTensor.one * 0.01, "D")
D = D_fl[0]
pos = nodes1.positions()
for i in range(nodes1.numInternalNodes):
    if abs(pos[i].x - x1) < breakRadius:
        D[i] = SymTensor.one

# For now we'll use a zero damage gradient.
# gradD_fl = db.newFluidVectorFieldList(Vector.zero, "grad D")

#-------------------------------------------------------------------------------
# Find the gradient of the damage using the simpler SPH definition.
#-------------------------------------------------------------------------------
Dt_fl = db.newFluidScalarFieldList(0.0, "damage trace")
Dt = Dt_fl[0]

db.updateConnectivityMap(True)
cm = db.connectivityMap()
position_fl = db.fluidPosition
H_fl = db.fluidHfield
mass_fl = db.fluidMass
rho_fl = db.fluidMassDensity
weight_fl = db.newFluidScalarFieldList(0.0, "weight")
weight = weight_fl[0]
mass = mass_fl[0]
rho = rho_fl[0]

for i in range(nodes1.numInternalNodes):
    weight[i] = mass[i]/rho[i]
    Dt[i] = D[i].Trace()

gradD_fl = gradient(Dt_fl, position_fl, weight_fl, mass_fl, rho_fl, H_fl, WT)
gradD_fl[0].name = "grad D"

#-------------------------------------------------------------------------------
# Initialize our field.
#-------------------------------------------------------------------------------
f = ScalarField("test field", nodes1)
for i in range(nodes1.numInternalNodes):
    x = nodes1.positions()[i].x
    if testCase == "linear":
        f[i] = y0 + m0*x
    elif testCase == "quadratic":
        f[i] = y2 + m2*x*x
    elif testCase == "step":
        if x < x1:
            f[i] = y0
        else:
            f[i] = 2*y0

#-------------------------------------------------------------------------------
# Prepare variables to accumulate the test values.
#-------------------------------------------------------------------------------
A_fl = db.newFluidScalarFieldList(0.0, "A")
B_fl = db.newFluidVectorFieldList(Vector.zero, "B")
gradA_fl = db.newFluidVectorFieldList(Vector.zero, "gradA")
gradB_fl = db.newFluidTensorFieldList(Tensor.zero, "gradB")

couple = DamagedNodeCoupling(D_fl, gradD_fl, H_fl)
computeCRKSPHCorrections(cm, WT, weight_fl, position_fl, H_fl, couple,
                         A_fl, B_fl, gradA_fl, gradB_fl)

# Extract the field state for the following calculations.
positions = position_fl[0]
weight = weight_fl[0]
H = H_fl[0]
A = A_fl[0]
B = B_fl[0]
gradA = gradA_fl[0]
gradB = gradB_fl[0]

#-------------------------------------------------------------------------------
# We also check the C++ interpolation and gradient methods.
#-------------------------------------------------------------------------------
f_fl = ScalarFieldList()
f_fl.appendField(f)
fCRKSPH_fl = interpolateCRKSPH(f_fl, position_fl, weight_fl, H_fl, A_fl, B_fl, 
                               cm, WT, couple)
dfCRKSPH_fl = gradientCRKSPH(f_fl, position_fl, weight_fl, H_fl,
                             A_fl, B_fl, gradA_fl, gradB_fl,
                             cm, WT, couple)
fCRKSPH = fCRKSPH_fl[0]
dfCRKSPH = dfCRKSPH_fl[0]

#-------------------------------------------------------------------------------
# Prepare the answer to check against.
#-------------------------------------------------------------------------------
xans = [positions[i].x for i in range(nodes1.numInternalNodes)]
yans = ScalarField("interpolation answer", nodes1)
dyans = ScalarField("derivative answer", nodes1)
for i in range(nodes1.numInternalNodes):
    if testCase == "linear":
        yans[i] = y0 + m0*xans[i]
        dyans[i] = m0
    elif testCase == "quadratic":
        yans[i] = y2 + m2*xans[i]*xans[i]
        dyans[i] = 2*m2*xans[i]
    elif testCase == "step":
        if i < nx1:
            yans[i] = y0
        else:
            yans[i] = 2*y0
        dyans[i] = 0.0

#-------------------------------------------------------------------------------
# Check our answers accuracy.
#-------------------------------------------------------------------------------
erryCRKSPH = ScalarField("CRKSPH interpolation error", nodes1)
errdyCRKSPH = ScalarField("CRKSPH derivative error", nodes1)
for i in range(nodes1.numInternalNodes):
    erryCRKSPH[i] =  fCRKSPH[i] - yans[i]
    errdyCRKSPH[i] = dfCRKSPH[i].x - dyans[i]

maxyCRKSPHerror = max([abs(x) for x in erryCRKSPH])
maxdyCRKSPHerror = max([abs(x) for x in errdyCRKSPH])

print("Maximum errors (interpolation): CRKSPH = %g" % (maxyCRKSPHerror))
print("Maximum errors   (derivatives): CRKSPH = %g" % (maxdyCRKSPHerror))

#-------------------------------------------------------------------------------
# Plot the things.
#-------------------------------------------------------------------------------
if graphics:
    from SpheralGnuPlotUtilities import *
    import Gnuplot
    xans = [positions[i].x for i in range(nodes1.numInternalNodes)]

    # Interpolated values.
    ansdata = Gnuplot.Data(xans, yans.internalValues(),
                           with_ = "lines",
                           title = "Answer",
                           inline = True)
    CRKSPHdata = Gnuplot.Data(xans, fCRKSPH.internalValues(),
                            with_ = "points",
                            title = "CRKSPH",
                            inline = True)
    errCRKSPHdata = Gnuplot.Data(xans, erryCRKSPH.internalValues(),
                               with_ = "points",
                               title = "CRKSPH",
                               inline = True)

    p1 = generateNewGnuPlot()
    p1.plot(ansdata)
    p1.replot(CRKSPHdata)
    p1("set key top left")
    p1.title("Interpolated values")
    p1.refresh()

    p2 = generateNewGnuPlot()
    p2.replot(errCRKSPHdata)
    p2.title("Error in interpolation")
    p2.refresh()

    # Derivative values.
    dansdata = Gnuplot.Data(xans, dyans.internalValues(),
                            with_ = "lines",
                            title = "Answer",
                            inline = True)
    dCRKSPHdata = Gnuplot.Data(xans, [x.x for x in dfCRKSPH.internalValues()],
                             with_ = "points",
                             title = "CRKSPH",
                             inline = True)
    errdCRKSPHdata = Gnuplot.Data(xans, errdyCRKSPH.internalValues(),
                                with_ = "points",
                                title = "CRKSPH",
                                inline = True)

    p3 = generateNewGnuPlot()
    p3.plot(dansdata)
    p3.replot(dCRKSPHdata)
    p3("set key top left")
    p3.title("Derivative values")
    p3.refresh()

    p4 = generateNewGnuPlot()
    p4.replot(errdCRKSPHdata)
    p4.title("Error in derivatives")
    p4.refresh()

    # If we're in 2D dump a silo file too.
    if testDim == "2d":
        from SpheralVoronoiSiloDump import SpheralVoronoiSiloDump
        dumper = SpheralVoronoiSiloDump("testInterpolation_%s_2d" % testCase,
                                        listOfFields = [fCRKSPH, dfCRKSPH,
                                                        yans, dyans,
                                                        erryCRKSPH, errdyCRKSPH],
                                        listOfFieldLists = [weight_fl, D_fl, gradD_fl,
                                                            A_fl, B_fl, gradA_fl, gradB_fl])
        dumper.dump(0.0, 0)
        # from siloPointmeshDump import siloPointmeshDump
        # siloPointmeshDump("testInterpolation_%s_2d" % testCase,
        #                   fields = [fSPH, fCRKSPH, dfSPH, dfCRKSPH,
        #                             yans, dyans,
        #                             errySPH, erryCRKSPH, errdySPH, errdyCRKSPH],
        #                   fieldLists = [weight_fl, 
        #                                 A_fl, B_fl, gradA_fl, gradB_fl,
        #                                 dfCRKSPH_fl])

#-------------------------------------------------------------------------------
# Check the maximum CRKSPH error and fail the test if it's out of bounds.
#-------------------------------------------------------------------------------
if maxyCRKSPHerror > interpolationTolerance:
    raise ValueError("CRKSPH interpolation error out of bounds: %g > %g" % (maxyCRKSPHerror, interpolationTolerance))

if maxdyCRKSPHerror > derivativeTolerance:
    raise ValueError("CRKSPH derivative error out of bounds: %g > %g" % (maxdyCRKSPHerror, derivativeTolerance))
