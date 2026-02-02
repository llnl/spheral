#ATS:test(SELF, "--graphics False --nx1 10 --nx2 10 --testDim 1d", label="CRKSPH sum mass density test -- 1D (serial)")
#ATS:test(SELF, "--graphics False --nx1 10 --nx2 10 --testDim 2d", label="CRKSPH sum mass density test -- 2D (serial)")
#ATS:test(SELF, "--graphics False --nx1 10 --nx2 10 --testDim 3d", label="CRKSPH sum mass density test -- 3D (serial)")
#-------------------------------------------------------------------------------
# Unit test of the CRKSPH sum density algorithm.
#-------------------------------------------------------------------------------
from Spheral import *
from SpheralTestUtilities import *

title("CRKSPH sum density tests")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(
    # Parameters for seeding nodes.
    nx1 = 50,     
    nx2 = 50,
    rho1 = 1.0,
    rho2 = 1.0,
    x0 = 0.0,
    x1 = 0.5,
    x2 = 1.0,
    nPerh = 1.25,
    hmin = 0.0001, 
    hmax = 10.0,

    # Should we randomly perturb the positions?
    ranfrac = 0.2,
    seed = 14892042,

    # What test problem are we doing?
    testDim = "1d",
    testCase = "linear",

    # The fields we're going to interpolate.
    # Linear coefficients: y = y0 + m0*x
    y0 = 1.0,
    m0 = 0.0,

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
    tolerance = 1.0e-8,

    graphics = True,
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
    gen1 = GenerateNodeDistribution2d(nx1, 2*nx1, rho1,
                                      distributionType = "lattice",
                                      xmin = (x0, x0),
                                      xmax = (x1, x2),
                                      nNodePerh = nPerh,
                                      SPH = True)
    gen2 = GenerateNodeDistribution2d(nx2, 2*nx2, rho2,
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
# Initialize the density field.
#-------------------------------------------------------------------------------
rho = nodes1.massDensity()
for i in range(nodes1.numInternalNodes):
    x = nodes1.positions()[i].x
    if testCase == "linear":
        rho[i] = y0 + m0*x
    elif testCase == "quadratic":
        rho[i] = y2 + m2*x*x
    elif testCase == "step":
        if x < x1:
            rho[i] = y0
        else:
            rho[i] = 2*y0

#-------------------------------------------------------------------------------
# Update the density.
#-------------------------------------------------------------------------------
nodes1.neighbor().updateNodes()
db.updateConnectivityMap(True)
cm = db.connectivityMap()
position_fl = db.fluidPosition
mass_fl = db.fluidMass
H_fl = db.fluidHfield
rho_fl = db.fluidMassDensity

# Compute the volumes to use as weighting.
vol_fl = db.newFluidScalarFieldList(1.0, "volume")
polyvol_fl = db.newFluidFacetedVolumeFieldList(FacetedVolume(), "polyvols")
#computeHullVolumes(cm, WT.kernelExtent, position_fl, H_fl, polyvol_fl, vol_fl)
#computeHVolumes(WT.kernelExtent, H_fl, vol_fl)
for i in range(nodes1.numInternalNodes):
    vol_fl[0][i] = mass_fl[0][i]/rho_fl[0][i]
boundaries = vector_of_Boundary()
#computeCRKSPHSumMassDensity(cm, WT, position_fl, mass_fl, vol_fl, H_fl, boundaries, rho_fl)
computeHullSumMassDensity(cm, WT, position_fl, mass_fl, H_fl, rho_fl)

#-------------------------------------------------------------------------------
# Prepare the answer to check against.
#-------------------------------------------------------------------------------
xans = [position_fl(0,i).x for i in range(nodes1.numInternalNodes)]
yans = ScalarField("Mass density answer", nodes1)
for i in range(nodes1.numInternalNodes):
    if testCase == "linear":
        yans[i] = y0 + m0*xans[i]
    elif testCase == "quadratic":
        yans[i] = y2 + m2*xans[i]*xans[i]
    elif testCase == "step":
        if i < nx1:
            yans[i] = y0
        else:
            yans[i] = 2*y0

#-------------------------------------------------------------------------------
# Check our answers accuracy.
#-------------------------------------------------------------------------------
erryCRKSPH =  ScalarField("CRKSPH mass density error", nodes1)
for i in range(nodes1.numInternalNodes):
    erryCRKSPH[i] =   rho[i] - yans[i]
maxyCRKSPHerror = erryCRKSPH.max()
print("Maximum error: CRKSPH = %g" % (maxyCRKSPHerror))

#-------------------------------------------------------------------------------
# Plot the things.
#-------------------------------------------------------------------------------
if graphics:
    from SpheralGnuPlotUtilities import *
    import Gnuplot
    xans = [position_fl(0,i).x for i in range(nodes1.numInternalNodes)]

    # Interpolated values.
    ansdata = Gnuplot.Data(xans, yans.internalValues(),
                           with_ = "lines",
                           title = "Answer",
                           inline = True)
    CRKSPHdata = Gnuplot.Data(xans, rho.internalValues(),
                            with_ = "points",
                            title = "CRKSPH",
                            inline = True)
    errCRKSPHdata = Gnuplot.Data(xans, erryCRKSPH.internalValues(),
                               with_ = "points",
                               title = "CRKSPH error",
                               inline = True)

    p1 = generateNewGnuPlot()
    p1.plot(ansdata)
    p1.replot(CRKSPHdata)
    p1("set key top left")
    p1.title("CRKSPH sum density")
    p1.refresh()

    p2 = generateNewGnuPlot()
    p2.replot(errCRKSPHdata)
    p2.title("Error in sum density")
    p2.refresh()

    # p3 = plotFieldList(vol_fl,
    #                    plotStyle = "lines",
    #                    lineStyle = "linetype 0",
    #                    winTitle = "volume")

    # If we're in 2D dump a silo file too.
    if testDim == "2d":
        from SpheralVoronoiSiloDump import SpheralVoronoiSiloDump
        dumper = SpheralVoronoiSiloDump("testCRKSPHSumDensity_%s_2d" % testCase,
                                        listOfFields = [rho, yans, erryCRKSPH],
                                        listOfFieldLists = [mass_fl, H_fl])
        dumper.dump(0.0, 0)

#-------------------------------------------------------------------------------
# Check the maximum CRKSPH error and fail the test if it's out of bounds.
#-------------------------------------------------------------------------------
if maxyCRKSPHerror > tolerance:
    raise ValueError("CRKSPH error out of bounds: %g > %g" % (maxyCRKSPHerror, tolerance))
