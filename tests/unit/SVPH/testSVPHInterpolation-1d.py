#ATS:test(SELF, "--linearConsistent True --graphics False", label="SVPH interpolation test -- 1-D (serial)", svph=True)
#-------------------------------------------------------------------------------
# A set of tests to compare how different meshless methods interpolate fields.
#-------------------------------------------------------------------------------
from Spheral1d import *
from SpheralTestUtilities import *
from generateMesh import *

title("SVPH Interpolation tests")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(
            # Parameters for seeding nodes.
            nx1 = 100,
            rho1 = 1.0,
            eps1 = 0.0,
            x0 = 0.0,
            x1 = 1.0,
            nPerh = 1.25,
            hmin = 0.0001, 
            hmax = 10.0,

            # Should we randomly perturb the positions?
            ranfrac = 0.2,
            seed = 14892042,

            # What test problem are we doing?
            testCase = "linear",

            # Should we use linearly consistent formalism?
            linearConsistent = True,
            
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
            maxHIterations = 100,
            Htolerance = 1.0e-4,

            # Parameters for passing the test
            interpolationTolerance = 5.0e-7,
            derivativeTolerance = 5.0e-5,

            graphics = True,
            plotKernels = False,
            )

assert testCase in ("constant", "linear", "quadratic")

#-------------------------------------------------------------------------------
# Which test case function are we doing?
#-------------------------------------------------------------------------------
if testCase == "constant":
    def func(x):
        return y0
    def dfunc(x):
        return 0.0
elif testCase == "linear":
    def func(x):
        return y0 + m0*x
    def dfunc(x):
        return m0
else:
    def func(x):
        return y2 + m2*x*x
    def dfunc(x):
        return 2*m2*x

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
from DistributeNodes import distributeNodesInRange1d
distributeNodesInRange1d([(nodes1, nx1, rho1, (x0, x1))], nPerh = nPerh)
output("nodes1.numNodes")

# Set node specific thermal energies
nodes1.specificThermalEnergy(ScalarField("tmp", nodes1, eps1))
nodes1.massDensity(ScalarField("tmp", nodes1, rho1))

#-------------------------------------------------------------------------------
# Optionally randomly jitter the node positions.
#-------------------------------------------------------------------------------
dx = (x1 - x0)/nx1
for i in range(nodes1.numInternalNodes):
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
# Construct some boundary conditions.
#-------------------------------------------------------------------------------
xbc0 = ReflectingBoundary(Plane(Vector(x0), Vector( 1.0)))
xbc1 = ReflectingBoundary(Plane(Vector(x1), Vector(-1.0)))
bounds = vector_of_Boundary([xbc0, xbc1])

#-------------------------------------------------------------------------------
# Iterate the h to convergence if requested.
#-------------------------------------------------------------------------------
if iterateH:
    method = SPHSmoothingScale(IdealH, WT)
    pkgs = vector_of_Physics([method])
    iterateIdealH(db,
                  pkgs,
                  vector_of_Boundary(),
                  maxHIterations,
                  Htolerance)

#-------------------------------------------------------------------------------
# Initialize our field.
#-------------------------------------------------------------------------------
f = ScalarField("test field", nodes1)
for i in range(nodes1.numInternalNodes):
    x = nodes1.positions()[i].x
    if testCase == "constant":
        f[i] = y0
    elif testCase == "linear":
        f[i] = y0 + m0*x
    elif testCase == "quadratic":
        f[i] = y2 + m2*x*x

#-------------------------------------------------------------------------------
# Build the tessellation.
#-------------------------------------------------------------------------------
mesh, void = generateLineMesh([nodes1],
                              bounds,
                              generateVoid = False,
                              removeBoundaryZones = False)

#-------------------------------------------------------------------------------
# Do that interpolation baby!
#-------------------------------------------------------------------------------
fl = ScalarFieldList()
fl.appendField(f)
db.updateConnectivityMap(True)
fSVPHfl = sampleFieldListSVPH(fl,
                              db.globalPosition,
                              db.globalHfield,
                              db.connectivityMap(),
                              WT,
                              mesh,
                              linearConsistent)
dfSVPHfl = gradientFieldListSVPH(fl,
                                 db.globalPosition,
                                 db.globalHfield,
                                 db.connectivityMap(),
                                 WT,
                                 mesh,
                                 linearConsistent)
fSVPH = fSVPHfl[0]
dfSVPH = dfSVPHfl[0]

#-------------------------------------------------------------------------------
# Prepare the answer to check against.
#-------------------------------------------------------------------------------
positions = nodes1.positions()
xans = [positions[i].x for i in range(nodes1.numInternalNodes)]
yans = [func(x) for x in xans]
dyans = [dfunc(x) for x in xans]

#-------------------------------------------------------------------------------
# Check our answers accuracy.
#-------------------------------------------------------------------------------
ySVPH = fSVPH.internalValues()
dySVPH = [x.x for x in dfSVPH.internalValues()]

errySVPH = [y - z for y, z in zip(ySVPH, yans)]
maxySVPHerror = max([abs(x) for x in errySVPH])

errdySVPH = [y - z for y, z in zip(dySVPH, dyans)]
maxdySVPHerror = max([abs(x) for x in errdySVPH])

# errdySPH = [y - z for y, z in zip(dySPH, dyans)]
# maxdySPHerror = max([abs(x) for x in errdySPH])

print("Maximum errors (interpolation, gradient): SVPH = (%g, %g)" % (maxySVPHerror,
                                                                     maxdySVPHerror))

#-------------------------------------------------------------------------------
# Plot the things.
#-------------------------------------------------------------------------------
if graphics:
    from SpheralMatplotlib import *

    p1 = newFigure()
    p1.plot(xans, yans, "k-", label = "Answer")
    p1.plot(xans, ySVPH, "ro", fillstyle = "none", label = "SVPH")
    p1.set_title("Interpolated values")
    p1.legend()

    p2 = newFigure()
    p2.plot(xans, errySVPH, "ro", fillstyle = "none", label = "SVPH")
    p2.set_title("Error in interpolation")

    p3 = newFigure()
    p3.plot(xans, dyans, "k-", label = "Answer")
    p3.plot(xans, dySVPH, "ro", fillstyle = "none", label = "SVPH")
    p3.set_title("Derivative values")
    p3.legend()

    p4 = newFigure()
    p4.plot(xans, errdySVPH, "ro", fillstyle = "none", label = "SVPH")
    p4.set_title("Error in derivatives")

#-------------------------------------------------------------------------------
# Check the maximum SVPH error and fail the test if it's out of bounds.
#-------------------------------------------------------------------------------
if maxySVPHerror > interpolationTolerance:
    raise ValueError("SVPH interpolation error out of bounds: %g > %g" % (maxySVPHerror, interpolationTolerance))
if maxdySVPHerror > interpolationTolerance:
    raise ValueError("SVPH gradient error out of bounds: %g > %g" % (maxdySVPHerror, interpolationTolerance))
