#ATS:test(SELF, "--linearConsistent True --graphics False", label="SVPH interpolation test -- 2-D (serial), svph=True")
#-------------------------------------------------------------------------------
# A set of tests to compare how different meshless methods interpolate fields.
#-------------------------------------------------------------------------------
from Spheral2d import *
from SpheralTestUtilities import *
from generateMesh import *
from GenerateNodeDistribution2d import *

title("SVPH Interpolation tests")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(
    # Parameters for seeding nodes.
    nx = 20,
    ny = 20,
    rho1 = 1.0,
    eps1 = 0.0,
    x0 = 0.0,
    x1 = 1.0,
    y0 = 0.0,
    y1 = 1.0,
    nPerh = 2.0,
    hmin = 0.0001, 
    hmax = 10.0,

    # Should we randomly perturb the positions?
    distribution = "random",
    seed = 14892042,

    # What test problem are we doing?
    testCase = "linear",

    # Should we use linearly consistent formalism?
    linearConsistent = True,
    
    # The fields we're going to interpolate.
    # Function coefficients: f(x,y) = A + B*x + C*y + D*x*y + E*x*x + F*y*y
    A = 1.0,
    B = 0.5,
    C = 2.0,
    D = 1.5,
    E = 0.25,
    F = 1.0,

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

    # Should we output 
    graphics = True,
    vizFileName = "SVPHInterpolation_test_viz",
    vizDir = ".",
)

assert testCase in ("constant", "linear", "quadratic", "step")
vizFileName += "_%s_linearConsistent=%s" % (testCase, linearConsistent)

#-------------------------------------------------------------------------------
# Which test case function are we doing?
#-------------------------------------------------------------------------------
if testCase == "constant":
    def func(pos):
        return A
    def dfunc(pos):
        return Vector(0.0, 0.0)
elif testCase == "linear":
    def func(pos):
        return A + B*pos.x + C*pos.y
    def dfunc(pos):
        return Vector(B, C)
elif testCase == "quadratic":
    def func(pos):
        return A + B*pos.x + C*pos.y + D*pos.x*pos.y + E*(pos.x)**2 + F*(pos.y)**2
    def dfunc(pos):
        return Vector(B + D*pos.y + 2.0*E*pos.x,
                      C + D*pos.x + 2.0*F*pos.y)
elif testCase == "step":
    def func(pos):
        plane = Plane(Vector(0.5, 0.5), Vector(1.0, 0.5).unitVector())
        if plane.compare(pos) > 0:
            return 1.0
        else:
            return 0.0
    def dfunc(pos):
        return Vector.zero

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
gen1 = GenerateNodeDistribution2d(nx, ny, rho1,
                                  distributionType = "lattice",
                                  xmin = (x0, y0),
                                  xmax = (x1, y1),
                                  nNodePerh = nPerh)

# If we're actually using a random distribution, replace the positions
# with new random sets.
if distribution == "random":
    nquant = 2**16
    usedset = set()
    i = 0
    j = random.randint(0, 2**32)
    dx, dy = (x1 - x0)/nquant, (y1 - y0)/nquant
    while i != nx*ny:
        while j in usedset:
            j = random.randint(0, 2**32)
        gen1.x[i] = x0 + (j % nquant + 0.5)*dx
        gen1.y[i] = y0 + (j // nquant + 0.5)*dy
        usedset.add(j)
        i += 1

if mpi.procs > 1:
    from VoronoiDistributeNodes import distributeNodes2d
else:
    from DistributeNodes import distributeNodes2d

distributeNodes2d((nodes1, gen1))
output("nodes1.numNodes")

# Set node specific thermal energies
nodes1.specificThermalEnergy(ScalarField("tmp", nodes1, eps1))
nodes1.massDensity(ScalarField("tmp", nodes1, rho1))

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
bounds = vector_of_Boundary()
xbc0 = ReflectingBoundary(Plane(Vector(x0, y0), Vector( 1.0,  0.0)))
xbc1 = ReflectingBoundary(Plane(Vector(x1, y1), Vector(-1.0,  0.0)))
ybc0 = ReflectingBoundary(Plane(Vector(x0, y0), Vector( 0.0,  1.0)))
ybc1 = ReflectingBoundary(Plane(Vector(x1, y1), Vector( 0.0, -1.0)))
bounds.append(xbc0)
bounds.append(xbc1)
bounds.append(ybc0)
bounds.append(ybc1)

#-------------------------------------------------------------------------------
# Iterate the h to convergence if requested.
#-------------------------------------------------------------------------------
if iterateH:
    method = ASPHSmoothingScale()
    emptyBounds = vector_of_Boundary()
    iterateIdealH(db,
                  emptyBounds,
                  WT,
                  method,
                  maxHIterations,
                  Htolerance)

#-------------------------------------------------------------------------------
# Initialize our field.
#-------------------------------------------------------------------------------
f = ScalarField("test field", nodes1)
positions = nodes1.positions()
for i in range(nodes1.numInternalNodes):
    f[i] = func(positions[i])

#-------------------------------------------------------------------------------
# Build the tessellation.
#-------------------------------------------------------------------------------
mesh, void = generatePolygonalMesh([nodes1],
                                   bounds,
                                   generateVoid = False,
                                   removeBoundaryZones = False)

#-------------------------------------------------------------------------------
# Do that interpolation baby!
#-------------------------------------------------------------------------------
fl = ScalarFieldList()
fl.appendField(f)
db.updateConnectivityMap(True)
mass = nodes1.mass()
rho = nodes1.massDensity()
weight = db.newFluidScalarFieldList(0.0, "weight")
for i in range(nodes1.numNodes):
    weight[0][i] = mass[i]/rho[i]
fSPHfl = smoothScalarFields(fl,
                            db.globalPosition,
                            weight,
                            db.globalMass,
                            db.fluidMassDensity,
                            db.globalHfield,
                            WT)
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
fSPH = fSPHfl[0]
fSVPH = fSVPHfl[0]
dfSVPH = dfSVPHfl[0]

#-------------------------------------------------------------------------------
# Prepare the answer to check against.
#-------------------------------------------------------------------------------
fans = [func(positions[i]) for i in range(nodes1.numInternalNodes)]
dfans = [dfunc(positions[i]) for i in range(nodes1.numInternalNodes)]

#-------------------------------------------------------------------------------
# Check our answers accuracy.
#-------------------------------------------------------------------------------
errfSPH = [y - z for y, z in zip(fSPH.internalValues(), fans)]
maxfSPHerror = max([abs(x) for x in errfSPH])

errfSVPH = [y - z for y, z in zip(fSVPH.internalValues(), fans)]
maxfSVPHerror = max([abs(x) for x in errfSVPH])

errdfSVPH = [y - z for y, z in zip(dfSVPH, dfans)]
maxdfSVPHerror = max([x.magnitude() for x in errdfSVPH])

print("Maximum errors (interpolation): SPH = %g" % maxfSVPHerror)
print("Maximum errors (interpolation, gradient): SVPH = (%g, %g)" % (maxfSVPHerror,
                                                                     maxdfSVPHerror))

#-------------------------------------------------------------------------------
# Plot the things.
#-------------------------------------------------------------------------------
if graphics:
    from SpheralVoronoiSiloDump import *
    fansField = ScalarField("answer interpolation", nodes1)
    dfansField = VectorField("answer gradient", nodes1)
    errfSPHfield = ScalarField("SPH interpolation error", nodes1)
    errfSVPHfield = ScalarField("SVPH interpolation error", nodes1)
    for i in range(nodes1.numInternalNodes):
        fansField[i] = fans[i]
        dfansField[i] = dfans[i]
        errfSPHfield[i] = errfSPH[i]
        errfSVPHfield[i] = errfSVPH[i]
    d = SpheralVoronoiSiloDump(baseFileName = vizFileName,
                               baseDirectory = vizDir,
                               listOfFields = [fSPH,
                                               fSVPH,
                                               dfSVPH,
                                               fansField,
                                               dfansField],
                               boundaries = bounds)
    d.dump(0.0, 0)

    # Write a file while we're at it.
    outfile = os.path.join(vizDir, "%s_tabular.txt" % vizFileName)
    f = open(outfile, "w")
    f.write(("#" + 7*"'%s'   " + "\n") % 
            ("x", "y", "fans", "fSPH", "fSVPH", "errSPH", "errSVPH"))
    for i in range(nodes1.numInternalNodes):
        f.write((7*"%16.13e " + "\n") %
                (positions[i].x, positions[i].y, 
                 fans[i], fSPH[i], fSVPH[i], 
                 errfSPH[i], errfSVPH[i]))
    f.close()

#-------------------------------------------------------------------------------
# 1D plotting
#-------------------------------------------------------------------------------
if graphics:
    from SpheralGnuPlotUtilities import *
    import Gnuplot

    xans = [positions[i].x for i in range(nodes1.numInternalNodes)]
    yans = [fansField[i]   for i in range(nodes1.numInternalNodes)]
    dyans = [dfansField[i].x for i in range(nodes1.numInternalNodes)]
    ySVPH = fSVPH.internalValues()
    dySVPH = [x.x for x in dfSVPH.internalValues()]
    errySVPH = [y - z for y, z in zip(ySVPH, yans)]
    errdySVPH = [y - z for y, z in zip(dySVPH, dyans)]

    # Interpolated values.
    ansdata = Gnuplot.Data(xans, yans,
                           with_ = "points",
                           title = "Answer",
                           inline = True)
    SVPHdata = Gnuplot.Data(xans, ySVPH,
                            with_ = "points",
                            title = "SVPH",
                            inline = True)
    errSVPHdata = Gnuplot.Data(xans, errySVPH,
                               with_ = "points",
                               title = "SVPH",
                               inline = True)

    p1 = generateNewGnuPlot()
    p1.plot(ansdata)
    p1.replot(SVPHdata)
    p1("set key top left")
    p1.title("Interpolated values")
    p1.refresh()

    p2 = generateNewGnuPlot()
    p2.replot(errSVPHdata)
    p2.title("Error in interpolation")
    p2.refresh()

    # Derivative values.
    dansdata = Gnuplot.Data(xans, dyans,
                            with_ = "points",
                            title = "Answer",
                            inline = True)
    dSVPHdata = Gnuplot.Data(xans, dySVPH,
                             with_ = "points",
                             title = "SVPH",
                             inline = True)
    errdSVPHdata = Gnuplot.Data(xans, errdySVPH,
                                with_ = "points",
                                title = "SVPH",
                                inline = True)

    p3 = generateNewGnuPlot()
    p3.plot(dansdata)
    p3.replot(dSVPHdata)
    p3("set key top left")
    p3.title("Derivative values")
    p3.refresh()

    p4 = generateNewGnuPlot()
    p4.replot(errdSVPHdata)
    p4.title("Error in derivatives")
    p4.refresh()



#-------------------------------------------------------------------------------
# Check the maximum SVPH error and fail the test if it's out of bounds.
#-------------------------------------------------------------------------------
if maxfSVPHerror > interpolationTolerance:
    raise ValueError("SVPH interpolation error out of bounds: %g > %g" % (maxfSVPHerror, interpolationTolerance))
# if maxdfSVPHerror > interpolationTolerance:
#     raise ValueError, "SVPH gradient error out of bounds: %g > %g" % (maxdfSVPHerror, interpolationTolerance)
