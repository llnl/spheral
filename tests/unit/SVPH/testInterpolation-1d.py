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
            linearConsistent = False,
            
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
bounds = vector_of_Boundary()
xbc0 = ReflectingBoundary(Plane(Vector(x0), Vector( 1.0)))
xbc1 = ReflectingBoundary(Plane(Vector(x1), Vector(-1.0)))
bounds.append(xbc0)
bounds.append(xbc1)

#-------------------------------------------------------------------------------
# Iterate the h to convergence if requested.
#-------------------------------------------------------------------------------
if iterateH:
    method = SPHSmoothingScale()
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
# Prepare variables to accumulate the test values.
#-------------------------------------------------------------------------------
fSPH = ScalarField("SPH interpolated values", nodes1)
dfSPH = VectorField("SPH derivative values", nodes1)
fSVPH = ScalarField("SVPH interpolated values", nodes1)
dfSVPH = VectorField("SVPH derivative values", nodes1)
norm = ScalarField("SVPH normalization", nodes1)

positions = nodes1.positions()
mass = nodes1.mass()
H = nodes1.Hfield()

# Prepare the connectivity
db.updateConnectivityMap(False)
cm = db.connectivityMap()

#-------------------------------------------------------------------------------
# If we're using linearly consistent FVSPH, we need the simple finite-difference
# estimate of the slope for each point.
#-------------------------------------------------------------------------------
G = ScalarField("test field FD linear gradient", nodes1)
if linearConsistent:
    for i in range(nodes1.numInternalNodes):
        ri = positions[i]
        zi = mesh.zone(i)
        normi = 0.0
        for fid in zi.faceIDs:
            if fid < 0:
                fid = ~fid
            face = mesh.face(fid)
            j = mesh.face(fid).oppositeZoneID(i)
            if j < 0:
                j = ~j
            if j >= 0 and j < nx1:
                rj = positions[j]
                xji = (rj - ri).x
                assert abs(xji) > 0.0
                normi += 0.5*abs(xji)
                G[i] += 0.5*abs(xji) * (f[j] - f[i])/xji
        assert normi > 0.0
        G[i] /= normi

#-------------------------------------------------------------------------------
# Measure the interpolated values and gradients.
#-------------------------------------------------------------------------------
for i in range(nodes1.numInternalNodes):
    ri = positions[i]
    Hi = H[i]
    Hdeti = H[i].Determinant()
    mi = mass[i]
    Vi = mesh.zone(i).volume
    fi = f[i]
    Gi = G[i]

    # Self contribution.
    W0 = WT.kernelValue(0.0, Hdeti)
    fSPH[i] = mi*W0 * fi
    fSVPH[i] = Vi*W0 * fi
    norm[i] = Vi*W0
    if linearConsistent:
        dfSVPH[i] += Vector(Vi*W0 * Gi)

    # Go over them neighbors.
    neighbors = cm.connectivityForNode(nodes1, i)
    assert len(neighbors) == 1
    for j in neighbors[0]:
        rj = positions[j]
        Hj = H[j]
        Hdetj = H[j].Determinant()
        mj = mass[j]
        Vj = mesh.zone(j).volume
        fj = f[j]
        Gj = G[j]

        # The standard SPH kernel and it's gradient.
        rij = ri - rj
        etaj = Hj*rij
        Wj = WT.kernelValue(etaj.magnitude(), Hdetj)
        gradWj = Hj*etaj.unitVector() * WT.gradValue(etaj.magnitude(), Hdetj)

        # Increment the result.
        fSPH[i] += mj*Wj * fj
        dfSPH[i] += mj*gradWj * fj
        norm[i] += Vj*Wj
        if linearConsistent:
            fSVPH[i] += Vj*Wj * (fj + Gj*rij.x)
            dfSVPH[i] += Vj*((fj - fi + Gj*rij.x)*gradWj + Vector(Gj*Wj))
        else:
            fSVPH[i] += Vj*Wj * fj
            dfSVPH[i] += Vj*gradWj * (fj - fi)

    # Finalize the SVPH values.
    assert norm[i] > 0.0
    fSVPH[i] /= norm[i]
    dfSVPH[i] /= norm[i]

#-------------------------------------------------------------------------------
# Prepare the answer to check against.
#-------------------------------------------------------------------------------
xans = [positions[i].x for i in range(nodes1.numInternalNodes)]
yans = [func(x) for x in xans]
dyans = [dfunc(x) for x in xans]

#-------------------------------------------------------------------------------
# Check our answers accuracy.
#-------------------------------------------------------------------------------
ySPH = fSPH.internalValues()
ySVPH = fSVPH.internalValues()

dySPH = [x.x for x in dfSPH.internalValues()]
dySVPH = [x.x for x in dfSVPH.internalValues()]

errySPH = [y - z for y, z in zip(ySPH, yans)]
errySVPH = [y - z for y, z in zip(ySVPH, yans)]
maxySPHerror = max([abs(x) for x in errySPH])
maxySVPHerror = max([abs(x) for x in errySVPH])

errdySPH = [y - z for y, z in zip(dySPH, dyans)]
errdySVPH = [y - z for y, z in zip(dySVPH, dyans)]
maxdySPHerror = max([abs(x) for x in errdySPH])
maxdySVPHerror = max([abs(x) for x in errdySVPH])

print("Maximum errors (interpolation): SPH = %g, SVPH = %g" % (maxySPHerror, maxySVPHerror))
print("Maximum errors   (derivatives): SPH = %g, SVPH = %g" % (maxdySPHerror, maxdySVPHerror))

#-------------------------------------------------------------------------------
# Plot the things.
#-------------------------------------------------------------------------------
if graphics:
    from SpheralGnuPlotUtilities import *
    import Gnuplot

    # Interpolated values.
    ansdata = Gnuplot.Data(xans, yans,
                           with_ = "lines",
                           title = "Answer",
                           inline = True)
    SPHdata = Gnuplot.Data(xans, ySPH,
                           with_ = "points",
                           title = "SPH",
                           inline = True)
    SVPHdata = Gnuplot.Data(xans, ySVPH,
                            with_ = "points",
                            title = "SVPH",
                            inline = True)
    errSPHdata = Gnuplot.Data(xans, errySPH,
                              with_ = "points",
                              title = "SPH",
                              inline = True)
    errSVPHdata = Gnuplot.Data(xans, errySVPH,
                               with_ = "points",
                               title = "SVPH",
                               inline = True)

    p1 = generateNewGnuPlot()
    p1.plot(ansdata)
    p1.replot(SPHdata)
    p1.replot(SVPHdata)
    p1("set key top left")
    p1.title("Interpolated values")
    p1.refresh()

    p2 = generateNewGnuPlot()
    p2.plot(errSPHdata)
    p2.replot(errSVPHdata)
    p2.title("Error in interpolation")
    p2.refresh()

    # Derivative values.
    dansdata = Gnuplot.Data(xans, dyans,
                            with_ = "lines",
                            title = "Answer",
                            inline = True)
    dSPHdata = Gnuplot.Data(xans, dySPH,
                            with_ = "points",
                            title = "SPH",
                            inline = True)
    dSVPHdata = Gnuplot.Data(xans, dySVPH,
                             with_ = "points",
                             title = "SVPH",
                             inline = True)
    errdSPHdata = Gnuplot.Data(xans, errdySPH,
                               with_ = "points",
                               title = "SPH",
                              inline = True)
    errdSVPHdata = Gnuplot.Data(xans, errdySVPH,
                                with_ = "points",
                                title = "SVPH",
                                inline = True)

    p3 = generateNewGnuPlot()
    p3.plot(dansdata)
    p3.replot(dSPHdata)
    p3.replot(dSVPHdata)
    p3("set key top left")
    p3.title("Derivative values")
    p3.refresh()

    p4 = generateNewGnuPlot()
    p4.plot(errdSPHdata)
    p4.replot(errdSVPHdata)
    p4.title("Error in derivatives")
    p4.refresh()

if plotKernels:
    import Gnuplot
    pk = generateNewGnuPlot()
    for i in range(nodes1.numInternalNodes):
        xi = positions[i].x
        Hi = H[i]
        Hdeti = Hi.Determinant()
        hi = 1.0/Hi.xx
        Ai = A[i]
        Bi = B[i]

        dx = 2.0*kernelExtent*hi/50
        x = [xi - kernelExtent*hi + (i + 0.5)*dx for i in range(50)]
        y = [Ai*(1.0 + Bi.x*(xi - xj))*WT.kernelValue(abs(xi - xj)/hi, Hdeti) for xj in x]
        d = Gnuplot.Data(x, y, with_="lines", inline=True)
        pk.replot(d)

#-------------------------------------------------------------------------------
# Check the maximum SVPH error and fail the test if it's out of bounds.
#-------------------------------------------------------------------------------
if maxySVPHerror > interpolationTolerance:
    raise ValueError("SVPH interpolation error out of bounds: %g > %g" % (maxySVPHerror, interpolationTolerance))

if maxdySVPHerror > derivativeTolerance:
    raise ValueError("SVPH derivative error out of bounds: %g > %g" % (maxdySVPHerror, derivativeTolerance))
