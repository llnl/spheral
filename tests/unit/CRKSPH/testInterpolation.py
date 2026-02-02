#ATS:test(SELF, "--graphics False --testSPH False --nx1 10 --nx2 100 --testDim 1d --testCase linear", label="RK linear interpolation test -- 1D (serial)")
#ATS:test(SELF, "--graphics False --testSPH False --nx1 10 --nx2 20  --testDim 2d --testCase linear", label="RK linear interpolation test -- 2D (serial)")
#ATS:test(SELF, "--graphics False --testSPH False --nx1 5  --nx2 10  --testDim 3d --testCase linear", label="RK linear interpolation test -- 3D (serial)")
#ATS:test(SELF, "--graphics False --testSPH False --nx1 10 --nx2 100 --testDim 1d --testCase quadratic --correctionOrder QuadraticOrder", label="RK quadratic interpolation test -- 1D (serial)")
#ATS:test(SELF, "--graphics False --testSPH False --nx1 10 --nx2 20  --testDim 2d --testCase quadratic --correctionOrder QuadraticOrder", label="RK quadratic interpolation test -- 2D (serial)")
#ATS:test(SELF, "--graphics False --testSPH False --nx1 5  --nx2 5   --testDim 3d --testCase quadratic --correctionOrder QuadraticOrder", label="RK quadratic interpolation test -- 3D (serial)")
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
    nPerh = 4.01,
    hmin = 0.0001, 
    hmax = 1000.0,

    # What order of reproducing kernel should we use (0,1,2)?
    correctionOrder = LinearOrder,
    
    # Should we randomly perturb the positions?
    ranfrac = 0.2,
    seed = 14892042,

    # What test problem are we doing?
    testDim = "1d",
    testCase = "linear",

    # Should we compare with SPH?
    testSPH = True,

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
    interpolationTolerance = 5.0e-7,
    derivativeTolerance = 5.0e-5,

    graphics = True,
    plotKernels = False,
    outputFile = None,
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
WT = TableKernel(WendlandC4Kernel(), 1000)
output("WT")
kernelExtent = WT.kernelExtent

#-------------------------------------------------------------------------------
# Make the NodeList.
#-------------------------------------------------------------------------------
nodes1 = makeFluidNodeList("nodes1", eos,
                           hmin = hmin,
                           hmax = hmax,
                           nPerh = nPerh)
nodes2 = makeFluidNodeList("nodes2", eos,
                           hmin = hmin,
                           hmax = hmax,
                           nPerh = nPerh)
nodeSet = [nodes1, nodes2]
for nodes in nodeSet:
    output("nodes")
    output("nodes.hmin")
    output("nodes.hmax")
    output("nodes.nodesPerSmoothingScale")

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
if testDim == "1d":
    from DistributeNodes import distributeNodesInRange1d
    distributeNodesInRange1d([(nodes1, [(nx1, rho1, (x0, x1))])], nPerh = nPerh)
    distributeNodesInRange1d([(nodes2, [(nx2, rho2, (x1, x2))])], nPerh = nPerh)
elif testDim == "2d":
    from DistributeNodes import distributeNodes2d
    from GenerateNodeDistribution2d import GenerateNodeDistribution2d
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
    distributeNodes2d((nodes1, gen1),
                      (nodes2, gen2))

elif testDim == "3d":
    from DistributeNodes import distributeNodes3d
    from GenerateNodeDistribution3d import GenerateNodeDistribution3d
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
    distributeNodes3d((nodes1, gen1),
                      (nodes2, gen2))

else:
    raise ValueError("Only tests cases for 1d,2d and 3d.") 

for nodes in nodeSet:
    output("nodes.name, nodes.numNodes")

# Set node properties.
for nodes, eps0 in ((nodes1, eps1),
                    (nodes2, eps2)):
    eps = nodes.specificThermalEnergy()
    for i in range(nodes.numInternalNodes):
        eps[i] = eps0

#-------------------------------------------------------------------------------
# Optionally randomly jitter the node positions.
#-------------------------------------------------------------------------------
dx1 = (x1 - x0)/nx1
dx2 = (x2 - x1)/nx2
dy = (x2 - x0)/(nx1 + nx2)
dz = (x2 - x0)/(nx1 + nx2)
for nodes, dx in ((nodes1, dx1),
                  (nodes2, dx2)):
    pos = nodes.positions()
    for i in range(nodes.numInternalNodes):
        if testDim == "1d":
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
for nodes in nodeSet:
    db.appendNodeList(nodes)
output("db")
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
# Initialize our field.
#-------------------------------------------------------------------------------
f = db.newFluidScalarFieldList(name="test field")
pos = db.fluidPosition
for iNodeList, nodes in enumerate(db.nodeLists):
    for i in range(nodes.numInternalNodes):
        x = pos(iNodeList, i).x
        if testCase == "linear":
            f[iNodeList][i] = y0 + m0*x
        elif testCase == "quadratic":
            f[iNodeList][i] = y2 + m2*x*x
        elif testCase == "step":
            if x < x1:
                f[iNodeList][i] = y0
            else:
                f[iNodeList][i] = 2*y0

#-------------------------------------------------------------------------------
# Prepare variables to accumulate the test values.
#-------------------------------------------------------------------------------
fSPH = db.newFluidScalarFieldList(name="SPH interpolated values")
dfSPH = db.newFluidVectorFieldList(name="SPH derivative values")

A = db.newFluidScalarFieldList(name="A")
B = db.newFluidVectorFieldList(name="B")
C = db.newFluidTensorFieldList(name="C")
gradA = db.newFluidVectorFieldList(name="gradA")
gradB = db.newFluidTensorFieldList(name="gradB")
gradC = db.newFluidThirdRankTensorFieldList(name="gradB")

M0 = db.newFluidScalarFieldList(name="M0")
M1 = db.newFluidVectorFieldList(name="M1")
M2 = db.newFluidSymTensorFieldList(name="M2")
M3 = db.newFluidThirdRankTensorFieldList(name="M3")
M4 = db.newFluidFourthRankTensorFieldList(name="M4")
gradM0 = db.newFluidVectorFieldList(name="grad M0")
gradM1 = db.newFluidTensorFieldList(name="grad M1")
gradM2 = db.newFluidThirdRankTensorFieldList(name="grad M2")
gradM3 = db.newFluidFourthRankTensorFieldList(name="grad M3")
gradM4 = db.newFluidFifthRankTensorFieldList(name="grad M4")

surfacePoint = db.newFluidIntFieldList(name="surface point")

db.updateConnectivityMap(True)
cm = db.connectivityMap()
position = db.fluidPosition
weight = db.fluidMass
weight /= db.fluidMassDensity
H = db.fluidHfield

# Compute the volumes to use as weighting.
#polyvol = db.newFluidFacetedVolumeFieldList(name=FacetedVolume(), "polyvols")
#weight = db.newFluidScalarFieldList(name=1.0, "volume")
#computeHullVolumes(cm, position, polyvol, weight)

computeCRKSPHMoments(cm, WT, weight, position, H, correctionOrder, NodeCoupling(),
                     M0, M1, M2, M3, M4, gradM0, gradM1, gradM2, gradM3, gradM4)
computeCRKSPHCorrections(M0, M1, M2, M3, M4, gradM0, gradM1, gradM2, gradM3, gradM4, H,
                         surfacePoint,
                         correctionOrder,
                         A, B, C, gradA, gradB, gradC)

#-------------------------------------------------------------------------------
# Measure the interpolated values and gradients.
#-------------------------------------------------------------------------------
if testSPH:
    for iNodeList, nodes in enumerate(db.nodeLists):
        for i in range(nodes.numInternalNodes):
            ri = position(iNodeList, i)
            Hi = H(iNodeList, i)
            Hdeti = Hi.Determinant()
            wi = weight(iNodeList, i)
            fi = f(iNodeList, i)

            # Self contribution.
            W0 = WT.kernelValue(0.0, Hdeti)
            fSPH[iNodeList][i] = wi*W0 * fi

            # Go over them neighbors.
            allneighbors = cm.connectivityForNode(iNodeList, i)
            for jNodeList, neighbors in enumerate(allneighbors):
                for j in neighbors:
                    rj = position(jNodeList, j)
                    Hj = H(jNodeList, j)
                    Hdetj = Hj.Determinant()
                    wj = weight(jNodeList, j)
                    fj = f(jNodeList, j)

                    # The standard SPH kernel and it's gradient.
                    rij = ri - rj
                    etai = Hi*rij
                    etaj = Hj*rij
                    Wj = WT.kernelValue(etaj.magnitude(), Hdetj)
                    gradWj = Hj*etaj.unitVector() * WT.gradValue(etaj.magnitude(), Hdetj)

                    # Increment our interpolated values.
                    fSPH[iNodeList][i] += fj * wj*Wj

                    # Increment the derivatives.
                    dfSPH[iNodeList][i] += fj * wj*gradWj

#-------------------------------------------------------------------------------
# Check the C++ interpolation and gradient methods.
#-------------------------------------------------------------------------------
fRK = interpolateCRKSPH(f, position, weight, H, A, B, C,
                            cm, correctionOrder, WT)
dfRK = gradientCRKSPH(f, position, weight, H,
                          A, B, C, gradA, gradB, gradC,
                          cm, correctionOrder, WT)

#-------------------------------------------------------------------------------
# Prepare the answer to check against.
#-------------------------------------------------------------------------------
yans = db.newFluidScalarFieldList(name="interpolation answer")
dyans = db.newFluidScalarFieldList(name="derivative answer")
for iNodeList in range(db.numNodeLists):
    n = yans[iNodeList].numInternalElements
    for i in range(n):
        xi = position(iNodeList, i).x
        if testCase == "linear":
            yans[iNodeList][i] = y0 + m0*xi
            dyans[iNodeList][i] = m0
        elif testCase == "quadratic":
            yans[iNodeList][i] = y2 + m2*xi*xi
            dyans[iNodeList][i] = 2*m2*xi
        elif testCase == "step":
            if iNodeList == 0:
                yans[iNodeList][i] = y0
            else:
                yans[iNodeList][i] = 2*y0
            dyans[iNodeList][i] = 0.0

#-------------------------------------------------------------------------------
# Check our answers accuracy.
#-------------------------------------------------------------------------------
def flattenFieldList(fl):
    result = []
    for f in fl:
        result += list(f.internalValues())
    return result

errySPH = flattenFieldList(fSPH - yans)
erryRK = flattenFieldList(fRK - yans)

errdySPH = []
errdyRK = []
for iNodeList in range(db.numNodeLists):
    n = fSPH[iNodeList].numInternalElements
    for i in range(n):
        errdySPH.append(dfSPH(iNodeList, i).x - dyans(iNodeList, i))
        errdyRK.append(dfRK(iNodeList, i).x - dyans(iNodeList, i))

maxySPHerror = max([abs(x) for x in errySPH])
maxdySPHerror = max([abs(x) for x in errdySPH])
maxyRKerror = max([abs(x) for x in erryRK])
maxdyRKerror = max([abs(x) for x in errdyRK])

print("Maximum errors (interpolation): SPH = %g, RK = %g" % (maxySPHerror, maxyRKerror))
print("Maximum errors   (derivatives): SPH = %g, RK = %g" % (maxdySPHerror, maxdyRKerror))

# Output timing tables.
Timer.TimerSummary()

#-------------------------------------------------------------------------------
# Plot the things.
#-------------------------------------------------------------------------------
if graphics:
    from SpheralMatplotlib import *

    xans = [x.x for x in flattenFieldList(position)]

    # Interpolated values.
    p1 = plotFieldList(fRK,
                       plotStyle = "g*",
                       lineTitle = "RK",
                       winTitle = "Interpolated values")
    if testSPH:
        plotFieldList(fSPH,
                      plotStyle = "r+",
                      lineTitle = "SPH",
                      plot = p1)
    plotFieldList(yans,
                  plotStyle = "k-",
                  lineTitle = "Answer",
                  plot = p1)

    # Interpolation error
    p2 = newFigure()
    p2.plot(xans, erryRK, "g*",
            label = "RK")
    if testSPH:
        p2.plot(xans, errySPH, "r+",
                label = "SPH")
    p2.axes.legend()
    plt.title("Error in interpolation")

    # Derivative values.
    p3 = plotFieldList(dfRK,
                       yFunction = "%s.x",
                       plotStyle = "g*",
                       lineTitle = "RK",
                       winTitle = "Derivative values")
    if testSPH:
        plotFieldList(dfSPH,
                      yFunction = "%s.x",
                      plotStyle = "r+",
                      lineTitle = "SPH",
                      plot = p3)
    plotFieldList(dyans,
                  plotStyle = "k-",
                  lineTitle = "Answer",
                  plot = p3)

    # Derivative error
    p4 = newFigure()
    p4.plot(xans, errdyRK, "g*",
            label = "RK")
    if testSPH:
        p4.plot(xans, errdySPH, "r+",
                label = "SPH")
    p4.axes.legend()
    plt.title("Error in derivatives")

    # Plot the kernel shapes as appropriate.
    if testDim == "1d":
        p7 = newFigure()
        j = -2 # int(nodes1.numInternalNodes/2)
        Hj = H[1][j]
        hj = 1.0/Hj.xx
        Hdetj = Hj.Determinant()
        Aj = A[1][j]
        Bj = B[1][j].x
        Cj = C[1][j].xx
        nsamp = 100
        dx = 4.0/nsamp
        xvals = [i*dx - 2.0 for i in range(nsamp)]
        W = [WT.kernelValue(abs(xi), Hdetj) for xi in xvals]
        WR = [Wi*Aj*(1.0 + Bj*(xi)*hj+Cj*(xi)*(xi)*hj*hj) for xi, Wi in zip(xvals, W)]
        p7.plot(xvals, W, "r-", label="SPH")
        p7.plot(xvals, WR, "g-", label="RK")
        p7.axes.legend()
        plt.title("Kernel")
        if outputFile:
            f = open("Kernel_" + outputFile, "w")
            f.write(("#" + 3*' "%20s"' + "\n") % ("eta", "Wj", "WRj"))
            for xi, Wi, WRi in zip(xvals, W, WR):
                f.write((3*" %20g" + "\n") % (xi, Wi, WRi))
            f.close()

    # We may want a gnu/pdv style text file.
    if outputFile and testDim == "2d":
        of = open(outputFile, "w")
        of.write(('#' + 7*' "%20s"' + '\n') % ("x", "interp answer", "grad answer", "interp SPH", "interp CRK", "grad SPH", "grad CRK"))
        for iNodeList, nodes in enumerate(db.nodeLists):
            for i in range(nodes.numInternalNodes):
                of.write((7*" %20g" + "\n") %
                         (position(iNodeList,i), yans(iNodeList,i), dyans(iNodeList,i), fSPH(iNodeList,i), fRK(iNodeList,i), dfSPH(iNodeList,i).x, dfRK(iNodeList,i).x))
        of.close()

    # If we're in 2D/3D dump a silo file too.
    if testDim != "1d":
        from siloPointmeshDump import siloPointmeshDump
        siloPointmeshDump("testInterpolation_%s_%s" % (testCase, testDim),
                          fieldLists = [fSPH, fRK, dfSPH, dfRK,
                                        yans, dyans,
                                        weight, H, A, B, gradA, gradB,
                                        M0, M1, M2])

if plotKernels:
    import Gnuplot
    pk = generateNewGnuPlot()
    for iNodeList, nodes in enumerate(db.nodeLists):
        for i in range(nodes.numInternalNodes):
            xi = positions(iNodeList,i).x
            Hi = H(iNodeList,i)
            Hdeti = Hi.Determinant()
            hi = 1.0/Hi.xx
            Ai = A(iNodeList,i)
            Bi = B(iNodeList,i)
            Ci = C(iNodeList,i)

            dx = 2.0*kernelExtent*hi/50
            x = [xi - kernelExtent*hi + (i + 0.5)*dx for i in range(50)]
            #y = [Ai*(1.0 + Bi.x*(xi - xj))*WT.kernelValue(abs(xi - xj)/hi, Hdeti) for xj in x]
            y = [Ai*(1.0 + Bi.x*(xi - xj)+Ci.xx*(xi-xj)*(xi-xj))*WT.kernelValue(abs(xi - xj)/hi, Hdeti) for xj in x]
            d = Gnuplot.Data(x, y, with_="lines", inline=True)
            pk.replot(d)

#-------------------------------------------------------------------------------
# Check the maximum RK error and fail the test if it's out of bounds.
#-------------------------------------------------------------------------------
if maxyRKerror > interpolationTolerance:
    raise ValueError("RK interpolation error out of bounds: %g > %g" % (maxyRKerror, interpolationTolerance))

if maxdyRKerror > derivativeTolerance:
    raise ValueError("RK derivative error out of bounds: %g > %g" % (maxdyRKerror, derivativeTolerance))
