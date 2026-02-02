#-------------------------------------------------------------------------------
# The Hydrostatic Equilibrium/Surface Tension Test
#-------------------------------------------------------------------------------
import shutil
from math import *
from Spheral2d import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from findLastRestart import *
from GenerateNodeDistribution2d import *
from CubicNodeGenerator import GenerateSquareNodeDistribution
from CentroidalVoronoiRelaxation import *

import mpi
import DistributeNodes

title("2-D integrated hydro test --  Hydrostatic Equilibrium/Surface Tension Test")

#-------------------------------------------------------------------------------
# A rejecter for making elliptical regions.
#-------------------------------------------------------------------------------
class EllipticalRejecter:
    def __init__(self, a, b, reverse):
        self.a = a
        self.b = b
        self.reverse = reverse
        return
    def __call__(self, x0, y0, m0, H0):
        x, y, m, H = [], [], [], []
        for i in range(len(x0)):
            xi = x0[i]
            yi = y0[i]
            test = (sqrt(((xi - 0.5)/a)**2 + ((yi - 0.5)/b)**2) <= 1.0)
            if self.reverse:
                test = not test
            if test:
                x.append(xi)
                y.append(yi)
                m.append(m0[i])
                H.append(H0[i])
        return x, y, m, H

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(
    # Should we do the box or ellipse?
    problem = "box",

    # Outer state.
    rho1 = 1.0,
    P1 = 1.0,
    gamma1 = 1.5,

    # Inner state
    rho2 = 4.0,
    P2 = 1.0,
    gamma2 = 1.5,

    # Geometry (box)
    x0 = 0.0,
    x1 = 0.25,
    x2 = 0.75,
    x3 = 1.0,
    y0 = 0.0,
    y1 = 0.25,
    y2 = 0.75,
    y3 = 1.0,

    # Geometry (ellipse)
    a = 0.25,   # axis length in x
    b = 0.125,  # axis length in y

    # Translation
    velx=0.0,
    vely=0.0,

    # Resolution and node seeding.
    nx1 = 100,
    ny1 = 100,

    nx2 = 50,
    ny2 = 50,

    # kernel
    HUpdate = IdealH,
    nPerh = 1.35,
    KernelConstructor = NBSplineKernel,
    order = 5,
    hmin = 1e-5,
    hmax = 0.5,
    hminratio = 0.1,

    # hydro 
    svph = False,
    crksph = False,
    fsisph = False,
    gsph = False,
    psph = False,

    # hydro parameters 
    asph = False,
    solid = False,
    filter = 0.0,  
    evolveTotalEnergy = False,
    XSPH = False,
    epsilonTensile = 0.0,
    nTensile = 8,
    densityUpdate = RigorousSumDensity,
    compatibleEnergy = True,
    gradhCorrection = True,
    correctVelocityGradient = True,

    # CRKSPH - specific parameters
    correctionOrder = LinearOrder,

    # fsisph-specific parameters
    fsiSurfaceCoefficient = 0.00,                # magnificaiton factor for interface face
    fsiRhoStabilizeCoeff = 0.1,                  # coeff for HLLC-pressure term style diffusion
    fsiXSPHCoeff=0.00,                           # xsph coeff
    fsiEpsDiffuseCoeff=0.1,                      # coeff second order artificial conduction
    fsiInterfaceMethod = HLLCInterface,          # how do we handle multimat? (HLLCInterface,ModulusInterface,NoInterface)
    fsiKernelMethod = AverageInterfaceKernels,   # should we avg the kernels? (AlwaysAverageKernels,NeverAverageKernels,AverageInterfaceKernels)
    fsiSumDensity = True,                        # apply sum density to air nodes

    # gsph-specific parameters
    gsphEpsDiffuseCoeff = 0.0,                   # 1st order artificial conduction coeff
    gsphLinearCorrect = True,                    # True - second order False - first order HLLC
    gsphGradientMethod = MixedMethodGradient,    # gradient def (RiemannGradient,HydroAccelerationGradient,SPHGradient,MixedMethodGradient)
    
    # SVPH-specific parameters
    linearConsistent = False,
    fcentroidal = 0.0,
    fcellPressure = 0.0,

    # Artificial Viscosity
    Cl = None,
    Cq = None,
    boolReduceViscosity = False,
    nh = 5.0,
    aMin = 0.1,
    aMax = 2.0,
    boolCullenViscosity = False,
    alphMax = 2.0,
    alphMin = 0.02,
    betaC = 0.7,
    betaD = 0.05,
    betaE = 1.0,
    fKern = 1.0/3.0,
    boolHopkinsCorrection = True,
    Qlimiter = None,
    balsaraCorrection = None,
    epsilon2 = None,
    
    # integrator
    cfl = 0.5,
    IntegratorConstructor = CheapSynchronousRK2Integrator,
    goalTime = 7.0,
    steps = None,
    vizCycle = None,
    vizTime = 0.1,
    dt = 0.0001,
    dtMin = 1.0e-5, 
    dtMax = 0.1,
    dtGrowth = 2.0,
    maxSteps = None,
    statsStep = 10,
    domainIndependent = False,
    rigorousBoundaries = False,
    dtverbose = False,
    redistributeStep = 1000,

    # output
    serialDump = False,
    clearDirectories = False,
    restoreCycle = -1,
    restartStep = 200,
    dataDir = "dumps-boxtension-xy",
    )

assert not svph 
assert not(boolReduceViscosity and boolCullenViscosity)
assert problem.lower() in ("box", "ellipse")

# Build our directory paths.
if svph:
    hydroname = "SVPH"
elif crksph:
    hydroname = "CRKSPH"
elif gsph:
    hydroname = "GSPH"
elif fsisph:
    hydroname = "FSISPH"
elif psph:
    hydroname = "PSPH"
else:
    hydroname = "SPH"
densityUpdateLabel = {IntegrateDensity : "IntegrateDensity",
                      SumDensity : "SumDensity",
                      RigorousSumDensity : "RigorousSumDensity",
                      SumVoronoiCellDensity : "SumVoronoiCellDensity"}
baseDir = os.path.join(dataDir,
                       problem,
                       hydroname,
                       densityUpdateLabel[densityUpdate],
                       "compatibleEnergy=%s" % compatibleEnergy,
                       "XSPH=%s" % XSPH,
                       "nPerh=%3.1f" % nPerh,
                       "fcentroidal=%1.3f" % fcentroidal,
                       "fcellPressure=%1.3f" % fcellPressure,
                       "%ix%i" % (nx1 + nx2, ny1 + ny2))
restartDir = os.path.join(baseDir, "restarts")
restartBaseName = os.path.join(restartDir, "boxtension-xy-%ix%i" % (nx1 + nx2, ny1 + ny2))

vizDir = os.path.join(baseDir, "visit")
if vizTime is None and vizCycle is None:
    vizBaseName = None
else:
    vizBaseName = "boxtension-xy-%ix%i" % (nx1 + nx2, ny1 + ny2)

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
import os, sys
if mpi.rank == 0:
    if clearDirectories and os.path.exists(baseDir):
        shutil.rmtree(baseDir)
    if not os.path.exists(restartDir):
        os.makedirs(restartDir)
    if not os.path.exists(vizDir):
        os.makedirs(vizDir)
mpi.barrier()

#-------------------------------------------------------------------------------
# Material properties.
#-------------------------------------------------------------------------------
mu = 1.0
eos1 = GammaLawGasMKS(gamma1, mu)
eos2 = GammaLawGasMKS(gamma2, mu)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
if KernelConstructor==NBSplineKernel:
  WBase = NBSplineKernel(order)
else:
  WBase = KernelConstructor()
WT = TableKernel(WBase,1000)
WTPi = WT
output("WT")
kernelExtent = WT.kernelExtent

#-------------------------------------------------------------------------------
# Make the NodeLists.
#-------------------------------------------------------------------------------
if solid:
    nodeListConstructor = makeSolidNodeList
else:
    nodeListConstructor = makeFluidNodeList
outerNodes = nodeListConstructor("outer", eos1,
                               hmin = hmin,
                               hmax = hmax,
                               hminratio = hminratio,
                               kernelExtent = kernelExtent,
                               nPerh = nPerh)
innerNodes = nodeListConstructor("inner", eos2,
                               hmin = hmin,
                               hmax = hmax,
                               hminratio = hminratio,
                               kernelExtent = kernelExtent,
                               nPerh = nPerh)
nodeSet = (outerNodes, innerNodes)
for nodes in nodeSet:
    output("nodes.name")
    output("    nodes.hmin")
    output("    nodes.hmax")
    output("    nodes.hminratio")
    output("    nodes.nodesPerSmoothingScale")
del nodes

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
if problem.lower() == "box":
    generatorOuter = GenerateNodeDistribution2d(nx1, ny1, rho1,
                                                distributionType = "lattice",
                                                xmin = (x0, y0),
                                                xmax = (x3, y3),
                                                xminreject = (x1, y1),
                                                xmaxreject = (x2, y2),
                                                nNodePerh = nPerh,
                                                SPH = not ASPH,
                                                reversereject = True)
    generatorInner = GenerateNodeDistribution2d(nx2, ny2, rho2,
                                                distributionType = "lattice",
                                                xmin = (x1, y1),
                                                xmax = (x2, y2),
                                                nNodePerh = nPerh,
                                                SPH = not ASPH)
else:
    generatorOuter = GenerateNodeDistribution2d(nx1, ny1, rho1,
                                                distributionType = "lattice",
                                                xmin = (x0, y0),
                                                xmax = (x3, y3),
                                                nNodePerh = nPerh,
                                                rejecter = EllipticalRejecter(a, b, True),
                                                SPH = not ASPH)
    generatorInner = GenerateNodeDistribution2d(nx2, ny2, rho2,
                                                distributionType = "lattice",
                                                xmin = (x1, y1),
                                                xmax = (x2, y2),
                                                nNodePerh = nPerh,
                                                rejecter = EllipticalRejecter(a, b, False),
                                                SPH = not ASPH)
if mpi.procs > 1:
    from VoronoiDistributeNodes import distributeNodes2d
else:
    from DistributeNodes import distributeNodes2d

distributeNodes2d((outerNodes, generatorOuter),
                  (innerNodes, generatorInner))
for nodes in nodeSet:
    print(nodes.name, ":")
    output("    mpi.reduce(nodes.numInternalNodes, mpi.MIN)")
    output("    mpi.reduce(nodes.numInternalNodes, mpi.MAX)")
    output("    mpi.reduce(nodes.numInternalNodes, mpi.SUM)")
del nodes

# Set node specific thermal energies
for (nodes, gamma, rho, P) in ((outerNodes, gamma1, rho1, P1),
                               (innerNodes, gamma2, rho2, P2)):
    eps0 = P/((gamma - 1.0)*rho)
    nodes.specificThermalEnergy(ScalarField("tmp", nodes, eps0))
del nodes

vel = outerNodes.velocity()
for i in range(outerNodes.numInternalNodes):
    vel[i]=Vector(velx,vely)
vel = innerNodes.velocity()
for i in range(innerNodes.numInternalNodes):
    vel[i]=Vector(velx,vely)

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node lists
#-------------------------------------------------------------------------------
db = DataBase()
output("db")
for nodes in nodeSet:
    db.appendNodeList(nodes)
del nodes
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
if svph:
    hydro = SVPH(dataBase = db,
                 W = WT, 
                 cfl = cfl,
                 compatibleEnergyEvolution = compatibleEnergy,
                 densityUpdate = densityUpdate,
                 XSVPH = XSPH,
                 linearConsistent = linearConsistent,
                 generateVoid = False,
                 HUpdate = HUpdate,
                 fcentroidal = fcentroidal,
                 fcellPressure = fcellPressure,
                 xmin = Vector(x0 - (x3 - x0), y0 - (y3 - y0)),
                 xmax = Vector(x3 + (x3 - x0), y3 + (y3 - y0)))
elif crksph:
    hydro = CRKSPH(dataBase = db,
                   filter = filter,
                   epsTensile = epsilonTensile,
                   nTensile = nTensile,
                   cfl = cfl,
                   compatibleEnergyEvolution = compatibleEnergy,
                   XSPH = XSPH,
                   order = correctionOrder,
                   densityUpdate = densityUpdate,
                   HUpdate = HUpdate)
elif psph:
    hydro = PSPH(dataBase = db,
                 W = WT,
                 filter = filter,
                 cfl = cfl,
                 compatibleEnergyEvolution = compatibleEnergy,
                 evolveTotalEnergy = evolveTotalEnergy,
                 correctVelocityGradient = correctVelocityGradient,
                 densityUpdate = densityUpdate,
                 HUpdate = HUpdate,
                 XSPH = XSPH)
elif fsisph: 
    sumDensityNodeLists = []
    slidePairs = []
    if fsiSumDensity:
        sumDensityNodeLists = [outerNodes,innerNodes]
    hydro = FSISPH(dataBase = db,
                W = WT,
                cfl = cfl,
                surfaceForceCoefficient = fsiSurfaceCoefficient,                   
                densityStabilizationCoefficient = fsiRhoStabilizeCoeff,         
                specificThermalEnergyDiffusionCoefficient = fsiEpsDiffuseCoeff,  
                xsphCoefficient = fsiXSPHCoeff,
                interfaceMethod = fsiInterfaceMethod,  
                kernelAveragingMethod = fsiKernelMethod,
                sumDensityNodeLists = sumDensityNodeLists,
                linearCorrectGradients = correctVelocityGradient,
                compatibleEnergyEvolution = compatibleEnergy,
                evolveTotalEnergy = evolveTotalEnergy,
                HUpdate = HUpdate,
                ASPH = asph,
                epsTensile = epsilonTensile,
                nTensile = nTensile)

elif gsph:
    limiter = VanLeerLimiter()
    waveSpeed = DavisWaveSpeed()
    solver = HLLC(limiter,waveSpeed,gsphLinearCorrect,gsphGradientMethod)
    hydro = GSPH(dataBase = db,
                riemannSolver = solver,
                W = WT,
                cfl=cfl,
                specificThermalEnergyDiffusionCoefficient = gsphEpsDiffuseCoeff,
                compatibleEnergyEvolution = compatibleEnergy,
                correctVelocityGradient= correctVelocityGradient,
                evolveTotalEnergy = evolveTotalEnergy,
                HUpdate = HUpdate,
                XSPH = XSPH,
                ASPH = asph,
                epsTensile = epsilonTensile,
                nTensile = nTensile)         

else:
    hydro = SPH(dataBase = db,
                W = WT,
                cfl = cfl,
                compatibleEnergyEvolution = compatibleEnergy,
                gradhCorrection = gradhCorrection,
                correctVelocityGradient = correctVelocityGradient,
                XSPH = XSPH,
                densityUpdate = densityUpdate,
                HUpdate = HUpdate,
                epsTensile = epsilonTensile,
                nTensile = nTensile)
output("hydro")
output("hydro.cfl")
output("hydro.compatibleEnergyEvolution")
output("hydro.densityUpdate")
output("hydro.HEvolution")

packages = [hydro]

#-------------------------------------------------------------------------------
# Tweak the artificial viscosity.
#-------------------------------------------------------------------------------
if not gsph:
    q = hydro.Q
    if not Cl is None:
        q.Cl = Cl
    if not Cq is None:
        q.Cq = Cq
    if not Qlimiter is None:
        q.limiter = Qlimiter
    if not epsilon2 is None:
        q.epsilon2 = epsilon2
    output("q")
    output("q.Cl")
    output("q.Cq")
    output("q.limiter")
    output("q.epsilon2")
    output("q.linearInExpansion")
    output("q.quadraticInExpansion")

    #-------------------------------------------------------------------------------
    # Construct the MMRV physics object.
    #-------------------------------------------------------------------------------
    if boolReduceViscosity:
        evolveReducingViscosityMultiplier = MorrisMonaghanReducingViscosity(nh,aMin,aMax)
        packages.append(evolveReducingViscosityMultiplier)
    elif boolCullenViscosity:
        evolveCullenViscosityMultiplier = CullenDehnenViscosity(WT,alphMax,alphMin,betaC,betaD,betaE,fKern,boolHopkinsCorrection)
        packages.append(evolveCullenViscosityMultiplier)

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
xPlane0 = Plane(Vector(x0, y0), Vector( 1.0,  0.0))
xPlane1 = Plane(Vector(x3, y0), Vector(-1.0,  0.0))
yPlane0 = Plane(Vector(x0, y0), Vector( 0.0,  1.0))
yPlane1 = Plane(Vector(x0, y3), Vector( 0.0, -1.0))

xbc = PeriodicBoundary(xPlane0, xPlane1)
ybc = PeriodicBoundary(yPlane0, yPlane1)

bcSet = [xbc, ybc]

for p in packages:
    for bc in bcSet:
        p.appendBoundary(bc)

#-------------------------------------------------------------------------------
# Construct a time integrator, and add the physics packages.
#-------------------------------------------------------------------------------
integrator = IntegratorConstructor(db)
for p in packages:
    integrator.appendPhysicsPackage(p)
integrator.cullGhostNodes = False
integrator.lastDt = dt
integrator.dtMin = dtMin
integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
integrator.domainDecompositionIndependent = domainIndependent
integrator.verbose = dtverbose
integrator.rigorousBoundaries = rigorousBoundaries
output("integrator")
output("integrator.havePhysicsPackage(hydro)")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.dtGrowth")
output("integrator.domainDecompositionIndependent")
output("integrator.rigorousBoundaries")
output("integrator.verbose")

#-------------------------------------------------------------------------------
# Make the problem controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            initializeDerivatives = True,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle,
                            redistributeStep = redistributeStep,
                            vizBaseName = vizBaseName,
                            vizDir = vizDir,
                            vizStep = vizCycle,
                            vizTime = vizTime,
                            SPH = (not ASPH))
output("control")

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
if not steps is None:
    control.step(steps)

else:
    control.advance(goalTime, maxSteps)
    control.updateViz(control.totalSteps, integrator.currentTime, 0.0)
    control.dropRestartFile()

if serialDump:
  procs = mpi.procs
  rank = mpi.rank
  serialData = []
  i,j = 0,0
  for i in range(procs):
    for nodeL in nodeSet:
      if rank == i:
        for j in range(nodeL.numInternalNodes):
          serialData.append([nodeL.positions()[j],3.0/(nodeL.Hfield()[j].Trace()),nodeL.mass()[j],nodeL.massDensity()[j],nodeL.specificThermalEnergy()[j]])
  serialData = mpi.reduce(serialData,mpi.SUM)
  if rank == 0:
    f = open(dataDir + "/serialDump.ascii",'w')
    for i in range(len(serialData)):
      f.write("{0} {1} {2} {3} {4} {5} {6} {7}\n".format(i,serialData[i][0][0],serialData[i][0][1],0.0,serialData[i][1],serialData[i][2],serialData[i][3],serialData[i][4]))
    f.close()

