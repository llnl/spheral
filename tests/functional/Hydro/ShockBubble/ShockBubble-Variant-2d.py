#-------------------------------------------------------------------------------
# Spheral++ script to generate the shock-bubble mesh.
#-------------------------------------------------------------------------------
import shutil
import mpi
import numpy
from math import *
from Spheral2d import *
from SpheralTestUtilities import *
from findLastRestart import *
from GenerateNodeDistribution2d import *
from CompositeNodeDistribution import *
from VoronoiDistributeNodes import distributeNodes2d as distributeNodes

title("2-D hydro test -- cylindrical He-Air bubble shock problem")

# Conversions.
m2cm = 1.0e2
cubm2cubc = m2cm**3
kg2g = 1.0e3
rhoconv = kg2g/cubm2cubc
Econv = kg2g/(m2cm)
Fconv = kg2g*m2cm
Pconv = Fconv/(m2cm**2)
vconv = m2cm

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(air2He1 = 2.0,            # Ratio of zone lengths in Air/He
            air2He2 = 4.0,
            nrHe = 25,

            # material properties
            airDensity   = 1.0 * rhoconv,
            airPressure  = 1.0 * Pconv,
            airHighPressure = 10.0 * Pconv,
            airMu = 28.96,
            airGamma     = 1.4,

            HeDensity   = 0.182 * rhoconv,
            HePressure  = 1.0 * Pconv,
            HeMu = 5.269,
            HeGamma     = 1.648,
            HeSeed = "constantDTheta",

            # Hydro parameters.
            nPerh = 1.51,
            SVPH = False,
            CRKSPH = False,
            ASPH = False,
            SPH = True,   # This just chooses the H algorithm -- you can use this with CRKSPH for instance.
            filter = 0.0,   # CRKSPH filtering
            Qconstructor = MonaghanGingoldViscosity,
            #Qconstructor = TensorMonaghanGingoldViscosity,
            Cl = 1.0, 
            Cq = 0.75,
            Qlimiter = False,
            balsaraCorrection = False,
            epsilon2 = 1e-2,
            hmin = 1.0e-10,
            hmax = 10.0,
            hminratio = 0.1,
            cfl = 0.5,
            XSPH = True,
            epsilonTensile = 0.0,
            nTensile = 8,
            hourglass = None,
            hourglassOrder = 0,
            hourglassLimiter = 0,
            hourglassFraction = 0.5,
            densityUpdate = RigorousSumDensity,
            compatibleEnergy = True,
            gradhCorrection = False,

            # Timestepping/advancement
            IntegratorConstructor = CheapSynchronousRK2Integrator,
            goalTime = 0.003,
            steps = None,
            vizCycle = 20,
            vizTime = 1e-4,
            dt = 1.0e-5,
            dtMin = 1.0e-8, 
            dtMax = 1.0e-2,
            dtGrowth = 2.0,
            maxSteps = None,
            statsStep = 10,
            smoothIters = 0,
            HUpdate = IdealH,
            domainIndependent = False,
            rigorousBoundaries = False,
            dtverbose = False,

            useVoronoiOutput = False,
            clearDirectories = False,
            restoreCycle = None,
            restartStep = 50,
            checkRestart = False,
            dataDir = "dumps-bubbleShock-variant-2d",
            vizName = "ShockBubble-variant-2d",
            outputFile = None,
            )

airEnergy = airPressure/((airGamma - 1.0)*airDensity)
airHighEnergy = airHighPressure/((airGamma - 1.0)*airDensity)
HeEnergy = HePressure/((HeGamma - 1.0)*HeDensity)

# Overall box dimensions.
x0, x1, x2 = 0.0, 1.0, 1.1
y0, y1 = 0.0, 0.4

# He bubble
x0He, y0He, rHe = 0.8, 0.2, 0.1
drHe = rHe/nrHe

# Air region 1.
x0Air1, x1Air1 = x0, x1
nxAir1 = int((x1Air1 - x0Air1)/(air2He1*drHe) + 0.5)
dxAir1 = (x1Air1 - x0Air1)/nxAir1
nyAir1 = int((y1 - y0)/dxAir1 + 0.5)
dyAir1 = (y1 - y0)/nyAir1

# Air region 2.
x0Air2, x1Air2 = x1, x2
nxAir2 = int((x1Air2 - x0Air2)/(air2He2*drHe) + 0.5)
dxAir2 = (x1Air2 - x0Air2)/nxAir2
nyAir2 = int((y1 - y0)/dxAir2 + 0.5)
dyAir2 = (y1 - y0)/nyAir2

# We need to create an extra layer of nodes around the bubble to give us the
# the expected smooth surface.
nintHe = 1
r1He = rHe + nintHe*drHe

# Decide on our hydro algorithm.
if SVPH:
    if ASPH:
        HydroConstructor = ASVPHFacetedHydro
    else:
        HydroConstructor = SVPHFacetedHydro
elif CRKSPH:
    if ASPH:
        HydroConstructor = ACRKSPHHydro
    else:
        HydroConstructor = CRKSPHHydro
    Qconstructor = LimitedMonaghanGingoldViscosity
else:
    if ASPH:
        HydroConstructor = ASPHHydro
    else:
        HydroConstructor = SPHHydro

dataDir = os.path.join(dataDir,
                       "nrHe=%i" % nrHe,
                       str(HydroConstructor).split("'")[1].split(".")[-1],
                       "densityUpdate=%s" % densityUpdate,
                       "XSPH=%s" % XSPH,
                       "%s-Cl=%g-Cq=%g" % (str(Qconstructor).split("'")[1].split(".")[-1], Cl, Cq),
                       "nPerh=%g" % nPerh)
restartDir = os.path.join(dataDir, "restarts")
vizDir = os.path.join(dataDir, "visit")
restartBaseName = os.path.join(restartDir, "BubbleShock-cylindrical-2d-%i" % (nrHe))

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
import os, sys
if mpi.rank == 0:
    for DIR in (restartDir, vizDir):
        if clearDirectories and os.path.exists(DIR):
            shutil.rmtree(DIR)
        if not os.path.exists(DIR):
            os.makedirs(DIR)
mpi.barrier()

#-------------------------------------------------------------------------------
# If we're restarting, find the set of most recent restart files.
#-------------------------------------------------------------------------------
if restoreCycle is None:
    restoreCycle = findLastRestart(restartBaseName)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
WT = TableKernel(BSplineKernel(), 100000)
output("WT")

#===============================================================================
# Material properties
#===============================================================================
airEOS = GammaLawGasCGS(airGamma, airMu)
HeEOS = GammaLawGasCGS(HeGamma, HeMu)

#===============================================================================
# Create the NodeLists.
#===============================================================================
nodesAir1 = makeFluidNodeList("air", airEOS,
                              hmin = hmin,
                              hmax = hmax,
                              hminratio = hminratio,
                              nPerh = nPerh)
nodesAir2 = makeFluidNodeList("high pressure air", airEOS,
                              hmin = hmin,
                              hmax = hmax,
                              hminratio = hminratio,
                              nPerh = nPerh)
nodesHe = makeFluidNodeList("bubble", HeEOS,
                            hmin = hmin,
                            hmax = hmax,
                            hminratio = hminratio,
                            nPerh = nPerh)

nodeSet = (nodesAir1, nodesAir2, nodesHe)

#===============================================================================
# Create the generators.
#===============================================================================
if restoreCycle is None:
    if HeSeed == "lattice":
        generatorHe = GenerateNodeDistribution2d(2*nrHe, 2*nrHe, HeDensity,
                                                 distributionType = HeSeed,
                                                 xmin = (x0He - rHe, 0.5*y1 - rHe),
                                                 xmax = (x0He + rHe, 0.5*y1 + rHe),
                                                 rreject = rHe,
                                                 originreject = (x0He, y0He),
                                                 theta = 2*pi,
                                                 nNodePerh = nPerh,
                                                 SPH = (HydroConstructor == SPHHydro))
        generatorAir = GenerateNodeDistribution2d(nxAir1, nyAir1, airDensity,
                                                  distributionType = "lattice",
                                                  xmin = (x0Air1, y0),
                                                  xmax = (x1Air1, y1),
                                                  rreject = rHe,
                                                  originreject = (x0He, y0He),
                                                  reversereject = True,
                                                  nNodePerh = nPerh,
                                                  SPH = (HydroConstructor == SPHHydro))
    else:
        generatorHe = GenerateNodeDistribution2d(nrHe, nrHe, HeDensity,
                                                 distributionType = HeSeed,
                                                 rmin = 0.0,
                                                 rmax = rHe,
                                                 theta = 2*pi,
                                                 offset = (x0He, y0He),
                                                 nNodePerh = nPerh,
                                                 SPH = (HydroConstructor == SPHHydro))
        generatorAir1 = GenerateNodeDistribution2d(nxAir1, nyAir1, airDensity,
                                                   distributionType = "lattice",
                                                   xmin = (x0Air1, y0),
                                                   xmax = (x1Air1, y1),
                                                   rreject = r1He,
                                                   originreject = (x0He, y0He),
                                                   reversereject = True,
                                                   nNodePerh = nPerh,
                                                   SPH = (HydroConstructor == ASPHHydro))
        # This is a litle hacky -- in order to get a smooth interface on the bubble we
        # want exactly the same number of nodes in a ring for the inner boundary of the
        # air as are in the outer ring of nodes in the He.
        radiiHe = [(Vector2d(xi, yi) - Vector2d(x0He, y0He)).magnitude() for xi, yi in zip(generatorHe.x, generatorHe.y)]
        nThetaAir1 = mpi.allreduce(len([x for x in radiiHe if x > rHe - drHe]), mpi.SUM)
        generatorAir1_interface = GenerateNodeDistribution2d(nintHe, nThetaAir1, airDensity,
                                                             distributionType = "constantDTheta",
                                                             rmin = rHe,
                                                             rmax = r1He,
                                                             theta = 2*pi,
                                                             offset = (x0He, y0He),
                                                             nNodePerh = nPerh,
                                                             SPH = (HydroConstructor == SPHHydro))
        generatorAir = CompositeNodeDistribution(generatorAir1, generatorAir1_interface)
    generatorAir2 = GenerateNodeDistribution2d(nxAir2, nyAir2, airDensity,
                                               distributionType = "lattice",
                                               xmin = (x0Air2, y0),
                                               xmax = (x1Air2, y1),
                                               nNodePerh = nPerh,
                                               SPH = (HydroConstructor == SPHHydro))

    distributeNodes((nodesAir1, generatorAir),
                    (nodesHe, generatorHe),
                    (nodesAir2, generatorAir2))
    for nodes in nodeSet:
        print("Num internal nodes for ", nodes.name, " : ", mpi.allreduce(nodes.numInternalNodes, mpi.SUM))

    # Set initial conditions.
    nodesAir1.specificThermalEnergy(ScalarField("eps", nodesAir1, airEnergy))
    nodesAir2.specificThermalEnergy(ScalarField("eps", nodesAir2, airHighEnergy))
    nodesHe.specificThermalEnergy(ScalarField("eps", nodesHe, HeEnergy))

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase()
output("db")
for nodes in nodeSet:
    output("db.appendNodeList(nodes)")
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Construct the artificial viscosity.
#-------------------------------------------------------------------------------
q = Qconstructor(Cl, Cq)
q.epsilon2 = epsilon2
q.limiter = Qlimiter
q.balsaraShearCorrection = balsaraCorrection
output("q")
output("q.Cl")
output("q.Cq")
output("q.epsilon2")
output("q.limiter")
output("q.balsaraShearCorrection")

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
if SVPH:
    hydro = HydroConstructor(W = WT,
                             Q = q,
                             cfl = cfl,
                             compatibleEnergyEvolution = compatibleEnergy,
                             densityUpdate = densityUpdate,
                             XSVPH = XSPH,
                             linearConsistent = linearConsistent,
                             generateVoid = False,
                             HUpdate = HUpdate,
                             fcentroidal = fcentroidal,
                             fcellPressure = fcellPressure,
                             xmin = Vector(-2.0, -2.0),
                             xmax = Vector(3.0, 3.0))
                             # xmin = Vector(x0 - 0.5*(x2 - x0), y0 - 0.5*(y2 - y0)),
                             # xmax = Vector(x2 + 0.5*(x2 - x0), y2 + 0.5*(y2 - y0)))
elif CRKSPH:
    hydro = HydroConstructor(W = WT,
                             Q = q,
                             filter = filter,
                             cfl = cfl,
                             compatibleEnergyEvolution = compatibleEnergy,
                             XSPH = XSPH,
                             densityUpdate = densityUpdate,
                             HUpdate = HUpdate)
else:
    hydro = HydroConstructor(W = WT,
                             Q = q,
                             cfl = cfl,
                             compatibleEnergyEvolution = compatibleEnergy,
                             gradhCorrection = gradhCorrection,
                             XSPH = XSPH,
                             densityUpdate = densityUpdate,
                             HUpdate = HUpdate,
                             epsTensile = epsilonTensile,
                             nTensile = nTensile)

output("hydro")
output("hydro.kernel()")
output("hydro.cfl")
output("hydro.compatibleEnergyEvolution")
output("hydro.densityUpdate")
output("hydro.HEvolution")

packages = [hydro]

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
xPlane0 = Plane(Vector(x0, y0), Vector( 1.0,  0.0))
xPlane1 = Plane(Vector(x2, y0), Vector(-1.0,  0.0))
yPlane0 = Plane(Vector(x0, y0), Vector( 0.0,  1.0))
yPlane1 = Plane(Vector(x0, y1), Vector( 0.0, -1.0))
xbc0 = ReflectingBoundary(xPlane0)
xbc1 = ReflectingBoundary(xPlane1)
ybc0 = ReflectingBoundary(yPlane0)
ybc1 = ReflectingBoundary(yPlane1)
for p in packages:
    for bc in (xbc0, xbc1, ybc0, ybc1):
        p.appendBoundary(bc)

#-------------------------------------------------------------------------------
# Construct a time integrator, and add the physics packages.
#-------------------------------------------------------------------------------
integrator = IntegratorConstructor(db)
for p in packages:
    integrator.appendPhysicsPackage(p)
integrator.lastDt = dt
integrator.dtMin = dtMin
integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
integrator.domainDecompositionIndependent = domainIndependent
integrator.verbose = dtverbose
integrator.rigorousBoundaries = rigorousBoundaries
output("integrator")
output("integrator.havePhysicsPackage(hydro)")
if hourglass:
    output("integrator.havePhysicsPackage(hg)")
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
if useVoronoiOutput:
    import SpheralVoronoiSiloDump
    vizMethod = SpheralVoronoiSiloDump.dumpPhysicsState
else:
    import SpheralPointmeshSiloDump
    vizMethod = SpheralPointmeshSiloDump.dumpPhysicsState
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle,
                            vizMethod = vizMethod,
                            vizBaseName = vizName,
                            vizDir = vizDir,
                            vizStep = vizCycle,
                            vizTime = vizTime)
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

