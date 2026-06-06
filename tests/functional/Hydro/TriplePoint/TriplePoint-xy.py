#-------------------------------------------------------------------------------
# The Infamous triple point shock test case.
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

title("2-D integrated hydro test -- triple point problem")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(
    # Left state.
    rho1 = 1.0,
    P1 = 1.0,
    gamma1 = 1.5,

    # Top state
    rho2 = 0.125,
    P2 = 0.1,
    gamma2 = 1.5,

    # Bottom state
    rho3 = 1.0,
    P3 = 0.1,
    gamma3 = 1.4,

    # Geometry 
    x0 = 0.0,
    x1 = 1.0,
    x2 = 7.0,
    y0 = 0.0,
    y1 = 1.5,
    y2 = 3.0,

    # Resolution and node seeding.
    nx1 = 20,
    ny1 = 60,

    nx2 = 40,
    ny2 = 10,

    nx3 = 120,
    ny3 = 30,

    # Optionally set the initial density with the sum definition, and recompute
    # the energy to re-establish the the initial pressures.
    relaxInitialDensity = False,

    nPerh = 1.51,
    KernelConstructor = NBSplineKernel,
    order = 5,

    svph = False,
    crksph = False,
    psph = False,
    asph = False,        # This just chooses the H algorithm -- you can use this with CRKSPH for instance.
    solid = False,       # If true, use the fluid limit of the solid hydro option
    filter = 0.0,        # For CRKSPH
    Qconstructor = None,
    boolReduceViscosity = False,
    nhQ = 5.0,
    nhL = 10.0,
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
    linearConsistent = False,
    fcentroidal = 0.0,
    fcellPressure = 0.0,
    Cl = None,
    Cq = None,
    linearInExpansion = None,
    quadraticInExpansion = None,
    Qlimiter = None,
    balsaraCorrection = None,
    epsilon2 = None,
    hmin = 1e-5,
    hmax = 0.5,
    hminratio = 0.1,
    cfl = 0.5,
    XSPH = False,
    epsilonTensile = 0.0,
    nTensile = 8,

    IntegratorConstructor = CheapSynchronousRK2Integrator,
    goalTime = 7.0,
    steps = None,
    vizCycle = 50,
    vizTime = 0.1,
    dt = 0.0001,
    dtMin = 1.0e-5, 
    dtMax = 0.1,
    dtGrowth = 2.0,
    maxSteps = None,
    statsStep = 10,
    HUpdate = IdealH,
    correctionOrder = LinearOrder,
    QcorrectionOrder = LinearOrder,
    volumeType = RKSumVolume,
    domainIndependent = False,
    rigorousBoundaries = False,
    dtverbose = False,

    densityUpdate = RigorousSumDensity, # VolumeScaledDensity,
    compatibleEnergy = True,
    gradhCorrection = True,
    correctVelocityGradient = True,
    HopkinsConductivity = False,     # For PSPH
    evolveTotalEnergy = False,       # Only for SPH variants -- evolve total rather than specific energy

    useVoronoiOutput = False,
    clearDirectories = False,
    restoreCycle = None,
    restartStep = 200,
    dataDir = "dumps-triplepoint-xy",
    serialDump = False, #whether to dump a serial ascii file at the end for viz
    )

assert not(boolReduceViscosity and boolCullenViscosity)

# Decide on our hydro algorithm.
if svph:
    hydroname = "SVPH"
elif crksph:
    hydroname = "CRKSPH"
elif psph:
    hydroname = "PSPH"
else:
    hydroname = "SPH"
if asph:
    hydroname = "A" + hydroname
if solid:
    hydroname = "Solid" + hydroname

# Build our directory paths.
baseDir = os.path.join(dataDir,
                       hydroname,
                       "densityUpdate=%s" % densityUpdate,
                       "linearConsistent=%s" % linearConsistent,
                       "XSPH=%s" % XSPH,
                       "Cullen=%s" % boolCullenViscosity,
                       "nPerh=%3.1f" % nPerh,
                       "fcentroidal=%1.3f" % fcentroidal,
                       "fcellPressure=%1.3f" % fcellPressure,
                       "filter=%f" % filter,
                       "%ix%i" % (nx1 + nx2, ny1 + ny2))
if Cl and Cq:
    baseDir = os.path.join(baseDir,
                           "Cl=%3.1f_Cq=%3.1f" % (Cl,Cq))
if crksph:
    baseDir = os.path.join(baseDir, 
                           "correctionOrder=%s" % correctionOrder,
                           "QcorrectionOrder=%s" % QcorrectionOrder)
restartDir = os.path.join(baseDir, "restarts")
restartBaseName = os.path.join(restartDir, "triplepoint-xy-%ix%i" % (nx1 + nx2, ny1 + ny2))

vizDir = os.path.join(baseDir, "visit")
if vizTime is None and vizCycle is None:
    vizBaseName = None
else:
    vizBaseName = "triplepoint-xy-%ix%i" % (nx1 + nx2, ny1 + ny2)

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
# If we're restarting, find the set of most recent restart files.
#-------------------------------------------------------------------------------
if restoreCycle is None:
    restoreCycle = findLastRestart(restartBaseName)

#-------------------------------------------------------------------------------
# Material properties.
#-------------------------------------------------------------------------------
units = CGS()
mu = 1.0
eos1 = GammaLawGas(gamma1, mu, minimumPressure = 0.0, constants = units)
eos2 = GammaLawGas(gamma2, mu, minimumPressure = 0.0, constants = units)
eos3 = GammaLawGas(gamma3, mu, minimumPressure = 0.0, constants = units)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
if KernelConstructor==NBSplineKernel:
  WT = TableKernel(NBSplineKernel(order), 1000)
else:
  WT = TableKernel(KernelConstructor(), 1000)
output("WT")
kernelExtent = WT.kernelExtent

#-------------------------------------------------------------------------------
# Make the NodeLists.
#-------------------------------------------------------------------------------
if solid:
    NodeListConstructor = makeSolidNodeList
else:
    NodeListConstructor = makeFluidNodeList

leftNodes = NodeListConstructor("Left", eos1,
                                hmin = hmin,
                                hmax = hmax,
                                hminratio = hminratio,
                                nPerh = nPerh,
                                kernelExtent = kernelExtent)
topNodes = NodeListConstructor("Top", eos2,
                               hmin = hmin,
                               hmax = hmax,
                               hminratio = hminratio,
                               nPerh = nPerh,
                               kernelExtent = kernelExtent)
bottomNodes = NodeListConstructor("Bottom", eos3,
                                  hmin = hmin,
                                  hmax = hmax,
                                  hminratio = hminratio,
                                  nPerh = nPerh,
                                  kernelExtent = kernelExtent)
nodeSet = (leftNodes, topNodes, bottomNodes)
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
if restoreCycle is None:
    generatorLeft = GenerateNodeDistribution2d(nx1, ny1, rho1,
                                               distributionType = "lattice",
                                               xmin = (x0, y0),
                                               xmax = (x1, y2),
                                               nNodePerh = nPerh,
                                               SPH = SPH)
    generatorTop = GenerateNodeDistribution2d(nx2, ny2, rho2,
                                              distributionType = "lattice",
                                              xmin = (x1, y1),
                                              xmax = (x2, y2),
                                              nNodePerh = nPerh,
                                              SPH = SPH)
    generatorBottom = GenerateNodeDistribution2d(nx3, ny3, rho3,
                                                 distributionType = "lattice",
                                                 xmin = (x1, y0),
                                                 xmax = (x2, y1),
                                                 nNodePerh = nPerh,
                                                 SPH = SPH)

    if mpi.procs > 1:
        from PeanoHilbertDistributeNodes import distributeNodes2d
    else:
        from DistributeNodes import distributeNodes2d

    distributeNodes2d((leftNodes, generatorLeft),
                      (topNodes,generatorTop),
                      (bottomNodes, generatorBottom))
    for nodes in nodeSet:
        print(nodes.name, ":")
        output("    mpi.reduce(nodes.numInternalNodes, mpi.MIN)")
        output("    mpi.reduce(nodes.numInternalNodes, mpi.MAX)")
        output("    mpi.reduce(nodes.numInternalNodes, mpi.SUM)")
    del nodes

    # Set node specific thermal energies
    for (nodes, gamma, rho, P) in ((leftNodes, gamma1, rho1, P1),
                                   (topNodes, gamma2, rho2, P2),
                                   (bottomNodes, gamma3, rho3, P3)):
        eps0 = P/((gamma - 1.0)*rho)
        nodes.specificThermalEnergy(ScalarField("tmp", nodes, eps0))
    del nodes

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
                 xmin = Vector(x0 - (x2 - x0), y0 - (y2 - y0)),
                 xmax = Vector(x2 + (x2 - x0), y2 + (y2 - y0)),
                 ASPH = asph)
elif crksph:
    hydro = CRKSPH(dataBase = db,
                   W = WT,
                   filter = filter,
                   cfl = cfl,
                   compatibleEnergyEvolution = compatibleEnergy,
                   XSPH = XSPH,
                   correctionOrder = correctionOrder,
                   volumeType = volumeType,
                   densityUpdate = densityUpdate,
                   HUpdate = HUpdate,
                   ASPH = asph)
elif psph:
    hydro = PSPH(dataBase = db,
                 W = WT,
                 filter = filter,
                 cfl = cfl,
                 compatibleEnergyEvolution = compatibleEnergy,
                 evolveTotalEnergy = evolveTotalEnergy,
                 HopkinsConductivity = HopkinsConductivity,
                 correctVelocityGradient = correctVelocityGradient,
                 densityUpdate = densityUpdate,
                 HUpdate = HUpdate,
                 XSPH = XSPH,
                 ASPH = asph)
else:
    hydro = SPH(dataBase = db,
                W = WT,
                filter = filter,
                cfl = cfl,
                compatibleEnergyEvolution = compatibleEnergy,
                evolveTotalEnergy = evolveTotalEnergy,
                gradhCorrection = gradhCorrection,
                correctVelocityGradient = correctVelocityGradient,
                densityUpdate = densityUpdate,
                HUpdate = HUpdate,
                XSPH = XSPH,
                epsTensile = epsilonTensile,
                nTensile = nTensile,
                ASPH = asph)
output("hydro")
output("hydro.kernel")
output("hydro.PiKernel")
output("hydro.cfl")
output("hydro.compatibleEnergyEvolution")
output("hydro.densityUpdate")
output("hydro.HEvolution")

packages = [hydro]

#-------------------------------------------------------------------------------
# Set the artificial viscosity parameters.
#-------------------------------------------------------------------------------
q = hydro.Q
if not Cl is None:
    q.Cl = Cl
if not Cq is None:
    q.Cq = Cq
if not epsilon2 is None:
    q.epsilon2 = epsilon2
if not Qlimiter is None:
    q.limiter = Qlimiter
if not balsaraCorrection is None:
    q.balsaraShearCorrection = balsaraCorrection
if not quadraticInExpansion is None:
    q.quadraticInExpansion = quadraticInExpansion
output("q")
output("q.Cl")
output("q.Cq")
output("q.epsilon2")
output("q.limiter")
output("q.balsaraShearCorrection")
try:
    output("q.linearInExpansion")
    output("q.quadraticInExpansion")
except:
    pass

#-------------------------------------------------------------------------------
# Construct the MMRV physics object.
#-------------------------------------------------------------------------------
if boolReduceViscosity:
    evolveReducingViscosityMultiplier = MorrisMonaghanReducingViscosity(nhQ,nhL,aMin,aMax)
    packages.append(evolveReducingViscosityMultiplier)
elif boolCullenViscosity:
    evolveCullenViscosityMultiplier = CullenDehnenViscosity(WT,alphMax,alphMin,betaC,betaD,betaE,fKern,boolHopkinsCorrection)
    packages.append(evolveCullenViscosityMultiplier)

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
xPlane0 = Plane(Vector(x0, y0), Vector( 1.0,  0.0))
xPlane1 = Plane(Vector(x2, y0), Vector(-1.0,  0.0))
yPlane0 = Plane(Vector(x0, y0), Vector( 0.0,  1.0))
yPlane1 = Plane(Vector(x0, y2), Vector( 0.0, -1.0))

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
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle,
                            vizBaseName = vizBaseName,
                            vizDir = vizDir,
                            vizStep = vizCycle,
                            vizTime = vizTime,
                            skipInitialPeriodicWork = svph,
                            SPH = not asph)
output("control")

#-------------------------------------------------------------------------------
# Optionally compute the sum density and recompute the initial conditions.
#-------------------------------------------------------------------------------
if relaxInitialDensity and control.totalSteps == 0:
    state = State(db, integrator.physicsPackages())
    derivs = StateDerivatives(db, integrator.physicsPackages())
    integrator.preStepInitialize(state, derivs)
    integrator.initializeDerivatives(0.0, 0.0, state, derivs)
    cm = db.connectivityMap()
    mass = state.scalarFields(HydroFieldNames.mass)
    position = state.vectorFields(HydroFieldNames.position)
    H = state.symTensorFields(HydroFieldNames.H)
    rho = state.scalarFields(HydroFieldNames.massDensity)
    if CRKSPH:
        vol = state.scalarFields(HydroFieldNames.volume)
        A = state.scalarFields(HydroFieldNames.A_CRKSPH)
        B = state.vectorFields(HydroFieldNames.B_CRKSPH)
        C = state.tensorFields(HydroFieldNames.C_CRKSPH)
        computeCRKSPHSumMassDensity(cm, WT, position, mass, vol, H, A, B, C, hydro.correctionOrder, rho)
    else:
        computeSPHSumMassDensity(cm, WT, True, position, mass, H, rho)
    for (nodes, P0) in ((leftNodes, P1),
                        (topNodes, P2),
                        (bottomNodes, P3)):
        rhof = nodes.massDensity()
        epsf = nodes.specificThermalEnergy()
        for i in range(nodes.numInternalNodes):
            epsf[i] = P0/((gamma1 - 1.0)*rhof[i])

    # Force another visit dumps so we can see what changed.
    control.dropViz(0, 0.0, 0.0)

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
if not steps is None:
    control.step(steps)

else:
    control.advance(goalTime, maxSteps)
    control.updateViz(control.totalSteps, integrator.currentTime, 0.0)
    control.dropRestartFile()

Eerror = (control.conserve.EHistory[-1] - control.conserve.EHistory[0])/max(1.0e-30, control.conserve.EHistory[0])
print("Total energy error: %g" % Eerror)

if serialDump:
  procs = mpi.procs
  rank = mpi.rank
  serialData = []
  i,j,k = 0,0,0
  for i in range(procs):
    if rank == i:
        k = 0
        for nodeL in nodeSet:
            for j in range(nodeL.numInternalNodes):
                serialData.append([nodeL.positions()[j],3.0/(nodeL.Hfield()[j].Trace()),nodeL.mass()[j],nodeL.massDensity()[j],nodeL.specificThermalEnergy()[j],k])
            k = k + 1
  serialData = mpi.reduce(serialData,mpi.SUM)
  if rank == 0:
    f = open(baseDir + "/serialDump.ascii",'w')
    for i in range(len(serialData)):
      f.write("{0} {1} {2} {3} {4} {5} {6} {7}\n".format(i,serialData[i][0][0],serialData[i][0][1],0.0,serialData[i][1],serialData[i][2],serialData[i][5],serialData[i][4]))
    f.close()
