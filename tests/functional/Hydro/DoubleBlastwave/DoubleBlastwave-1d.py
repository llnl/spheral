import os, sys
import shutil
from Spheral1d import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *

title("Woodward-Colella double blastwave test.")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(

    # Resolution
    nx = 400,
    nPerh = 1.25,
    hmin = 1e-10,
    hmax = 1.0,
    hsmooth = 0.0,    # Optionally smooth the initial discontinuities
            
    # Hydro algorithm.
    SVPH = False,
    CRKSPH = False,
    PSPH = False,
    cfl = 0.5,
    XSPH = False,
    epsilonTensile = 0.0,
    nTensile = 8,
    filter = 0.00,
    HUpdate = IdealH,
    volumeType = RKSumVolume,
    densityUpdate = RigorousSumDensity,
    correctionOrder = LinearOrder,
    compatibleEnergy = True,
    gradhCorrection = True,
    correctVelocityGradient = True,
    linearConsistent = False,
    KernelConstructor = BSplineKernel,
    order = 5,
    HopkinsConductivity = False,     # For PSPH
    evolveTotalEnergy = False,       # Only for SPH variants -- evolve total rather than specific energy

    # Q
    Qconstructor = MonaghanGingoldViscosity,
    Cl = 1.0,
    Cq = 1.0,
    linearInExpansion = False,
    Qlimiter = False,
    epsilon2 = 1e-4,
    etaCritFrac = 1.0,
    etaFoldFrac = 0.2,
  
    boolReduceViscosity = False,
    nh = 5.0,
    aMin = 0.1,
    aMax = 2.0,
    Qhmult = 1.0,
    boolCullenViscosity = False,
    alphMax = 2.0,
    alphMin = 0.02,
    betaC = 0.7,
    betaD = 0.05,
    betaE = 1.0,
    fKern = 1.0/3.0,
    boolHopkinsCorrection = True,

    # Integrator.
    IntegratorConstructor = CheapSynchronousRK2Integrator,
    dtverbose = False,
    rigorousBoundaries = False,
    dt = 1e-6,
    dtMin = 1.0e-7,
    dtMax = 0.1,
    dtGrowth = 2.0,

    # Times.
    goalTimes = (0.01, 0.016, 0.026, 0.028, 0.030, 0.032, 0.034, 0.038), 
    steps = None,
    maxSteps = None,
    statsStep = 10,
    numHIterationsBetweenCycles = 0,

    # Output
    clearDirectories = False,
    restoreCycle = None,
    restartStep = 200,
    dataDirBase = "dumps-DoubleBlastwave-1d",
    restartBaseName = "DoubleBlastwave-restart",
    graphics = True,
)
assert not(boolReduceViscosity and boolCullenViscosity)

if SVPH:
    HydroConstructor = SVPHFacetedHydro
elif CRKSPH:
    HydroConstructor = CRKSPHHydro
    Qconstructor = LimitedMonaghanGingoldViscosity
elif PSPH:
    HydroConstructor = PSPHHydro
else:
    HydroConstructor = SPHHydro

dataDir = os.path.join(dataDirBase, 
                       HydroConstructor.__name__,
                       Qconstructor.__name__,
                       "Cl=%g_Cq=%g" % (Cl, Cq),
                       "%i" % (nx),
                       "nPerh=%s" % nPerh)
restartDir = os.path.join(dataDir, "restarts")
restartBaseName = os.path.join(restartDir, "DoubleBlastwave-1d-%i")
plotDir = os.path.join(dataDir, "plots")

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
if mpi.rank == 0:
    if clearDirectories and os.path.exists(dataDir):
        shutil.rmtree(dataDir)
    if not os.path.exists(restartDir):
        os.makedirs(restartDir)
    if not os.path.exists(plotDir):
        os.makedirs(plotDir)
mpi.barrier()

#-------------------------------------------------------------------------------
# Material properties.
#-------------------------------------------------------------------------------
eos = GammaLawGasMKS(1.4, 2.0)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
if KernelConstructor==NBSplineKernel:
  WT = TableKernel(NBSplineKernel(order), 1000)
else:
  WT = TableKernel(KernelConstructor(), 1000)
kernelExtent = WT.kernelExtent
output("WT")

#-------------------------------------------------------------------------------
# Make the NodeLists.
#-------------------------------------------------------------------------------
nodes = makeFluidNodeList("nodes", eos, 
                           hmin = hmin,
                           hmax = hmax,
			   kernelExtent = kernelExtent,
                           nPerh = nPerh)
nodeSet = [nodes]

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
from DistributeNodes import distributeNodesInRange1d
rho0 = 1.0
distributeNodesInRange1d([(nodes, [(nx, rho0, (0.0, 1.0))])])
output("nodes.numNodes")

# Set the initial conditions.
eps1, eps2, eps3 = 1000.0/0.4, 0.01/0.4, 100.0/0.4
hfold = hsmooth*1.0/nx
def epsfunc(x):
    if x < 0.1 - hfold:
        return eps1
    elif x < 0.1 + hfold:
        f = 0.5*(sin(0.5*pi*(x - 0.1)/hfold) + 1.0)
        return (1.0 - f)*eps1 + f*eps2
    elif x < 0.9 - hfold:
        return eps2
    elif x < 0.9 + hfold:
        f = 0.5*(sin(0.5*pi*(x - 0.9)/hfold) + 1.0)
        return (1.0 - f)*eps2 + f*eps3
    else:
        return eps3
pos = nodes.positions()
eps = nodes.specificThermalEnergy()
for i in range(nodes.numInternalNodes):
    eps[i] = epsfunc(pos[i].x)

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
try:
    q = Qconstructor(Cl, Cq, linearInExpansion)
except:
    q = Qconstructor(Cl, Cq)
q.limiter = Qlimiter
q.epsilon2 = epsilon2
output("q")
output("q.Cl")
output("q.Cq")
output("q.limiter")
output("q.epsilon2")
try:
    output("q.linearInExpansion")
    output("q.quadraticInExpansion")
except:
    pass

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
if SVPH:
    hydro = HydroConstructor(W = WT, 
                             Q = q,
                             cfl = cfl,
                             compatibleEnergyEvolution = compatibleEnergy,
                             XSVPH = XSPH,
                             linearConsistent = linearConsistent,
                             generateVoid = False,
                             densityUpdate = densityUpdate,
                             HUpdate = HUpdate,
                             xmin = Vector(-100.0),
                             xmax = Vector( 100.0))
elif CRKSPH:
    hydro = HydroConstructor(W = WT, 
                             Q = q,
                             filter = filter,
                             cfl = cfl,
                             volumeType = volumeType,
                             compatibleEnergyEvolution = compatibleEnergy,
                             evolveTotalEnergy = evolveTotalEnergy,
                             XSPH = XSPH,
                             correctionOrder = correctionOrder,
                             densityUpdate = densityUpdate,
                             HUpdate = HUpdate)
    q.etaCritFrac = etaCritFrac
    q.etaFoldFrac = etaFoldFrac
    q.QcorrectionOrder = correctionOrder
elif PSPH:
    hydro = HydroConstructor(W = WT,
                             Q = q,
                             filter = filter,
                             cfl = cfl,
                             compatibleEnergyEvolution = compatibleEnergy,
                             evolveTotalEnergy = evolveTotalEnergy,
                             HopkinsConductivity = HopkinsConductivity,
                             densityUpdate = densityUpdate,
                             correctVelocityGradient = correctVelocityGradient,
                             HUpdate = HUpdate,
                             XSPH = XSPH)
else:
    hydro = HydroConstructor(W = WT,
                             Q = q,
                             cfl = cfl,
                             compatibleEnergyEvolution = compatibleEnergy,
                             evolveTotalEnergy = evolveTotalEnergy,
                             gradhCorrection = gradhCorrection,
                             densityUpdate = densityUpdate,
                             HUpdate = HUpdate,
                             XSPH = XSPH,
                             correctVelocityGradient = correctVelocityGradient,
                             epsTensile = epsilonTensile,
                             nTensile = nTensile)
output("hydro")

packages = [hydro]

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
xPlane0 = Plane(Vector(0.0), Vector( 1.0))
xPlane1 = Plane(Vector(1.0), Vector(-1.0))
xbc0 = ReflectingBoundary(xPlane0)
xbc1 = ReflectingBoundary(xPlane1)
for p in packages:
    p.appendBoundary(xbc0)
    p.appendBoundary(xbc1)

#-------------------------------------------------------------------------------
# Construct an integrator, and add the one physics package.
#-------------------------------------------------------------------------------
integrator = IntegratorConstructor(db)
for p in packages:
    integrator.appendPhysicsPackage(p)
integrator.lastDt = dt
integrator.dtGrowth = dtGrowth
if dtMin:
    integrator.dtMin = dtMin
if dtMax:
    integrator.dtMax = dtMax
integrator.rigorousBoundaries = rigorousBoundaries
integrator.verbose = dtverbose
output("integrator")
output("integrator.havePhysicsPackage(hydro)")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.rigorousBoundaries")

#-------------------------------------------------------------------------------
# Make the problem controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle,
                            numHIterationsBetweenCycles = numHIterationsBetweenCycles)
output("control")

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
if not steps is None:
    control.step(steps)

else:
    rhoPlots = []
    for goalTime in goalTimes:
        control.advance(goalTime, maxSteps)
        
        # Write the current state.
        xprof = mpi.reduce([x.x for x in nodes.positions().internalValues()], mpi.SUM)
        rhoprof = mpi.reduce(nodes.massDensity().internalValues(), mpi.SUM)
        P = ScalarField("pressure", nodes)
        nodes.pressure(P)
        Pprof = mpi.reduce(P.internalValues(), mpi.SUM)
        vprof = mpi.reduce([v.x for v in nodes.velocity().internalValues()], mpi.SUM)
        epsprof = mpi.reduce(nodes.specificThermalEnergy().internalValues(), mpi.SUM)
        hprof = mpi.reduce([1.0/H.xx for H in nodes.Hfield().internalValues()], mpi.SUM)
        if mpi.rank == 0:
            multiSort(xprof, rhoprof, Pprof, vprof, epsprof, hprof)
            outputFile = os.path.join(dataDir, "DoubleBlastwave_state_t=%4.3f.ascii" % control.time())
            f = open(outputFile, "w")
            f.write(("#  " + 6*"'%s' " + "\n") % ("x", "rho", "P", "v", "eps", "h"))
            for (xi, rhoi, Pi, vi, epsi, hi) in zip(xprof, rhoprof, Pprof, vprof, epsprof, hprof):
                f.write((6*"%16.12e " + "\n") % (xi, rhoi, Pi, vi, epsi, hi))
            f.close()

        #---------------------------------------------------------------------------
        # Interactively plot the density if requested.
        #---------------------------------------------------------------------------
        if graphics:
            rhoPlots.append(plotFieldList(db.fluidMassDensity,
                                          plotStyle = "lines",
                                          winTitle = "time=%4.3f" % control.time()))
            rhoPlots[-1].hardcopy(os.path.join(plotDir, "DB_%i_rho_plot_t=%4.3f.png" % (nx, control.time())), terminal="png")
