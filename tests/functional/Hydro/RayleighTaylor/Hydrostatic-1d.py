#-------------------------------------------------------------------------------
# This is a test of maintaining hydrostatic equilibrium with constant boundaries,
# akin to how we handle the Rayleigh-Taylor problem.
#-------------------------------------------------------------------------------
import shutil
from math import *
from Spheral1d import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from findLastRestart import *

import mpi
import DistributeNodes

class ExponentialDensity:
    def __init__(self,
                 rho1,
                 rho2,
                 delta):
        self.rho1 = rho1
        self.rho2 = rho2
        self.delta = delta
        return
    def __call__(self, r):
        return self.rho1+(self.rho2-self.rho1)/(1+exp(-(r-0.5)/delta))

title("Rayleigh-Taylor hydrostatic equilibrium test problem in 1D")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(nx1 = 128,
            rhoT = 2.0,
            rhoB = 1.0,
            x0 = 0.0,
            x1 = 1.0,
            gval = -0.5,
	    w0  = 0.025,
            delta = 0.025, 
            gamma = 1.4,
            mu = 1.0,
            
            nPerh = 2.01,
            
            SVPH = False,
            CRKSPH = False,
            SPH = True,   # This just chooses the H algorithm -- you can use this with CRKSPH for instance.
            filter = 0.0,   # CRKSPH filtering
            Qconstructor = MonaghanGingoldViscosity,
            #Qconstructor = TensorMonaghanGingoldViscosity,
            linearConsistent = False,
            fcentroidal = 0.0,
            fcellPressure = 0.0,
            nh = 5.0,
            aMin = 0.1,
            aMax = 2.0,
            Qhmult = 1.0,
            Cl = 1.0,
            Cq = 1.0,
            linearInExpansion = False,
            Qlimiter = False,
            balsaraCorrection = False,
            epsilon2 = 1e-2,
            hmin = 0.0001,
            hmax = 0.5,
            hminratio = 0.1,
            cfl = 0.5,
            useVelocityMagnitudeForDt = False,
            XSPH = False,
            epsilonTensile = 0.0,
            nTensile = 8,
            
            IntegratorConstructor = CheapSynchronousRK2Integrator,
            goalTime = 20.0,
            steps = None,
            vizCycle = None,
            vizTime = 0.01,
            dt = 0.0001,
            dtMin = 1.0e-8,
            dtMax = 0.1,
            dtGrowth = 2.0,
            maxSteps = None,
            statsStep = 10,
            smoothIters = 0,
            HUpdate = IdealH,
            domainIndependent = False,
            rigorousBoundaries = False,
            dtverbose = False,
            
            densityUpdate = RigorousSumDensity, # VolumeScaledDensity,
            compatibleEnergy = True,            # <--- Important!  rigorousBoundaries does not work with the compatibleEnergy algorithm currently.
            gradhCorrection = False,
            
            clearDirectories = False,
            restoreCycle = -1,
            restartStep = 100,
            redistributeStep = 500,
            checkRestart = False,
            dataDir = "dumps-Rayleigh-Taylor-1d_hopkins",

            serialDump = False, #whether to dump a serial ascii file at the end for viz
            graphics = True,
            
            bArtificialConduction = False,
            arCondAlpha = 0.5,
            )

# Decide on our hydro algorithm.
if SVPH:
    HydroConstructor = SVPHFacetedHydro
elif CRKSPH:
    Qconstructor = LimitedMonaghanGingoldViscosity
    HydroConstructor = CRKSPHHydro
else:
    HydroConstructor = SPHHydro

dataDir = os.path.join(dataDir,
                       "gval=%g" % (gval),
                       str(HydroConstructor).split("'")[1].split(".")[-1],
                       "densityUpdate=%s" % (densityUpdate),
                       "XSPH=%s" % XSPH,
                       "filter=%s" % filter,
                       "compatible=%s" % compatibleEnergy,
                       "%s-Cl=%g-Cq=%g" % (str(Qconstructor).split("'")[1].split(".")[-1], Cl, Cq),
                       "%i" % nx1,
                       "nPerh=%g-Qhmult=%g" % (nPerh, Qhmult))
restartDir = os.path.join(dataDir, "restarts")
vizDir = os.path.join(dataDir, "visit")
restartBaseName = os.path.join(restartDir, "Rayleigh-Taylor-1d")
vizBaseName = "Rayleigh-Taylor-1d"

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
import os, sys
if mpi.rank == 0:
    if clearDirectories and os.path.exists(dataDir):
        shutil.rmtree(dataDir)
    if not os.path.exists(restartDir):
        os.makedirs(restartDir)
    if not os.path.exists(vizDir):
        os.makedirs(vizDir)
mpi.barrier()

#-------------------------------------------------------------------------------
# Material properties.
#-------------------------------------------------------------------------------
eos = GammaLawGasMKS(gamma, mu)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
WT = TableKernel(BSplineKernel(), 10000)
output("WT")
kernelExtent = WT.kernelExtent

#-------------------------------------------------------------------------------
# Make the NodeList.
#-------------------------------------------------------------------------------
nodes = makeFluidNodeList("gas", eos,
                           hmin = hmin,
                           hmax = hmax,
                           hminratio = hminratio,
                           nPerh = nPerh)
output("nodes.name")
output("nodes.hmin")
output("nodes.hmax")
output("nodes.hminratio")
output("nodes.nodesPerSmoothingScale")

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
# Add some points above and below the problem to represent the infinite atmosphere.
nxbound = 20
dx = (x1 - x0)/nx1
from DistributeNodes import distributeNodesInRange1d
distributeNodesInRange1d([(nodes, nx1 + 2*nxbound, rhoT,
                           (x0 - nxbound*dx,
                            x1 + nxbound*dx))],
                         nPerh = nPerh)

#Set IC
eps = nodes.specificThermalEnergy()
pos = nodes.positions()
rho = nodes.massDensity()
mass = nodes.mass()
rhoFunc = ExponentialDensity(rhoB,
                             rhoT,
                             delta)
for i in range(nodes.numInternalNodes):
    xi = pos[i].x
    P0 = rhoT/gamma
    rho[i] = rhoFunc(xi)
    mass[i] = dx*rho[i]
    Pi = P0 + gval*rho[i]*(xi-0.5)
    eps0 = Pi/((gamma - 1.0)*rho[i])
    eps[i] = eps0

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase()
output("db")
db.appendNodeList(nodes)
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Construct the artificial viscosity.
#-------------------------------------------------------------------------------
q = Qconstructor(Cl, Cq, linearInExpansion)
q.epsilon2 = epsilon2
q.limiter = Qlimiter
q.balsaraShearCorrection = balsaraCorrection
output("q")
output("q.Cl")
output("q.Cq")
output("q.epsilon2")
output("q.limiter")
output("q.balsaraShearCorrection")
output("q.linearInExpansion")
output("q.quadraticInExpansion")

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
if SVPH:
    hydro = HydroConstructor(W = WT,
                             Q = q,
                             cfl = cfl,
                             useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
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
elif CRKSPH:
    hydro = HydroConstructor(W = WT,
                             Q = q,
                             filter = filter,
                             cfl = cfl,
                             useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                             compatibleEnergyEvolution = compatibleEnergy,
                             XSPH = XSPH,
                             densityUpdate = densityUpdate,
                             HUpdate = HUpdate)
else:
    hydro = HydroConstructor(W = WT,
                             Q = q,
                             cfl = cfl,
                             useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                             compatibleEnergyEvolution = compatibleEnergy,
                             gradhCorrection = gradhCorrection,
                             XSPH = XSPH,
                             densityUpdate = densityUpdate,
                             HUpdate = HUpdate,
                             epsTensile = epsilonTensile,
                             nTensile = nTensile)
output("hydro")
output("hydro.kernel()")
output("hydro.PiKernel()")
output("hydro.cfl")
output("hydro.compatibleEnergyEvolution")
output("hydro.densityUpdate")
output("hydro.HEvolution")

packages = [hydro]

#-------------------------------------------------------------------------------
# Construct the Artificial Conduction physics object.
#-------------------------------------------------------------------------------
if bArtificialConduction:
    #q.reducingViscosityCorrection = True
    ArtyCond = ArtificialConduction(WT,arCondAlpha)
    
    packages.append(ArtyCond)

#-------------------------------------------------------------------------------
# Construct the gravitational acceleration object.
#-------------------------------------------------------------------------------
pos = nodes.positions()
nodeIndices = vector_of_int()
for i in range(nodes.numInternalNodes):
    if pos[i].x > x0 and pos[i].x < x1:
        nodeIndices.append(i)

gravity = ConstantAcceleration1d(Vector1d(gval),
                                  nodes,
                                  nodeIndices)

packages.append(gravity)

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
xp1 = Plane(Vector(x0), Vector( 1.0))
xp2 = Plane(Vector(x1), Vector(-1.0))

# The x boundary will be a snapshot of the state of the points above and below
# the x-cutoffs.
pos = nodes.positions()
xlow, xhigh = vector_of_int(), vector_of_int()
for i in range(nodes.numInternalNodes):
    if pos[i].x < x0:
        xlow.append(i)
    elif pos[i].x > x1:
        xhigh.append(i)
xbc1 = ConstantBoundary(nodes, xlow, xp1)
xbc2 = ConstantBoundary(nodes, xhigh, xp2)

bcSet = [xbc1, xbc2]

for bc in bcSet:
    for p in packages:
        p.appendBoundary(bc)
del bc

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
                            initializeDerivatives = True,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle,
                            redistributeStep = None,
                            SPH = SPH)
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

if graphics:
    from SpheralGnuPlotUtilities import *
    rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotState(db, plotGhosts=True)

if serialDump:
    procs = mpi.procs
    rank = mpi.rank
    serialData = []
    i,j = 0,0
    for i in range(procs):
        if rank == i:
            for j in range(nodes.numInternalNodes):
                serialData.append([nodes.positions()[j],3.0/(nodes.Hfield()[j].Trace()),nodes.mass()[j],nodes.massDensity()[j],nodes.specificThermalEnergy()[j]])
    serialData = mpi.reduce(serialData,mpi.SUM)
    if rank == 0:
        f = open(dataDir + "/serialDump.ascii",'w')
        for i in range(len(serialData)):
            f.write("{0} {1} {2} {3} {4} {5} {6} {7}\n".format(i,serialData[i][0][0],serialData[i][0][1],0.0,serialData[i][1],serialData[i][2],serialData[i][3],serialData[i][4]))
        f.close()
