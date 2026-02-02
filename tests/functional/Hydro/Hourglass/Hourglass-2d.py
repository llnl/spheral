#ATS:test(SELF, "--graphics False", label="Hourglass test problem -- 2-D (serial)")
#-------------------------------------------------------------------------------
# A made up 2-D problem to test the anti-hourglassing algorithms.
#-------------------------------------------------------------------------------
import shutil
from Spheral2d import *
from SpheralTestUtilities import *
from NodeHistory import *

title("2-D hourglassing test")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(nx = 10,
            ny = 10,
            rho1 = 1.0,
            eps1 = 1.0,
            xmin = (0.0, 0.0),
            xmax = (1.0, 1.0),
            nPerh = 4.01,
            
            gamma = 5.0/3.0,
            mu = 1.0,

            wavelength = 0.05,
            amplitude = 0.25,

            hydroType = "SPH",
            fhourglass = 0.0,
            xhourglass = 0.0,
            filter = 0.0,

            hmin = 0.0001, 
            hmax = 100.0,
            cfl = 0.25,
            XSPH = False,
            epsilonTensile = 0.0,
            nTensile = 8,

            goalTime = 1.0,
            steps = None,
            dt = 0.0001,
            dtMin = 1.0e-5, 
            dtMax = None,
            dtGrowth = 2.0,
            dtverbose = False,
            maxSteps = None,
            statsStep = 1,
            smoothIters = 0,
            HUpdate = IdealH,
            useVelocityMagnitudeForDt = False,
            evolveTotalEnergy = False,         # Only for SPH variants -- evolve total rather than specific energy
            densityUpdate = RigorousSumDensity, # VolumeScaledDensity,
            compatibleEnergy = True,
            gradhCorrection = False,
            domainIndependent = False,
            clearDirectories = False,

            restoreCycle = None,
            restartStep = 10000,
            restartBaseName = "Hourglass-2d",

            vizTime = None,
            vizCycle = 1,
            vizDerivs = False,
            dataDir = "dumps-HourGlass-2d",

            graphics = True,
            )

dataDir = os.path.join(dataDir,
                       f"amplitude={amplitude}",
                       f"nPerh={nPerh}",
                       f"compatibleEnergy={compatibleEnergy}",
                       f"filter={filter}",
                       f"fhourglass={fhourglass}",
                       f"nx={nx}_ny={ny}")
restartDir = os.path.join(dataDir, "restarts")
restartBaseName = os.path.join(restartDir, f"Hourglass-2d-{nx}x{ny}")

vizDir = os.path.join(dataDir, "visit")
if vizTime is None and vizCycle is None:
    vizBaseName = None
else:
    vizBaseName = f"Hourglass-2d-{nx}x{ny}"

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
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
WT = TableKernel(WendlandC4Kernel(), 100)
output("WT")

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
from DistributeNodes import distributeNodes2d
from GenerateNodeDistribution2d import GenerateNodeDistribution2d
generator = GenerateNodeDistribution2d(nx, ny, rho1, "lattice",
                                       xmin = xmin,
                                       xmax = xmax,
                                       nNodePerh = nPerh)
distributeNodes2d((nodes1, generator))
output("nodes1.numNodes")

# Set node specific thermal energies
nodes1.specificThermalEnergy(ScalarField("tmp", nodes1, eps1))

# Displace the nodes in a pattern that looks like hourglassing.
dx = (xmax[0] - xmin[0])/nx
dy = (xmax[1] - xmin[1])/ny
pos = nodes1.positions()
for i in range(nodes1.numInternalNodes):
    ix = int((pos[i].x - xmin[0])/dx)
    iy = int((pos[i].y - xmin[1])/dy)
    delta = amplitude * Vector((-1)**(ix % 2) * dx,
                               (-1)**(iy % 2) * dy)
    pos[i] += delta

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase()
output("db")
output("db.appendNodeList(nodes1)")
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
if hydroType == "SVPH":
    hydro = SVPH(dataBase = db,
                 W = WT,
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
                 xmin = Vector(-100.0),
                 xmax = Vector( 100.0))
elif hydroType == "CRKSPH":
    hydro = CRKSPH(dataBase = db,
                   W = WT,
                   order = correctionOrder,
                   filter = filter,
                   cfl = cfl,
                   useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                   compatibleEnergyEvolution = compatibleEnergy,
                   evolveTotalEnergy = evolveTotalEnergy,
                   XSPH = XSPH,
                   densityUpdate = densityUpdate,
                   HUpdate = HUpdate,
                   crktype = crktype)
elif hydroType == "PSPH":
    hydro = PSPH(dataBase = db,
                 W = WT,
                 filter = filter,
                 cfl = cfl,
                 useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                 compatibleEnergyEvolution = compatibleEnergy,
                 evolveTotalEnergy = evolveTotalEnergy,
                 densityUpdate = densityUpdate,
                 HUpdate = HUpdate,
                 XSPH = XSPH)

elif hydroType == "FSISPH":
    hydro = FSISPH(dataBase = db,
                   W = WT,
                   cfl = cfl,
                   interfaceMethod = HLLCInterface,
                   sumDensityNodeLists=[nodes1],                       
                   densityStabilizationCoefficient = 0.00,
                   useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                   compatibleEnergyEvolution = compatibleEnergy,
                   evolveTotalEnergy = evolveTotalEnergy,
                   HUpdate = HUpdate)
elif hydroType == "GSPH":
    limiter = VanLeerLimiter()
    waveSpeed = DavisWaveSpeed()
    solver = HLLC(limiter,
                  waveSpeed,
                  True)
    hydro = GSPH(dataBase = db,
                riemannSolver = solver,
                W = WT,
                cfl=cfl,
                useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                compatibleEnergyEvolution = compatibleEnergy,
                evolveTotalEnergy = evolveTotalEnergy,
                XSPH = XSPH,
                gradientType = gsphReconstructionGradient,
                densityUpdate=densityUpdate,
                HUpdate = IdealH,
                epsTensile = epsilonTensile,
                nTensile = nTensile)
elif hydroType == "MFM":
    limiter = VanLeerLimiter()
    waveSpeed = DavisWaveSpeed()
    solver = HLLC(limiter,
                  waveSpeed,
                  True)
    hydro = MFM(dataBase = db,
                riemannSolver = solver,
                W = WT,
                cfl=cfl,
                useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                compatibleEnergyEvolution = compatibleEnergy,
                evolveTotalEnergy = evolveTotalEnergy,
                XSPH = XSPH,
                gradientType = gsphReconstructionGradient,
                densityUpdate=densityUpdate,
                HUpdate = IdealH,
                epsTensile = epsilonTensile,
                nTensile = nTensile)
else:
    assert hydroType == "SPH"
    hydro = SPH(dataBase = db,
                W = WT,
                filter = filter,
                cfl = cfl,
                useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                compatibleEnergyEvolution = compatibleEnergy,
                evolveTotalEnergy = evolveTotalEnergy,
                gradhCorrection = gradhCorrection,
                densityUpdate = densityUpdate,
                HUpdate = HUpdate,
                XSPH = XSPH,
                epsTensile = epsilonTensile,
                nTensile = nTensile)
output("hydro")
try:
    output("hydro.kernel")
    output("hydro.PiKernel")
    output("hydro.XSPH")
except:
    pass
output("hydro.cfl")
output("hydro.compatibleEnergyEvolution")
output("hydro.densityUpdate")

packages = [hydro]

#-------------------------------------------------------------------------------
# Optionally construct an hourglass control object.
#-------------------------------------------------------------------------------
hg = SubPointPressureHourglassControl(fhourglass, xhourglass)
output("hg")
output("hg.fHG")
output("hg.xfilter")
packages.append(hg)

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
xbc0 = ReflectingBoundary2d(Plane2d(Vector2d(*xmin), Vector2d(1.0, 0.0)))
xbc1 = ReflectingBoundary2d(Plane2d(Vector2d(*xmax), Vector2d(-1.0, 0.0)))
ybc0 = ReflectingBoundary2d(Plane2d(Vector2d(*xmin), Vector2d(0.0, 1.0)))
ybc1 = ReflectingBoundary2d(Plane2d(Vector2d(*xmax), Vector2d(0.0, -1.0)))

for bc in [xbc0, xbc1, ybc0, ybc1]:
    for p in packages:
        p.appendBoundary(bc)

#-------------------------------------------------------------------------------
# Construct a predictor corrector integrator, and add the one physics package.
#-------------------------------------------------------------------------------
integrator = CheapSynchronousRK2Integrator(db)
for p in packages:
    integrator.appendPhysicsPackage(p)
integrator.lastDt = dt
if dtMin:
    integrator.dtMin = dtMin
if dtMax:
    integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
integrator.verbose = dtverbose
output("integrator")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.dtGrowth")
output("integrator.verbose")

#-------------------------------------------------------------------------------
# Track the history of the motion of select points
#-------------------------------------------------------------------------------
def samplefunc(nodes, indices):
    i = indices[0]
    pos = nodes.positions()
    vel = nodes.velocity()
    DvDt_hydro = hydro.DvDt
    DvDt_hg = hg.DvDt
    return pos[i].x, pos[i].y, vel[i].x, vel[i].y, vel[i].magnitude(), DvDt_hydro(0,i).x, DvDt_hydro(0,i).y, DvDt_hydro(0,i).magnitude(), DvDt_hg(0,i).x, DvDt_hg(0,i).y, DvDt_hg(0,i).magnitude()

histories = [NodeHistory(nodes1, [i], samplefunc, os.path.join(dataDir, f"NodeHistory{i}")) for i in range(nx)]

#-------------------------------------------------------------------------------
# Make the problem controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            vizBaseName = vizBaseName,
                            vizDir = vizDir,
                            vizStep = vizCycle,
                            vizTime = vizTime,
                            vizDerivs = vizDerivs,
                            periodicWork = [(hist, 1) for hist in histories])
output("control")

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
if steps is None:
    if control.time() < goalTime:
        control.step(5)
        control.advance(goalTime, maxSteps)
else:
    control.step(steps)

Eerror = (control.conserve.EHistory[-1] - control.conserve.EHistory[0])/control.conserve.EHistory[0]
print("Total energy error: %g" % Eerror)

# #-------------------------------------------------------------------------------
# # Plot the final state.
# #-------------------------------------------------------------------------------
# import Gnuplot
# from SpheralGnuPlotUtilities import *
# rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotState(db)
# Eplot = plotEHistory(control.conserve)
# xplot = plotNodePositions2d(db,
#                             title = "Positions")
