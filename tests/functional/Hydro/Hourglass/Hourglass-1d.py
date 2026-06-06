#ATS:test(SELF, "--graphics False", label="Planar Hourglass test problem -- 1-D (serial)")
#-------------------------------------------------------------------------------
# A made up 1-D problem to test the anti-hourglassing algorithms.
#-------------------------------------------------------------------------------
from Spheral1d import *
from SpheralTestUtilities import *
from NodeHistory import *

title("1-D planar hourglassing test")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(nx1 = 20,
            rho1 = 1.0,
            eps1 = 1.0,
            x0 = 0.0,
            x1 = 1.0,
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

            restoreCycle = None,
            restartStep = 10000,
            restartBaseName = "Hourglass-1d",

            graphics = True,
            )

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
from DistributeNodes import distributeNodesInRange1d
distributeNodesInRange1d([(nodes1, nx1, rho1, (x0, x1))])
output("nodes1.numNodes")

# Set node specific thermal energies
nodes1.specificThermalEnergy(ScalarField("tmp", nodes1, eps1))

# Displace the nodes in a pattern that looks like the tensile instability clumping
dx = (x1 - x0)/nx1
pos = nodes1.positions()
for i in range(nodes1.numInternalNodes):
    delta = amplitude*((-1.0)**(i % 2))*dx # amplitude*sin(2.0*pi*nodes1.positions()[i].x/wavelength)
    pos[i].x += delta

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
xPlane0 = Plane(Vector(x0), Vector( 1.0))
xPlane1 = Plane(Vector(x1), Vector(-1.0))
xbc0 = ReflectingBoundary(xPlane0)
xbc1 = ReflectingBoundary(xPlane1)
for p in packages:
    p.appendBoundary(xbc0)
    p.appendBoundary(xbc1)

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
    return pos[i].x, vel[i].x, DvDt_hydro(0,i).x, DvDt_hg(0,i).x

histories = [NodeHistory(nodes1, [i], samplefunc, "/dev/null") for i in range(nx1)]

#-------------------------------------------------------------------------------
# Make the problem controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
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

#-------------------------------------------------------------------------------
# Plot the final state.
#-------------------------------------------------------------------------------
from SpheralMatplotlib import *
rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotState(db)
EPlot = plotEHistory(control.conserve)
a = hydro.DvDt
aplot = plotFieldList(a,
                      yFunction = "%s.x",
                      winTitle = "Acceleration (Hydro)")
ahg = hg.DvDt
ahgplot = plotFieldList(ahg,
                        yFunction = "%s.x",
                        winTitle = "Acceleration (Anti-hourglass)")

def computeNearestNeighborDistance():
    result = ScalarField("nearest neighbor distance", nodes1, 1e10)
    db.updateConnectivityMap()
    cm = db.connectivityMap()
    pairs = cm.nodePairList
    pos = nodes1.positions()
    for pair in pairs:
        i, j = pair.i_node, pair.j_node
        rij = (pos(i) - pos(j)).magnitude()
        result[i] = min(result[i], rij)
        result[j] = min(result[j], rij)
    return result

def plotit(x, y,
           style = "ro",
           title = None,
           xlabel = None,
           ylabel = None):
    fig = newFigure()
    plt.plot(x, y, style)
    if title:  plt.title(title)
    if xlabel: plt.xlabel(xlabel)
    if ylabel: plt.ylabel(ylabel)
    return fig

nearestPlot = plotField(computeNearestNeighborDistance(),
                        winTitle = "Nearest neighbor distance",
                        xlabel = "x",
                        ylabel = "min(d)")

x0hist = plotit(histories[0].timeHistory, [s[0] for s in histories[0].sampleHistory],
                title = "Node 0 position",
                xlabel = "time",
                ylabel = "x")
v0hist = plotit(histories[0].timeHistory, [s[1] for s in histories[0].sampleHistory],
                title = "Node 0 velocity",
                xlabel = "time",
                ylabel = "vel")
a0hist = plotit(histories[0].timeHistory, [s[2] for s in histories[0].sampleHistory],
                title = "Node 0 acceleration (hydro)",
                xlabel = "time",
                ylabel = "accel")
hg0hist = plotit(histories[0].timeHistory, [s[3] for s in histories[0].sampleHistory],
                 title = "Node 0 acceleration (anti-HG)",
                 xlabel = "time",
                 ylabel = "accel")

plots = [(rhoPlot, "Hourglass-1d-rho.png"),
         (velPlot, "Hourglass-1d-vel.png"),
         (epsPlot, "Hourglass-1d-eps.png"),
         (PPlot, "Hourglass-1d-P.png"),
         (HPlot, "Hourglass-1d-h.png"),
         (aplot, "Hourglass-1d-accel.png")]
