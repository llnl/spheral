#ATS:t0 = test(      SELF, "--graphics None --clearDirectories True  --checkError True   --restartStep 20", label="Planar Noh problem -- 1-D (serial)")
#ATS:t1 = testif(t0, SELF, "--graphics None --clearDirectories False --checkError False  --restartStep 20 --restoreCycle 20 --steps 20 --checkRestart True", label="Planar Noh problem -- 1-D (serial) RESTART CHECK")
#ATS:t2 = test(      SELF, "--graphics None --clearDirectories True  --checkError True  --dataDir 'dumps-planar-restartcheck' --restartStep 20", np=2, label="Planar Noh problem -- 1-D (parallel)")
#ATS:t3 = testif(t2, SELF, "--graphics None --clearDirectories False --checkError False --dataDir 'dumps-planar-restartcheck' --restartStep 20 --restoreCycle 20 --steps 20 --checkRestart True", np=2, label="Planar Noh problem -- 1-D (parallel) RESTART CHECK")
#ATS:t4 = test(      SELF, "--graphics None --clearDirectories True  --checkError True  --dataDir 'dumps-planar-reproducing' --domainIndependent True", label="Planar Noh problem -- 1-D (serial reproducing test setup)")
#ATS:t5 = testif(t4, SELF, "--graphics None --clearDirectories False  --checkError True  --dataDir 'dumps-planar-reproducing' --domainIndependent True", np=4, label="Planar Noh problem -- 1-D (4 proc reproducing test)")
#-------------------------------------------------------------------------------
# The Planar Noh test case run in 1-D.
#
# W.F. Noh 1987, JCP, 72, 78-120.
#-------------------------------------------------------------------------------
import os, shutil
from Spheral2d import *
from SpheralTestUtilities import *
from GenerateNodeDistribution2d import *
import DistributeNodes

title("1-D integrated hydro test -- discontinuities")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(KernelConstructor = BSplineKernel,
            nx1 = 100,
            ny1 = 100,
            nPerh = 1.25,
            
            eps1 = 1.0,
            rho1 = 1.0,
            x0 = 0.0,
            x1 = 1.0,
            y0 = 0.0,
            y1 = 1.0,
            xwall = 1.0,
            
            dfx1 = 0,   #positions of delta functions
            dfx2 = 0.5,

            gamma = 5.0/3.0,
            mu = 1.0,
            
            SVPH = False,
            CRKSPH = False,
            SPH = True,
            Qconstructor = MonaghanGingoldViscosity,
            #Qconstructor = TensorMonaghanGingoldViscosity,
            
            linearConsistent = False,  #what does this do?
            fcentroidal = 0.0,
            fcellPressure = 0.0,
            Qhmult = 1.0,
            Cl = 1.0,
            Cq = 1.0,
            Qlimiter = False,
            epsilon2 = 1e-2,
            hmin = 0.0001,
            hmax = 0.1,
            cfl = 0.5,
            XSPH = True,
            epsilonTensile = 0.3,
            nTensile = 4.0,
            hourglass = None,
            hourglassOrder = 0,
            hourglassLimiter = 0,
            hourglassFraction = 0.5,
            filter = 0.0,
            
            IntegratorConstructor = CheapSynchronousRK2Integrator,
            goalTime = 0.6,
            steps = 100,
            vizCycle = None,
            vizTime = 0.1,
            dt = 0.0001,
            dtMin = 1.0e-5,
            dtMax = 0.1,
            dtGrowth = 2.0,
            rigorousBoundaries = False,
            updateBoundaryFrequency = 1,
            maxSteps = None,
            statsStep = 1,
            smoothIters = 0,
            HUpdate = IdealH,
            #densityUpdate = RigorousSumDensity, # VolumeScaledDensity,
            densityUpdate = IntegrateDensity,
            compatibleEnergy = True,
            gradhCorrection = True,
            domainIndependent = True,
            cullGhostNodes = True,
            
            bArtificialConduction = True,
            arCondAlpha = 0.5,
            
            useVoronoiOutput = False,
            clearDirectories = True,
            checkError = True,
            checkRestart = False,
            checkEnergy = True,
            restoreCycle = None,
            restartStep = 10000,
            redistributeStep = 200,
            dataDir = "dumps-planar",
            restartBaseName = "Noh-planar-1d",

            # Parameters for the test
            scalePressure = 5.0,
            scaleEnergy = 2.0,
            showDecay = False,
            zerovpkg = False,
            
            graphics = "gnu",
            periodic = True
            )

restartDir = os.path.join(dataDir, "restarts")
restartBaseName = os.path.join(restartDir, "discontinuity-1d-%i" % nx1)
vizBaseName = "DiscontinuityTest2d"
vizDir = os.path.join(dataDir, "visit")

dx = (x1 - x0)/nx1

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
import os, sys
if mpi.rank == 0:
    if clearDirectories and os.path.exists(dataDir):
        shutil.rmtree(dataDir)
    if not os.path.exists(restartDir):
        os.makedirs(restartDir)
mpi.barrier()

#-------------------------------------------------------------------------------
# Material properties.
#-------------------------------------------------------------------------------
eos = GammaLawGasMKS(gamma, mu)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
WT = TableKernel(KernelConstructor(), 1000)
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
if restoreCycle is None:
    generator1 = GenerateNodeDistribution2d(nx1, ny1,
                                            rho = rho1,
                                            distributionType = "lattice",
                                            xmin = (0.0,  0.0),
                                            xmax = (1.0,  1.0),
                                            nNodePerh = nPerh,
                                            SPH = SPH)
    if mpi.procs > 1:
        from VoronoiDistributeNodes import distributeNodes2d
    else:
        from DistributeNodes import distributeNodes2d

    distributeNodes2d((nodes1,generator1))


# Set node specific thermal energies
nodes1.specificThermalEnergy(ScalarField("tmp", nodes1, eps1))
#nodes1.massDensity(ScalarField("tmp", nodes1, rho1))

nodeSet = [nodes1]

# Set node specific thermal energies
xeps = x1/nx1
for nodes in nodeSet:
    pos = nodes.positions()
    eps = nodes.specificThermalEnergy()
    rho = nodes.massDensity()
    vel = nodes.velocity()
    for i in range(nodes.numInternalNodes):
        x = pos[i].x
        y = pos[i].y
        vel[i][0] = sin(2*3.14159*y)
        if ((abs(x-dfx2) < xeps and abs(y-dfx2) < xeps)):
            eps[i] = eps1 + scaleEnergy
            #vel[i][1] = 2.0
        else:
            eps[i] = eps1
            
#rho[i] = specificEnergy(pos[i].x)


#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase()
output("db")
output("db.appendNodeList(nodes1)")
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Construct the artificial viscosity.
#-------------------------------------------------------------------------------
q = Qconstructor(Cl, Cq)
q.epsilon2 = epsilon2
q.limiter = Qlimiter
output("q")
output("q.Cl")
output("q.Cq")
output("q.epsilon2")
output("q.limiter")

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
if SVPH:
    hydro = SVPHFacetedHydro(W = WT, 
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
                             xmin = Vector(-100.0),
                             xmax = Vector( 100.0))
elif CRKSPH:
    hydro = CRKSPHHydro(W = WT, 
                        Q = q,
                        filter = filter,
                        cfl = cfl,
                        compatibleEnergyEvolution = compatibleEnergy,
                        XSPH = XSPH,
                        densityUpdate = densityUpdate,
                        HUpdate = HUpdate)
else:
    hydro = SPHHydro(W = WT, 
                     Q = q,
                     cfl = cfl,
                     compatibleEnergyEvolution = compatibleEnergy,
                     gradhCorrection = gradhCorrection,
                     densityUpdate = densityUpdate,
                     HUpdate = HUpdate,
                     XSPH = XSPH,
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
# zero velocity package
#-------------------------------------------------------------------------------


if zerovpkg:
    class zeroV_pkg(Physics):
        def __init__(self):
            Physics.__init__(self)
            return
        
        def evaluateDerivatives(self, t, dt, db, state, derivs):
            DepsDt = derivs.scalarFields("delta " + HydroFieldNames.specificThermalEnergy)
            DvDt = derivs.vectorFields("delta " + HydroFieldNames.velocity)
            DepsDt.Zero()
            DvDt.Zero()
            return
        
        def dt(self, db, state, derivs, t):
            return pair_double_string(1e100, "No vote")
        
        def registerState(self, dt, state):
            return
        
        def registerDerivatives(self, db, derivs):
            return
        
        def label(self):
            return "zeroV package"
        
        def initialize(self, t, dt, db, state, derivs):
            return

    zv = zeroV_pkg()

    packages.append(zv)

#-------------------------------------------------------------------------------
# Construct the Artificial Conduction physics object.
#-------------------------------------------------------------------------------

if bArtificialConduction:
    #q.reducingViscosityCorrection = True
    ArtyCond = ArtificialConduction(WT,arCondAlpha)
    
    packages.append(ArtyCond)

#-------------------------------------------------------------------------------
# Optionally construct an hourglass control object.
#-------------------------------------------------------------------------------
if hourglass:
    mask = db.newFluidIntFieldList(1, "mask")
    pos = nodes1.positions()
    for i in range(nodes1.numInternalNodes):
        if pos[i].x > (x1 - dx):
            mask[0][i] = 0
    hg = hourglass(WT,
                   order = hourglassOrder,
                   limiter = hourglassLimiter,
                   fraction = hourglassFraction,
                   mask = mask)
    output("hg")
    output("hg.order")
    output("hg.limiter")
    output("hg.fraction")
    packages.append(hg)

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
xp1 = Plane(Vector(0.0, 0.0), Vector( 1.0, 0.0))
xp2 = Plane(Vector(1.0, 0.0), Vector(-1.0, 0.0))
yp1 = Plane(Vector(0.0, 0.0), Vector(0.0,  1.0))
yp2 = Plane(Vector(0.0, 1.0), Vector(0.0, -1.0))
xbc = PeriodicBoundary(xp1, xp2)
ybc = PeriodicBoundary(yp1, yp2)
ybc1 = ReflectingBoundary(yp1)
ybc2 = ReflectingBoundary(yp2)
if periodic:
    bcSet = [xbc,ybc]
else:
    bcSet = [xbc, ybc1, ybc2]

for p in packages:
    for bc in bcSet:
        p.appendBoundary(bc)

#-------------------------------------------------------------------------------
# Construct an integrator.
#-------------------------------------------------------------------------------
integrator = IntegratorConstructor(db)
for p in packages:
    integrator.appendPhysicsPackage(p)
del p
integrator.lastDt = dt
integrator.dtMin = dtMin
integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
integrator.rigorousBoundaries = rigorousBoundaries
integrator.updateBoundaryFrequency = updateBoundaryFrequency
integrator.domainDecompositionIndependent = domainIndependent
integrator.cullGhostNodes = cullGhostNodes
output("integrator")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.dtGrowth")
output("integrator.rigorousBoundaries")
output("integrator.updateBoundaryFrequency")
output("integrator.domainDecompositionIndependent")
output("integrator.cullGhostNodes")

#-------------------------------------------------------------------------------
# Make the problem controller.
#-------------------------------------------------------------------------------
if useVoronoiOutput:
    import SpheralVoronoiSiloDump
    vizMethod = SpheralVoronoiSiloDump.dumpPhysicsState
else:
    import SpheralVisitDump
    vizMethod = SpheralVisitDump.dumpPhysicsState
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle,
                            redistributeStep = redistributeStep,
                            vizMethod = vizMethod,
                            vizBaseName = vizBaseName,
                            vizDir = vizDir,
                            vizStep = vizCycle,
                            vizTime = vizTime,
                            SPH = SPH)
output("control")

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
if not steps is None:
    control.step(steps)
else:
    control.advance(goalTime,maxSteps)
    control.updateViz(control.totalSteps, integrator.currentTime, 0.0)
    control.dropRestartFile()
#control.advance(goalTime, maxSteps)


#-------------------------------------------------------------------------------
# Plot the final state.
#-------------------------------------------------------------------------------
if graphics == "matplot":
    import pylab
    from SpheralMatplotlibUtilities import *
    rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotState(db)
    plotAnswer(answer, control.time(), rhoPlot, velPlot, epsPlot, PPlot, HPlot)
    plotEHistory(control.conserve)

elif graphics == "gnu":
    from SpheralGnuPlotUtilities import *



    #EPlot = plotEHistory(control.conserve)
    
    
    if showDecay:
        dudtPlot = generateNewGnuPlot()
        dudtData = Gnuplot.Data(timeArray,eps50,with_ = "points", title="eps50")
        dudtPlot.plot(dudtData)
        #dudtPlot('set logscale y')
        dudtPlot.refresh()
    rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotState(db)
#plotAnswer(control.time(), rhoPlot, velPlot, epsPlot, PPlot, HPlot)


Eerror = (control.conserve.EHistory[-1] - control.conserve.EHistory[0])/control.conserve.EHistory[0]
print("Total energy error: %g" % Eerror)
if checkEnergy and abs(Eerror) > 1e-13:
    raise ValueError("Energy error outside allowed bounds.")
