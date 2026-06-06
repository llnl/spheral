import os, sys
import shutil
from SolidSpheral2d import *
from SpheralTestUtilities import *
from SodAnalyticSolution import SodSolutionGasGas as SodSolution

title("2-D integrated hydro test -- planar Sod problem")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(nx1 = 200,
            ny1 = 80,

            nx2 = 100,
            ny2 = 40,

            initialRotation = 0.0,     # Degrees, optionally rotate and clip the initial positions

            rho1 = 1.0,
            rho2 = 0.25,
            P1 = 1.0,
            P2 = 0.1795,
            gamma1 = 5.0/3.0,
            gamma2 = 5.0/3.0,

            numNodeLists = 1,

            x0 = -0.5,
            x1 = 0.0,
            x2 = 0.5,

            y0 = 0.0,
            y1 = 0.2,
            yplotbuf = 0.0,

            hsmooth = 0.5,             # Optionally smooth initial discontinuity
            sumInitialDensity = False, # Optionally sum the initial density before setting the pressure and such

            nPerh = 1.25,

            mu = 1.0,
            
            svph = False,
            crksph = False,
            psph = False,
            fsisph = False,
            gsph= False,
            evolveTotalEnergy = False,  # Only for SPH variants -- evolve total rather than specific energy
            solid = False,    # If true, use the fluid limit of the solid hydro option
            boolReduceViscosity = False,
            nh = 5.0,
            aMin = 0.1,
            aMax = 2.0,
            Cl = 1.0,
            Cq = 1.5,
            etaCritFrac = 1.0,
            etaFoldFrac = 0.2,
            boolCullenViscosity = False,
            alphMax = 2.0,
            alphMin = 0.02,
            betaC = 0.7,
            betaD = 0.05,
            betaE = 1.0,
            fKern = 1.0/3.0,
            boolHopkinsCorrection = True,
            HopkinsConductivity = False,
            linearInExpansion = False,
            quadraticInExpansion = False,
            Qlimiter = False,
            epsilon2 = 1e-4,
            hmin = 1e-10,
            hmax = 1.0,
            cfl = 0.5,
            XSPH = False,
            epsilonTensile = 0.0,
            nTensile = 8,
            rhoMin = 0.01,
            filter = 0.00,
            KernelConstructor = BSplineKernel,
            order = 5,
            
            bArtificialConduction = False,
            arCondAlpha = 0.5,

            IntegratorConstructor = CheapSynchronousRK2Integrator,
            dtverbose = False,
            steps = None,
            goalTime = 0.15,
            dt = 1e-6,
            dtMin = 1.0e-6,
            dtMax = 0.1,
            dtGrowth = 2.0,
            rigorousBoundaries = False,
            maxSteps = None,
            statsStep = 10,
            HUpdate = IdealH,
            correctionOrder = LinearOrder,
            volumeType = RKSumVolume,
            densityUpdate = RigorousSumDensity,
            compatibleEnergy = True,
            correctVelocityGradient = True,
            gradhCorrection = True,
            linearConsistent = False,

            clearDirectories = False,
            restoreCycle = -1,
            restartStep = 100,
            dataDirBase = "dumps-Sod-planar-2d",
            restartBaseName = "Sod-planar-2d-restart",
            outputFile = None,
            checkRestart = False,

            vizCycle = None,
            vizTime = 0.01,
            graphics = True,
            )

assert not(boolReduceViscosity and boolCullenViscosity)
assert not (GSPH and (boolReduceViscosity or boolCullenViscosity))
assert numNodeLists in (1, 2)
assert not (fsisph and not solid)
assert not svph

if svph:
    hydroname = "SVPH"
elif crksph:
    hydroname = "CRKSPH"
elif psph:
    hydroname = "PSPH"
elif fsisph:
    hydroname = "FSISPH"
elif gsph:
    hydroname = "GSPH"
else:
    hydroname = "SPH"
if solid:
    hydroname = "Solid"+hydroname


dataDir = os.path.join(dataDirBase,
                       "rotation=%f" % initialRotation,
                       hydroname,
                       "nPerh=%f" % nPerh,
                       "compatibleEnergy=%s" % compatibleEnergy,
                       "correctionOrder=%s" % correctionOrder,
                       "Cullen=%s" % boolCullenViscosity,
                       "Condc=%s" % HopkinsConductivity,
                       "filter=%f" % filter,
                       "%i" % (nx1 + nx2))
restartDir = os.path.join(dataDir, "restarts")
restartBaseName = os.path.join(restartDir, "Sod-planar-2d-%i" % (nx1 + nx2))

vizDir = os.path.join(dataDir, "visit")
if vizTime is None and vizCycle is None:
    vizBaseName = None
else:
    vizBaseName = "Sod-planar-2d-%i" % (nx1 + nx2)

if numNodeLists == 1 and abs(gamma1-gamma2) > 0.001:
    gamma2 = gamma1
    print(" ")
    print("============================================================")
    print("warning: because we're running 1 nodelist, gamma2 had to be ")
    print("         set equal to gamma1")
    print("============================================================")
    print(" ")

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
eos1 = GammaLawGasMKS(gamma1, mu)
eos2 = GammaLawGasMKS(gamma2, mu)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
WT = TableKernel(NBSplineKernel(order), 1000)
kernelExtent = WT.kernelExtent
output("WT")

#-------------------------------------------------------------------------------
# Make the NodeLists.
#-------------------------------------------------------------------------------
if solid:
    makeNL = makeSolidNodeList
else:
    makeNL = makeFluidNodeList
nodes1 = makeNL("nodes1", eos1, 
                hmin = hmin,
                hmax = hmax,
                nPerh = nPerh,
                kernelExtent = kernelExtent,
                rhoMin = rhoMin)
nodes2 = makeNL("nodes2", eos2, 
                hmin = hmin,
                hmax = hmax,
                nPerh = nPerh,
                kernelExtent = kernelExtent,
                rhoMin = rhoMin)
nodeSet = [nodes1, nodes2]

#-------------------------------------------------------------------------------
# Functions to specify the initial density and specific thermal energy with
# optional smoothing.
#-------------------------------------------------------------------------------
dx1 = (x1 - x0)/nx1
dx2 = (x2 - x1)/nx2
hfold = hsmooth*max(dx1, dx2)
def rho_initial(xi):
    if hfold > 0.0:
        return rho1 + (rho2 - rho1)/(1.0 + exp(-(xi - x1)/hfold))
    elif xi <= x1:
        return rho1
    else:
        return rho2

def specificEnergy(xi, rhoi, gammai):
    if hfold > 0.0:
        Pi = P1 + (P2 - P1)/(1.0 + exp((-xi - x1)/hfold))
    elif xi <= x1:
        Pi = P1
    else:
        Pi = P2
    return Pi/((gammai - 1.0)*rhoi)

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
from GenerateNodeDistribution2d import GenerateNodeDistribution2d
from GenerateNodeProfile import GeneratePlanarNodeProfile2d
if initialRotation == 0.0:
    gen1 = GeneratePlanarNodeProfile2d(nx = nx1,
                                       ny = ny1,
                                       rho = rho_initial,
                                       xmin = (x0, y0),
                                       xmax = (x1, y1),
                                       nNodePerh = nPerh,
                                       SPH = True)
    gen2 = GeneratePlanarNodeProfile2d(nx = nx2,
                                       ny = ny2,
                                       rho = rho_initial,
                                       xmin = (x1, y0),
                                       xmax = (x2, y1),
                                       nNodePerh = nPerh,
                                       SPH = True)
else:
    if x1 - x0 > y1 - y0:
        dl = 0.5*(x1 - x0)
        n = nx1
    else:
        dl = 0.5*(y1 - x0)
        n = ny1
    gen1 = GenerateNodeDistribution2d(n, n, rho1, "lattice",
                                      xmin = (-dl, -dl),
                                      xmax = (dl, dl),
                                      rotation = initialRotation * pi/180.0,
                                      #offset = (0.5*(x0 + x1), 0.5*(y0 + y1) + 0.5*dl/n),
                                      offset = (0.5*(x0 + x1), 0.5*(y0 + y1)),
                                      xminreject = (x0, y0),
                                      #xmaxreject = (x1, y1),
                                      xmaxreject = (x1 - 0.001*dl/n, y1 - 0.1*dl/n),
                                      nNodePerh = nPerh,
                                      SPH = True)
    if x2 - x1 > y1 - y0:
        dl = 0.5*(x2 - x1)
        n = nx2
    else:
        dl = 0.5*(y1 - x0)
        n = ny2
    gen2 = GenerateNodeDistribution2d(n, n, rho2, "lattice",
                                      xmin = (-dl, -dl),
                                      xmax = (dl, dl),
                                      rotation = initialRotation * pi/180.0,
                                      offset = (0.5*(x1 + x2), 0.5*(y0 + y1)),
                                      #offset = (0.5*(x1 + x2), 0.5*(y0 + y1) + 0.5*dl/n),
                                      xminreject = (x1, y0),
                                      #xmaxreject = (x2, y1),
                                      xmaxreject = (x2 - 0.001*dl/n, y1 - 0.1*dl/n),
                                      nNodePerh = nPerh,
                                      SPH = True)

if mpi.procs > 1:
    from VoronoiDistributeNodes import distributeNodes2d
else:
    from DistributeNodes import distributeNodes2d

if numNodeLists == 1:
    from CompositeNodeDistribution import CompositeNodeDistribution
    gen = CompositeNodeDistribution(gen1, gen2)
    distributeNodes2d((nodes1, gen))
else:
    distributeNodes2d((nodes1, gen1),
                      (nodes2, gen2))

for n in nodeSet:
    output("n.name")
    output("   mpi.reduce(n.numInternalNodes, mpi.MIN)")
    output("   mpi.reduce(n.numInternalNodes, mpi.MAX)")
    output("   mpi.reduce(n.numInternalNodes, mpi.SUM)")
del n

# Set node specific thermal energies
for nodes in nodeSet:

    pos = nodes.positions()
    eps = nodes.specificThermalEnergy()
    rho = nodes.massDensity()
    gammai = nodes.eos.gamma

    for i in range(nodes.numInternalNodes):
        eps[i] = specificEnergy(pos[i].x, rho[i], gammai)

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
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
if svph:
    hydro = SVPH(dataBase = db,
                 W = WT,
                 cfl = cfl,
                 compatibleEnergyEvolution = compatibleEnergy,
                 XSVPH = XSPH,
                 linearConsistent = linearConsistent,
                 generateVoid = False,
                 densityUpdate = densityUpdate,
                 HUpdate = HUpdate,
                 xmin = Vector(-100.0),
                 xmax = Vector( 100.0))
elif crksph:
    hydro = CRKSPH(dataBase = db,
                   order = correctionOrder,
                   filter = filter,
                   cfl = cfl,
                   compatibleEnergyEvolution = compatibleEnergy,
                   evolveTotalEnergy = evolveTotalEnergy,
                   XSPH = XSPH,
                   densityUpdate = densityUpdate,
                   HUpdate = HUpdate)
elif psph:
    hydro = PSPH(dataBase = db,
                 W = WT,
                 cfl = cfl,
                 compatibleEnergyEvolution = compatibleEnergy,
                 evolveTotalEnergy = evolveTotalEnergy,
                 densityUpdate = densityUpdate,
                 HUpdate = HUpdate,
                 XSPH = XSPH,
                 correctVelocityGradient = correctVelocityGradient)
elif fsisph:
    sumDensityNodeLists = [nodes1]
    if numNodeLists == 2:
        sumDensityNodeLists += [nodes2]
    hydro = FSISPH(dataBase = db,
                   W = WT,
                   cfl = cfl,
                   sumDensityNodeLists=sumDensityNodeLists,                       
                   densityStabilizationCoefficient = 0.1,
                   specificThermalEnergyDiffusionCoefficient = 0.1,
                   interfaceMethod = HLLCInterface,
                   compatibleEnergyEvolution = compatibleEnergy,
                   evolveTotalEnergy = evolveTotalEnergy,
                   linearCorrectGradients = correctVelocityGradient,
                   HUpdate = HUpdate)
elif gsph:
    limiter = VanLeerLimiter()
    waveSpeed = DavisWaveSpeed()
    solver = HLLC(limiter,waveSpeed,True)
    hydro = GSPH(dataBase = db,
                riemannSolver = solver,
                W = WT,
                cfl=cfl,
                compatibleEnergyEvolution = compatibleEnergy,
                correctVelocityGradient=correctVelocityGradient,
                evolveTotalEnergy = evolveTotalEnergy,
                XSPH = XSPH,
                densityUpdate=densityUpdate,
                HUpdate = IdealH,
                epsTensile = epsilonTensile,
                nTensile = nTensile)
else:
    hydro = SPH(dataBase = db,
                W = WT,
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
output("hydro.compatibleEnergyEvolution")
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
    if not linearInExpansion is None:
        q.linearInExpansion = linearInExpansion
    if not quadraticInExpansion is None:
        q.quadraticInExpansion = quadraticInExpansion
    if not Qlimiter is None:
        q.limiter = Qlimiter
    if not epsilon2 is None:
        q.epsilon2 = epsilon2
    if not etaCritFrac is None:
        try:
            q.etaCritFrac = etaCritFrac
        except:
            pass
    if not etaFoldFrac is None:
        try:
            q.etaFoldFrac = etaFoldFrac
        except:
            pass
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
        evolveReducingViscosityMultiplier = MorrisMonaghanReducingViscosity(nh,nh,aMin,aMax)
        packages.append(evolveReducingViscosityMultiplier)
    elif boolCullenViscosity:
        evolveCullenViscosityMultiplier = CullenDehnenViscosity(WT,alphMax,alphMin,betaC,betaD,betaE,fKern,boolHopkinsCorrection)
        packages.append(evolveCullenViscosityMultiplier)

#-------------------------------------------------------------------------------
# Construct the Artificial Conduction physics object.
#-------------------------------------------------------------------------------
if bArtificialConduction:
    #q.reducingViscosityCorrection = True
    ArtyCond = ArtificialConduction(WT,arCondAlpha)
    
    packages.append(ArtyCond)

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
xPlane0 = Plane(Vector(x0, y0), Vector( 1.0,  0.0))
yPlane0 = Plane(Vector(x0, y0), Vector( 0.0,  1.0))
xPlane1 = Plane(Vector(x2, y1), Vector(-1.0,  0.0))
yPlane1 = Plane(Vector(x2, y1), Vector( 0.0, -1.0))

xbc0 = ReflectingBoundary(xPlane0)
xbc1 = ReflectingBoundary(xPlane1)
ybc = PeriodicBoundary(yPlane0, yPlane1)

for p in packages:
    for bc in [xbc0, xbc1, ybc]:
        p.appendBoundary(bc)
del p, bc

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
                            vizBaseName = vizBaseName,
                            vizDir = vizDir,
                            vizStep = vizCycle,
                            vizTime = vizTime,
                            SPH = True)
output("control")

#-------------------------------------------------------------------------------
# If requested, reset the intial energies and densities.
#-------------------------------------------------------------------------------
if sumInitialDensity and control.totalSteps == 0:
    packages = integrator.physicsPackages()
    state = State(db, packages)
    derivs = StateDerivatives(db, packages)
    integrator.setGhostNodes()
    integrator.preStepInitialize(state, derivs)
    cm = db.connectivityMap()
    pos = db.fluidPosition
    mass = db.fluidMass
    H = db.fluidHfield
    rho = db.fluidMassDensity
    computeSPHSumMassDensity(cm, WT, True, pos, mass, H, rho)
    for nodes in nodeSet:
        pos = nodes.positions()
        eps = nodes.specificThermalEnergy()
        rho = nodes.massDensity()
        gammai = nodes.eos.gamma
        for i in range(nodes.numInternalNodes):
            eps[i] = specificEnergy(pos[i].x, rho[i],gammai)

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
if not steps is None:
    control.step(steps)

    # Are we doing the restart test?
    if checkRestart:
        state0 = State(db, integrator.physicsPackages())
        state0.copyState()
        control.loadRestartFile(control.totalSteps)
        state1 = State(db, integrator.physicsPackages())
        if not state1 == state0:
            raise ValueError("The restarted state does not match!")
        else:
            print("Restart check PASSED.")

else:
    control.advance(goalTime, maxSteps)
    control.dropRestartFile()

#-------------------------------------------------------------------------------
# Plot the final state.
#-------------------------------------------------------------------------------
dx1 = (x1 - x0)/nx1
dx2 = (x2 - x1)/nx2
h1 = 1.0/(nPerh*dx1)
h2 = 1.0/(nPerh*dx2)
answer = SodSolution(nPoints=nx1 + nx2,
                     gamma1 = gamma1,
                     gamma2 = gamma2,
                     rho1 = rho1,
                     P1 = P1,
                     rho2 = rho2,
                     P2 = P2,
                     x0 = x0,
                     x1 = x1,
                     x2 = x2,
                     h1 = 1.0/h1,
                     h2 = 1.0/h2)

cs = hydro.soundSpeed

# Provide a function to select the points we want to plot.
def plotFilter(pos):
    return (pos.y >= y0 + yplotbuf and pos.y <= y1 - yplotbuf)

# Make a flat list from a FieldList
pos = db.fluidPosition
def createList(x):
    result = []
    for i in range(len(x)):
        for j in range(x[i].numInternalElements):
            if plotFilter(pos(i,j)):
                result.append(x(i,j))
    return mpi.allreduce(result)

gamma = db.newFluidScalarFieldList(0.0, "specific heat ratio")
db.fluidGamma(gamma)
gammaList = createList(gamma)

# Compute the simulated specific entropy.
rho = createList(db.fluidMassDensity)
P = createList(hydro.pressure)
A = [Pi/rhoi**gammai for (Pi, rhoi, gammai) in zip(P, rho, gammaList)]

# The analytic solution for the simulated entropy.
xprof = [x.x for x in createList(db.fluidPosition)]
vprof = [v.x for v in createList(db.fluidVelocity)]
epsprof = createList(db.fluidSpecificThermalEnergy)
hprof = [1.0/Hi.xx for Hi in createList(db.fluidHfield)]
multiSort(xprof, rho, P, A, vprof, epsprof, hprof)
xans, vans, uans, rhoans, Pans, hans, gammaans = answer.solution(control.time(), xprof)
Aans = [Pi/rhoi**gammai for (Pi, rhoi, gammai) in zip(Pans,  rhoans, gammaans)]
csAns = [sqrt(gammai*Pi/rhoi) for (Pi, rhoi, gammai) in zip(Pans,  rhoans, gammaans)]

if graphics:
    import Gnuplot
    from SpheralGnuPlotUtilities import *

    rPlot = plotNodePositions2d(db, colorNodeLists=0, colorDomains=1)
    rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotState(db, plotStyle="points", filterFunc=plotFilter)
    plotAnswer(answer, control.time(),
               rhoPlot, velPlot, epsPlot, PPlot, HPlot)
    pE = plotEHistory(control.conserve)

    csPlot = plotFieldList(cs, winTitle="Sound speed", plotStyle="points", filterFunc=plotFilter)
    csAnsData = Gnuplot.Data(xans, csAns, 
                             with_ = "lines",
                             title = "Analytic")
    csPlot.replot(csAnsData)

    Aplot = generateNewGnuPlot()
    Adata = Gnuplot.Data(xprof, A,
                         with_ = "points",
                         title = "P/rho^\gamma",
                         inline = True)
    AansData = Gnuplot.Data(xprof, Aans,
                         with_ = "points",
                         title = "Solution",
                         inline = True)
    Aplot.replot(Adata)
    Aplot.replot(AansData)

    plots = [(rhoPlot, "Sod-planar-rho.png"),
             (velPlot, "Sod-planar-vel.png"),
             (epsPlot, "Sod-planar-eps.png"),
             (PPlot, "Sod-planar-P.png"),
             (HPlot, "Sod-planar-h.png"),
             (csPlot, "Sod-planar-cs.png"),
             (Aplot, "Sod-planar-entropy.png")]
    if not gsph:
        viscPlot = plotFieldList(hydro.maxViscousPressure,
                             winTitle = "max(rho^2 Piij)",
                             plotStyle = "points",
                             colorNodeLists = False,
                             filterFunc = plotFilter)
        plots.append((viscPlot, "Sod-planar-viscosity.png"))
    if crksph:
        volPlot = plotFieldList(control.RKCorrections.volume, 
                                winTitle = "volume",
                                colorNodeLists = False, plotGhosts = False)
        splot = plotFieldList(control.RKCorrections.surfacePoint,
                              winTitle = "surface point",
                              colorNodeLists = False)
        plots += [(volPlot, "Sod-planar-vol.png"),
                   (splot, "Sod-planar-surfacePoint.png")]

    if boolCullenViscosity:
        cullAlphaPlot = plotFieldList(q.ClMultiplier,
                                      plotStyle = "points",
                                      winTitle = "Cullen alpha",
                                      filterFunc = plotFilter)
        cullDalphaPlot = plotFieldList(evolveCullenViscosityMultiplier.DalphaDt(),
                                       plotStyle = "points",
                                       winTitle = "Cullen DalphaDt",
                                       filterFunc = plotFilter)
        plots += [(cullAlphaPlot, "Sod-planar-Cullen-alpha.png"),
                  (cullDalphaPlot, "Sod-planar-Cullen-DalphaDt.png")]

    if boolReduceViscosity:
        alphaPlot = plotFieldList(q.ClMultiplier,
                                  winTitle = "rvAlpha",
                                  plotStyle = "points",
                                  colorNodeLists = False,
                                  filterFunc = plotFilter)

    # Make hardcopies of the plots.
    for p, filename in plots:
        p.hardcopy(os.path.join(dataDir, filename), terminal="png")

print("Energy conservation: original=%g, final=%g, error=%g" % (control.conserve.EHistory[0],
                                                                control.conserve.EHistory[-1],
                                                                (control.conserve.EHistory[-1] - control.conserve.EHistory[0])/control.conserve.EHistory[0]))

#-------------------------------------------------------------------------------
# If requested, write out the state in a global ordering to a file.
#-------------------------------------------------------------------------------
rmin = x0
rmax = x2
if mpi.rank == 0:
    if outputFile:
        outputFile = os.path.join(dataDir, outputFile)
        f = open(outputFile, "w")
        f.write(("#" + 12*" '%s'" + "\n") % ("x", "rho", "P", "v", "eps", "A", "h", 
                                             "rhoans", "Pans", "vans", "Aans", "hans"))
        for (xi, rhoi, Pi, vi, epsi, Ai, hi, 
             rhoansi, Pansi, vansi, Aansi, hansi) in zip(xprof, rho, P, vprof, epsprof, A, hprof, 
                                                         rhoans, Pans, vans, Aans, hans):
            f.write((12*" %16.12e" + '\n') % 
                    (xi, rhoi, Pi, vi, epsi, Ai, hi, 
                     rhoansi, Pansi, vansi, Aansi, hansi))
        f.close()

    import Pnorm
    print("\tQuantity \t\tL1 \t\t\tL2 \t\t\tLinf")
    failure = False
    hD = []
    for (name, data, ans) in [("Mass Density", rho, rhoans),
                              ("Pressure", P, Pans),
                              ("Velocity", vprof, vans),
                              ("Thermal E", epsprof, uans),
                              ("Entropy", A, Aans),
                              ("h       ", hprof, hans)]:
        assert len(data) == len(ans)
        error = [data[i] - ans[i] for i in range(len(data))]
        Pn = Pnorm.Pnorm(error, xprof)
        L1 = Pn.gridpnorm(1, rmin, rmax)
        L2 = Pn.gridpnorm(2, rmin, rmax)
        Linf = Pn.gridpnorm("inf", rmin, rmax)
        print("\t%s \t\t%g \t\t%g \t\t%g" % (name, L1, L2, Linf))
        #f.write(("\t\t%g") % (L1))
        hD.append([L1,L2,Linf])
    #f.write("\n")

    print("%d\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t" % (nx1+nx2,hD[0][0],hD[1][0],hD[2][0],hD[3][0],
                                                                                hD[0][1],hD[1][1],hD[2][1],hD[3][1],
                                                                                hD[0][2],hD[1][2],hD[2][2],hD[3][2]))
