import os, sys
import shutil
from SolidSpheral2d import *
from SpheralTestUtilities import *
from SodAnalyticSolution import SodSolutionStiffGasStiffGas as SodSolution

from GenerateNodeDistribution2d import GenerateNodeDistribution2d
from GenerateNodeProfile import GeneratePlanarNodeProfile2d
if mpi.procs > 1:
    from VoronoiDistributeNodes import distributeNodes2d
else:
    from DistributeNodes import distributeNodes2d

title("2-D integrated hydro test -- planar Sod problem Water-Gas")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(initialRotation = 0.0,     # Degrees, optionally rotate and clip the initial positions

            # resolution
            nx1 = 320,      # number of nodes
            nx2 = 80,
            ny1 = 128,
            ny2 = 32,

            # initial conditions
            rho1 = 1000.0,  # density
            rho2 = 50.0,
            P1 = 1.0e9,     # Pressure
            P2 = 1.0e5,
            Po1 = 6.6e8,    # stiffening pressure
            Po2 = 0.0,
            gamma1 = 4.4,   # specific heat ratio
            gamma2 = 1.4,
            Cv1 = 4184.0,   # specific heat (irrevlant on this prob need for initialization)
            Cv2 = 50.0,   

            numNodeLists = 2,

            # domain
            x0 = -0.5,
            x1 = 0.2,
            x2 = 0.5,

            y0 = 0.0,
            y1 = 0.2,
            yplotbuf = 0.0,

            # Kernel Things
            nPerh = 3.01,
            KernelConstructor = WendlandC2Kernel,
            order = 5,
            HUpdate = IdealH,
            gradhCorrection = True,

            # which hydro
            crksph = False,
            psph = False,
            fsisph = False,

            # FSI parameters
            fsiRhoStabilizeCoeff = 0.1,          # diffusion operating through the vel-gradient
            fsiEpsDiffuseCoeff = 0.1,            # diffusion coeff for specific thermal energy
            fsiXSPHCoeff = 0.0,                  # ramps xsph up
            fsiInterfaceMethod = HLLCInterface,  # (HLLCInterface, ModulusInterface, NoInterface)

            # CRKSPH inputs
            correctionOrder = LinearOrder,
            volumeType = RKSumVolume,

            # hydro parameters
            densityUpdate = RigorousSumDensity,
            correctVelocityGradient = True,
            compatibleEnergy = True,
            evolveTotalEnergy = False,  # Only for SPH variants -- evolve total rather than specific energy
            solid = False,    # If true, use the fluid limit of the solid hydro option
            
            # artificial viscosity parameters
            boolReduceViscosity = False,
            nh = 5.0,
            aMin = 0.1,
            aMax = 2.0,
            Cl = None,
            Cq = None,
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
            linearInExpansion = False,
            quadraticInExpansion = False,
            Qlimiter = False,
            epsilon2 = 1e-4,
            hmin = 1e-10,
            hmax = 1.0,
            XSPH = False,
            epsilonTensile = 0.0,
            nTensile = 8,
            rhoMin = 0.01,
            filter = 0.00,
            
            # artificial conduction
            bArtificialConduction = False,
            arCondAlpha = 0.5,

            # integrator parameters
            IntegratorConstructor = CheapSynchronousRK2Integrator,
            cfl = 0.25,
            goalTime = 240e-6,
            dt = 1.0e-12,
            dtMin = 1.0e-12,
            dtMax = 0.1,
            dtGrowth = 2.0,
            dtverbose = False,
            rigorousBoundaries = False,
            steps = None,
            maxSteps = None,
            statsStep = 10,

            # IO settings
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
assert not (fsisph and not solid)

if crksph:
    hydroname = "CRKSPH"
elif psph:
    hydroname = "PSPH"
elif fsisph:
    hydroname = "FSISPH"
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
                       "filter=%f" % filter,
                       "%i" % (nx1 + nx2))
restartDir = os.path.join(dataDir, "restarts")
restartBaseName = os.path.join(restartDir, "Sod-planar-2d-WG-%i" % (nx1 + nx2))

vizDir = os.path.join(dataDir, "visit")
if vizTime is None and vizCycle is None:
    vizBaseName = None
else:
    vizBaseName = "Sod-planar-2d-WG-%i" % (nx1 + nx2)

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
units = MKS()
eos1 = StiffenedGas(gamma=gamma1, 
                    P0=Po1, 
                    Cv=Cv1,
                    constants=units)
eos2 = StiffenedGas(gamma=gamma2, 
                    P0=Po2, 
                    Cv=Cv2,
                    constants=units)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
if KernelConstructor==NBSplineKernel:
    assert order in (3,5,7)
    WT = TableKernel(KernelConstructor(order), 1000)
else:
    WT = TableKernel(KernelConstructor(), 1000) 
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

def rho_initial(xi):
    if xi < x1:
        return rho1
    else:
        return rho2

def specificEnergy(xi, rhoi):
    if xi < x1:
        Pi = P1
        Poi = Po1
        gammai = gamma1
    else:
        Pi = P2
        Poi = Po2
        gammai = gamma2
    return (Pi+gammai*Poi)/((gammai - 1.0)*rhoi)


#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------

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

    for i in range(nodes.numInternalNodes):
        eps[i] = specificEnergy(pos[i].x, rho[i])

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
if crksph:
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
    hydro = FSISPH(dataBase = db,
                   W = WT,
                   cfl = cfl,
                   sumDensityNodeLists=[],                       
                   densityStabilizationCoefficient = fsiRhoStabilizeCoeff,
                   specificThermalEnergyDiffusionCoefficient = fsiEpsDiffuseCoeff,
                   xsphCoefficient = fsiXSPHCoeff,
                   interfaceMethod = fsiInterfaceMethod,
                   compatibleEnergyEvolution = compatibleEnergy,
                   evolveTotalEnergy = evolveTotalEnergy,
                   linearCorrectGradients = correctVelocityGradient,
                   HUpdate = HUpdate)
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
                     Po1 = Po1,
                     rho2 = rho2,
                     P2 = P2,
                     Po2 = Po2,
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
# A = [Pi/rhoi**gammai for (Pi, rhoi, gammai) in zip(P, rho, gammaList)]

# The analytic solution for the simulated entropy.
xprof = [x.x for x in createList(db.fluidPosition)]
vprof = [v.x for v in createList(db.fluidVelocity)]
epsprof = createList(db.fluidSpecificThermalEnergy)
hprof = [1.0/Hi.xx for Hi in createList(db.fluidHfield)]
multiSort(xprof, rho, P, vprof, epsprof, hprof)
xans, vans, uans, rhoans, Pans, hans, gammaans = answer.solution(control.time(), xprof)
# Aans = [Pi/rhoi**gammai for (Pi, rhoi, gammai) in zip(Pans,  rhoans, gammaans)]
# csAns = [sqrt(gammai*Pi/rhoi) for (Pi, rhoi, gammai) in zip(Pans,  rhoans, gammaans)]

if graphics:
    import Gnuplot
    from SpheralGnuPlotUtilities import *

    #rPlot = plotNodePositions2d(db, colorNodeLists=0, colorDomains=1)
    rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotState(db, plotStyle="points", filterFunc=plotFilter)
    plotAnswer(answer, control.time(),
               rhoPlot, velPlot, epsPlot, PPlot, HPlot)
    pE = plotEHistory(control.conserve)

    #csPlot = plotFieldList(cs, winTitle="Sound speed", plotStyle="points", filterFunc=plotFilter)
    #csAnsData = Gnuplot.Data(xans, csAns, 
    #                         with_ = "lines",
    #                         title = "Analytic")
    #csPlot.replot(csAnsData)

    #Aplot = generateNewGnuPlot()
    #Adata = Gnuplot.Data(xprof, A,
    #                     with_ = "points",
    #                     title = "P/rho^\gamma",
    #                     inline = True)
    #AansData = Gnuplot.Data(xprof, Aans,
    #                     with_ = "points",
    #                     title = "Solution",
    #                     inline = True)
    #Aplot.replot(Adata)
    #Aplot.replot(AansData)

    plots = [(rhoPlot, "Sod-planar-rho.png"),
             (velPlot, "Sod-planar-vel.png"),
             (epsPlot, "Sod-planar-eps.png"),
             (PPlot, "Sod-planar-P.png"),
             (HPlot, "Sod-planar-h.png")]
    
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
        f.write(("#" + 10*" '%s'" + "\n") % ("x", "rho", "P", "v", "eps", "h", 
                                             "rhoans", "Pans", "vans",  "hans"))
        for (xi, rhoi, Pi, vi, epsi, hi, 
             rhoansi, Pansi, vansi, hansi) in zip(xprof, rho, P, vprof, epsprof, hprof, 
                                                         rhoans, Pans, vans, hans):
            f.write((10*" %16.12e" + '\n') % 
                    (xi, rhoi, Pi, vi, epsi, hi, 
                     rhoansi, Pansi, vansi, hansi))
        f.close()

    import Pnorm
    print("\tQuantity \t\tL1 \t\t\tL2 \t\t\tLinf")
    failure = False
    hD = []
    for (name, data, ans) in [("Mass Density", rho, rhoans),
                              ("Pressure", P, Pans),
                              ("Velocity", vprof, vans),
                              ("Thermal E", epsprof, uans),
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
