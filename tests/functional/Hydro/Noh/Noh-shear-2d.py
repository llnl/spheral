#-------------------------------------------------------------------------------
# This is our shearing variant of the planar Noh problem.
#-------------------------------------------------------------------------------
from math import *
import shutil
import mpi
from Spheral2d import *
from SpheralTestUtilities import *
from SpheralMatplotlib import *
from findLastRestart import *
from SpheralVisitDump import dumpPhysicsState

from GenerateNodeDistribution2d import *
from CompositeNodeDistribution import *
from CentroidalVoronoiRelaxation import *

title("2-D integrated hydro test -- Shearing Planar Noh problem")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(seed = "lattice",

            nx = 20,
            ny = 100,
            nPerh = 4.01,
            KernelConstructor = WendlandC4Kernel,

            x0 = 0.0,
            x1 = 0.2,
            y0 = 0.0,
            y1 = 1.0,

            rho1 = 1.0,
            eps1 = 0.0,
            P1 = 0.0,
            vshear = 1.0,
            vy = -1.0,

            gamma = 5.0/3.0,
            mu = 1.0,

            hydroType = "SPH",                 # one of (SPH, SVPH, CRKSPH, PSPH, FSISPH, GSPH, MFM)
            ASPH = False,   # This just chooses the H algorithm -- you can use this with CRKSPH for instance.
            Qconstructor = MonaghanGingoldViscosity,
            #Qconstructor = TensorMonaghanGingoldViscosity,
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
            Cl = 1.0, 
            Cq = 0.75,
            Qlimiter = False,
            linearInExpansion = False,
            balsaraCorrection = False,
            epsilon2 = 1e-2,
            fslice = 0.5,
            hmin = 0.0001, 
            hmax = 0.5,
            hminratio = 0.1,
            cfl = 0.5,
            useVelocityMagnitudeForDt = False,
            XSPH = False,
            epsilonTensile = 0.0,
            nTensile = 8,
            filter = 0.0,
            hourglass = None,
            hourglassOrder = 0,
            hourglassLimiter = 0,
            hourglassFraction = 0.5,
            HopkinsConductivity = False,     # For PSPH
            evolveTotalEnergy = False,       # Only for SPH variants -- evolve total rather than specific energy

            IntegratorConstructor = CheapSynchronousRK2Integrator,
            goalTime = 0.6,
            steps = None,
            vizCycle = None,
            vizTime = 0.05,
            dt = 0.0001,
            dtMin = 1.0e-5, 
            dtMax = 0.1,
            dtGrowth = 2.0,
            maxSteps = None,
            statsStep = 10,
            smoothIters = 0,
            HUpdate = IdealH,
            HEvolution = IdealH,
            domainIndependent = False,
            rigorousBoundaries = False,
            dtverbose = False,

            correctionOrder = LinearOrder,
            densityUpdate = RigorousSumDensity, # VolumeScaledDensity,
            compatibleEnergy = True,
            gradhCorrection = True,
            correctVelocityGradient = True,

            clearDirectories = False,
            vizDerivs = False,
            vizGhosts = False,
            checkRestart = False,
            redistributeStep = 500,
            restoreCycle = None,
            restartStep = 1000,
            dataRoot = "dumps-shearingNoh-2d",

            graphics = True,
            outputFile = None,
            comparisonFile = None,

            )
assert not(boolReduceViscosity and boolCullenViscosity)

hydroType = hydroType.upper()

if dataRoot:
    dataDir = os.path.join(dataRoot,
                           hydroType,
                           Qconstructor.__name__,
                           "basaraShearCorrection=%s_Qlimiter=%s" % (balsaraCorrection, Qlimiter),
                           "nperh=%4.2f" % nPerh,
                           "XSPH=%s" % XSPH,
                           "densityUpdate=%s" % densityUpdate,
                           "compatibleEnergy=%s" % compatibleEnergy,
                           "Cullen=%s" % boolCullenViscosity,
                           "gradhCorrection=%s" % gradhCorrection,
                           "nx=%i_ny=%i" % (nx, ny))
    restartDir = os.path.join(dataDir, "restarts")
    vizDir = os.path.join(dataDir, "visit")
    restartBaseName = os.path.join(restartDir, "Noh-shear-2d-%ix%i" % (nx, ny))
else:
    restartDir = None
    vizDir = None
    restartBaseName = None
if vizTime is None and vizCycle is None:
    vizBaseName = None
else:
    vizBaseName = "Noh-shear-2d-%ix%i" % (nx, ny)

xmin = (x0, y0)
xmax = (x1, y1)

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
# If we're restarting, find the set of most recent restart files.
#-------------------------------------------------------------------------------
if restoreCycle is None:
    restoreCycle = findLastRestart(restartBaseName)

#-------------------------------------------------------------------------------
# Material properties.
#-------------------------------------------------------------------------------
eos = GammaLawGasMKS(gamma, mu)

#-------------------------------------------------------------------------------
# Create our interpolation kernels -- one for normal hydro interactions, and
# one for use with the artificial viscosity
#-------------------------------------------------------------------------------
if KernelConstructor==NBSplineKernel:
  WT = TableKernel(NBSplineKernel(order))
else:
  WT = TableKernel(KernelConstructor())
output("WT")
kernelExtent = WT.kernelExtent

#-------------------------------------------------------------------------------
# Create a NodeList and associated Neighbor object.
#-------------------------------------------------------------------------------
nodes1 = makeFluidNodeList("nodes1", eos, 
                           hmin = hmin,
                           hmax = hmax,
                           kernelExtent = kernelExtent,
                           nPerh = nPerh)
output("nodes1")
output("nodes1.hmin")
output("nodes1.hmax")
output("nodes1.hminratio")
output("nodes1.nodesPerSmoothingScale")
#-------------------------------------------------------------------------------
# Set node properties.
#-------------------------------------------------------------------------------
if restoreCycle is None:
    from DistributeNodes import distributeNodes2d
    print("Generating node distribution.")
    generator1 = GenerateNodeDistribution2d(nx, ny, rho1, seed,
                                            xmin = xmin,
                                            xmax = xmax,
                                            nNodePerh = nPerh,
                                            SPH = True)
    n1 = generator1.globalNumNodes()

    if mpi.procs > 1:
        #from VoronoiDistributeNodes import distributeNodes2d
        from PeanoHilbertDistributeNodes import distributeNodes2d
    else:
        from DistributeNodes import distributeNodes2d
    distributeNodes2d((nodes1, generator1))
    output("mpi.reduce(nodes1.numInternalNodes, mpi.MIN)")
    output("mpi.reduce(nodes1.numInternalNodes, mpi.MAX)")
    output("mpi.reduce(nodes1.numInternalNodes, mpi.SUM)")

    # Set node specific thermal energies
    nodes1.specificThermalEnergy(ScalarField("tmp", nodes1, P1/((gamma - 1.0)*rho1)))

    # Set node velocities
    pos = nodes1.positions()
    vel = nodes1.velocity()
    for i in range(nodes1.numInternalNodes):
        vel[i] = Vector(vshear*cos(2.0*pi*pos[i].y), vy)

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
if Qconstructor is TensorSVPHViscosity:
    q = Qconstructor(Cl, Cq, fslice)
else:
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

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
if hydroType == "SVPH":
    hydro = SVPH(dataBase = db,
                 W = WT, 
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
                 xmin = Vector(-1.1, -1.1),
                 xmax = Vector( 1.1,  1.1),
                 ASPH = ASPH)
elif hydroType == "CRKSPH":
    hydro = CRKSPH(dataBase = db,
                   W = WT,
                   Q = q,
                   filter = filter,
                   cfl = cfl,
                   compatibleEnergyEvolution = compatibleEnergy,
                   XSPH = XSPH,
                   correctionOrder = correctionOrder,
                   densityUpdate = densityUpdate,
                   HUpdate = HUpdate,
                   ASPH = ASPH)
elif hydroType == "PSPH":
   hydro = PSPH(dataBase = db,
                W = WT,
                Q = q,
                filter = filter,
                cfl = cfl,
                compatibleEnergyEvolution = compatibleEnergy,
                evolveTotalEnergy = evolveTotalEnergy,
                HopkinsConductivity = HopkinsConductivity,
                correctVelocityGradient = correctVelocityGradient,
                densityUpdate = densityUpdate,
                HUpdate = HUpdate,
                XSPH = XSPH,
                ASPH = ASPH)
else:
    hydro = SPH(dataBase = db,
                W = WT, 
                Q = q,
                filter = filter,
                cfl = cfl,
                useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                compatibleEnergyEvolution = compatibleEnergy,
                evolveTotalEnergy = evolveTotalEnergy,
                gradhCorrection = gradhCorrection,
                correctVelocityGradient = correctVelocityGradient,
                densityUpdate = densityUpdate,
                HUpdate = HUpdate,
                XSPH = XSPH,
                epsTensile = epsilonTensile,
                nTensile = nTensile,
                ASPH = ASPH)

output("hydro")
output("hydro.kernel")
output("hydro.PiKernel")
output("hydro.cfl")
output("hydro.compatibleEnergyEvolution")
output("hydro.densityUpdate")
output("hydro.HEvolution")

packages = [hydro]

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
xPlane0 = Plane(Vector(*xmin), Vector( 1.0, 0.0))
xPlane1 = Plane(Vector(*xmax), Vector(-1.0, 0.0))
yPlane0 = Plane(Vector(*xmin), Vector( 0.0, 1.0))
yPlane1 = Plane(Vector(*xmax), Vector( 0.0,-1.0))
xbc = PeriodicBoundary(xPlane0, xPlane1)
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
integrator.cullGhostNodes = False
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
if vizGhosts:
   from SpheralPointmeshSiloDump import dumpPhysicsState
else:
   dumpPhysicsState = None
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle,
                            redistributeStep = redistributeStep,
                            vizMethod = dumpPhysicsState,
                            vizBaseName = vizBaseName,
                            vizDir = vizDir,
                            vizStep = vizCycle,
                            vizTime = vizTime,
                            vizDerivs = vizDerivs,
                            vizGhosts = vizGhosts,
                            skipInitialPeriodicWork = (hydroType == "SVPH"))
output("control")

# Smooth the initial conditions.
if restoreCycle is not None:
    control.loadRestartFile(restoreCycle)
else:
    control.smoothState(smoothIters)
    if densityUpdate in (VoronoiCellDensity, SumVoronoiCellDensity):
        print("Reinitializing node masses.")
        control.voronoiInitializeMass()
    control.dropRestartFile()
    control.dropViz()

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
    control.updateViz(control.totalSteps, integrator.currentTime, 0.0)
    control.dropRestartFile()
#-------------------------------------------------------------------------------
# Plot the results.
#-------------------------------------------------------------------------------
if graphics:

    # Plot the node positions.
    rPlot = plotNodePositions2d(db, colorNodeLists=0, colorDomains=1)

    # Plot the final state.
    rhoPlot, vrPlot, epsPlot, PPlot, HPlot = plotState(db,
                                                       xFunction = "%s.y",
                                                       vecyFunction = "%s.y",
                                                       tenyFunction = "%s.yy ** -1")

    # Overplot the analytic solution.
    import NohAnalyticSolution
    answer = NohAnalyticSolution.NohSolution(1,
                                             h0 = nPerh*y1/ny)
    plotAnswer(answer, control.time(),
               rhoPlot = rhoPlot,
               velPlot = vrPlot,
               epsPlot = epsPlot,
               PPlot = PPlot,
               HPlot = HPlot)

    plots = [(rhoPlot, "Noh-shear-rho.png"),
             (vrPlot, "Noh-shear-vel.png"),
             (epsPlot, "Noh-shear-eps.png"),
             (PPlot, "Noh-shear-P.png"),
             (HPlot, "Noh-shear-h.png")]

    # Make hardcopies of the plots.
    for p, filename in plots:
        p.figure.savefig(os.path.join(dataDir, filename))

    # Report the error norms.
    rmin, rmax = 0.05, 0.35
    r = mpi.allreduce([x.y for x in nodes1.positions().internalValues()], mpi.SUM)
    rho = mpi.allreduce(list(nodes1.massDensity().internalValues()), mpi.SUM)
    v = mpi.allreduce([x.y for x in nodes1.velocity().internalValues()], mpi.SUM)
    eps = mpi.allreduce(list(nodes1.specificThermalEnergy().internalValues()), mpi.SUM)
    Pf = ScalarField("pressure", nodes1)
    nodes1.pressure(Pf)
    P = mpi.allreduce(list(Pf.internalValues()), mpi.SUM)
    if mpi.rank == 0:
        from SpheralTestUtilities import multiSort
        import Pnorm
        multiSort(r, rho, v, eps, P)
        rans, vans, epsans, rhoans, Pans, hans = answer.solution(control.time(), r)
        print("\tQuantity \t\tL1 \t\t\tL2 \t\t\tLinf")
        for (name, data, ans) in [("Mass Density", rho, rhoans),
                                  ("Pressure", P, Pans),
                                  ("Velocity", v, vans),
                                  ("Thermal E", eps, epsans)]:
            assert len(data) == len(ans)
            error = [data[i] - ans[i] for i in range(len(data))]
            Pn = Pnorm.Pnorm(error, r)
            L1 = Pn.gridpnorm(1, rmin, rmax)
            L2 = Pn.gridpnorm(2, rmin, rmax)
            Linf = Pn.gridpnorm("inf", rmin, rmax)
            print("\t%s \t\t%g \t\t%g \t\t%g" % (name, L1, L2, Linf))

#-------------------------------------------------------------------------------
# If requested, write out the state in a global ordering to a file.
#-------------------------------------------------------------------------------
if outputFile:
    outputFile = os.path.join(dataDir, outputFile)
    from SpheralTestUtilities import multiSort
    P = ScalarField("pressure", nodes1)
    nodes1.pressure(P)
    xprof = mpi.reduce([x.x for x in nodes1.positions().internalValues()], mpi.SUM)
    yprof = mpi.reduce([x.y for x in nodes1.positions().internalValues()], mpi.SUM)
    rprof = mpi.reduce([x.y for x in nodes1.positions().internalValues()], mpi.SUM)
    rhoprof = mpi.reduce(nodes1.massDensity().internalValues(), mpi.SUM)
    Pprof = mpi.reduce(P.internalValues(), mpi.SUM)
    vprof = mpi.reduce([v.x for v in nodes1.velocity().internalValues()], mpi.SUM)
    epsprof = mpi.reduce(nodes1.specificThermalEnergy().internalValues(), mpi.SUM)
    hprof = mpi.reduce([1.0/sqrt(H.Determinant()) for H in nodes1.Hfield().internalValues()], mpi.SUM)
    mof = mortonOrderIndices(db)
    mo = mpi.reduce(mof[0].internalValues(), mpi.SUM)
    if mpi.rank == 0:
        import NohAnalyticSolution
        answer = NohAnalyticSolution.NohSolution(1,
                                             h0 = nPerh*y1/ny)
        multiSort(rprof, mo, xprof, yprof, rhoprof, Pprof, vprof, epsprof, hprof)
        rans, vans, epsans, rhoans, Pans, hans = answer.solution(control.time(), rprof)
        f = open(outputFile, "w")
        f.write(("# " + 21*"%15s " + "\n") % ("r", "x", "y", "rho", "P", "v", "eps", "h", "mortonOrder", "QWork",
                                              "rhoans", "Pans", "vans", "epsans",
                                              "x_uu", "y_uu", "rho_uu", "P_uu", "v_uu", "eps_uu", "h_uu"))
        for (ri, xi, yi, rhoi, Pi, vi, epsi, hi, mi,
             rhoansi, Pansi, vansi, epsansi)  in zip(rprof, xprof, yprof, rhoprof, Pprof, vprof, epsprof, hprof, mo,
                                                     rhoans, Pans, vans, epsans):
            f.write((8*"%16.12e " + "%i " + 4*"%16.12e " + 7*"%i " + "\n") % (ri, xi, yi, rhoi, Pi, vi, epsi, hi, mi,
                                                                              rhoansi, Pansi, vansi, epsansi,
                                                                              unpackElementUL(packElementDouble(xi)),
                                                                              unpackElementUL(packElementDouble(yi)),
                                                                              unpackElementUL(packElementDouble(rhoi)),
                                                                              unpackElementUL(packElementDouble(Pi)),
                                                                              unpackElementUL(packElementDouble(vi)),
                                                                              unpackElementUL(packElementDouble(epsi)),
                                                                              unpackElementUL(packElementDouble(hi))))
        f.close()

        #---------------------------------------------------------------------------
        # Also we can optionally compare the current results with another file.
        #---------------------------------------------------------------------------
        if comparisonFile:
            comparisonFile = os.path.join(dataDir, comparisonFile)
            import filecmp
            assert filecmp.cmp(outputFile, comparisonFile)
Eerror = (control.conserve.EHistory[-1] - control.conserve.EHistory[0])/max(1.0e-30, control.conserve.EHistory[0])
print("Total energy error: %g" % Eerror)
if compatibleEnergy and abs(Eerror) > 1e-13:
    raise ValueError("Energy error outside allowed bounds.")
