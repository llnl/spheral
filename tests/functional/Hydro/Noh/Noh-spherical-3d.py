#-------------------------------------------------------------------------------
# The Spherical Noh test case run in 3-D.
#
# W.F. Noh 1987, JCP, 72, 78-120.
#-------------------------------------------------------------------------------

import os, shutil, sys
from math import *
from SolidSpheral3d import *
from SpheralTestUtilities import *
from GenerateNodeDistribution3d import *

import mpi

if mpi.procs > 1:
    from PeanoHilbertDistributeNodes import distributeNodes3d
else:
    from DistributeNodes import distributeNodes3d

title("3-D integrated hydro test -- spherical Noh problem")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(order = 5,
            seed = "lattice",

            nx = 50,
            ny = 50,
            nz = 50,

            x0 = 0.0,
            y0 = 0.0,
            z0 = 0.0,
            x1 = 1.0,
            y1 = 1.0,
            z1 = 1.0,
            rmin = 0.0,
            rmax = 1.0,
            nPerh = 2.01,
            rho0 = 1.0,
            eps0 = 0.0,
            smallPressure = False,

            vr0 = -1.0, 

            gamma = 5.0/3.0,
            mu = 1.0,

            solid = False,                   # If true, use the fluid limit of the solid hydro option

            svph = False,
            crksph = False,
            psph = False,
            fsisph = False,
            gsph = False,
            mfm = False,
            mfv = False,

            asph = False,                    # This just chooses the H algorithm -- you can use this with CRKSPH for instance.
            radialOnly = False,              # Force ASPH tensors to be aligned and evolve radially
            boolReduceViscosity = False,
            HopkinsConductivity = False,     # For PSPH
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

            Cl = None, 
            Cq = None,
            linearInExpansion = None,
            Qlimiter = None,
            balsaraCorrection = None,
            epsilon2 = None,
            hmin = 0.0001, 
            hmax = 0.5,
            hminratio = 0.1,
            cfl = 0.5,
            XSPH = False,
            epsilonTensile = 0.0,
            nTensile = 8,
            xfilter = 0.0,

            IntegratorConstructor = CheapSynchronousRK2Integrator,
            goalTime = 0.6,
            steps = None,
            vizCycle = None,
            vizTime = 0.1,
            dt = 0.0001,
            dtMin = 1.0e-5, 
            dtMax = 0.1,
            dtGrowth = 2.0,
            maxSteps = None,
            statsStep = 10,
            smoothIters = 0,
            HUpdate = IdealH,
            correctionOrder = LinearOrder,
            volumeType = RKSumVolume,
            domainIndependent = False,
            rigorousBoundaries = False,
            dtverbose = False,

            densityUpdate = RigorousSumDensity, # VolumeScaledDensity,
            evolveTotalEnergy = False,  # Only for SPH variants -- evolve total rather than specific energy
            compatibleEnergy = True,
            gradhCorrection = True,
            correctVelocityGradient = True,

            useVoronoiOutput = True,
            clearDirectories = False,
            vizDerivs = False,
            restoreCycle = -1,
            restartStep = 1000,
            checkRestart = False,
            dataDir = "dumps-spherical-Noh",
            outputFile = "Noh_spherical_profiles.gnu",
            comparisonFile = None,
            doCompare = True,

            graphics = True,
            )

assert not(boolReduceViscosity and boolCullenViscosity)
assert not((gsph or mfm) and (boolReduceViscosity and boolCullenViscosity))
assert not(fsisph and not solid)

if smallPressure:
   P0 = 1.0e-6
   eps0 = P0/((gamma - 1.0)*rho0)

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
elif mfm:
    hydroname = "MFM"
elif mfv:
    hydroname = "MFV"
else:
    hydroname = "SPH"
if asph:
    hydroname = "A" + hydroname

if solid:
    hydroname = "Solid" + hydroname

if dataDir:
    dataDir = os.path.join(dataDir,
                           hydroname,
                           "nPerh=%f" % nPerh,
                           "compatibleEnergy=%s" % compatibleEnergy,
                           "Cullen=%s" % boolCullenViscosity,
                           "xfilter=%f" % xfilter,
                           "nx=%i_ny=%i_nz=%i" % (nx, ny, nz))
    restartDir = os.path.join(dataDir, "restarts")
    restartBaseName = os.path.join(restartDir, "Noh-spherical-3d-%ix%ix%i" % (nx, ny, nz))
    vizDir = os.path.join(dataDir, "visit")
else:
    restartBaseName = None
    vizDir = None
if vizTime is None and vizCycle is None:
    vizBaseName = None
else:
    vizBaseName = "Noh-spherical-3d-%ix%ix%i" % (nx, ny, nz)

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
if mpi.rank == 0 and dataDir:
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
WT = TableKernel(NBSplineKernel(order), 10)
output("WT")
kernelExtent = WT.kernelExtent

#-------------------------------------------------------------------------------
# Make the NodeList.
#-------------------------------------------------------------------------------
if solid:
    nodes1 = makeSolidNodeList("nodes1", eos,
                               hmin = hmin,
                               hmax = hmax,
                               kernelExtent = kernelExtent,
                               hminratio = hminratio,
                               nPerh = nPerh)
else:
    nodes1 = makeFluidNodeList("nodes1", eos,
                               hmin = hmin,
                               hmax = hmax,
                               kernelExtent = kernelExtent,
                               hminratio = hminratio,
                               nPerh = nPerh)
output("nodes1")
output("nodes1.hmin")
output("nodes1.hmax")
output("nodes1.hminratio")
output("nodes1.nodesPerSmoothingScale")

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
pos = nodes1.positions()
vel = nodes1.velocity()
generator = GenerateNodeDistribution3d(nx, ny, nz, rho0, seed,
                                       rmin = rmin,
                                       rmax = rmax,
                                       xmin = (x0, y0, z0),
                                       xmax = (x1, y1, z1),
                                       nNodePerh = nPerh,
                                       SPH = True)

distributeNodes3d((nodes1, generator))
output("mpi.reduce(nodes1.numInternalNodes, mpi.MIN)")
output("mpi.reduce(nodes1.numInternalNodes, mpi.MAX)")
output("mpi.reduce(nodes1.numInternalNodes, mpi.SUM)")

# Set node specific thermal energies
nodes1.specificThermalEnergy(ScalarField("tmp", nodes1, eps0))

# Set node velocities
for nodeID in range(nodes1.numNodes):
    vel[nodeID] = pos[nodeID].unitVector()*vr0

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
                 xmin = Vector(-1.1, -1.1, -1.1),
                 xmax = Vector( 1.1,  1.1, -1.1),
                 ASPH = asph)
elif crksph:
    hydro = CRKSPH(dataBase = db,
                   W = WT,
                   order = correctionOrder,
                   filter = xfilter,
                   cfl = cfl,
                   compatibleEnergyEvolution = compatibleEnergy,
                   XSPH = XSPH,
                   densityUpdate = densityUpdate,
                   HUpdate = HUpdate,
                   ASPH = asph)
elif fsisph:
    hydro = FSISPH(dataBase = db,
                   W = WT,
                   cfl = cfl,
                   interfaceMethod = HLLCInterface,
                   sumDensityNodeLists=[nodes1],                       
                   densityStabilizationCoefficient = 0.00,
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
                gradientType = RiemannGradient,
                densityUpdate=densityUpdate,
                HUpdate = IdealH,
                epsTensile = epsilonTensile,
                nTensile = nTensile)
elif mfm:
    limiter = VanLeerLimiter()
    waveSpeed = DavisWaveSpeed()
    solver = HLLC(limiter,waveSpeed,True)
    hydro = MFM(dataBase = db,
                riemannSolver = solver,
                W = WT,
                cfl=cfl,
                compatibleEnergyEvolution = compatibleEnergy,
                correctVelocityGradient=correctVelocityGradient,
                evolveTotalEnergy = evolveTotalEnergy,
                XSPH = XSPH,
                gradientType = RiemannGradient,
                densityUpdate=densityUpdate,
                HUpdate = IdealH,
                epsTensile = epsilonTensile,
                nTensile = nTensile)
elif mfv:
    limiter = VanLeerLimiter()
    waveSpeed = DavisWaveSpeed()
    solver = HLLC(limiter,waveSpeed,True)
    hydro = MFV(dataBase = db,
                riemannSolver = solver,
                W = WT,
                cfl=cfl,
                compatibleEnergyEvolution = compatibleEnergy,
                correctVelocityGradient=correctVelocityGradient,
                evolveTotalEnergy = evolveTotalEnergy,
                XSPH = XSPH,
                gradientType = RiemannGradient,
                densityUpdate=densityUpdate,
                HUpdate = IdealH,
                epsTensile = epsilonTensile,
                nTensile = nTensile)
elif psph:
    hydro = PSPH(dataBase = db,
                 W = WT,
                 filter = xfilter,
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
                filter = xfilter,
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
output("hydro.cfl")
output("hydro.compatibleEnergyEvolution")
try:
    output("hydro.PiKernel")
except:
    pass
if not fsisph:
    output("hydro.densityUpdate")
output("hydro._smoothingScaleMethod.HEvolution")
if radialOnly:
    assert asph
    hydro._smoothingScaleMethod.radialOnly = True
    output("hydro._smoothingScaleMethod.radialOnly")

packages = [hydro]

#-------------------------------------------------------------------------------
# Set the artificial viscosity parameters.
#-------------------------------------------------------------------------------
if not (gsph or mfm or mfv):
    q = hydro.Q
    if Cl:
        q.Cl = Cl
    if Cq:
        q.Cq = Cq
    if epsilon2:
        q.epsilon2 = epsilon2
    if Qlimiter:
        q.limiter = Qlimiter
    if balsaraCorrection:
        q.balsaraShearCorrection = balsaraCorrection
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
bcs = []
if x0 == 0.0:
    xPlane0 = Plane(Vector(0.0, 0.0, 0.0), Vector(1.0, 0.0, 0.0))
    bcs.append(ReflectingBoundary(xPlane0))
if y0 == 0.0:
    yPlane0 = Plane(Vector(0.0, 0.0, 0.0), Vector(0.0, 1.0, 0.0))
    bcs.append(ReflectingBoundary(yPlane0))
if z0 == 0.0:
    zPlane0 = Plane(Vector(0.0, 0.0, 0.0), Vector(0.0, 0.0, 1.0))
    bcs.append(ReflectingBoundary(zPlane0))

for p in packages:
    for bc in bcs:
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
if useVoronoiOutput:
    import SpheralVoronoiSiloDump
    vizMethod = SpheralVoronoiSiloDump.dumpPhysicsState
else:
    import SpheralPointmeshSiloDump
    vizMethod = SpheralPointmeshSiloDump.dumpPhysicsState
    #import SpheralVisitDump
    #vizMethod = SpheralVisitDump.dumpPhysicsState

control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle,
                            vizMethod = vizMethod,
                            vizBaseName = vizBaseName,
                            vizDir = vizDir,
                            vizStep = vizCycle,
                            vizTime = vizTime,
                            vizDerivs = vizDerivs,
                            #skipInitialPeriodicWork = SVPH,
                            SPH = True        # Only for iterating H
                            )
output("control")

# Do some startup stuff (unless we're restarting).
if restoreCycle is None:
    control.smoothState(smoothIters)
    if densityUpdate in (VoronoiCellDensity, SumVoronoiCellDensity):
        print("Reinitializing node masses.")
        control.voronoiInitializeMass()
    control.dropRestartFile()

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

if not doCompare:
    sys.exit(0)

#-------------------------------------------------------------------------------
# Plot the results.
#-------------------------------------------------------------------------------
import NohAnalyticSolution
answer = NohAnalyticSolution.NohSolution(3,
                                         h0 = nPerh*rmax/nx)

r = mpi.allreduce([x.magnitude() for x in nodes1.positions().internalValues()], mpi.SUM)
rho = mpi.allreduce(list(nodes1.massDensity().internalValues()), mpi.SUM)
rans, vans, epsans, rhoans, Pans, hans = answer.solution(control.time(), r)
if mpi.rank == 0:
        L1 = 0.0
        for i in range(len(rho)):
          L1 = L1 + abs(rho[i]-rhoans[i])
        L1_tot = L1 / len(rho)
        print("L1=",L1_tot,"\n")
        with open("Converge.txt", "a") as myfile:
          myfile.write("%s %s %s %s %s\n" % (nx, ny,nz,nx+ny+nz, L1_tot))
if graphics:

    # Plot the node positions.
    from SpheralGnuPlotUtilities import *

    # Plot the final state.
    rhoPlot, vrPlot, epsPlot, PPlot, HPlot = plotRadialState(db)
    del HPlot
    Hinverse = db.newFluidSymTensorFieldList()
    db.fluidHinverse(Hinverse)
    hr = db.newFluidScalarFieldList()
    ht = db.newFluidScalarFieldList()
    for Hfield, hrfield, htfield in zip(Hinverse,
                                        hr,
                                        ht):
        n = Hfield.numElements
        assert hrfield.numElements == n
        assert htfield.numElements == n
        positions = Hfield.nodeList().positions()
        for i in range(n):
            runit = positions[i].unitVector()
            tunit = Vector(-(positions[i].y), positions[i].x).unitVector()
            hrfield[i] = (Hfield[i]*runit).magnitude()
            htfield[i] = (Hfield[i]*tunit).magnitude()
    hrPlot = plotFieldList(hr, xFunction="%s.magnitude()", plotStyle="points", winTitle="h_r")
    htPlot = plotFieldList(ht, xFunction="%s.magnitude()", plotStyle="points", winTitle="h_t")

    # Overplot the analytic solution.
    plotAnswer(answer, control.time(),
               rhoPlot = rhoPlot,
               velPlot = vrPlot,
               epsPlot = epsPlot,
               PPlot = PPlot,
               HPlot = hrPlot)

    if mpi.rank == 0:
        r, hrans, htans = answer.hrtsolution(control.time())
        htData = Gnuplot.Data(r, htans,
                              title = "Solution",
                              with_ = "lines",
                              inline = "true")
        htPlot.replot(htData)
    plots = [(rhoPlot, "Noh-spherical-rho.png"),
             (vrPlot, "Noh-spherical-vel.png"),
             (epsPlot, "Noh-spherical-eps.png"),
             (PPlot, "Noh-spherical-P.png"),
             (hrPlot, "Noh-spherical-hr.png"),
             (htPlot, "Noh-spherical-ht.png")]

    # Make hardcopies of the plots.
    for p, filename in plots:
        p.hardcopy(os.path.join(dataDir, filename), terminal="png")

    # Report the error norms.
    rmin, rmax = 0.05, 0.35
    r = mpi.allreduce([x.magnitude() for x in nodes1.positions().internalValues()], mpi.SUM)
    rho = mpi.allreduce(list(nodes1.massDensity().internalValues()), mpi.SUM)
    v = mpi.allreduce([x.magnitude() for x in nodes1.velocity().internalValues()], mpi.SUM)
    eps = mpi.allreduce(list(nodes1.specificThermalEnergy().internalValues()), mpi.SUM)
    Pf = ScalarField("pressure", nodes1)
    nodes1.pressure(Pf)
    P = mpi.allreduce(list(Pf.internalValues()), mpi.SUM)
    if mpi.rank == 0:
        from SpheralTestUtilities import multiSort
        import Pnorm
        multiSort(r, rho, v, eps, P)
        rans, vans, epsans, rhoans, Pans, hans = answer.solution(control.time(), r)
        with open("Converge.txt", "a") as myfile:
          myfile.write("%s %s %s %s " % (nx, ny,nz,nx+ny+nz))
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
            myfile.write("\t\t%g \t\t%g \t\t%g" % (L1, L2, Linf))
          myfile.write("\n")

#-------------------------------------------------------------------------------
# If requested, write out the state in a global ordering to a file.
#-------------------------------------------------------------------------------
rmaxnorm = 0.35
rminnorm = 0.05

if outputFile:
    outputFile = os.path.join(dataDir, outputFile)
    from SpheralTestUtilities import multiSort
    P = ScalarField("pressure", nodes1)
    nodes1.pressure(P)
    xprof = mpi.reduce([x.x for x in nodes1.positions().internalValues()], mpi.SUM)
    yprof = mpi.reduce([x.y for x in nodes1.positions().internalValues()], mpi.SUM)
    zprof = mpi.reduce([x.z for x in nodes1.positions().internalValues()], mpi.SUM)
    rprof = mpi.reduce([ri.magnitude() for ri in nodes1.positions().internalValues()],mpi.SUM)
    rhoprof = mpi.reduce(nodes1.massDensity().internalValues(), mpi.SUM)
    Pprof = mpi.reduce(P.internalValues(), mpi.SUM)
    vprof = mpi.reduce([vi.dot(ri.unitVector()) for ri,vi in zip(nodes1.positions().internalValues(),nodes1.velocity().internalValues())],mpi.SUM)
    #vprof = mpi.reduce([v.magnitude() for v in nodes1.velocity().internalValues()], mpi.SUM)
    epsprof = mpi.reduce(nodes1.specificThermalEnergy().internalValues(), mpi.SUM)
    hprof = mpi.reduce([1.0/sqrt(H.Determinant()) for H in nodes1.Hfield().internalValues()], mpi.SUM)
    mof = mortonOrderIndices(db)
    mo = mpi.reduce(mof[0].internalValues(), mpi.SUM)
    if mpi.rank == 0:
        from Pnorm import Pnorm
        multiSort(rprof, mo, xprof, yprof, zprof, rhoprof, Pprof, vprof, epsprof, hprof)
        rans, vans, epsans, rhoans, Pans, hans = answer.solution(control.time(), rprof)
        velans = vans
        L1rho = Pnorm(rhoprof, rprof, rhoans).pnorm(1, rmin=rminnorm, rmax=rmaxnorm) 
        L2rho = Pnorm(rhoprof, rprof, rhoans).pnorm(2, rmin=rminnorm, rmax=rmaxnorm)
        Linfrho = Pnorm(rhoprof, rprof, rhoans).pnorm("inf", rmin=rminnorm, rmax=rmaxnorm)
        L1eps = Pnorm(epsprof, rprof, epsans).pnorm(1, rmin=rminnorm, rmax=rmaxnorm)
        L2eps = Pnorm(epsprof, rprof, epsans).pnorm(2, rmin=rminnorm, rmax=rmaxnorm)
        Linfeps = Pnorm(epsprof, rprof, epsans).pnorm("inf", rmin=rminnorm, rmax=rmaxnorm)
        L1vel = Pnorm(vprof, rprof, velans).pnorm(1, rmin=rminnorm, rmax=rmaxnorm)
        L2vel = Pnorm(vprof, rprof, velans).pnorm(2, rmin=rminnorm, rmax=rmaxnorm)
        Linfvel = Pnorm(vprof, rprof, velans).pnorm("inf", rmin=rminnorm, rmax=rmaxnorm)
        L1P = Pnorm(Pprof, rprof, Pans).pnorm(1, rmin=rminnorm, rmax=rmaxnorm)
        L2P = Pnorm(Pprof, rprof, Pans).pnorm(2, rmin=rminnorm, rmax=rmaxnorm)
        LinfP = Pnorm(Pprof, rprof, velans).pnorm("inf", rmin=rminnorm, rmax=rmaxnorm)
        with open("convergeNoh3d-CRK-%s-cullen-%s-PSPH-%s.txt" % (crksph,boolCullenViscosity,psph), "a") as myfile:
          myfile.write(("#" + 14*"%16s\t " + "%16s\n") % ("N", "L1rho", "L1eps", "L1vel", "L2rho", "L2eps", "L2vel", "Linfrho", "Linfeps", "Linfvel", "L1P", "L2P", "LinfP", "cycles", "runtime"))
          myfile.write((14*"%16s\t " + "%16s\n") % (nx, L1rho, L1eps, L1vel, L2rho, L2eps, L2vel, Linfrho, Linfeps, Linfvel, L1P, L2P, LinfP, control.totalSteps, control.stepTimer.elapsedTime))
        

        f = open(outputFile, "w")
        f.write(("# " + 22*"%15s " + "\n") % ("r", "x", "y", "z", "rho", "P", "v", "eps", "h", "mortonOrder",
                                              "rhoans", "Pans", "vans", "epsans",
                                              "x_uu", "y_uu", "z_uu", "rho_uu", "P_uu", "v_uu", "eps_uu", "h_uu"))
        for (ri, xi, yi, zi, rhoi, Pi, vi, epsi, hi, mi, 
             rhoansi, Pansi, vansi, epsansi)  in zip(rprof, xprof, yprof, zprof, rhoprof, Pprof, vprof, epsprof, hprof, mo,
                                                     rhoans, Pans, vans, epsans):
            f.write((9*"%16.12e " + "%i " + 4*"%16.12e " + 8*"%i " + "\n") % (ri, xi, yi, zi, rhoi, Pi, vi, epsi, hi, mi,
                                                                              rhoansi, Pansi, vansi, epsansi,
                                                                              unpackElementUL(packElementDouble(xi)),
                                                                              unpackElementUL(packElementDouble(yi)),
                                                                              unpackElementUL(packElementDouble(zi)),
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
Eerror = (control.conserve.EHistory[-1] - control.conserve.EHistory[0])/control.conserve.EHistory[0]
print("Total energy error: %g" % Eerror)
if compatibleEnergy and abs(Eerror) > 1e-10:
    raise ValueError("Energy error outside allowed bounds.")

