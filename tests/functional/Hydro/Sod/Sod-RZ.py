#-------------------------------------------------------------------------------
# The Sod test case run in RZ symmetry.
#-------------------------------------------------------------------------------
import os, shutil, mpi
from SolidSpheralRZ import *
from SpheralTestUtilities import *

from GenerateNodeDistribution2d import *
if mpi.procs > 1:
    from VoronoiDistributeNodes import distributeNodes2d
else:
    from DistributeNodes import distributeNodes2d

title("RZ hydro test -- Sod problem")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(problem = "planar",     # one of (planar, cylindrical, spherical)
            KernelConstructor = NBSplineKernel,
            order = 5,

            rho1 = 1.0,
            rho2 = 0.25,
            eps1 = 1.0,
            eps2 = 1.0,

            n1 = 100,
            n2 = 50,

            nPerh = 1.35,

            gammaGas = 1.4,
            mu = 1.0,

            solid = False,    # If true, use the fluid limit of the solid hydro option

            CRKSPH = False,
            PSPH = False,
            SPH = True,       # Choose the H advancement
            evolveTotalEnergy = False,  # Only for SPH variants -- evolve total rather than specific energy
            Qconstructor = MonaghanGingoldViscosity,
            boolReduceViscosity = False,
            HopkinsConductivity = False,     # For PSPH
            nhQ = 5.0,
            nhL = 10.0,
            aMin = 0.1,
            aMax = 2.0,
            boolCullenViscosity = False,
            cullenUseHydroDerivatives = True,  # Reuse the hydro calculation of DvDx.
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
            Qhmult = 1.0,
            Cl = 1.0, 
            Cq = 1.0,
            etaCritFrac = 1.0,
            etaFoldFrac = 0.2,
            Qlimiter = False,
            balsaraCorrection = False,
            epsilon2 = 1e-2,
            hmin = 0.0001, 
            hmax = 0.1,
            hminratio = 0.1,
            cfl = 0.5,
            useVelocityMagnitudeForDt = False,
            XSPH = False,
            epsilonTensile = 0.0,
            nTensile = 4.0,
            hourglass = None,
            hourglassOrder = 0,
            hourglassLimiter = 0,
            hourglassFraction = 0.5,
            filter = 0.0,

            IntegratorConstructor = CheapSynchronousRK2Integrator,
            goalTime = 0.1,
            steps = None,
            dt = 1.0e-4,
            dtMin = 1.0e-6, 
            dtMax = 0.1,
            dtGrowth = 2.0,
            dtverbose = False,
            rigorousBoundaries = False,
            maxSteps = None,
            statsStep = 1,
            vizCycle = None,
            vizTime = 0.1,
            vizDerivs = False,
            HUpdate = IdealH,
            correctionOrder = LinearOrder,
            QcorrectionOrder = LinearOrder,
            volumeType = RKSumVolume,
            densityUpdate = RigorousSumDensity, # VolumeScaledDensity,
            compatibleEnergy = True,
            gradhCorrection = False,
            correctVelocityGradient = True,
            domainIndependent = False,
            cullGhostNodes = True,
            
            bArtificialConduction = False,
            arCondAlpha = 0.5,

            clearDirectories = True,
            checkError = False,
            checkRestart = False,
            checkEnergy = False,
            restoreCycle = -1,
            restartStep = 100,
            outputFile = None,
            comparisonFile = None,
            normOutputFile = None,
            writeOutputLabel = True,

            graphics = True,
            )

assert not(boolReduceViscosity and boolCullenViscosity)
assert problem in ("planar", "cylindrical", "spherical")

if CRKSPH:
   if solid:
      if SPH:
         HydroConstructor = SolidCRKSPHHydro
      else:
         HydroConstructor = SolidACRKSPHHydro
   else:
      if SPH:
         HydroConstructor = CRKSPHHydro
      else:
         HydroConstructor = ACRKSPHHydro
      Qconstructor = LimitedMonaghanGingoldViscosity
      gradhCorrection = False
else:
   if solid:
      if SPH:
         HydroConstructor = SolidSPHHydro
      else:
         HydroConstructor = SolidASPHHydro
   else:
      if SPH:
         HydroConstructor = SPHHydro
      else:
         HydroConstructor = ASPHHydro

dataDir = os.path.join("dumps-%s-Sod-RZ" % problem,
                       HydroConstructor.__name__,
                       Qconstructor.__name__,
                       "nPerh=%f" % nPerh,
                       "compatibleEnergy=%s" % compatibleEnergy,
                       "Cullen=%s" % boolCullenViscosity,
                       "filter=%f" % filter)
restartDir = os.path.join(dataDir, "restarts")
restartBaseName = os.path.join(restartDir, "Sod-%s-RZ" % problem)

vizDir = os.path.join(dataDir, "visit")
if vizTime is None and vizCycle is None:
    vizBaseName = None
else:
    vizBaseName = "Sod-%s-RZ" % problem

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
eos = GammaLawGasMKS(gammaGas, mu)
strength = NullStrength()

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
if KernelConstructor==NBSplineKernel:
    Wbase = NBSplineKernel(order)
else:
    Wbase = KernelConstructor()
WT = TableKernel(Wbase, 1000)
kernelExtent = WT.kernelExtent
output("WT")

#-------------------------------------------------------------------------------
# Make the NodeList.
#-------------------------------------------------------------------------------
if solid:
    nodes1 = makeSolidNodeList("nodes1", eos, strength,
                               hmin = hmin,
                               hmax = hmax,
                               hminratio = hminratio,
                               nPerh = nPerh,
                               xmin = (-10.0, -10.0),
                               xmax = ( 10.0,  10.0),
                               kernelExtent = kernelExtent)
    nodes2 = makeSolidNodeList("nodes2", eos, strength,
                               hmin = hmin,
                               hmax = hmax,
                               hminratio = hminratio,
                               nPerh = nPerh,
                               xmin = (-10.0, -10.0),
                               xmax = ( 10.0,  10.0),
                               kernelExtent = kernelExtent)
else:
    nodes1 = makeFluidNodeList("nodes1", eos, 
                               hmin = hmin,
                               hmax = hmax,
                               hminratio = hminratio,
                               nPerh = nPerh,
                               xmin = (-10.0, -10.0),
                               xmax = ( 10.0,  10.0),
                               kernelExtent = kernelExtent)
    nodes2 = makeFluidNodeList("nodes2", eos, 
                               hmin = hmin,
                               hmax = hmax,
                               hminratio = hminratio,
                               nPerh = nPerh,
                               xmin = (-10.0, -10.0),
                               xmax = ( 10.0,  10.0),
                               kernelExtent = kernelExtent)
    
nodeSet = (nodes1, nodes2)
for n in nodeSet:
    output("n")
    output("n.name")
    output("n.hmin")
    output("n.hmax")
    output("n.nodesPerSmoothingScale")
del n

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
if problem == "planar":
    r0, r1 = 4.0, 4.2
    z0, z1, z2 = 0.0, 0.5, 1.0
    nmatch2 = int(float(n2)*(r1 - r0)/(z2 - z1) + 0.5)
    if SPH:
        nmatch1 = int(float(n1)*(r1 - r0)/(z1 - z0) + 0.5)
    else:
        nmatch1 = nmatch2
    gen1 = RZGenerator(GenerateNodeDistribution2d(n1, nmatch1, rho1, "lattice",
                                                  xmin = (z0, r0),
                                                  xmax = (z1, r1),
                                                  nNodePerh = nPerh,
                                                  SPH = SPH))
    gen2 = RZGenerator(GenerateNodeDistribution2d(n2, nmatch2, rho2, "lattice",
                                                  xmin = (z1, r0),
                                                  xmax = (z2, r1),
                                                  nNodePerh = nPerh,
                                                  SPH = SPH))
elif problem == "cylindrical":
    r0, r1, r2 = 4.0, 5.0, 6.0,
    z0, z1 = 0.0, 0.2
    nmatch1 = int(float(n1)*(z1 - z0)/(r1 - r0) + 0.5)
    nmatch2 = int(float(n2)*(z1 - z0)/(r2 - r1) + 0.5)
    gen1 = RZGenerator(GenerateNodeDistribution2d(nmatch1, n1, rho1, "lattice",
                                                  xmin = (z0, r0),
                                                  xmax = (z1, r1),
                                                  nNodePerh = nPerh,
                                                  SPH = SPH))
    gen2 = RZGenerator(GenerateNodeDistribution2d(nmatch2, n2, rho2, "lattice",
                                                  xmin = (z0, r1),
                                                  xmax = (z1, r2),
                                                  nNodePerh = nPerh,
                                                  SPH = SPH))
else:
    assert problem == "spherical"
    rmax1 = 1.0
    rmax2 = 2.0
    dr2 = (rmax2 - rmax1)/n2
    nrind = 10
    gen1 = RZGenerator(GenerateNodeDistribution2d(n1, n1, rho1, "constantDTheta",
                                                  xmin = (0.0, 0.0),
                                                  xmax = (rmax1, rmax1),
                                                  rmin = 0.0,
                                                  rmax = rmax1,
                                                  nNodePerh = nPerh,
                                                  SPH = SPH))
    gen2 = RZGenerator(GenerateNodeDistribution2d(n2 + nrind, n2, rho2, "constantDTheta",
                                                  xmin = (0.0, 0.0),
                                                  xmax = (rmax2, rmax2),
                                                  rmin = rmax1,
                                                  rmax = rmax2 + nrind*dr2,
                                                  nNodePerh = nPerh,
                                                  SPH = SPH))

distributeNodes2d((nodes1, gen1),
                  (nodes2, gen2))
for n in nodeSet:
    output("mpi.reduce(n.numInternalNodes, mpi.MIN)")
    output("mpi.reduce(n.numInternalNodes, mpi.MAX)")
    output("mpi.reduce(n.numInternalNodes, mpi.SUM)")
del n

# Set node specific thermal energies
for n, eps in ((nodes1, eps1),
               (nodes2, eps2)):
    n.specificThermalEnergy(ScalarField("tmp", n, eps))
del n

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase()
output("db")
for n in nodeSet:
    db.appendNodeList(n)
del n
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Construct the artificial viscosity.
#-------------------------------------------------------------------------------
q = Qconstructor(Cl, Cq)
q.epsilon2 = epsilon2
q.limiter = Qlimiter
q.balsaraShearCorrection = balsaraCorrection
q.QcorrectionOrder = QcorrectionOrder
output("q")
output("q.Cl")
output("q.Cq")
output("q.epsilon2")
output("q.limiter")
output("q.balsaraShearCorrection")

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
if CRKSPH:
    hydro = HydroConstructor(W = WT,
                             Q = q,
                             filter = filter,
                             cfl = cfl,
                             useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                             compatibleEnergyEvolution = compatibleEnergy,
                             evolveTotalEnergy = evolveTotalEnergy,
                             XSPH = XSPH,
                             correctionOrder = correctionOrder,
                             volumeType = volumeType,
                             densityUpdate = densityUpdate,
                             HUpdate = HUpdate)
    q.etaCritFrac = etaCritFrac
    q.etaFoldFrac = etaFoldFrac
else:
    hydro = HydroConstructor(W = WT,
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
# Construct the MMRV physics object.
#-------------------------------------------------------------------------------
if boolReduceViscosity:
    evolveReducingViscosityMultiplier = MorrisMonaghanReducingViscosity(nhQ,nhL,aMin,aMax)
    packages.append(evolveReducingViscosityMultiplier)
elif boolCullenViscosity:
    evolveCullenViscosityMultiplier = CullenDehnenViscosity(WT,alphMax,alphMin,betaC,betaD,betaE,fKern,boolHopkinsCorrection,cullenUseHydroDerivatives)
    packages.append(evolveCullenViscosityMultiplier)

#-------------------------------------------------------------------------------
# Construct the Artificial Conduction physics object.
#-------------------------------------------------------------------------------
if bArtificialConduction:
    ArtyCond = ArtificialConduction(WT,arCondAlpha)
    packages.append(ArtyCond)

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
if problem == "planar":
    bcs = [ReflectingBoundary(Plane(Vector(z0, r0), Vector( 1.0,  0.0))),
           ReflectingBoundary(Plane(Vector(z2, r0), Vector(-1.0,  0.0))),
           ReflectingBoundary(Plane(Vector(z2, r1), Vector( 0.0, -1.0)))]
    if r0 > 0.0:
        bcs.append(ReflectingBoundary(Plane(Vector(z0, r0), Vector(0.0, 1.0))))
elif problem == "cylindrical":
    bcs = [ReflectingBoundary(Plane(Vector(z0, r0), Vector( 1.0,  0.0))),
           ReflectingBoundary(Plane(Vector(z1, r0), Vector(-1.0,  0.0))),
           ReflectingBoundary(Plane(Vector(z0, r2), Vector( 0.0, -1.0)))]
    if r0 > 0.0:
        bcs.append(ReflectingBoundary(Plane(Vector(z0, r0), Vector(0.0, 1.0))))
else:
    assert problem == "spherical"
    boundNodes = vector_of_int()
    pos = nodes2.positions()
    for i in range(nodes2.numInternalNodes):
        if pos[i].magnitude() > rmax2:
            boundNodes.append(i)
    print("Selected %i boundary nodes" % mpi.allreduce(len(boundNodes), mpi.SUM))
    denialPlane = Plane(Vector(-2.0*rmax2, 0.0), Vector(1.0, 0.0))  # A fake denial plane since we're working in circles.
    bcs = [ReflectingBoundary(Plane(Vector(0.0, 0.0), Vector( 1.0,  0.0))),
           ConstantBoundary(nodes2, boundNodes, denialPlane)]

for bc in bcs:
    for p in packages:
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
integrator.domainDecompositionIndependent = domainIndependent
integrator.verbose = dtverbose
output("integrator")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.dtGrowth")
output("integrator.rigorousBoundaries")
output("integrator.domainDecompositionIndependent")
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
                            vizDerivs = vizDerivs,
                            SPH = SPH)
output("control")

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
if not steps is None:
    control.step(steps)
else:
   control.advance(goalTime, maxSteps)
   control.dropRestartFile()

#-------------------------------------------------------------------------------
# Plot the final state.
#-------------------------------------------------------------------------------
if graphics:
    from SpheralGnuPlotUtilities import *
    if problem == "planar":
        rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotState(db, xFunction="%s.x", vecyFunction="%s.x", tenyFunction="1.0/%s.xx")
    elif problem == "cylindrical":
        rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotState(db, xFunction="%s.y", vecyFunction="%s.y", tenyFunction="1.0/%s.yy")
    else:
        rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotRadialState(db)
    EPlot = plotEHistory(control.conserve)
    plots = [(rhoPlot, "Sod-%s-rho-RZ.png" % problem),
             (velPlot, "Sod-%s-vel-RZ.png" % problem),
             (epsPlot, "Sod-%s-eps-RZ.png" % problem),
             (PPlot, "Sod-%s-P-RZ.png" % problem),
             (HPlot, "Sod-%s-h-RZ.png" % problem)]

    # If this is the planar problem we can compare with the solution.
    if problem == "planar":
        from SodAnalyticSolution import *
        dz1 = (z1 - z0)/n1
        dz2 = (z2 - z1)/n2
        h1 = 1.0/(nPerh*dz1)
        h2 = 1.0/(nPerh*dz2)
        answer = SodSolution(nPoints = n1 + n2,
                             gamma = gammaGas,
                             rho1 = rho1,
                             P1 = (gammaGas - 1.0)*rho1*eps1,
                             rho2 = rho2,
                             P2 = (gammaGas - 1.0)*rho2*eps2,
                             x0 = z0,
                             x1 = z1,
                             x2 = z2,
                             h1 = 1.0/h1,
                             h2 = 1.0/h2)
        plotAnswer(answer, control.time(),
                   rhoPlot, velPlot, epsPlot, PPlot, HPlot)

Eerror = (control.conserve.EHistory[-1] - control.conserve.EHistory[0])/control.conserve.EHistory[0]
print("Total energy error: %g" % Eerror)
if checkEnergy and abs(Eerror) > 1e-13:
    raise ValueError("Energy error outside allowed bounds.")
