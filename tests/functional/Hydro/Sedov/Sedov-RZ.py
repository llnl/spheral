#-------------------------------------------------------------------------------
# The Sedov test case run in RZ symmetry.
#-------------------------------------------------------------------------------
import os, shutil, mpi
from SpheralRZ import *
from SpheralTestUtilities import *

from GenerateNodeDistribution2d import *
if mpi.procs > 1:
    from VoronoiDistributeNodes import distributeNodes2d
else:
    from DistributeNodes import distributeNodes2d

title("RZ hydro test -- Sedov problem")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(problem = "planar",     # one of (planar, cylindrical, spherical)
            KernelConstructor = NBSplineKernel,
            order = 5,

            n1 = 100,
            n2 = 20,

            nPerh = 1.35,

            gamma = 5.0/3.0,
            mu = 1.0,

            rho0 = 1.0,
            eps0 = 0.0,
            Espike = 1.0,
            smoothSpikeScale = 0.25,   # How much to smooth the spike in eta space

            solid = False,             # If true, use the fluid limit of the solid hydro option

            crksph = False,
            psph = False,
            asph = False,              # Choose the H advancement
            compatibleEnergy = False,
            evolveTotalEnergy = True,  # Only for SPH variants -- evolve total rather than specific energy
            Cl = 1.0, 
            Cq = 1.0,
            etaCritFrac = 1.0,
            etaFoldFrac = 0.2,
            Qlimiter = False,
            balsaraCorrection = False,
            epsilon2 = 1e-2,
            hmin = 1e-10,
            hmax = 1e10,
            hminratio = 0.02,
            cfl = 0.5,
            useVelocityMagnitudeForDt = False,
            xsph = False,
            etaMinAxis = 0.01,        # r at which we start to modify hydro significantly

            IntegratorConstructor = VerletIntegrator,
            goalTime = None,
            goalRadius = 0.8,
            steps = None,
            dt = 1e-8,
            dtMin = 1.0e-8,
            dtMax = 0.1,
            dtGrowth = 1.1,
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
            gradhCorrection = False,
            correctVelocityGradient = True,
            domainIndependent = False,
            cullGhostNodes = True,
            
            bArtificialConduction = False,
            arCondAlpha = 0.5,

            clearDirectories = False,
            checkError = False,
            checkRestart = False,
            checkEnergy = False,
            restoreCycle = -1,
            restartStep = 10000,
            comparisonFile = None,
            normOutputFile = None,
            writeOutputLabel = True,

            graphics = True,
            )

outputFile = "Sedov-%s-RZ.gnu" % problem

assert problem in ("planar", "cylindrical", "spherical")
assert not (compatibleEnergy and evolveTotalEnergy)

hydroname = ""
if solid:
    hydroname += "Solid"
if asph:
    hydroname += "A"
if crksph:
    hydroname += "CRKSPH"
else:
    hydroname += "SPH"

dataDir = os.path.join("dumps-%s-Sedov-RZ" % problem,
                       hydroname,
                       "nPerh=%f" % nPerh)
if compatibleEnergy:
    dataDir = os.path.join(dataDir, "compatibleEnergy")
elif evolveTotalEnergy:
    dataDir = os.path.join(dataDir, "evolveTotalEnergy")
else:
    dataDir = os.path.join(dataDir, "nonconservative")
restartDir = os.path.join(dataDir, "restarts")
restartBaseName = os.path.join(restartDir, "Sedov-%s-RZ" % problem)

vizDir = os.path.join(dataDir, "visit")
if vizTime is None and vizCycle is None:
    vizBaseName = None
else:
    vizBaseName = "Sedov-%s-RZ" % problem

# Figure out what our goal time should be.
import SedovAnalyticSolution
h0 = 1.0/n1*nPerh
if problem == "planar":
    ndim = 1
elif problem == "cylindrical":
    ndim = 2
else:
    ndim = 3
answer = SedovAnalyticSolution.SedovSolution(nDim = ndim,
                                             gamma = gamma,
                                             rho0 = rho0,
                                             E0 = Espike,
                                             h0 = h0)
if goalTime is None:
    assert not goalRadius is None
    nu1 = 1.0/(answer.nu + 2.0)
    nu2 = 2.0*nu1
    goalTime = (goalRadius*(answer.alpha*rho0/Espike)**nu1)**(1.0/nu2)
vs, r2, v2, rho2, P2 = answer.shockState(goalTime)
print("Predicted shock position %g at goal time %g." % (r2, goalTime))

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
eos = GammaLawGasMKS(gamma, mu, minimumPressure=0.0)
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
                               xmin = Vector(-10.0, -10.0),
                               xmax = Vector( 10.0,  10.0),
                               kernelExtent = kernelExtent)
else:
    nodes1 = makeFluidNodeList("nodes1", eos, 
                               hmin = hmin,
                               hmax = hmax,
                               hminratio = hminratio,
                               nPerh = nPerh,
                               xmin = Vector(-10.0, -10.0),
                               xmax = Vector( 10.0,  10.0),
                               kernelExtent = kernelExtent)
    
output("nodes1")
output("nodes1.hmin")
output("nodes1.hmax")
output("nodes1.nodesPerSmoothingScale")

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
if problem == "planar":
    nz = n1
    nr = n2
    z0, z1 = 0.0, 1.0
    r0, r1 = 0.0, 0.2
    rmin, rmax = None, None
elif problem == "cylindrical":
    nz = n2
    nr = n1
    z0, z1 = 0.0, 0.2
    r0, r1 = 0.0, 1.0
    rmin, rmax = None, None
else:
    assert problem == "spherical"
    nz = n1
    nr = n1
    rmin, rmax = 0.0, 1.0
    z0, z1 = 0.0, 1.0
    r0, r1 = 0.0, 1.0

generator = RZGenerator(GenerateNodeDistribution2d(nz, nr, rho0, "lattice",
                                                   xmin = (z0, r0),
                                                   xmax = (z1, r1),
                                                   rmin = rmin,
                                                   rmax = rmax,
                                                   nNodePerh = nPerh,
                                                   SPH = not asph))

distributeNodes2d((nodes1, generator))
output("mpi.reduce(nodes1.numInternalNodes, mpi.MIN)")
output("mpi.reduce(nodes1.numInternalNodes, mpi.MAX)")
output("mpi.reduce(nodes1.numInternalNodes, mpi.SUM)")

#-------------------------------------------------------------------------------
# Set the point source of energy.
#-------------------------------------------------------------------------------
pos = nodes1.positions()
vel = nodes1.velocity()
mass = nodes1.mass()
eps = nodes1.specificThermalEnergy()
H = nodes1.Hfield()
Esum = 0.0
dr = (r1 - r0)/nr
dz = (z1 - z0)/nz
msum = 0.0
if problem == "planar":
    epsi = 0.5*Espike/(rho0*dz)
    for i in range(nodes1.numInternalNodes):
        if pos[i].x < z0 + dz:
            eps[i] += epsi
            Esum += mass[i]*epsi
elif problem == "cylindrical":
    epsi = Espike/(rho0*pi*dr*dr)
    for i in range(nodes1.numInternalNodes):
        if pos[i].y < r0 + dr:
            eps[i] += epsi
            Esum += mass[i]*epsi
else:
    Wsum = 0.0
    for i in range(nodes1.numInternalNodes):
        Hi = H[i]
        etaij = (Hi*pos[i]).magnitude()
        Wi = WT.kernelValue(etaij/smoothSpikeScale, 1.0) * pos[i].y
        Ei = Wi*0.25*Espike
        eps[i] = Ei
        Wsum += Wi
    Wsum = mpi.allreduce(Wsum, mpi.SUM)
    assert Wsum > 0.0
    for i in range(nodes1.numInternalNodes):
        eps[i] = eps[i]/(Wsum*mass[i])
        Esum += eps[i]*mass[i]
        eps[i] += eps0
Eglobal = mpi.allreduce(Esum, mpi.SUM)
if problem == "planar":
    Eexpect = 0.5*Espike*pi*(r1*r1 - r0*r0)
elif problem == "cylindrical":
    Eexpect = Espike*(z1 - z0)
else:
    Eexpect = 0.25*Espike
print("Initialized a total energy of", Eglobal, Eexpect, Eglobal/Eexpect)
assert fuzzyEqual(Eglobal, Eexpect)

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
if crksph:
    hydro = CRKSPHRZ(dataBase = db,
                     cfl = cfl,
                     useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                     compatibleEnergyEvolution = compatibleEnergy,
                     evolveTotalEnergy = evolveTotalEnergy,
                     XSPH = xsph,
                     correctionOrder = correctionOrder,
                     densityUpdate = densityUpdate,
                     HUpdate = HUpdate,
                     etaMinAxis = etaMinAxis)
    hydro.Q.etaCritFrac = etaCritFrac
    hydro.Q.etaFoldFrac = etaFoldFrac
else:
    hydro = SPHRZ(dataBase = db,
                  W = WT,
                  cfl = cfl,
                  useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                  compatibleEnergyEvolution = compatibleEnergy,
                  evolveTotalEnergy = evolveTotalEnergy,
                  gradhCorrection = gradhCorrection,
                  correctVelocityGradient = correctVelocityGradient,
                  densityUpdate = densityUpdate,
                  HUpdate = HUpdate,
                  XSPH = xsph,
                  etaMinAxis = etaMinAxis)
output("hydro")
output("hydro.cfl")
output("hydro.compatibleEnergyEvolution")
output("hydro.densityUpdate")
output("hydro.HEvolution")

packages = [hydro]

#-------------------------------------------------------------------------------
# Construct the artificial viscosity.
#-------------------------------------------------------------------------------
q = hydro.Q
q.Cl = Cl
q.Cq = Cq
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
# Construct the Artificial Conduction physics object.
#-------------------------------------------------------------------------------
if bArtificialConduction:
    #q.reducingViscosityCorrection = True
    ArtyCond = ArtificialConduction(WT,arCondAlpha)
    
    packages.append(ArtyCond)

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
if problem == "planar":
    bcs = [ReflectingBoundary(Plane(Vector(z0, r0), Vector( 1.0,  0.0))),
           ReflectingBoundary(Plane(Vector(z0, r1), Vector( 0.0, -1.0)))]
    if r0 != 0.0:
        bcs.append(ReflectingBoundary(Plane(Vector(z0, r0), Vector( 0.0, 1.0))))
elif problem == "cylindrical":
    bcs = [ReflectingBoundary(Plane(Vector(z0, r0), Vector( 1.0,  0.0))),
           ReflectingBoundary(Plane(Vector(z1, r0), Vector(-1.0,  0.0)))]
else:
    assert problem == "spherical"
    bcs = [ReflectingBoundary(Plane(Vector(z0, r0), Vector( 1.0,  0.0)))]

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
integrator.allowDtCheck = True
output("integrator")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.dtGrowth")
output("integrator.rigorousBoundaries")
output("integrator.domainDecompositionIndependent")
output("integrator.verbose")
output("integrator.allowDtCheck")

#-------------------------------------------------------------------------------
# Make the problem controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            volumeType = volumeType,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle,
                            vizBaseName = vizBaseName,
                            vizDir = vizDir,
                            vizStep = vizCycle,
                            vizTime = vizTime,
                            vizDerivs = vizDerivs,
                            SPH = not asph)
output("control")

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
if not steps is None:
    control.step(steps)

else:
   control.advance(goalTime, maxSteps)

#-------------------------------------------------------------------------------
# Compute the analytic answer.
#-------------------------------------------------------------------------------
import mpi
import SedovAnalyticSolution
if problem == "planar":
    xprof = mpi.allreduce([x.x for x in nodes1.positions().internalValues()], mpi.SUM)
elif problem == "cylindrical":
    xprof = mpi.allreduce([x.y for x in nodes1.positions().internalValues()], mpi.SUM)
else:
    xprof = mpi.allreduce([x.magnitude() for x in nodes1.positions().internalValues()], mpi.SUM)

# Compute the simulated specific entropy.
rho = mpi.allreduce(nodes1.massDensity().internalValues(), mpi.SUM)
Pf = ScalarField("pressure", nodes1)
nodes1.pressure(Pf)
P = mpi.allreduce(Pf.internalValues(), mpi.SUM)
A = [Pi/rhoi**gamma for (Pi, rhoi) in zip(P, rho)]

# Solution profiles.
xans, vans, uans, rhoans, Pans, Aans, hans = answer.solution(control.time(), xprof)
L1 = 0.0
for i in range(len(rho)):
    L1 = L1 + abs(rho[i]-rhoans[i])
L1_tot = L1 / len(rho)
# if mpi.rank == 0 and outputFile:
#     print "L1=",L1_tot,"\n"
#     with open("Converge.txt", "a") as myfile:
#         myfile.write("%s %s\n" % (nz, L1_tot))

#-------------------------------------------------------------------------------
# Plot the final state.
#-------------------------------------------------------------------------------
if graphics:
    from SpheralMatplotlib import *
    if problem == "planar":
        rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotState(db, xFunction="%s.x", vecyFunction="%s.x", tenyFunction="1.0/%s.xx")
    elif problem == "cylindrical":
        rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotState(db, xFunction="%s.y", vecyFunction="%s.y", tenyFunction="1.0/%s.yy")
    else:
        rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotRadialState(db)
    APlot = newFigure()
    APlot.plot(xprof, A, marker='o', label="Simulation")
    plotAnswer(answer, control.time(), rhoPlot, velPlot, epsPlot, PPlot, APlot, HPlot)
    EPlot = plotEHistory(control.conserve)
    plots = [(rhoPlot, "Sedov-%s-rho-RZ.png" % problem),
             (velPlot, "Sedov-%s-vel-RZ.png" % problem),
             (epsPlot, "Sedov-%s-eps-RZ.png" % problem),
             (PPlot, "Sedov-%s-P-RZ.png" % problem),
             (APlot, "Sedov-planar-A.png"),
             (HPlot, "Sedov-%s-h-RZ.png" % problem)]

    if crksph:
        volPlot = plotFieldList(hydro.volume, 
                                winTitle = "volume",
                                colorNodeLists = False, plotGhosts = False)
        plots.append((volPlot, "Sedov-%s-vol.png" % problem))

    # Make hardcopies of the plots.
    for p, filename in plots:
        p.figure.savefig(os.path.join(dataDir, filename))

#-------------------------------------------------------------------------------
# Measure the difference between the simulation and analytic answer.
#-------------------------------------------------------------------------------
rmin, rmax = 0.05, 0.35   # Throw away anything with r < rwall to avoid wall heating.
rhoprof = mpi.reduce(nodes1.massDensity().internalValues(), mpi.SUM)
P = ScalarField("pressure", nodes1)
nodes1.pressure(P)
Pprof = mpi.reduce(P.internalValues(), mpi.SUM)
vprof = mpi.reduce([v.x for v in nodes1.velocity().internalValues()], mpi.SUM)
epsprof = mpi.reduce(nodes1.specificThermalEnergy().internalValues(), mpi.SUM)
hprof = mpi.reduce([1.0/H.xx for H in nodes1.Hfield().internalValues()], mpi.SUM)
xprof = mpi.reduce([x.magnitude() for x in nodes1.positions().internalValues()], mpi.SUM)

#-------------------------------------------------------------------------------
# If requested, write out the state in a global ordering to a file.
#-------------------------------------------------------------------------------
if outputFile:
    outputFile = os.path.join(dataDir, outputFile)
    from SpheralTestUtilities import multiSort
    mof = mortonOrderIndices(db)
    mo = mpi.reduce(mof[0].internalValues(), mpi.SUM)
    mprof = mpi.reduce(nodes1.mass().internalValues(), mpi.SUM)
    rhoprof = mpi.reduce(nodes1.massDensity().internalValues(), mpi.SUM)
    P = ScalarField("pressure", nodes1)
    nodes1.pressure(P)
    Pprof = mpi.reduce(P.internalValues(), mpi.SUM)
    vprof = mpi.reduce([v.x for v in nodes1.velocity().internalValues()], mpi.SUM)
    epsprof = mpi.reduce(nodes1.specificThermalEnergy().internalValues(), mpi.SUM)
    hprof = mpi.reduce([1.0/H.xx for H in nodes1.Hfield().internalValues()], mpi.SUM)
    if mpi.rank == 0:
        multiSort(xprof, rhoprof, Pprof, vprof, epsprof, hprof, mo,
                  rhoans, Pans, vans, uans, hans)
        f = open(outputFile, "w")
        f.write(("#  " + 20*"'%s' " + "\n") % ("x", "m", "rho", "P", "v", "eps", "h", "mo",
                                               "rhoans", "Pans", "vans", "epsans", "hans",
                                               "x_UU", "m_UU", "rho_UU", "P_UU", "v_UU", "eps_UU", "h_UU"))
        for (xi, mi, rhoi, Pi, vi, epsi, hi, moi,
             rhoansi, Pansi, vansi, uansi, hansi) in zip(xprof, mprof, rhoprof, Pprof, vprof, epsprof, hprof, mo,
                                                         rhoans, Pans, vans, uans, hans):
            f.write((7*"%16.12e " + "%i " + 5*"%16.12e " + 7*"%i " + '\n') % 
                    (xi, mi, rhoi, Pi, vi, epsi, hi, moi,
                     rhoansi, Pansi, vansi, uansi, hansi,
                     unpackElementUL(packElementDouble(xi)),
                     unpackElementUL(packElementDouble(mi)),
                     unpackElementUL(packElementDouble(rhoi)),
                     unpackElementUL(packElementDouble(Pi)),
                     unpackElementUL(packElementDouble(vi)),
                     unpackElementUL(packElementDouble(epsi)),
                     unpackElementUL(packElementDouble(hi))))
        f.close()

        # #---------------------------------------------------------------------------
        # # Also we can optionally compare the current results with another file.
        # #---------------------------------------------------------------------------
        # if comparisonFile:
        #     comparisonFile = os.path.join(dataDir, comparisonFile)
        #     import filecmp
        #     assert filecmp.cmp(outputFile, comparisonFile)

Eerror = (control.conserve.EHistory[-1] - control.conserve.EHistory[0])/control.conserve.EHistory[0]
print("Total energy error: %g" % Eerror)
if checkEnergy and abs(Eerror) > 1e-13:
    raise ValueError("Energy error outside allowed bounds.")
