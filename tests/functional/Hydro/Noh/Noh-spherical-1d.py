#-------------------------------------------------------------------------------
# The Spherical Noh test case run in 1D spherical coordinates
#
# W.F. Noh 1987, JCP, 72, 78-120.
#-------------------------------------------------------------------------------
#
# Ordinary SPH
#
#ATS:t0 = test(      SELF, "--graphics None --clearDirectories True  --checkError True   --restartStep 20", label="Spherical Noh problem -- 1-D (serial)")
#ATS:t1 = testif(t0, SELF, "--graphics None --clearDirectories False --checkError False  --restartStep 20 --restoreCycle 20 --steps 20 --checkRestart True", label="Spherical Noh problem -- 1-D (serial) RESTART CHECK")
#ATS:t2 = test(      SELF, "--graphics None --clearDirectories True  --checkError True  --dataDirBase 'dumps-spherical-restartcheck' --restartStep 20", np=2, label="Spherical Noh problem -- 1-D (parallel)")
#ATS:t3 = testif(t2, SELF, "--graphics None --clearDirectories False --checkError False --dataDirBase 'dumps-spherical-restartcheck' --restartStep 20 --restoreCycle 20 --steps 20 --checkRestart True", np=2, label="Spherical Noh problem -- 1-D (parallel) RESTART CHECK")
#ATS:t4 = test(      SELF, "--graphics None --clearDirectories True  --checkError True  --dataDirBase 'dumps-spherical-reproducing' --domainIndependent True --outputFile 'Noh-spherical-1proc-reproducing.txt'", label="Spherical Noh problem -- 1-D (serial reproducing test setup)")
#ATS:t5 = testif(t4, SELF, "--graphics None --clearDirectories False  --checkError True  --dataDirBase 'dumps-spherical-reproducing' --domainIndependent True --outputFile 'Noh-spherical-4proc-reproducing.txt' --comparisonFile 'Noh-spherical-1proc-reproducing.txt'", np=4, label="Spherical Noh problem -- 1-D (4 proc reproducing test)")
#
# Ordinary solid SPH
#
#ATS:t100 = test(        SELF, "--solid True --graphics None --clearDirectories True  --checkError True   --restartStep 20", label="Spherical Noh problem with solid SPH -- 1-D (serial)")
#ATS:t101 = testif(t100, SELF, "--solid True --graphics None --clearDirectories False --checkError False  --restartStep 20 --restoreCycle 20 --steps 20 --checkRestart True", label="Spherical Noh problem with solid SPH -- 1-D (serial) RESTART CHECK")
#ATS:t102 = test(        SELF, "--solid True --graphics None --clearDirectories True  --checkError True  --dataDirBase 'dumps-spherical-restartcheck' --restartStep 20", np=2, label="Spherical Noh problem with solid SPH -- 1-D (parallel)")
#ATS:t103 = testif(t102, SELF, "--solid True --graphics None --clearDirectories False --checkError False --dataDirBase 'dumps-spherical-restartcheck' --restartStep 20 --restoreCycle 20 --steps 20 --checkRestart True", np=2, label="Spherical Noh problem with solid SPH -- 1-D (parallel) RESTART CHECK")
#ATS:t104 = test(        SELF, "--solid True --graphics None --clearDirectories True  --checkError True  --dataDirBase 'dumps-spherical-reproducing' --domainIndependent True --outputFile 'Noh-spherical-1proc-reproducing.txt'", label="Spherical Noh problem with solid SPH -- 1-D (serial reproducing test setup)")
#ATS:t105 = testif(t104, SELF, "--solid True --graphics None --clearDirectories False  --checkError True  --dataDirBase 'dumps-spherical-reproducing' --domainIndependent True --outputFile 'Noh-spherical-4proc-reproducing.txt' --comparisonFile 'Noh-spherical-1proc-reproducing.txt'", np=4, label="Spherical Noh  problem with solid SPH -- 1-D (4 proc reproducing test)")

import os, sys, shutil
from SphericalSpheral import *
from GenerateSphericalNodeProfile1d import *
from SortAndDivideRedistributeNodes import distributeNodes1d
import mpi
import NohAnalyticSolution
from SpheralTestUtilities import *

title("Spherical Noh test run in spherical coordinates")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(KernelConstructor = WendlandC4Kernel3d,
            order = 5,

            nr = 100,
            rho1 = 1.0,
            eps1 = 0.0,
            smallPressure = False, #If set to True eps is not zero but small. 
            x0 = 0.0,
            x1 = 1.0,
            nPerh = 4.01,

            vr0 = -1.0, 
            vrSlope = 0.0,

            gamma = 5.0/3.0,
            mu = 1.0,

            solid = False,    # If true, use the fluid limit of the solid hydro option
            inflow = False,   # Should we impose inflow boundaries?

            crksph = False,
            psph = False,
            fsisph = False,
            gsph = False,
            crktype = "default",        # one of ("default", "variant")
            gsphReconstructionGradient = RiemannGradient, #one of (RiemannGradient, HydroAccelerationGradient, SPHGradient, MixedGradient, OnlyDvDxGradient)
            evolveTotalEnergy = False,  # Only for SPH variants -- evolve total rather than specific energy
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
            Cl = None, 
            Cq = None,
            Qself = 2.0,
            etaCritFrac = None,
            linearInExpansion = None,
            quadraticInExpansion = None,
            Qlimiter = None,
            balsaraCorrection = None,
            epsilon2 = None,
            QcorrectionOrder = None,
            hmin = 0.0001, 
            hmax = 0.1,
            cfl = 0.5,
            useVelocityMagnitudeForDt = False,
            XSPH = False,
            epsilonTensile = 0.0,
            nTensile = 4.0,

            IntegratorConstructor = VerletIntegrator,
            goalTime = 0.6,
            steps = None,
            dt = 0.0001,
            dtMin = 1.0e-10, 
            dtMax = 0.1,
            dtGrowth = 2.0,
            dtverbose = False,
            rigorousBoundaries = False,
            updateBoundaryFrequency = 1,
            maxSteps = None,
            statsStep = 1,
            smoothIters = 0,
            HUpdate = IdealH,
            correctionOrder = LinearOrder,
            volumeType = RKSumVolume,
            densityUpdate = RigorousSumDensity, # VolumeScaledDensity,
            compatibleEnergy = True,
            gradhCorrection = False,
            correctVelocityGradient = True,
            domainIndependent = True,
            cullGhostNodes = True,
            
            bArtificialConduction = False,
            arCondAlpha = 0.5,

            clearDirectories = True,
            checkError = False,
            checkRestart = False,
            restoreCycle = None,
            restartStep = 10000,
            dataDirBase = "dumps-spherical-Noh",
            restartBaseName = "Noh-spherical-1d",
            outputFile = None,
            comparisonFile = None,
            normOutputFile = None,
            writeOutputLabel = True,

            # Parameters for the test acceptance.,
            L1rho =   2.69219,   
            L2rho =   0.281965,  
            Linfrho = 30.5929,   
                                 
            L1P =     0.278906,  
            L2P =     0.0707854, 
            LinfP =   10.0544,   
                                 
            L1v =     0.0242732, 
            L2v =     0.00819671,
            Linfv =   0.91712,   
                                 
            L1eps =   0.0211726, 
            L2eps =   0.00273052,
            Linfeps = 0.325865,  
                                 
            L1h =     0.0013171, 
            L2h =     0.00036826,
            Linfh =   0.0267045, 

            tol = 1.0e-5,

            graphics = True,
            )

assert not(boolReduceViscosity and boolCullenViscosity)
assert not(gsph and (boolReduceViscosity or boolCullenViscosity))
assert not(fsisph and not solid)
if smallPressure:
    P0 = 1.0e-6
    eps1 = P0/((gamma - 1.0)*rho1)
   
if crksph:
    hydroname = os.path.join("CRKSPH",
                             str(volumeType),
                             str(correctionOrder))
elif fsisph:
    hydroname = "FSISPH"
elif gsph:
    hydroname = os.path.join("GSPH",str(gsphReconstructionGradient))
elif psph:
    hydroname = "PSPH"
else:
    hydroname = "SPH"
if solid:
    hydroname = "Solid" + hydroname

dataDir = os.path.join(dataDirBase,
                       hydroname,
                       "nPerh=%f" % nPerh,
                       "compatibleEnergy=%s" % compatibleEnergy,
                       "Cullen=%s" % boolCullenViscosity)
restartDir = os.path.join(dataDir, "restarts")
restartBaseName = os.path.join(restartDir, "Noh-spherical-1d-%i" % nr)

dx = (x1 - x0)/nr

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
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
# Note since this is in spherical coordinates, we need to build the 3D kernel
#-------------------------------------------------------------------------------
if KernelConstructor==NBSplineKernel3d:
    Wbase = NBSplineKernel3d(order)
else:
    Wbase = KernelConstructor()
WT = TableKernel3d(Wbase, 1000)
kernelExtent = WT.kernelExtent
output("WT")

#-------------------------------------------------------------------------------
# Make the NodeList.
#-------------------------------------------------------------------------------
if solid:
    nodes1 = makeSolidNodeList("nodes1", eos, 
                               hmin = hmin,
                               hmax = hmax,
                               nPerh = nPerh,
                               kernelExtent = kernelExtent)
else:
    nodes1 = makeFluidNodeList("nodes1", eos, 
                               hmin = hmin,
                               hmax = hmax,
                               nPerh = nPerh,
                               kernelExtent = kernelExtent)
    
output("nodes1")
output("nodes1.hmin")
output("nodes1.hmax")
output("nodes1.nodesPerSmoothingScale")

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
gen = GenerateSphericalNodeProfile1d(nr = nr,
                                     rho = rho1,
                                     rmin = x0,
                                     rmax = x1,
                                     nNodePerh = nPerh)
distributeNodes1d((nodes1, gen))
output("nodes1.numNodes")

# Set node specific thermal energies
nodes1.specificThermalEnergy(ScalarField("tmp", nodes1, eps1))
nodes1.massDensity(ScalarField("tmp", nodes1, rho1))

# Set node velocities
# Set the velocity of the node closest to the origin to zero to kick things off
pos = nodes1.positions()
vel = nodes1.velocity()
for ix in range(nodes1.numNodes):
    vel[ix].x = vr0 + vrSlope*pos[ix].x

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
    hydro = CRKSPH(dataBase = db,
                   order = correctionOrder,
                   cfl = cfl,
                   useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                   compatibleEnergyEvolution = compatibleEnergy,
                   evolveTotalEnergy = evolveTotalEnergy,
                   XSPH = XSPH,
                   densityUpdate = densityUpdate,
                   HUpdate = HUpdate,
                   crktype = crktype)
elif psph:
    hydro = PSPH(dataBase = db,
                 W = WT,
                 cfl = cfl,
                 useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                 compatibleEnergyEvolution = compatibleEnergy,
                 evolveTotalEnergy = evolveTotalEnergy,
                 correctVelocityGradient = correctVelocityGradient,
                 densityUpdate = densityUpdate,
                 HUpdate = HUpdate,
                 XSPH = XSPH)

elif fsisph:
    hydro = FSISPH(dataBase = db,
                   W = WT,
                   cfl = cfl,
                   interfaceMethod = ModulusInterface,
                   sumDensityNodeLists=[nodes1],                       
                   densityStabilizationCoefficient = 0.00,
                   useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                   compatibleEnergyEvolution = compatibleEnergy,
                   evolveTotalEnergy = evolveTotalEnergy,
                   correctVelocityGradient = correctVelocityGradient,
                   HUpdate = HUpdate)
elif gsph:
    limiter = VanLeerLimiter()
    waveSpeed = DavisWaveSpeed()
    solver = HLLC(limiter,
                  waveSpeed,
                  True,                           # False - first order , True - second order
                  gsphReconstructionGradient)     # what gradient are we using in reconstruction
    hydro = GSPH(dataBase = db,
                riemannSolver = solver,
                W = WT,
                cfl=cfl,
                useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
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
    hydro.Qself = Qself
    output("hydro.Qself")
output("hydro")
try:
    output("hydro.kernel")
    output("hydro.PiKernel")
except:
    pass
output("hydro.cfl")
output("hydro.compatibleEnergyEvolution")
output("hydro.densityUpdate")
output("hydro.XSPH")

packages = [hydro]

#-------------------------------------------------------------------------------
# Set the artificial viscosity parameters.
#-------------------------------------------------------------------------------
if not gsph:
    q = hydro.Q
    if not Cl is None:
        q.Cl = Cl
    if not Cq is None:
        q.Cq = Cq
    if not epsilon2 is None:
        q.epsilon2 = epsilon2
    if not Qlimiter is None:
        q.limiter = Qlimiter
    if not balsaraCorrection is None:
        q.balsaraShearCorrection = balsaraCorrection
    if not QcorrectionOrder is None:
        q.QcorrectionOrder = QcorrectionOrder
    output("q")
    output("q.Cl")
    output("q.Cq")
    output("q.epsilon2")
    output("q.limiter")
    output("q.balsaraShearCorrection")
    if hasattr(q, "linearInExpansion") and not linearInExpansion is None:
        q.linearInExpansion = linearInExpansion
        output("q.linearInExpansion")
    if hasattr(q, "quadraticInExpansion") and not quadraticInExpansion is None:
        q.quadraticInExpansion = quadraticInExpansion
        output("q.quadraticInExpansion")
    if hasattr(q, "etaCritFrac") and not etaCritFrac is None:
        q.etaCritFrac = etaCritFrac
        output("q.etaCritFrac")

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
bcs = []
if inflow:
    right_inflow = InflowOutflowBoundary(db, Plane(Vector(x1), Vector(-1)))
    bcs.append(right_inflow)
    packages.append(right_inflow)

for p in packages:
    for bc in bcs:
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
integrator.verbose = dtverbose
output("integrator")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.dtGrowth")
output("integrator.rigorousBoundaries")
output("integrator.updateBoundaryFrequency")
output("integrator.domainDecompositionIndependent")
output("integrator.cullGhostNodes")
output("integrator.verbose")

#-------------------------------------------------------------------------------
# Make the problem controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator,
                            kernel = hydro.kernel.baseKernel1d,
                            volumeType = volumeType,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle
                            )
output("control")

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
if not steps is None:
    if checkRestart:
        control.setRestartBaseName(restartBaseName + "_CHECK")
    control.step(steps)
    if checkRestart:
        control.setRestartBaseName(restartBaseName)

    # Are we doing the restart test?
    if checkRestart:
        state0 = State(db, integrator.physicsPackages())
        state0.copyState()
        print(control.totalSteps)
        control.loadRestartFile(control.totalSteps)
        state1 = State(db, integrator.physicsPackages())
        if not state1 == state0:
            raise ValueError("The restarted state does not match!")
        else:
            print("Restart check PASSED.")

else:
    if control.time() < goalTime:
        control.step(5)
        control.advance(goalTime, maxSteps)

#-------------------------------------------------------------------------------
# Compute the analytic answer.
#-------------------------------------------------------------------------------
rlocal = [pos.x for pos in nodes1.positions().internalValues()]
r = mpi.reduce(rlocal, mpi.SUM)
h1 = 1.0/(nPerh*dx)
answer = NohAnalyticSolution.NohSolution(3,
                                         r = r,
                                         v0 = -1.0,
                                         h0 = 1.0/h1)

# Compute the simulated specific entropy.
rho = mpi.allreduce(nodes1.massDensity().internalValues(), mpi.SUM)
Pf = ScalarField("pressure", nodes1)
nodes1.pressure(Pf)
P = mpi.allreduce(Pf.internalValues(), mpi.SUM)
A = [Pi/rhoi**gamma for (Pi, rhoi) in zip(P, rho)]

# The analytic solution for the simulated entropy.
xprof = mpi.allreduce([x.x for x in nodes1.positions().internalValues()], mpi.SUM)
xans, vans, uans, rhoans, Pans, hans = answer.solution(control.time(), xprof)
Aans = [Pi/rhoi**gamma for (Pi, rhoi) in zip(Pans,  rhoans)]
L1 = 0.0
for i in range(len(rho)):
  L1 = L1 + abs(rho[i]-rhoans[i])
L1_tot = L1 / len(rho)
if mpi.rank == 0 and outputFile:
 print("L1=",L1_tot,"\n")
 with open("Converge.txt", "a") as myfile:
    myfile.write("%s %s\n" % (nr, L1_tot))

#-------------------------------------------------------------------------------
# Plot the final state.
#-------------------------------------------------------------------------------
if graphics:
    from SpheralMatplotlib import *
    rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotState(db)
    plotAnswer(answer, control.time(), rhoPlot, velPlot, epsPlot, PPlot, HPlot = HPlot)
    EPlot = plotEHistory(control.conserve)
    plots = [(rhoPlot, "Noh-spherical-rho.png"),
             (velPlot, "Noh-spherical-vel.png"),
             (epsPlot, "Noh-spherical-eps.png"),
             (PPlot, "Noh-spherical-P.png"),
             (HPlot, "Noh-spherical-h.png")]

    # Plot the specific entropy.
    Aplot = newFigure()
    Aplot.plot(xprof, A, "ro", label="Simulation")
    Aplot.plot(xprof, Aans, "b-", label="Solution")
    Aplot.set_title("Specific entropy")
    plots.append((Aplot, "Noh-spherical-A.png"))
    
    if crksph:
        volPlot = plotFieldList(control.RKCorrections.volume, 
                                winTitle = "volume",
                                colorNodeLists = False, plotGhosts = False)
        splot = plotFieldList(control.RKCorrections.surfacePoint,
                              winTitle = "surface point",
                              colorNodeLists = False)
        plots += [(volPlot, "Noh-spherical-vol.png"),
                  (splot, "Noh-spherical-surfacePoint.png")]

    if boolCullenViscosity:
        cullAlphaPlot = plotFieldList(q.ClMultiplier(),
                                      winTitle = "Cullen alpha")
        cullDalphaPlot = plotFieldList(evolveCullenViscosityMultiplier.DalphaDt(),
                                       winTitle = "Cullen DalphaDt")
        plots += [(cullAlphaPlot, "Noh-spherical-Cullen-alpha.png"),
                  (cullDalphaPlot, "Noh-spherical-Cullen-DalphaDt.png")]

    if boolReduceViscosity:
        alphaPlotQ = plotFieldList(q.reducingViscosityMultiplierQ(),
                                  winTitle = "rvAlphaQ",
                                  colorNodeLists = False, plotGhosts = False)
        alphaPlotL = plotFieldList(q.reducingViscosityMultiplierL(),
                                   winTitle = "rvAlphaL",
                                   colorNodeLists = False, plotGhosts = False)

    # # Plot the grad h correction term (omega)
    # omegaPlot = plotFieldList(hydro.omegaGradh(),
    #                           winTitle = "grad h correction",
    #                           colorNodeLists = False)

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
xprof = mpi.reduce([x.x for x in nodes1.positions().internalValues()], mpi.SUM)

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
        multiSort(mo, xprof, rhoprof, Pprof, vprof, epsprof, hprof)
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

        #---------------------------------------------------------------------------
        # Also we can optionally compare the current results with another file.
        #---------------------------------------------------------------------------
        if comparisonFile:
            comparisonFile = os.path.join(dataDir, comparisonFile)
            import filecmp
            assert filecmp.cmp(outputFile, comparisonFile)

#------------------------------------------------------------------------------
# Compute the error.
#------------------------------------------------------------------------------
if mpi.rank == 0 :
    xans, vans, epsans, rhoans, Pans, hans = answer.solution(control.time(), xprof)
    import Pnorm
    print("\tQuantity \t\tL1 \t\t\tL2 \t\t\tLinf")
    failure = False
    hD = []

    if normOutputFile:
       f = open(normOutputFile, "a")
       if writeOutputLabel:
          f.write(("#" + 13*"%17s " + "\n") % ('"nx"',
                                               '"rho L1"', '"rho L2"', '"rho Linf"',
                                               '"P L1"',   '"P L2"',   '"P Linf"',
                                               '"vel L1"', '"vel L2"', '"vel Linf"',
                                               '"E L1"', '"E L2"', '"E Linf"',
                                               '"h L1"',   '"h L2"',   '"h Linf"'))
       f.write("%5i " % nr)
    for (name, data, ans,
         L1expect, L2expect, Linfexpect) in [("Mass Density", rhoprof, rhoans, L1rho, L2rho, Linfrho),
                                             ("Pressure", Pprof, Pans, L1P, L2P, LinfP),
                                             ("Velocity", vprof, vans, L1v, L2v, Linfv),
                                             ("Thermal E", epsprof, epsans, L1eps, L2eps, Linfeps),
                                             ("h       ", hprof, hans, L1h, L2h, Linfh)]:
        assert len(data) == len(ans)
        error = [data[i] - ans[i] for i in range(len(data))]
        Pn = Pnorm.Pnorm(error, xprof)
        L1 = Pn.gridpnorm(1, rmin, rmax)
        L2 = Pn.gridpnorm(2, rmin, rmax)
        Linf = Pn.gridpnorm("inf", rmin, rmax)
        print("\t%s \t\t%g \t\t%g \t\t%g" % (name, L1, L2, Linf))
        if normOutputFile:
           f.write((3*"%16.12e ") % (L1, L2, Linf))
        hD.append([L1,L2,Linf])

        if checkError:
            if not crksph and not psph and not fsisph and not gsph: # if sph use the known error norms
                if not fuzzyEqual(L1, L1expect, tol):
                    print("L1 error estimate for %s outside expected bounds: %g != %g" % (name,
                                                                                          L1,
                                                                                          L1expect))
                    failure = True
                if not fuzzyEqual(L2, L2expect, tol):
                    print("L2 error estimate for %s outside expected bounds: %g != %g" % (name,
                                                                                          L2,
                                                                                          L2expect))
                    failure = True
                if not fuzzyEqual(Linf, Linfexpect, tol):
                    print("Linf error estimate for %s outside expected bounds: %g != %g" % (name,
                                                                                            Linf,
                                                                                            Linfexpect))
                    failure = True

            if fsisph or gsph: # for fsi check if the norms are order of mag same as sph 
            
                if L1 > 2.0*L1expect:
                    print("L1 error estimate for %s outside expected bounds: %g != %g" % (name,
                                                                                          L1,
                                                                                          L1expect))
                    failure = True
                if L2 > 2.0*L2expect:
                    print("L2 error estimate for %s outside expected bounds: %g != %g" % (name,
                                                                                          L2,
                                                                                          L2expect))
                    failure = True
                if Linf > 2.0 * Linfexpect:
                    print("Linf error estimate for %s outside expected bounds: %g != %g" % (name,
                                                                                            Linf,
                                                                                            Linfexpect))
                    failure = True
    if normOutputFile:
       f.write("\n")
                                             
    if failure:
        raise ValueError("Error bounds violated.")

Eerror = (control.conserve.EHistory[-1] - control.conserve.EHistory[0])/control.conserve.EHistory[0]
print("Total energy error: %g" % Eerror)
if compatibleEnergy and abs(Eerror) > 1e-13:
    raise ValueError("Energy error outside allowed bounds.")
