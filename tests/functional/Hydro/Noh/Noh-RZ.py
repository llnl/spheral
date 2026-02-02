#-------------------------------------------------------------------------------
# The Noh test case run in RZ symmetry.
#
# W.F. Noh 1987, JCP, 72, 78-120.
#-------------------------------------------------------------------------------
#
# SPH
#
#ATS:sph0 = test(        SELF, "--hydroType SPH --goalTime 0.3 --graphics None --clearDirectories True  --checkError True   --restartStep 20", label="Planar RZ Noh problem (SPH serial)")
#ATS:sph1 = testif(sph0, SELF, "--hydroType SPH --goalTime 0.3 --graphics None --clearDirectories False --checkError False  --restartStep 20 --restoreCycle 20 --steps 20 --checkRestart True", label="Planar RZ Noh problem (SPH serial) RESTART CHECK")
#ATS:sph2 = test(        SELF, "--hydroType SPH --goalTime 0.3 --graphics None --clearDirectories True  --checkError True  --dataDirBase 'dumps-rz-planar-restartcheck' --restartStep 20", np=8, label="Planar Noh RZ problem (SPH parallel)")
#ATS:sph3 = testif(sph2, SELF, "--hydroType SPH --goalTime 0.3 --graphics None --clearDirectories False --checkError False --dataDirBase 'dumps-rz-planar-restartcheck' --restartStep 20 --restoreCycle 20 --steps 20 --checkRestart True", np=8, label="Planar RZ Noh problem -- (SPH parallel) RESTART CHECK")
#
# ASPH
#
#ATS:asph2 = test(         SELF, "--hydroType SPH --goalTime 0.3 --asph True --graphics None --clearDirectories True  --checkError True  --dataDirBase 'dumps-rz-planar-restartcheck' --restartStep 20", np=8, label="Planar Noh RZ problem (ASPH parallel)")
#ATS:asph3 = testif(asph2, SELF, "--hydroType SPH --goalTime 0.3 --asph True --graphics None --clearDirectories False --checkError False --dataDirBase 'dumps-rz-planar-restartcheck' --restartStep 20 --restoreCycle 20 --steps 20 --checkRestart True", np=8, label="Planar RZ Noh problem -- (ASPH parallel) RESTART CHECK")
#
# ASPHClassic
#
#ATS:acsph2 = test(          SELF, "--hydroType SPH --goalTime 0.3 --asph Classic --graphics None --clearDirectories True  --checkError True  --dataDirBase 'dumps-rz-planar-restartcheck' --restartStep 20", np=8, label="Planar Noh RZ problem (ASPH Classic parallel)")
#ATS:acsph3 = testif(acsph2, SELF, "--hydroType SPH --goalTime 0.3 --asph Classic --graphics None --clearDirectories False --checkError False --dataDirBase 'dumps-rz-planar-restartcheck' --restartStep 20 --restoreCycle 20 --steps 20 --checkRestart True", np=8, label="Planar RZ Noh problem -- (ASPH Classic parallel) RESTART CHECK")
#
# CRKSPH
#
#ATS:crk2 = test(        SELF, "--hydroType CRKSPH --goalTime 0.3 --graphics None --clearDirectories True  --checkError True  --dataDirBase 'dumps-rz-planar-restartcheck' --restartStep 20 --tol 5e-3", np=8, label="Planar Noh RZ problem (CRKSPH parallel)")               # Only need tolerance override for BlueOS
#ATS:crk3 = testif(crk2, SELF, "--hydroType CRKSPH --goalTime 0.3 --graphics None --clearDirectories False --checkError False --dataDirBase 'dumps-rz-planar-restartcheck' --restartStep 20 --restoreCycle 20 --steps 20 --checkRestart True", np=8, label="Planar RZ Noh problem -- (CRKSPH parallel) RESTART CHECK")

import os, sys, shutil, mpi
import numpy as np
from SpheralRZ import *
from SpheralTestUtilities import *

from GenerateNodeDistribution2d import *
if mpi.procs > 1:
    from VoronoiDistributeNodes import distributeNodes2d
else:
    from DistributeNodes import distributeNodes2d

title("RZ hydro test -- Noh problem")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(problem = "planar",     # one of (planar, cylindrical, spherical)
            KernelConstructor = WendlandC4Kernel,

            n1 = 100,
            n2 = 20,

            nPerh = 4.01,

            gamma = 5.0/3.0,
            mu = 1.0,

            solid = False,                     # If true, use the fluid limit of the solid hydro option
            asph = False,                      # This just chooses the H algorithm -- you can use this with CRKSPH for instance.

            hydroType = "SPH",                 # one of (SPH, CRKSPH)
            evolveTotalEnergy = False,         # Only for SPH variants -- evolve total rather than specific energy
            boolReduceViscosity = False,
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
            Qlimiter = False,
            balsaraCorrection = False,
            epsilon2 = 1e-2,
            hmin = 0.0001, 
            hmax = 0.1,
            hminratio = 0.1,
            cfl = 0.5,
            useVelocityMagnitudeForDt = False,
            XSPH = True,
            epsilonTensile = 0.0,
            nTensile = 4.0,

            IntegratorConstructor = CheapSynchronousRK2Integrator,
            goalTime = 0.6,
            steps = None,
            dt = 0.0001,
            dtMin = 1.0e-5, 
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
            checkRestart = False,
            checkEnergy = False,
            restoreCycle = -1,
            restartStep = 10000,
            dataDirBase = "dumps-rz-Noh",
            outputFile = "Noh-RZ.gnu",
            comparisonFile = None,
            normOutputFile = None,
            writeOutputLabel = True,

            # Parameters for the test acceptance.,
            tol = 1.0e-5,
            checkError = False,

            graphics = True,
            )

hydroType = hydroType.upper()

assert not(boolReduceViscosity and boolCullenViscosity)
   
assert problem in ("planar", "cylindrical", "spherical")
rho0 = 1.0
eps0 = 0.0

if hydroType == "CRKSPH":
    gradhCorrection = False

hydroPath = (("Solid" if solid else "") +
             ("AC" if asph == "Classic" else "A" if asph else "") +
             hydroType)
dataDir = os.path.join(dataDirBase,
                       hydroPath,
                       "nPerh=%f" % nPerh,
                       "compatibleEnergy=%s" % compatibleEnergy,
                       "Cullen=%s" % boolCullenViscosity)
restartDir = os.path.join(dataDir, "restarts")
restartBaseName = os.path.join(restartDir, "Noh-%s-RZ" % problem)

vizDir = os.path.join(dataDir, "visit")
if vizTime is None and vizCycle is None:
    vizBaseName = None
else:
    vizBaseName = "Noh-%s-RZ" % problem

#-------------------------------------------------------------------------------
# The reference values for error norms checking for pass/fail
#-------------------------------------------------------------------------------
LnormRef = {"SPH": {"Mass density" : {"L1"   : 0.927051,   
                                      "L2"   : 0.0257097,  
                                      "Linf" : 3.09951},         
                    "Pressure    " : {"L1"   : 0.417473,   
                                      "L2"   : 0.0117583,  
                                      "Linf" : 1.97256},   
                    "Velocity    " : {"L1"   : 0.323763,   
                                      "L2"   : 0.00903822, 
                                      "Linf" : 1.02008},   
                    "Spec Therm E" : {"L1"   : 0.168691,   
                                      "L2"   : 0.00466419, 
                                      "Linf" : 0.987794},  
                    "h           " : {"L1"   : 0.0082967,  
                                      "L2"   : 0.000186032,
                                      "Linf" : 0.0190145}},
            "ASPH": {"Mass density" : {"L1"   : 0.931764,    
                                       "L2"   : 0.0259673,   
                                       "Linf" : 3.07466},     
                     "Pressure    " : {"L1"   : 0.415028,    
                                       "L2"   : 0.0116435,   
                                       "Linf" : 1.52774},    
                     "Velocity    " : {"L1"   : 0.322906,    
                                       "L2"   : 0.00897806,  
                                       "Linf" : 1.01521},    
                     "Spec Therm E" : {"L1"   : 0.165872,    
                                       "L2"   : 0.00453749,  
                                       "Linf" : 0.86993},    
                     "h           " : {"L1"   : 0.0100794,  
                                       "L2"   : 0.000266552, 
                                       "Linf" : 0.0300606}},
            "ACSPH": {"Mass density" : {"L1"   : 0.910927,   
                                        "L2"   : 0.0257385,  
                                        "Linf" : 3.06255},    
                      "Pressure    " : {"L1"   : 0.40629,    
                                        "L2"   : 0.0115813,  
                                        "Linf" : 1.59406},   
                      "Velocity    " : {"L1"   : 0.31841,    
                                        "L2"   : 0.00894237, 
                                        "Linf" : 1.01271},   
                      "Spec Therm E" : {"L1"   : 0.162679,   
                                        "L2"   : 0.00452432, 
                                        "Linf" : 0.899942},  
                      "h           " : {"L1"   : 0.00941208, 
                                        "L2"   : 0.000237763,
                                        "Linf" : 0.0268561}},
            "CRKSPH": {"Mass density" : {"L1"   : 0.918847,    
                                         "L2"   : 0.0251823,   
                                         "Linf" : 3.29814},     
                       "Pressure    " : {"L1"   : 0.403824,    
                                         "L2"   : 0.0113766,   
                                         "Linf" : 1.76299},    
                       "Velocity    " : {"L1"   : 0.315356,    
                                         "L2"   : 0.00890281,  
                                         "Linf" : 1.06319},    
                       "Spec Therm E" : {"L1"   : 0.160593,    
                                         "L2"   : 0.00451318,  
                                         "Linf" : 0.861714},   
                       "h           " : {"L1"   : 0.00799764, 
                                         "L2"   : 0.000181642, 
                                         "Linf" : 0.0196302}},
}

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
eos = GammaLawGasMKS(gamma, mu, minimumPressure=0.0)
strength = NullStrength()

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
WT = TableKernel(KernelConstructor(), 200)
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
                               kernelExtent = kernelExtent)
else:
    nodes1 = makeFluidNodeList("nodes1", eos, 
                               hmin = hmin,
                               hmax = hmax,
                               hminratio = hminratio,
                               nPerh = nPerh,
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
    vz0 = -1.0
    vr0 = 0.0
elif problem == "cylindrical":
    nz = n2
    nr = n1
    z0, z1 = 0.0, 0.2
    r0, r1 = 0.0, 1.0
    rmin, rmax = None, None
    vz0 = 0.0
    vr0 = -1.0
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

# Set node specific thermal energies
nodes1.specificThermalEnergy(ScalarField("tmp", nodes1, eps0))
nodes1.massDensity(ScalarField("tmp", nodes1, rho0))

# Set node velocities
pos = nodes1.positions()
vel = nodes1.velocity()
if problem == "spherical":
    for i in range(nodes1.numNodes):
        vel[i] = -1.0 * pos[i].unitVector()
else:
    for i in range(nodes1.numNodes):
        vel[i] = Vector(vz0, vr0)

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
if hydroType == "CRKSPH":
    hydro = CRKSPH(dataBase = db,
                   W = WT,
                   cfl = cfl,
                   useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                   compatibleEnergyEvolution = compatibleEnergy,
                   evolveTotalEnergy = evolveTotalEnergy,
                   XSPH = XSPH,
                   order = correctionOrder,
                   densityUpdate = densityUpdate,
                   HUpdate = HUpdate,
                   ASPH = asph)
else:
    assert hydroType == "SPH"
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
                ASPH = asph,
                XSPH = XSPH,
                epsTensile = epsilonTensile,
                nTensile = nTensile)
output("hydro")
output("hydro.cfl")
output("hydro.compatibleEnergyEvolution")
output("hydro.densityUpdate")
output("hydro._smoothingScaleMethod.HEvolution")

packages = [hydro]

#-------------------------------------------------------------------------------
# Set the artificial viscosity parameters.
#-------------------------------------------------------------------------------
q = hydro.Q
q.Cl = Cl
q.Cq = Cq
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
    #q.reducingViscosityCorrection = True
    ArtyCond = ArtificialConduction(WT,arCondAlpha)
    
    packages.append(ArtyCond)

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
if problem == "planar":
    bcs = [ReflectingBoundary(Plane(Vector(z0, r0), Vector( 1.0,  0.0))),
           ReflectingBoundary(Plane(Vector(z1, r0), Vector(-1.0,  0.0))),
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
import NohAnalyticSolution
if problem == "planar":
    xprof = mpi.allreduce([x.x for x in nodes1.positions().internalValues()], mpi.SUM)
    h1 = 1.0/(nPerh/n1)
    answer = NohAnalyticSolution.NohSolution(1,
                                             r = xprof,
                                             v0 = -1.0,
                                             h0 = 1.0/h1)
elif problem == "cylindrical":
    xprof = mpi.allreduce([x.y for x in nodes1.positions().internalValues()], mpi.SUM)
    h1 = 1.0/(nPerh/n1)
    answer = NohAnalyticSolution.NohSolution(2,
                                             r = xprof,
                                             v0 = -1.0,
                                             h0 = 1.0/h1)
else:
    xprof = mpi.allreduce([x.magnitude() for x in nodes1.positions().internalValues()], mpi.SUM)
    h1 = 1.0/(nPerh/n1)
    answer = NohAnalyticSolution.NohSolution(3,
                                             r = xprof,
                                             v0 = -1.0,
                                             h0 = 1.0/h1)

# Compute the simulated specific entropy.
rho = mpi.allreduce(nodes1.massDensity().internalValues(), mpi.SUM)
Pf = ScalarField("pressure", nodes1)
nodes1.pressure(Pf)
P = mpi.allreduce(Pf.internalValues(), mpi.SUM)
A = [Pi/rhoi**gamma for (Pi, rhoi) in zip(P, rho)]

# Solution profiles.
xans, vans, uans, rhoans, Pans, hans = answer.solution(control.time(), xprof)
Aans = [Pi/rhoi**gamma for (Pi, rhoi) in zip(Pans,  rhoans)]
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
    plotAnswer(answer, control.time(), rhoPlot=rhoPlot, velPlot=velPlot, epsPlot=epsPlot, PPlot=PPlot, HPlot=HPlot,
               plotStyle = "kx")
    EPlot = plotEHistory(control.conserve)
    plots = [(rhoPlot, "Noh-%s-rho-RZ.png" % problem),
             (velPlot, "Noh-%s-vel-RZ.png" % problem),
             (epsPlot, "Noh-%s-eps-RZ.png" % problem),
             (PPlot, "Noh-%s-P-RZ.png" % problem),
             (HPlot, "Noh-%s-h-RZ.png" % problem)]

    # Plot the specific entropy.
    Aplot = newFigure()
    Aplot.plot(xprof, A, "ro")
    Aplot.plot(xprof, Aans, "kx")
    plt.title("Specific entropy")
    plots.append((Aplot, "Noh-%s-A.png" % problem))
    
    # Throw the positions out there too.
    posPlot = plotNodePositions2d(db)
    plt.xlabel("z")
    plt.ylabel("r")
    plt.title("Node positions @ t=%g" % control.time())
    plots.append((posPlot, "Noh-%s-positions.png" % problem))

    if hydroType == "CRKSPH":
        volPlot = plotFieldList(control.RKCorrections.volume,
                                xFunction = "%s.y",
                                winTitle = "volume",
                                colorNodeLists = False, plotGhosts = False)
        plots.append((volPlot, "Noh-%s-vol.png" % problem))

    if boolCullenViscosity:
        cullAlphaPlot = plotFieldList(q.ClMultiplier(),
                                      xFunction = "%s.y",
                                      winTitle = "Cullen alpha")
        cullDalphaPlot = plotFieldList(evolveCullenViscosityMultiplier.DalphaDt(),
                                       xFunction = "%s.y",
                                       winTitle = "Cullen DalphaDt")
        plots += [(cullAlphaPlot, "Noh-%s-Cullen-alpha.png" % problem),
                  (cullDalphaPlot, "Noh-%s-Cullen-DalphaDt.png" % problem)]

    if boolReduceViscosity:
        alphaPlotQ = plotFieldList(q.reducingViscosityMultiplierQ(),
                                   xFunction = "%s.y",
                                  winTitle = "rvAlphaQ",
                                  colorNodeLists = False, plotGhosts = False)
        alphaPlotL = plotFieldList(q.reducingViscosityMultiplierL(),
                                   xFunction = "%s.y",
                                   winTitle = "rvAlphaL",
                                   colorNodeLists = False, plotGhosts = False)

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
    mprof = mpi.reduce(nodes1.mass().internalValues(), mpi.SUM)
    rhoprof = mpi.reduce(nodes1.massDensity().internalValues(), mpi.SUM)
    P = ScalarField("pressure", nodes1)
    nodes1.pressure(P)
    Pprof = mpi.reduce(P.internalValues(), mpi.SUM)
    vprof = mpi.reduce([v.x for v in nodes1.velocity().internalValues()], mpi.SUM)
    epsprof = mpi.reduce(nodes1.specificThermalEnergy().internalValues(), mpi.SUM)
    hprof = mpi.reduce([1.0/H.xx for H in nodes1.Hfield().internalValues()], mpi.SUM)
    if mpi.rank == 0:
        multiSort(xprof, rhoprof, Pprof, vprof, epsprof, hprof,
                  rhoans, Pans, vans, uans, hans)
        f = open(outputFile, "w")
        f.write(("#  " + 12*"'%s' " + "\n") % ("x", "m", "rho", "P", "v", "eps", "h",
                                               "rhoans", "Pans", "vans", "epsans", "hans"))
        for (xi, mi, rhoi, Pi, vi, epsi, hi, 
             rhoansi, Pansi, vansi, uansi, hansi) in zip(xprof, mprof, rhoprof, Pprof, vprof, epsprof, hprof, 
                                                         rhoans, Pans, vans, uans, hans):
            f.write((12*"%16.12e " + '\n') % 
                    (xi, mi, rhoi, Pi, vi, epsi, hi, 
                     rhoansi, Pansi, vansi, uansi, hansi))
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
if mpi.rank == 0:
    xans, vans, epsans, rhoans, Pans, hans = answer.solution(control.time(), xprof)
    import Pnorm
    print("\tQuantity \t\tL1 \t\t\tL2 \t\t\tLinf")
    failure = False

    if normOutputFile:
       f = open(normOutputFile, "a")
       if writeOutputLabel:
          f.write(("#" + 13*"%17s " + "\n") % ('"n"',
                                               '"rho L1"', '"rho L2"', '"rho Linf"',
                                               '"P L1"',   '"P L2"',   '"P Linf"',
                                               '"vel L1"', '"vel L2"', '"vel Linf"',
                                               '"E L1"', '"E L2"', '"E Linf"',
                                               '"h L1"',   '"h L2"',   '"h Linf"'))
       f.write("%5i " % nz)
    for (name, data, ans) in [("Mass density", rhoprof, rhoans),
                              ("Pressure    ", Pprof, Pans),
                              ("Velocity    ", vprof, vans),
                              ("Spec Therm E", epsprof, epsans),
                              ("h           ", hprof, hans)]:
        assert len(data) == len(ans)
        error = [data[i] - ans[i] for i in range(len(data))]
        Pn = Pnorm.Pnorm(error, xprof)
        L1 = Pn.gridpnorm(1, rmin, rmax)
        L2 = Pn.gridpnorm(2, rmin, rmax)
        Linf = Pn.gridpnorm("inf", rmin, rmax)
        print("\t%s \t\t%g \t\t%g \t\t%g" % (name, L1, L2, Linf))
        if normOutputFile:
           f.write((3*"%16.12e ") % (L1, L2, Linf))

        if checkError and not (np.allclose(L1, LnormRef[hydroPath][name]["L1"], tol, tol) and
                               np.allclose(L2, LnormRef[hydroPath][name]["L2"], tol, tol) and
                               np.allclose(Linf, LnormRef[hydroPath][name]["Linf"], tol, tol)):
            print("Failing Lnorm tolerance for ", name, (L1, L2, Linf), LnormRef[hydroPath][name])
            failure = True
  
    if normOutputFile:
       f.write("\n")

    if checkError and failure:
        raise ValueError("Error bounds violated.")

Eerror = (control.conserve.EHistory[-1] - control.conserve.EHistory[0])/control.conserve.EHistory[0]
print("Total energy error: %g" % Eerror)
if checkEnergy and abs(Eerror) > 1e-13:
    raise ValueError("Energy error outside allowed bounds.")
