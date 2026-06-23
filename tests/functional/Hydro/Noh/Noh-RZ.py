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
# CRKSPH   # Suspeneded until updated to new RZ formalism
#
# #ATS:crk2 = test(        SELF, "--hydroType CRKSPH --goalTime 0.3 --graphics None --clearDirectories True  --checkError True  --dataDirBase 'dumps-rz-planar-restartcheck' --restartStep 20 --tol 5e-3", np=8, label="Planar Noh RZ problem (CRKSPH parallel)")               # Only need tolerance override for BlueOS
# #ATS:crk3 = testif(crk2, SELF, "--hydroType CRKSPH --goalTime 0.3 --graphics None --clearDirectories False --checkError False --dataDirBase 'dumps-rz-planar-restartcheck' --restartStep 20 --restoreCycle 20 --steps 20 --checkRestart True", np=8, label="Planar RZ Noh problem -- (CRKSPH parallel) RESTART CHECK")

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

            seed = "lattice",
            nPerh = 6.01,

            gamma = 5.0/3.0,
            mu = 1.0,

            vetaramp = 1.0,                    # Optionally ramp velocity to zero as we approach the origin (0 => no ramp)

            solid = False,                     # If true, use the fluid limit of the solid hydro option
            asph = False,                      # This just chooses the H algorithm -- you can use this with CRKSPH for instance.

            hydroType = "SPH",                 # one of (SPH, CRKSPH)
            Q = None,                          # optionally override viscosity choice
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
            Cl = None, #1.5,
            Cq = None, #1.0,
            linearInExpansion = None,
            quadraticInExpansion = None,
            balsaraCorrection = None,
            epsilon2 = None,
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
            dtMin = 1.0e-8, 
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
                       problem,
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
LnormRef = {"SPH": {"Mass density" : {"L1"   : 0.203217,   
                                      "L2"   : 0.0254493,  
                                      "Linf" : 4.66596},         
                    "Pressure    " : {"L1"   : 0.0395288,  
                                      "L2"   : 0.00473042, 
                                      "Linf" : 0.689009},  
                    "Velocity    " : {"L1"   : 0.0254315,  
                                      "L2"   : 0.00535197, 
                                      "Linf" : 0.918963},  
                    "Spec Therm E" : {"L1"   : 0.0227909,  
                                      "L2"   : 0.00309336, 
                                      "Linf" : 0.381547},  
                    "h           " : {"L1"   : 0.00508196, 
                                      "L2"   : 0.000938975,
                                      "Linf" : 0.0399}},   
            "ASPH": {"Mass density" : {"L1"   : 0.197549,   
                                       "L2"   : 0.0306595,  
                                       "Linf" : 15.8962},    
                     "Pressure    " : {"L1"   : 0.0336372,  
                                       "L2"   : 0.00574558, 
                                       "Linf" : 1.08521},   
                     "Velocity    " : {"L1"   : 0.0214818,  
                                       "L2"   : 0.00529686, 
                                       "Linf" : 0.945891},  
                     "Spec Therm E" : {"L1"   : 0.0217109,  
                                       "L2"   : 0.00362169, 
                                       "Linf" : 0.47045},   
                     "h           " : {"L1"   : 0.00518618,
                                       "L2"   : 0.00140424, 
                                       "Linf" : 0.0399}},  
            "ACSPH": {"Mass density" : {"L1"   : 0.194015,   
                                        "L2"   : 0.0286555,  
                                        "Linf" : 7.66739},    
                      "Pressure    " : {"L1"   : 0.0366169,  
                                        "L2"   : 0.00541615, 
                                        "Linf" : 0.859485},  
                      "Velocity    " : {"L1"   : 0.0252846,  
                                        "L2"   : 0.00550918, 
                                        "Linf" : 0.943121},  
                      "Spec Therm E" : {"L1"   : 0.0232088,  
                                        "L2"   : 0.00350656, 
                                        "Linf" : 0.43536},   
                      "h           " : {"L1"   : 0.00476212, 
                                        "L2"   : 0.00114821, 
                                        "Linf" : 0.0399}},   
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
    theta = pi/2.0
elif problem == "cylindrical":
    nz = n2
    nr = n1
    z0, z1 = 0.0, 0.2
    r0, r1 = 0.0, 1.0
    rmin, rmax = None, None
    vz0 = 0.0
    vr0 = -1.0
    theta = pi/2.0
else:
    assert problem == "spherical"
    nz = n1
    nr = n1
    rmin, rmax = 0.0, 1.0
    z0, z1 = 0.0, 1.0
    r0, r1 = 0.0, 1.0
    theta = pi/2.0

generator = GenerateNodeDistribution2d(nz, nr, rho0, seed,
                                       xmin = (z0, r0),
                                       xmax = (z1, r1),
                                       rmin = rmin,
                                       rmax = rmax,
                                       theta = theta,
                                       nNodePerh = nPerh,
                                       SPH = not asph)

distributeNodes2d((nodes1, generator))
output("mpi.reduce(nodes1.numInternalNodes, mpi.MIN)")
output("mpi.reduce(nodes1.numInternalNodes, mpi.MAX)")
output("mpi.reduce(nodes1.numInternalNodes, mpi.SUM)")

# Set node specific thermal energies
nodes1.specificThermalEnergy(ScalarField("tmp", nodes1, eps0))
nodes1.massDensity(ScalarField("tmp", nodes1, rho0))

def velRamp(posi, Hi):
    def etai(eta):
        if problem == "planar":
            return eta.x
        elif problem == "cylindrical":
            return eta.y
        else:
            return eta.magnitude()
    if vetaramp > 0.0:
        xeta = etai(Hi*posi)
        return min(1.0, xeta/vetaramp)
    else:
        return 1.0

# Set node velocities
pos = nodes1.positions()
vel = nodes1.velocity()
H = nodes1.Hfield()
if problem == "spherical":
    for i in range(nodes1.numNodes):
        vel[i] = -1.0 * pos[i].unitVector() * velRamp(pos[i], H[i])
else:
    for i in range(nodes1.numNodes):
        vel[i] = Vector(vz0, vr0) * velRamp(pos[i], H[i])

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
if not Q is None:
    Q = Q(Cl, Cq, WT)

if hydroType == "CRKSPH":
    hydro = CRKSPH(dataBase = db,
                   W = WT,
                   Q = Q,
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
                Q = Q,
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
output("hydro.XSPH")
output("hydro._smoothingScaleMethod.HEvolution")

packages = [hydro]

#-------------------------------------------------------------------------------
# Set the artificial viscosity parameters.
#-------------------------------------------------------------------------------
q = hydro.Q
def setOptionalParam(thing, param, val):
    if not val is None:
        exec(f"thing.{param} = {val}")
    return
setOptionalParam(q, "Cl", Cl)
setOptionalParam(q, "Cq", Cq)
setOptionalParam(q, "epsilon2", epsilon2)
setOptionalParam(q, "linearInExpansion", linearInExpansion)
setOptionalParam(q, "quadraticInExpansion", quadraticInExpansion)
setOptionalParam(q, "balsaraShearCorrection", balsaraCorrection)
output("q")
output("q.Cl")
output("q.Cq")
output("q.epsilon2")
output("q.balsaraShearCorrection")
output("q.linearInExpansion")
output("q.quadraticInExpansion")

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
bcs = []
if problem == "planar":
    bcs += [ReflectingBoundary(Plane(Vector(z0, r0), Vector( 1.0,  0.0))),
            ReflectingBoundary(Plane(Vector(z1, r0), Vector(-1.0,  0.0))),
            ReflectingBoundary(Plane(Vector(z0, r1), Vector( 0.0, -1.0)))]
    if r0 != 0.0:
        bcs.append(ReflectingBoundary(Plane(Vector(z0, r0), Vector( 0.0, 1.0))))
elif problem == "cylindrical":
    bcs += [ReflectingBoundary(Plane(Vector(z0, r0), Vector( 1.0,  0.0))),
            ReflectingBoundary(Plane(Vector(z1, r0), Vector(-1.0,  0.0)))]
else:
    assert problem == "spherical"
    if theta == pi/2.0:
        bcs += [ReflectingBoundary(Plane(Vector(z0, r0), Vector( 1.0,  0.0)))]

for bc in bcs:
    for p in packages:
        p.appendBoundary(bc)

#-------------------------------------------------------------------------------
# Construct an integrator.
#-------------------------------------------------------------------------------
integrator = IntegratorConstructor(db, packages)
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
        rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotState(db, xFunction="%s.x", vecyFunction="%s.x", tenyFunction="1.0/%s.xx", plotAverage=True)
        xfunc = "%s.x"
    elif problem == "cylindrical":
        rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotState(db, xFunction="%s.y", vecyFunction="%s.y", tenyFunction="1.0/%s.yy", plotAverage=True)
        xfunc = "%s.y"
    else:
        rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotRadialState(db, plotAverage=True)
        xfunc = "%s.magnitude()"
    plotAnswer(answer, control.time(), rhoPlot=rhoPlot, velPlot=velPlot, epsPlot=epsPlot, PPlot=PPlot, HPlot=HPlot,
               plotStyle = "k-",
               x = np.linspace(0.0, max(xprof), 500))
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

    Qplot = plotFieldList(hydro.Q.maxViscousPressure,
                          xFunction = xfunc,
                          winTitle = "Max Q pressure")

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
    if problem == "planar":
        vprof = mpi.reduce([v.x for v in nodes1.velocity().internalValues()], mpi.SUM)
    elif problem == "cylindrical":
        vprof = mpi.reduce([v.y for v in nodes1.velocity().internalValues()], mpi.SUM)
    else:
        assert problem == "spherical"
        vprof = mpi.reduce([v.dot(x.unitVector()) for v, x in zip(nodes1.velocity().internalValues(),
                                                                  nodes1.positions().internalValues())], mpi.SUM)
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
