#-------------------------------------------------------------------------------
# A disk of stainless steel undergoing tensile strain.  This is intended as a
# test of the cracking/failure models.
#-------------------------------------------------------------------------------
from SolidSpheral import *
from SpheralTestUtilities import *
from findLastRestart import *
from SpheralVisitDump import dumpPhysicsState
from identifyFragments import identifyFragments, fragmentProperties
from math import *

from WeibullDamage import *

# Load the mpi module if we're parallel.
import loadmpi
mpi, procID, numProcs = loadmpi.loadmpi()

#-------------------------------------------------------------------------------
# Identify ourselves!
#-------------------------------------------------------------------------------
title("3-D Tensile disk strength/damage model test")

#-------------------------------------------------------------------------------
# Generic problem parameters
# All CGS units.
#-------------------------------------------------------------------------------
commandLine(SolidNodeListConstructor = SphSolidNodeList3d,
            FluidNodeListConstructor = SphNodeList3d,

            # How much of the 2 Pi geometry are we doing?
            phiFactor = 2.0,

            radius = 10.0,
            thickness = 0.3,
            nx = 5,
            ny = 100,
            nPerh = 2.01,

            rho0 = 7.9,

            # Material specific bounds on the mass density.
            etamin = 0.5,
            etamax = 1.5,

            # Parameters for the time dependent strain and cracking.
            DamageModelConstructor = GradyKippTensorDamage3d, # GradyKippScalarDamage3d, # WeibullTensorDamage3d, # GradyKippVectorDamage3d, # 
            v0 = 5.0e4,
            kWeibullFactor = 1.0,
            mWeibullFactor = 1.0,
            randomSeed = 548928513,
            strainType = TensorDamageModel3d.TensorStrainAlgorithm.StrainHistory,
            useDamageGradient = False,
            emulateScalarDamage = False,
            cullToWeakestFlaws = False,

            Qconstructor = MonaghanGingoldViscosity3d,
            Cl = 1.0,
            Cq = 1.0,
            Qlimiter = False,
            balsaraCorrection = False,
            epsilon2 = 1e-2,
            negligibleSoundSpeed = 1e-5,
            csMultiplier = 1e-4,
            hmin = 1e-5,
            hmax = 20.0,
            cfl = 0.5,
            useVelocityMagnitudeForDt = False,
            XSPH = True,
            epsilonTensile = 0.0,
            nTensile = 4,
            hybridMassDensityThreshold = 0.01,

            goalTime = 100.0e-6,
            dt = 1e-10,
            dtMin = 1e-12,
            dtMax = 1e-5,
            dtGrowth = 2.0,
            dumpFrac = 0.005,
            maxSteps = None,
            statsStep = 10,
            smoothIters = 0,
            HEvolution = Hydro3d.HEvolutionType.IdealH,
            sumForMassDensity = Hydro3d.MassDensityType.IntegrateDensity, # HybridDensity, # CorrectedSumDensity, # VolumeScaledDensity, #
            compatibleEnergyEvolution = True,

            restoreCycle = None,
            restartStep = 1000,

            basedir = ".",

            )

phi = pi * phiFactor

# Assuming strain ~ 4.5e3 sec^-1, L ~ 1.4cm
kWeibull = 8.8e4 * kWeibullFactor
mWeibull = 2.63  * mWeibullFactor

dataDir = "%s/dumps-TensileDisk-3d-%ix%i/%s/%s/gradDamage=%s/k=%4.2f_m=%4.2f" % (basedir, nx, ny,
                                                                                 str(DamageModelConstructor).split("'")[1],
                                                                                 str(SolidNodeListConstructor).split("'")[1],
                                                                                 useDamageGradient,
                                                                                 kWeibull,
                                                                                 mWeibull)
restartDir = dataDir + "/restarts/proc-%04i" % procID
visitDir = dataDir + "/visit"
restartBaseName = restartDir + "/TensileDisk-%ix%i" % (nx, ny)

xmin = (0.0, 0.0)
xmax = (thickness, radius)

volume = pi*radius**2 * thickness
dx = thickness/nx
dy = radius/ny

dtSample = dumpFrac*goalTime

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
import os, sys
if mpi.rank == 0:
    if not os.path.exists(restartDir):
        os.makedirs(restartDir)
    if not os.path.exists(visitDir):
        os.makedirs(visitDir)
mpi.barrier()
if not os.path.exists(restartDir):
    os.makedirs(restartDir)
mpi.barrier()

#-------------------------------------------------------------------------------
# If we're restarting, find the set of most recent restart files.
#-------------------------------------------------------------------------------
if restoreCycle is None:
    restoreCycle = findLastRestart(restartBaseName)

#-------------------------------------------------------------------------------
# Stainless steel material properties.
#-------------------------------------------------------------------------------
eos = GruneisenEquationOfStateCGS3d(7.90,    # reference density  
                                    etamin,  # etamin             
                                    etamax,  # etamax             
                                    0.457e6, # C0                 
                                    1.49,    # S1                 
                                    0.0,     # S2                 
                                    0.0,     # S3                 
                                    1.93,    # gamma0             
                                    0.5,     # b                  
                                    55.350)  # atomic weight
coldFit = NinthOrderPolynomialFit(-1.06797724e10,
                                  -2.06872020e10,
                                   8.24893246e11,
                                  -2.39505843e10,
                                  -2.44522017e10,
                                   5.38030101e10,
                                   0.0,
                                   0.0,
                                   0.0,
                                   0.0)
meltFit = NinthOrderPolynomialFit(7.40464217e10,
                                  2.49802214e11,
                                  1.00445029e12,
                                 -1.36451475e11,
                                  7.72897829e9,
                                  5.06390305e10,
                                  0.0,
                                  0.0,
                                  0.0,
                                  0.0)
strengthModel = SteinbergGuinanStrengthCGS3d(eos,
                                             7.700000e11,        # G0
                                             2.2600e-12,         # A
                                             4.5500e-04,          # B
                                             3.4000e9,           # Y0
                                             2.5e10,             # Ymax
                                             1.0e-3,             # Yp
                                             43.0000,            # beta
                                             0.0,                # gamma0
                                             0.35,               # nhard
                                             coldFit,
                                             meltFit)

# Construct another equation of state for the damaged material.
eosDamaged = GruneisenEquationOfStateCGS3d(7.90,    # reference density  
                                           etamin,  # etamin             
                                           etamax,  # etamax             
                                           0.457e6, # C0                 
                                           1.49,    # S1                 
                                           0.0,     # S2                 
                                           0.0,     # S3                 
                                           1.93,    # gamma0             
                                           0.5,     # b                  
                                           55.350,  # atomic weight
                                           0.0,     # external pressure,
                                           0.0)     # minimum pressure

#-------------------------------------------------------------------------------
# Create our interpolation kernels -- one for normal hydro interactions, and
# one for use with the artificial viscosity
#-------------------------------------------------------------------------------
WT = TableKernel3d(BSplineKernel3d(), 1000)
WTPi = TableKernel3d(BSplineKernel3d(), 1000)
output("WT")
output("WTPi")
kernelExtent = WT.kernelExtent()

#-------------------------------------------------------------------------------
# Create the NodeLists.
#-------------------------------------------------------------------------------
nodes = SolidNodeListConstructor("Stainless steel", eos, strengthModel, WT, WTPi)
nodesDamaged = FluidNodeListConstructor("Damaged stainless steel", eosDamaged, WT, WTPi)
nodeSet = [nodes, nodesDamaged]
for n in nodeSet:
    n.nodesPerSmoothingScale = nPerh
    n.epsilonTensile = epsilonTensile
    n.nTensile = nTensile
    n.XSPH = XSPH
    n.hmin = hmin
    n.hmax = hmax
    n.rhoMin = etamin*rho0
    n.rhoMax = etamax*rho0
    output("n.name()")
    output("  n.nodesPerSmoothingScale")
    output("  n.epsilonTensile")
    output("  n.nTensile")
    output("  n.XSPH")
    output("  n.hmin")
    output("  n.hmax")
    output("  n.rhoMin")
    output("  n.rhoMax")
del n

#-------------------------------------------------------------------------------
# Construct the neighbor objects and associate them with the node lists.
#-------------------------------------------------------------------------------
neighborTimer = SpheralTimer("Neighbor initialization.")
neighborTimer.start()
cache = []
for n in nodeSet:
    neighbor = TreeNeighbor3d(n,
                              kernelExtent = kernelExtent)
    n.registerNeighbor(neighbor)
    cache.append(neighbor)
neighborTimer.stop()
neighborTimer.printStatus()
del n

#-------------------------------------------------------------------------------
# Set node properties (positions, masses, H's, etc.)
#-------------------------------------------------------------------------------
if restoreCycle is None:
    print("Generating node distribution.")
    from GenerateNodeDistribution3d import *
    from ParMETISDistributeNodes import distributeNodes3d
    generator = GenerateCylindricalNodeDistribution3d(nx,
                                                      ny,
                                                      rho0,
                                                      "lattice",
                                                      xmin = xmin,
                                                      xmax = xmax,
                                                      nNodePerh = nPerh,
                                                      phi = phi,
                                                      SPH = not isinstance(nodes,
                                                                           AsphSolidNodeList3d),
                                                      )
    distributeNodes3d((nodes, generator))
    output("mpi.reduce(nodes.numInternalNodes, mpi.MIN)")
    output("mpi.reduce(nodes.numInternalNodes, mpi.MAX)")
    output("mpi.reduce(nodes.numInternalNodes, mpi.SUM)")

    # Set node velocites.
    for i in range(nodes.numInternalNodes):
        ri = nodes.positions()[i]
        vunit = Vector3d(0.0, ri.y, ri.z).unitVector()
        vr = sqrt(ri.y*ri.y + ri.z*ri.z)/radius * v0
        nodes.velocity()[i] = vunit*vr

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase3d()
for n in nodeSet:
    db.appendNodeList(n)
del n
output("db")
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Construct constant velocity boundary conditions to be applied to the outside
# of the disk.
#-------------------------------------------------------------------------------
## rNodes = vector_of_ULL()
## dr = radius/ny
## rNodes.extend([i for i in xrange(nodes.numInternalNodes)
##                if sqrt(nodes.positions()[i].y**2 +
##                        nodes.positions()[i].z**2) > radius - 2*dr])
## print "Selected %i constant velocity nodes." % mpi.allreduce(len(rNodes), mpi.SUM)
## rbc = ConstantVelocityBoundary3d(nodes, rNodes)

#-------------------------------------------------------------------------------
# Also create reflecting boundaries on the x=0 and y=0 lines.
#-------------------------------------------------------------------------------
bcs = [] # [rbc]
if phiFactor < 2.0:
    yPlane = Plane3d(Vector3d(0.0, 0.0, 0.0), Vector3d(0.0, 1.0, 0.0))
    zPlane = Plane3d(Vector3d(0.0, 0.0, 0.0), Vector3d(0.0, 0.0, 1.0))
    ybc = ReflectingBoundary3d(yPlane)
    zbc = ReflectingBoundary3d(zPlane)
    bcs.extend([ybc, zbc])

#-------------------------------------------------------------------------------
# Construct the artificial viscosities for the problem.
#-------------------------------------------------------------------------------
q = Qconstructor(Cl, Cq)
q.limiter = Qlimiter
q.balsaraShearCorrection = balsaraCorrection
q.epsilon2 = epsilon2
q.negligibleSoundSpeed = negligibleSoundSpeed
q.csMultiplier = csMultiplier
output("q")
output("q.Cl")
output("q.Cq")
output("q.limiter")
output("q.epsilon2")
output("q.negligibleSoundSpeed")
output("q.csMultiplier")
output("q.balsaraShearCorrection")

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
hydro = Hydro3d(WT, WTPi, q, compatibleEnergyEvolution)
hydro.cfl = cfl
hydro.useVelocityMagnitudeForDt = True
hydro.HEvolution = HEvolution
hydro.sumForMassDensity = sumForMassDensity
hydro.HsmoothMin = hmin
hydro.HsmoothMax = hmax
hydro.hybridMassDensityThreshold = hybridMassDensityThreshold
output("hydro")
output("hydro.cfl")
output("hydro.compatibleEnergyEvolution")
output("hydro.useVelocityMagnitudeForDt")
output("hydro.HEvolution")
output("hydro.sumForMassDensity")
output("hydro.HsmoothMin")
output("hydro.HsmoothMax")
output("hydro.kernel()")
output("hydro.PiKernel()")
output("hydro.valid()")
output("hydro.hybridMassDensityThreshold")

#-------------------------------------------------------------------------------
# Construct a strength physics object.
#-------------------------------------------------------------------------------
strength = Strength3d()
output("strength")

#-------------------------------------------------------------------------------
# Construct a damage model.
#-------------------------------------------------------------------------------
nflaw = 10
if (DamageModelConstructor is GradyKippTensorDamage3d or
    DamageModelConstructor is WeibullTensorDamage3d):
    damageModel = DamageModelConstructor(nodes,
                                         kWeibull,
                                         mWeibull,
                                         volume,
                                         WT,
                                         randomSeed,
                                         strainType,
                                         useDamageGradient,
                                         0.4,
                                         1,
                                         nflaw)
    damageModel.emulateScalarDamage = emulateScalarDamage
    output("damageModel")
    output("damageModel.useDamageGradient")
    output("damageModel.emulateScalarDamage")

elif DamageModelConstructor is GradyKippScalarDamage3d:
    damageModel = DamageModelConstructor(nodes,
                                         nodesDamaged,
                                         kWeibull,
                                         mWeibull,
                                         volume,
                                         WT,
                                         randomSeed,
                                         0.4,
                                         1,
                                         nflaw)
    output("damageModel")

if cullToWeakestFlaws:
    damageModel.cullToWeakestFlaws()
# damageModel.excludeNodes = rNodes

#-------------------------------------------------------------------------------
# Construct a predictor corrector integrator.
#-------------------------------------------------------------------------------
integrator = SynchronousRK2Integrator3d(db)
integrator.appendPhysicsPackage(hydro)
integrator.appendPhysicsPackage(strength)
integrator.appendPhysicsPackage(damageModel)
integrator.lastDt = dt
if dtMin:
    integrator.dtMin = dtMin
if dtMax:
    integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
output("integrator")
output("integrator.havePhysicsPackage(hydro)")
output("integrator.havePhysicsPackage(strength)")
output("integrator.havePhysicsPackage(damageModel)")
output("integrator.valid()")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.dtGrowth")

#-------------------------------------------------------------------------------
# Add the boundary conditions.
#-------------------------------------------------------------------------------
for package in integrator.physicsPackages():
    for bc in bcs:
        package.appendBoundary(bc)

#-------------------------------------------------------------------------------
# Build the controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            initializeMassDensity = False)
output("control")

#-------------------------------------------------------------------------------
# Monitor the evolution of the mass averaged strain.
#-------------------------------------------------------------------------------
sys.path.append("../Utilities")
from AverageStrain import AverageStrain
strainHistory = AverageStrain(damageModel,
                              dataDir + "/strainhistory.txt")
control.appendPeriodicWork(strainHistory.sample, 5)

#-------------------------------------------------------------------------------
# Helper method to dump viz files.
#-------------------------------------------------------------------------------
def viz(nodes, nodesDamaged, damageModel):
    if isinstance(damageModel, TensorDamageModel3d):
        tdamage = damageModel.damage()
        tstrain = damageModel.strain()
        sdamage = ScalarField3d("damage magnitude", nodes)
        sdamagemin = ScalarField3d("damage magnitude min", nodes)
        sdamagemax = ScalarField3d("damage magnitude max", nodes)
        sstrain = ScalarField3d("strain average", nodes)
        sstrainmin = ScalarField3d("strain min", nodes)
        sstrainmax = ScalarField3d("strain max", nodes)
        for i in range(nodes.numInternalNodes):
            sdamage[i] = tdamage[i].Trace()
            sdamagemin[i] = tdamage[i].eigenValues().minElement()
            sdamagemax[i] = tdamage[i].eigenValues().maxElement()
            sstrain[i] = tstrain[i].Trace()/2.0
            sstrainmin[i] = tstrain[i].eigenValues().minElement()
            sstrainmax[i] = tstrain[i].eigenValues().maxElement()
        dumpPhysicsState(integrator,
                         "TensileDisk-3d-visit",
                         visitDir,
                         fields = [damageModel.sumActivationEnergiesPerNode(),
                                   damageModel.numFlawsPerNode(),
                                   sdamage,
                                   sdamagemin,
                                   sdamagemax,
                                   sstrain,
                                   sstrainmin,
                                   sstrainmax,
                                   #ufragIndex
                                   ]
                         )

    else:
        raise "We need to add support for your damage model to the viz() function."

#-------------------------------------------------------------------------------
# Smooth the initial conditions/restore state.
#-------------------------------------------------------------------------------
if restoreCycle is not None:
    control.loadRestartFile(restoreCycle)
    strainHistory.flushHistory()

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
while control.time() < goalTime:
    nextGoalTime = min(control.time() + dtSample, goalTime)
    control.advance(nextGoalTime, maxSteps)
    control.dropRestartFile()
    viz(nodes, nodesDamaged, damageModel)

#-------------------------------------------------------------------------------
# Now we can collect some info on the fragment population.
#-------------------------------------------------------------------------------
from identifyFragments import *
linkRadius = 1.0
damageThreshold = 1.0
fragIDs = identifyFragments(nodes,
                            linkRadius,
                            damageModel.damage(),
                            damageThreshold,
                            True)
fragProps = fragmentProperties(nodes, fragIDs, damageModel.strain())
f = open("frag-props-t=%gusec-lr=%f-dt=%f" % (control.time(),
                                              linkRadius, damageThreshold), "w")
import pickle
pickle.dump(fragProps, f)
f.close()
