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

# Load the mpi module if we're parallel.
import loadmpi
mpi, procID, numProcs = loadmpi.loadmpi()

#-------------------------------------------------------------------------------
# Identify ourselves!
#-------------------------------------------------------------------------------
title("2-D Tensile disk strength/damage model test")

#-------------------------------------------------------------------------------
# Generic problem parameters
# All CGS units.
#-------------------------------------------------------------------------------
commandLine(SolidNodeListConstructor = AsphSolidNodeList2d,
            FluidNodeListConstructor = AsphNodeList2d,

            seed = "constantDTheta",

            radius = 10.0,
            thickness = 0.5,
            thetaFactor = 2.0,
            nr = 100,
            ntheta = 100,
            nPerh = 2.01,

            neff = None,

            rho0 = 7.9,

            # Material specific bounds on the mass density.
            etamin = 0.5,
            etamax = 1.5,

            # Parameters for the time dependent strain and cracking.
            DamageModelConstructor = GradyKippTensorDamage2d, # GradyKippScalarDamage2d, # WeibullTensorDamage2d, # GradyKippVectorDamage2d, # 
            v0 = 5e4,
            kWeibullFactor = 1.0,
            mWeibullFactor = 1.0,
            randomSeed = 548928513,
            strainType = TensorDamageModel2d.TensorStrainAlgorithm.StrainHistory,
            useDamageGradient = True,
            cullToWeakestFlaws = False,

            Qconstructor = MonaghanGingoldViscosity2d,
            Cl = 1.0,
            Cq = 1.0,
            Qlimiter = False,
            balsaraCorrection = False,
            epsilon2 = 1e-2,
            negligibleSoundSpeed = 1e-5,
            csMultiplier = 1e-4,
            hmin = 1e-5,
            hmax = 20.0,
            hminratio = 0.05,
            cfl = 0.5,
            useVelocityMagnitudeForDt = False,
            XSPH = True,
            epsilonTensile = 0.3,
            nTensile = 4,
            hybridMassDensityThreshold = 0.01,

            goalTime = 100.0e-6,
            dtSample = 1.0e-6,
            steps = None,
            dt = 1e-10,
            dtMin = 1e-9,
            dtMax = 1e-5,
            dtGrowth = 2.0,
            dtverbose = False,
            maxSteps = None,
            statsStep = 10,
            smoothIters = 0,
            HEvolution = Hydro2d.HEvolutionType.IdealH,
            sumForMassDensity = Hydro2d.MassDensityType.IntegrateDensity, # HybridDensity, # CorrectedSumDensity, # VolumeScaledDensity, #
            compatibleEnergyEvolution = True,
            gradhCorrection = True,
            postIterateHCycle = 5,
            criticalNodesPerSmoothingScale = 0.99,

            restoreCycle = None,
            restartStep = 500,

            basedir = ".",

            )

theta = thetaFactor*pi

# Assuming strain ~ 4.5e3 sec^-1, L ~ 1.4cm
kWeibull = 8.8e4 * kWeibullFactor
mWeibull = 2.63  * mWeibullFactor

dataDir = ("%s/dumps-TensileDisk-2d/%ix%i/nPerh=%4.2f/%s/%s/XSPH=%s/%s/gradDamage=%s/k=%4.2f_m=%4.2f" %
           (basedir,
            nr, ntheta,
            nPerh,
            str(DamageModelConstructor).split("'")[1],
            str(SolidNodeListConstructor).split("'")[1],
            XSPH,
            str(strainType),
            useDamageGradient,
            kWeibull,
            mWeibull))
restartDir = dataDir + "/restarts/proc-%04i" % procID
visitDir = dataDir + "/visit"
restartBaseName = restartDir + "/TensileDisk-%ix%i" % (nr, ntheta)

xmin = (0.0, 0.0)
xmax = (thickness, radius)

volume = pi*radius**2 * thickness

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
eos = GruneisenEquationOfStateCGS2d(7.90,    # reference density  
                                    etamin,  # etamin             
                                    etamax,  # etamax             
                                    0.457e6, # C0                 
                                    1.49,    # S1                 
                                    0.0,     # S2                 
                                    0.0,     # S3                 
                                    1.93,    # gamma0             
                                    0.5,     # b                  
                                    55.350,  # atomic weight
                                    0.0,     # external pressure
                                   -1.0e50,  # minimum pressure
                                    1.0e11)  # maximum pressure
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
strengthModel = SteinbergGuinanStrengthCGS2d(eos,
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
eosDamaged = GruneisenEquationOfStateCGS2d(7.90,    # reference density  
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
WT = TableKernel2d(BSplineKernel2d(), 1000)
output("WT")
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
    n.hminratio = hminratio
    n.rhoMin = etamin*rho0
    n.rhoMax = etamax*rho0
    output("n.name()")
    output("  n.nodesPerSmoothingScale")
    output("  n.epsilonTensile")
    output("  n.nTensile")
    output("  n.XSPH")
    output("  n.hmin")
    output("  n.hmax")
    output("  n.hminratio")
    output("  n.rhoMin")
    output("  n.rhoMax")
del n

#-------------------------------------------------------------------------------
# Construct the neighbor objects and associate them with the node lists.
#-------------------------------------------------------------------------------
cache = []
for n in nodeSet:
    neighbor = TreeNeighbor2d(n,
                              kernelExtent = kernelExtent)
    n.registerNeighbor(neighbor)
    cache.append(neighbor)
del n

#-------------------------------------------------------------------------------
# Set node properties (positions, masses, H's, etc.)
#-------------------------------------------------------------------------------
eps0 = 0.0
if restoreCycle is None:
    print("Generating node distribution.")
    from GenerateNodeDistribution2d import *
    from ParMETISDistributeNodes import distributeNodes2d
    generator = GenerateNodeDistribution2d(nr,
                                           ntheta,
                                           rho0,
                                           seed,
                                           rmin = 0.0,
                                           rmax = radius,
                                           theta = theta,
                                           nNodePerh = nPerh,
                                           SPH = not isinstance(nodes,
                                                                AsphSolidNodeList2d),
                                           )
    distributeNodes2d((nodes, generator))
    output("mpi.reduce(nodes.numInternalNodes, mpi.MIN)")
    output("mpi.reduce(nodes.numInternalNodes, mpi.MAX)")
    output("mpi.reduce(nodes.numInternalNodes, mpi.SUM)")

    # Set node specific thermal energies
    eps0 = eos.specificThermalEnergy(rho0, 300.0)
    nodes.specificThermalEnergy(ScalarField2d("tmp", nodes, eps0))

    # Set node velocites.
    for i in range(nodes.numInternalNodes):
        vr = nodes.positions()[i].magnitude()/radius*v0
        vunit = nodes.positions()[i].unitVector()
        nodes.velocity()[i] = vunit*vr

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase2d()
for n in nodeSet:
    db.appendNodeList(n)
del n
output("db")
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Create reflecting boundaries on the x=0 and y=0 lines.
#-------------------------------------------------------------------------------
bcs = []
if theta == 0.5*pi:
    xPlane = Plane2d(Vector2d(0.0, 0.0), Vector2d(1.0, 0.0))
    yPlane = Plane2d(Vector2d(0.0, 0.0), Vector2d(0.0, 1.0))
    xbc = ReflectingBoundary2d(xPlane)
    ybc = ReflectingBoundary2d(yPlane)
    bcs += [xbc, ybc]

#-------------------------------------------------------------------------------
# Construct a constant velocity boundary conditions to be applied to the outside
# of the disk.
#-------------------------------------------------------------------------------
rNodes = vector_of_ULL()
dr = radius/nr
rNodes.extend([i for i in range(nodes.numInternalNodes)
               if nodes.positions()[i].magnitude() > radius - 4*dr])
print("Selected %i constant velocity nodes." % mpi.allreduce(len(rNodes), mpi.SUM))
rbc = ConstantVelocityBoundary2d(nodes, rNodes)
# bcs += [rbc]

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
hydro = Hydro2d(W = WT,
                Q = q,
                compatibleEnergyEvolution = compatibleEnergyEvolution,
                gradhCorrection = gradhCorrection,
                densityUpdate = sumForMassDensity,
                HUpdate = HEvolution)
hydro.cfl = cfl
hydro.useVelocityMagnitudeForDt = useVelocityMagnitudeForDt
output("hydro")
output("hydro.cfl")
output("hydro.compatibleEnergyEvolution")
output("hydro.gradhCorrection")
output("hydro.useVelocityMagnitudeForDt")
output("hydro.HEvolution")
output("hydro.sumForMassDensity")
output("hydro.kernel()")
output("hydro.PiKernel()")
output("hydro.postIterateHCycle")
output("hydro.valid()")

#-------------------------------------------------------------------------------
# Construct a strength physics object.
#-------------------------------------------------------------------------------
strength = Strength2d(criticalNodesPerSmoothingScale)
output("strength")
output("strength.criticalNodesPerSmoothingScale")

#-------------------------------------------------------------------------------
# Construct a damage model.
#-------------------------------------------------------------------------------
if DamageModelConstructor is StochasticWeibullTensorDamageModel2d:
    damageModel = DamageModelConstructor(nodes,
                                         kWeibull,
                                         mWeibull,
                                         volume,
                                         WT,
                                         randomSeed,
                                         strainType,
                                         useDamageGradient,
                                         0.4,
                                         thickness,
                                         int(thickness/dy + 0.5))
    output("damageModel")
    output("damageModel.useDamageGradient")
    output("damageModel.nodalVolumeMultiplier")
    output("damageModel.maxNumFlawsPerNode")

else:
    # Compute the minimum number of flaws we want to seed.
    nthick = int(thickness/(radius/nr) + 0.5)
    if neff is None:
        neff = max(1, mpi.allreduce(nodes.numInternalNodes, mpi.SUM)*nthick)
    nflaws = int(neff * log(neff))
    print("nthick, nflaws = ", nthick, nflaws)
    damageModel = DamageModelConstructor(nodes,
                                         kWeibull,
                                         mWeibull,
                                         volume,
                                         WT,
                                         randomSeed,
                                         strainType,
                                         useDamageGradient,
                                         0.4,
                                         nthick,
                                         nflaws)
    output("damageModel")
    output("damageModel.useDamageGradient")

    if cullToWeakestFlaws:
        damageModel.cullToWeakestFlaws()

# damageModel.excludeNodes = rNodes

#-------------------------------------------------------------------------------
# Construct a time integrator.
#-------------------------------------------------------------------------------
integrator = SynchronousRK2Integrator2d(db)
integrator.appendPhysicsPackage(hydro)
integrator.appendPhysicsPackage(strength)
integrator.appendPhysicsPackage(damageModel)
integrator.lastDt = dt
if dtMin:
    integrator.dtMin = dtMin
if dtMax:
    integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
integrator.verbose = dtverbose
output("integrator")
output("integrator.havePhysicsPackage(hydro)")
output("integrator.havePhysicsPackage(strength)")
output("integrator.havePhysicsPackage(damageModel)")
output("integrator.valid()")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.dtGrowth")
output("integrator.verbose")

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
    if isinstance(damageModel, TensorDamageModel2d):
        tdamage = nodes.damage()
        tstrain = damageModel.effectiveStrain()
        vr = ScalarField2d("radial velocity", nodes)
        sdamage = ScalarField2d("damage magnitude", nodes)
        sdamagemin = ScalarField2d("damage magnitude min", nodes)
        sdamagemax = ScalarField2d("damage magnitude max", nodes)
        sstrain = ScalarField2d("strain average", nodes)
        sstrainmin = ScalarField2d("strain min", nodes)
        sstrainmax = ScalarField2d("strain max", nodes)
        for i in range(nodes.numInternalNodes):
            runit = nodes.positions()[i].unitVector()
            vr[i] = nodes.velocity()[i].dot(runit)
            sdamage[i] = tdamage[i].Trace()
            sdamagemin[i] = tdamage[i].eigenValues().minElement()
            sdamagemax[i] = tdamage[i].eigenValues().maxElement()
            sstrain[i] = tstrain[i].Trace()/2.0
            sstrainmin[i] = tstrain[i].eigenValues().minElement()
            sstrainmax[i] = tstrain[i].eigenValues().maxElement()
        dumpPhysicsState(integrator,
                         "TensileDisk-2d-visit",
                         visitDir,
                         fields = [vr,
                                   damageModel.sumActivationEnergiesPerNode(),
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
if not steps is None:
    control.step(steps)
else:
    while control.time() < goalTime:
        nextGoalTime = min(control.time() + dtSample, goalTime)
        control.advance(nextGoalTime, maxSteps)
        control.dropRestartFile()
        viz(nodes, nodesDamaged, damageModel)

    #---------------------------------------------------------------------------
    # Now we can collect some info on the fragment population.
    #---------------------------------------------------------------------------
    from identifyFragments import *
    linkRadius = 1.25
    damageThreshold = 1.0
    fragIDs = identifyFragments(nodes,
                                linkRadius,
                                nodes.damage(),
                                damageThreshold,
                                True)
    fragProps = fragmentProperties(nodes, fragIDs, damageModel.strain())
    f = open("frag-props-t=%gusec-lr=%f-dt=%f" % (control.time(),
                                                  linkRadius, damageThreshold), "w")
    import pickle
    pickle.dump(fragProps, f)
    f.close()
