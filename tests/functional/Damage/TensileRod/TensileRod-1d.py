# Grady-Kipp-Owen damage 
#ATS:t10 = test(SELF, "--DamageModelConstructor GradyKippTensorDamageOwen --graphics False --clearDirectories True --domainIndependent True --outputFile 'TensileRod-1d-1proc.gnu' --checkRef True", np=1, label="Tensile rod (GradyKippOwen damage) domain independence test SERIAL RUN")
#ATS:t11 = testif(t10, SELF, "--DamageModelConstructor GradyKippTensorDamageOwen --graphics False --clearDirectories False --domainIndependent True --outputFile 'TensileRod-1d-4proc.gnu' --comparisonFile 'TensileRod-1d-1proc.gnu'", np=4, label="Tensile rod (GradyKippOwen damage) domain independence test 4 DOMAIN RUN")
#ATS:t12 = testif(t11, SELF, "--DamageModelConstructor GradyKippTensorDamageOwen --graphics False --clearDirectories False --domainIndependent True --outputFile 'TensileRod-1d-1proc-restart.gnu' --comparisonFile 'TensileRod-1d-1proc.gnu' --restoreCycle 500", np=1, label="Tensile rod (GradyKippOwen damage) domain independence test SERIAL RESTART RUN")
#ATS:t13 = testif(t11, SELF, "--DamageModelConstructor GradyKippTensorDamageOwen --graphics False --clearDirectories False --domainIndependent True --outputFile 'TensileRod-1d-4proc-restart.gnu' --comparisonFile 'TensileRod-1d-1proc.gnu' --restoreCycle 500", np=4, label="Tensile rod (GradyKippOwen damage) domain independence test 4 DOMAIN RESTART RUN")
#
# Probabilistic damage
#ATS:t20 = test(SELF, "--DamageModelConstructor ProbabilisticDamageModel --graphics False --clearDirectories True --domainIndependent True --outputFile 'TensileRod-1d-1proc.gnu' --checkRef True", np=1, label="Tensile rod (probabilistic damage) domain independence test SERIAL RUN")
#ATS:t21 = testif(t20, SELF, "--DamageModelConstructor ProbabilisticDamageModel --graphics False --clearDirectories False --domainIndependent True --outputFile 'TensileRod-1d-4proc.gnu' --comparisonFile 'TensileRod-1d-1proc.gnu'", np=4, label="Tensile rod (probabilistic damage) domain independence test 4 DOMAIN RUN")
#ATS:t22 = testif(t21, SELF, "--DamageModelConstructor ProbabilisticDamageModel --graphics False --clearDirectories False --domainIndependent True --outputFile 'TensileRod-1d-1proc-restart.gnu' --comparisonFile 'TensileRod-1d-1proc.gnu' --restoreCycle 500", np=1, label="Tensile rod (probabilistic damage) domain independence test SERIAL RESTART RUN")
#ATS:t23 = testif(t21, SELF, "--DamageModelConstructor ProbabilisticDamageModel --graphics False --clearDirectories False --domainIndependent True --outputFile 'TensileRod-1d-4proc-restart.gnu' --comparisonFile 'TensileRod-1d-1proc.gnu' --restoreCycle 500", np=4, label="Tensile rod (probabilistic damage) domain independence test 4 DOMAIN RESTART RUN")

#-------------------------------------------------------------------------------
# A rod of stainless steel undergoing tensile strain.  This is intended as a
# test of the cracking/failure models.
#
# See Benz & Asphaug (1994), Icarus, 107, 98
#-------------------------------------------------------------------------------
from SolidSpheral1d import *
from SpheralTestUtilities import *
from identifyFragments import identifyFragments, fragmentProperties
from math import *
import os, shutil
import mpi

#-------------------------------------------------------------------------------
# Identify ourselves!
#-------------------------------------------------------------------------------
title("1-D Tensile rod strength/damage model test")

#-------------------------------------------------------------------------------
# Stupid little class to override the mass density evolution of the control
# boundary nodes.
#-------------------------------------------------------------------------------
class OverrideNodeProperties:
    def __init__(self,
                 nodeList,
                 rho0,
                 eps0,
                 controlNodeIDs):
        self.restart = RestartableObject(self)
        self.nodeList = nodeList
        self.rho0 = rho0
        self.eps0 = eps0
        self.controlNodeFlags = ScalarField("flag nodes",
                                              nodeList,
                                              0.0)
        for i in controlNodeIDs:
            assert i < nodeList.numInternalNodes
            self.controlNodeFlags[i] = 1.0
        return

    def override(self, cycle, t, dt):
        ids = [i for i in range(self.nodeList.numInternalNodes)
               if self.controlNodeFlags[i] > 0.5]
        for i in ids:
            self.nodeList.massDensity()[i] = self.rho0
            self.nodeList.specificThermalEnergy()[i] = self.eps0
        return

    def label(self):
        return "OverrideNodeProperties"

    def dumpState(self, file, path):
        file.writeObject(self.rho0, path + "/rho0")
        file.writeObject(self.eps0, path + "/eps0")
        file.write(self.controlNodeFlags, path + "/controlNodeFlags")
        return

    def restoreState(self, file, path):
        self.rho0 = file.readObject(path + "/rho0")
        self.eps0 = file.readObject(path + "/eps0")
        file.read(self.controlNodeFlags, path + "/controlNodeFlags")
        return

#-------------------------------------------------------------------------------
# Build our units.
#-------------------------------------------------------------------------------
units = PhysicalConstants(0.01,  # Unit length in m
                          0.001, # Unit mass in kg
                          1e-6)  # Unit length in sec

#-------------------------------------------------------------------------------
# Generic problem parameters
# All (cm, gm, usec) units.
#-------------------------------------------------------------------------------
commandLine(length = 3.0,
            radius = 0.5,
            nx = 100,

            rho0 = 7.9,

            # Material specific bounds on the mass density.
            etamin = 0.5,
            etamax = 1.5,

            # Parameters for the time dependent strain and cracking.
            useDamage = True,
            DamageModelConstructor = ProbabilisticDamageModel,
            volumeMultiplier = (3.0/100.0)**2,
            numFlawsPerNode = 1,
            v0 = 1e-2,
            kWeibullFactor = 1.0,
            mWeibullFactor = 1.0,
            randomSeed = 548928513,
            strainType = PseudoPlasticStrain,
            damageCoupling = PairMaxDamage,
            cullToWeakestFlaws = False,
            damageInCompression = False,

            # Johnson-Cook choices
            D1 = 0.0,
            D2 = 2.0,
            D3 = -1.5,
            D4 = 0.0,
            D5 = 0.0,
            epsilondot0 = 0.01,
            Tcrit = -2.0,
            sigmamax = -3.0,
            efailmin = 0.01,
            sigmaD1 = 0.0,  # Gaussian
            aD1 = 0.0,      # Weibull
            bD1 = 0.0,      # Weibull
            eps0D1 = 0.0,   # Weibull
            sigmaD2 = 0.05, # Gaussian
            aD2 = 0.065,    # Weibull
            bD2 = 2.0,      # Weibull
            eps0D2 = 0.165, # Weibull

            # Optionally we can initialize a break near the origin.
            initialBreakRadius = 0.0,
            
            hydroType = "SPH",                 # one of (SPH, CRKSPH, FSISPH)
            hmin = 1e-5,
            hmax = 1.0,
            cfl = 0.5,
            useVelocityMagnitudeForDt = False,
            XSPH = False,
            epsilonTensile = 0.3,
            nTensile = 4,
            hybridMassDensityThreshold = 0.01,
            filter = 0.0,
            volumeType = RKSumVolume,

            IntegratorConstructor = CheapSynchronousRK2Integrator,
            goalTime = 50.0,
            steps = None,
            dt = 1e-10,
            dtMin = 1e-6,
            dtMax = 10.0,
            dtGrowth = 2.0,
            dumpFrac = 0.005,
            maxSteps = None,
            statsStep = 10,
            smoothIters = 0,
            HUpdate = IdealH,
            densityUpdate = IntegrateDensity,
            compatibleEnergy = True,
            gradhCorrection = True,
            correctVelocityGradient = True,
            domainIndependent = True,
            dtverbose = False,

            restoreCycle = -1,
            restartStep = 500,

            graphics = True,
            vizCycle = None,
            vizTime = 10.0,

            testtol = 1.0e-4,
            clearDirectories = False,
            dataDirBase = "dumps-TensileRod-1d",
            checkRef = False,
            outputFile = None,
            comparisonFile = None,
            )

# On the IBM BlueOS machines we have some tolerance issues...
if "SYS_TYPE" in os.environ and os.environ["SYS_TYPE"] == "blueos_3_ppc64le_ib_p9":
    testtol *= 20.0

hydroType = hydroType.upper()
assert hydroType in ("SPH", "CRKSPH", "FSISPH")

nPerh = 1.51
order = 5

if DamageModelConstructor in (GradyKippTensorDamage, GradyKippTensorDamageOwen, ProbabilisticDamageModel):
    damageName = os.path.join(str(DamageModelConstructor.__name__), str(damageCoupling))
else:
    damageName = DamageModelConstructor.__name__

#kWeibull = 8.8e4 * kWeibullFactor
#kWeibull = 6.52e3 * kWeibullFactor
kWeibull = 6.52e5 * kWeibullFactor
mWeibull = 2.63   * mWeibullFactor

dataDir = os.path.join(dataDirBase,
                       hydroType,
                       damageName,
                       "nx=%i" % nx,
                       "k=%4.2f_m=%4.2f" % (kWeibull, mWeibull))
refDir = os.path.join("Reference",
                       hydroType,
                       damageName,
                       "nx=%i" % nx,
                       "k=%4.2f_m=%4.2f" % (kWeibull, mWeibull))

restartDir = os.path.join(dataDir, "restarts")
restartBaseName = os.path.join(restartDir, "TensileRod-%i" % nx)
vizDir = os.path.join(dataDir, "viz")
vizBaseName = "TensileRod-%i" % nx

xmin = -0.5*length
xmax =  0.5*length

volume = pi*radius**2 * length
dx = length/nx

dtSample = dumpFrac*goalTime

#-------------------------------------------------------------------------------
# Sampling function to measure the average strain in the volume of the rod.
#-------------------------------------------------------------------------------
class AverageStrain:
    def __init__(self, damageModel, filename):
        self.restart = RestartableObject(self)
        self.damageModel = damageModel
        self.filename = filename
        self.timeHistory = []
        self.strainHistory = []

        self.file = None
        if mpi.rank == 0:
            self.file = open(self.filename, "w")
            assert self.file is not None

        return

    def sample(self, cycle, atime, dt):
        nodes = self.damageModel.nodeList
        mass = nodes.mass()
        strain = self.damageModel.strain

        n = nodes.numInternalNodes
        result = (mpi.allreduce(sum([mass[i]*(strain[i].Trace()) for i in range(n)]), mpi.SUM)/
                  mpi.allreduce(sum(mass.internalValues()), mpi.SUM))

        self.timeHistory.append(atime)
        self.strainHistory.append(result)

        if mpi.rank == 0:
            self.file.write("%g   %g\n" % (atime, result))
            self.file.flush()

##         # Check for out of control nodes!
##         vel = [v.x for v in nodes.velocity().internalValues()]
##         for i in xrange(nodes.numInternalNodes):
##             if abs(vel[i]) > 1e6:
##                 print "Bad node: ", i, " vel = ", vel[i]

        return

    def flushHistory(self):
        if mpi.rank == 0:
            assert self.file is not None
            n = len(self.timeHistory)
            assert len(self.strainHistory) == n
            if mpi.rank == 0:
                for i in range(n):
                    self.file.write("%g   %g\n" % (self.timeHistory[i],
                                                   self.strainHistory[i]))
            self.file.flush()

    def label(self):
        return "AverageStrain"

    def dumpState(self, file, path):
        file.writeObject(self.timeHistory, path + "/timeHistory")
        file.writeObject(self.strainHistory, path + "/strainHistory")
        return

    def restoreState(self, file, path):
        self.timeHistory = file.readObject(path + "/timeHistory")
        self.strainHistory = file.readObject(path + "/strainHistory")
        return

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
# Stainless steel material properties.
#-------------------------------------------------------------------------------
eos = GruneisenEquationOfState(rho0,    # reference density  
                               etamin,  # etamin             
                               etamax,  # etamax             
                               0.457, # C0                 
                               1.49,    # S1                 
                               0.0,     # S2                 
                               0.0,     # S3                 
                               1.93,    # gamma0             
                               0.5,     # b                  
                               55.350,  # atomic weight
                               units)
coldFit = NinthOrderPolynomialFit(-1.06797724e-2,
                                  -2.06872020e-2,
                                   8.24893246e-1,
                                  -2.39505843e-2,
                                  -2.44522017e-2,
                                   5.38030101e-2,
                                   0.0,
                                   0.0,
                                   0.0,
                                   0.0)
meltFit = NinthOrderPolynomialFit(7.40464217e-2,
                                  2.49802214e-2,
                                  1.00445029e-2,
                                 -1.36451475e-1,
                                  7.72897829e-3,
                                  5.06390305e-2,
                                  0.0,
                                  0.0,
                                  0.0,
                                  0.0)
strengthModel = SteinbergGuinanStrength(eos,
                                        7.700000e-1,        # G0
                                        2.2600,             # A
                                        4.5500e-04,         # B
                                        3.4000e-3,          # Y0
                                        2.5e-2,             # Ymax
                                        1.0e-3,             # Yp
                                        43.0000,            # beta
                                        0.0,                # gamma0
                                        0.35,               # nhard
                                        coldFit,
                                        meltFit)

#-------------------------------------------------------------------------------
# Create our interpolation kernels -- one for normal hydro interactions, and
# one for use with the artificial viscosity
#-------------------------------------------------------------------------------
WT = TableKernel(NBSplineKernel(order), 1000)
output("WT")

#-------------------------------------------------------------------------------
# Create the NodeLists.
#-------------------------------------------------------------------------------
nodes = makeSolidNodeList("Stainless steel", eos, strengthModel,
                          nPerh = nPerh,
                          hmin = hmin,
                          hmax = hmax,
                          rhoMin = etamin*rho0,
                          rhoMax = etamax*rho0,
                          xmin = -100.0*Vector.one,
                          xmax =  100.0*Vector.one)
nodeSet = [nodes]

#-------------------------------------------------------------------------------
# Set node properties (positions, masses, H's, etc.)
#-------------------------------------------------------------------------------
eps0 = 0.0
print("Generating node distribution.")
from DistributeNodes import distributeNodesInRange1d
distributeNodesInRange1d([(nodes, nx, rho0, (xmin, xmax))])
output("mpi.reduce(nodes.numInternalNodes, mpi.MIN)")
output("mpi.reduce(nodes.numInternalNodes, mpi.MAX)")
output("mpi.reduce(nodes.numInternalNodes, mpi.SUM)")

# # Set node specific thermal energies
# eps0 = eos.specificThermalEnergy(rho0, 300.0)
# nodes.specificThermalEnergy(ScalarField("tmp", nodes, eps0))

# Set node velocites.
pos = nodes.positions()
vel = nodes.velocity()
for i in range(nodes.numInternalNodes):
    vel[i].x = pos[i].x/(0.5*length)*v0

# Set an initial damage if requested.
if initialBreakRadius > 0.0:
    D = nodes.damage()
    fragIDs = nodes.fragmentIDs()
    for i in range(nodes.numInternalNodes):
        if abs(pos[i].x) < initialBreakRadius:
            D[i] = SymTensor.one
        if pos[i].x < 0.0:
            fragIDs[i] = 1
        else:
            fragIDs[i] = 2

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase()
for n in nodeSet:
    db.appendNodeList(n)
del n
output("db")
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Construct constant velocity boundary conditions to be applied to the rod ends.
#-------------------------------------------------------------------------------
x0Nodes = vector_of_ULL()
x1Nodes = vector_of_ULL()
dummy = [x0Nodes.append(i) for i in range(nodes.numInternalNodes) if pos[i].x < -0.5*length + 5*dx]
dummy = [x1Nodes.append(i) for i in range(nodes.numInternalNodes) if pos[i].x >  0.5*length - 5*dx]
print("Selected %i constant velocity nodes." % (mpi.allreduce(len(x0Nodes) + len(x1Nodes), mpi.SUM)))

# Set the nodes we're going to control to one single velocity at each end.
v0 = mpi.allreduce(max([vel[i].x for i in x0Nodes] + [-100.0]), mpi.MAX)
v1 = mpi.allreduce(min([vel[i].x for i in x1Nodes] + [ 100.0]), mpi.MIN)
for i in x0Nodes:
    vel[i].x = v0
for i in x1Nodes:
    vel[i].x = v1

xbc0 = ConstantVelocityBoundary(nodes, x0Nodes)
xbc1 = ConstantVelocityBoundary(nodes, x1Nodes)

bcs = [xbc0, xbc1]

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
if hydroType == "CRKSPH":
    hydro = CRKSPH(dataBase = db,
                   W = WT,
                   cfl = cfl,
                   useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                   compatibleEnergyEvolution = compatibleEnergy,
                   XSPH = XSPH,
                   densityUpdate = densityUpdate,
                   HUpdate = HUpdate)

elif hydroType == "FSISPH":
    hydro = FSISPH(dataBase = db,
                   W = WT,
                   cfl = cfl,
                   interfaceMethod = HLLCInterface,
                   #decoupleDamagedMaterial = True,
                   densityStabilizationCoefficient = 0.00,
                   useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                   compatibleEnergyEvolution = compatibleEnergy,
                   linearCorrectGradients = correctVelocityGradient,
                   HUpdate = HUpdate)
else:
    assert hydroType == "SPH"
    hydro = SPH(dataBase = db,
                W = WT,
                cfl = cfl,
                useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                compatibleEnergyEvolution = compatibleEnergy,
                gradhCorrection = gradhCorrection,
                XSPH = XSPH,
                correctVelocityGradient = correctVelocityGradient,
                densityUpdate = densityUpdate,
                epsTensile = epsilonTensile,
                nTensile = nTensile,
                HUpdate = HUpdate)

output("hydro")
output("hydro.cfl")
output("hydro.useVelocityMagnitudeForDt")
output("hydro._smoothingScaleMethod.HEvolution")
output("hydro.densityUpdate")
output("hydro.compatibleEnergyEvolution")

packages = [hydro]

#-------------------------------------------------------------------------------
# Construct a damage model.
#-------------------------------------------------------------------------------
vizFields = []
if DamageModelConstructor is GradyKippTensorDamage:
    damageModel = DamageModelConstructor(nodeList = nodes,
                                         kWeibull = kWeibull,
                                         mWeibull = mWeibull,
                                         volume = volume,
                                         kernel = WT,
                                         seed = randomSeed,
                                         strainAlgorithm = strainType,
                                         damageCouplingAlgorithm = damageCoupling,
                                         damageInCompression = damageInCompression)

elif DamageModelConstructor is GradyKippTensorDamageOwen:
    damageModel = DamageModelConstructor(nodeList = nodes,
                                         kWeibull = kWeibull,
                                         mWeibull = mWeibull,
                                         kernel = WT,
                                         seed = randomSeed,
                                         volumeMultiplier = volumeMultiplier,
                                         strainAlgorithm = strainType,
                                         damageCouplingAlgorithm = damageCoupling,
                                         damageInCompression = damageInCompression)

elif DamageModelConstructor is JohnsonCookDamageWeibull:
    damageModel = DamageModelConstructor(nodes,
                                         D1 = D1,
                                         D2 = D2,
                                         D3 = D3,
                                         D4 = D4,
                                         D5 = D5,
                                         epsilondot0 = epsilondot0,
                                         Tcrit = Tcrit,
                                         sigmamax = sigmamax,
                                         efailmin = efailmin,
                                         aD1 = aD1,
                                         bD1 = bD1,
                                         eps0D1 = eps0D1,
                                         aD2 = aD2,
                                         bD2 = bD2,
                                         eps0D2 = eps0D2,
                                         seed = randomSeed,
                                         domainIndependent = domainIndependent)
    vizFields += [damageModel.D1(), damageModel.D2()]

elif DamageModelConstructor is JohnsonCookDamageGaussian:
    damageModel = DamageModelConstructor(nodes,
                                         D1 = D1,
                                         D2 = D2,
                                         D3 = D3,
                                         D4 = D4,
                                         D5 = D5,
                                         epsilondot0 = epsilondot0,
                                         Tcrit = Tcrit,
                                         sigmamax = sigmamax,
                                         efailmin = efailmin,
                                         sigmaD1 = sigmaD1,
                                         sigmaD2 = sigmaD2,
                                         seed = randomSeed,
                                         domainIndependent = domainIndependent)

elif DamageModelConstructor is ProbabilisticDamageModel:
    damageModel = DamageModelConstructor(nodeList = nodes,
                                         kernel = WT,
                                         kWeibull = kWeibull,
                                         mWeibull = mWeibull,
                                         seed = randomSeed,
                                         volumeMultiplier = volumeMultiplier,
                                         strainAlgorithm = strainType,
                                         damageCouplingAlgorithm = damageCoupling,
                                         damageInCompression = damageInCompression)

output("damageModel")
output("damageModel.strainAlgorithm")
output("damageModel.damageCouplingAlgorithm")

if DamageModelConstructor in (GradyKippTensorDamage, GradyKippTensorDamageOwen):
    if cullToWeakestFlaws:
        damageModel.cullToWeakestFlaws()

# damageModel.excludeNodes = xNodes

if useDamage:
    packages.append(damageModel)

#-------------------------------------------------------------------------------
# Construct a time integrator.
#-------------------------------------------------------------------------------
integrator = IntegratorConstructor(db, packages)
integrator.lastDt = dt
if dtMin:
    integrator.dtMin = dtMin
if dtMax:
    integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
integrator.domainDecompositionIndependent = domainIndependent
integrator.verbose = dtverbose
output("integrator")
output("integrator.havePhysicsPackage(hydro)")
output("integrator.havePhysicsPackage(damageModel)")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.dtGrowth")
output("integrator.domainDecompositionIndependent")

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
                            volumeType = volumeType,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle,
                            vizBaseName = vizBaseName,
                            vizDir = vizDir,
                            vizStep = vizCycle,
                            vizTime = vizTime,
                            vizFields = vizFields)
output("control")

#-------------------------------------------------------------------------------
# Since we are creating a forced velocity gradient on the control nodes, we have
# to override their mass density evolution or they falsely wander away from the
# reference density.
#-------------------------------------------------------------------------------
## override = OverrideNodeProperties(nodes,
##                                   rho0,
##                                   eps0,
##                                   xNodes)
## control.appendPeriodicWork(override.override, 1)

#-------------------------------------------------------------------------------
# Monitor the evolution of the mass averaged strain.
#-------------------------------------------------------------------------------
if DamageModelConstructor in (GradyKippTensorDamageBenzAsphaug,
                              GradyKippTensorDamageOwen,
                              ProbabilisticDamageModel):
    strainHistory = AverageStrain(damageModel,
                                  os.path.join(dataDir, "strainhistory.gnu"))
    control.appendPeriodicWork(strainHistory.sample, 1)

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
if not steps is None:
    control.step(steps)
else:
    control.advance(goalTime)
    control.dropRestartFile()

#-------------------------------------------------------------------------------
# Plot the state.
#-------------------------------------------------------------------------------
if graphics:
    from SpheralMatplotlib import *
    state = State(db, integrator.physicsPackages())
    H = state.symTensorFields("H")
    h = db.newFluidScalarFieldList(0.0, "h")
    for i in range(nodes.numInternalNodes):
        h[0][i] = 1.0/H[0][i].xx
    rhoPlot = plotFieldList(state.scalarFields(HydroFieldNames.massDensity),
                            plotStyle="o-",
                            winTitle="rho @ %g %i" % (control.time(), mpi.procs))
    velPlot = plotFieldList(state.vectorFields(HydroFieldNames.velocity),
                            yFunction = "%s.x",
                            plotStyle="o-",
                            winTitle="vel @ %g %i" % (control.time(), mpi.procs))
    mPlot = plotFieldList(state.scalarFields(HydroFieldNames.mass),
                          plotStyle="o-",
                          winTitle="mass @ %g %i" % (control.time(), mpi.procs))
    PPlot = plotFieldList(state.scalarFields(HydroFieldNames.pressure),
                          plotStyle="o-",
                          winTitle="pressure @ %g %i" % (control.time(), mpi.procs))
    uPlot = plotFieldList(state.scalarFields(HydroFieldNames.specificThermalEnergy),
                          plotStyle="o-",
                          winTitle="specific thermal energy @ %g %i" % (control.time(), mpi.procs))
    SPlot = plotFieldList(state.symTensorFields(SolidFieldNames.deviatoricStress),
                          yFunction = "%s.xx",
                          plotStyle="o-",
                          winTitle="$S_{xx}$ @ %g %i" % (control.time(), mpi.procs))
    hPlot = plotFieldList(h,
                          plotStyle="o-",
                          winTitle="h @ %g %i" % (control.time(), mpi.procs))

    d = state.symTensorFields(SolidFieldNames.tensorDamage)
    dPlot = plotFieldList(d,
                          yFunction = "%s.xx",
                          winTitle="damage @ %g %i" % (control.time(), mpi.procs),
                          plotStyle="o-")
    plots = [(rhoPlot, "rho.png"),
             (velPlot, "vel.png"),
             (mPlot, "mass.png"),
             (PPlot, "pressure.png"),
             (SPlot, "devstress.png"),
             (uPlot, "u.png"),
             (hPlot, "h.png"),
             (dPlot, "damage.png")]

    if hydroType == "FSISPH":
        DPPlot = plotFieldList(hydro.damagedPressure,
                               plotStyle="o-",
                               winTitle="damaged pressure @ %g %i" % (control.time(), mpi.procs))
        plots += [(DPPlot, "pressure_damaged.png")]

    if DamageModelConstructor in (GradyKippTensorDamage, GradyKippTensorDamageOwen):
        ts = damageModel.strain
        s = ScalarField("strain", nodes)
        for i in range(nodes.numInternalNodes):
            s[i] = ts[i].xx
        sl = ScalarFieldList()
        sl.appendField(s)
        sPlot = plotFieldList(sl, winTitle="strain @ %g %i" % (control.time(), mpi.procs),
                              plotStyle="o-")
        eps = damageModel.sumActivationEnergiesPerNode
        nflaws = damageModel.numFlawsPerNode
        for i in range(nodes.numInternalNodes):
            assert nflaws[i] > 0
            eps[i] /= nflaws[i]
        epsl = ScalarFieldList()
        epsl.appendField(eps)
        epsPlot = plotFieldList(epsl, winTitle="Flaw activation strains",
                                plotStyle="o-")
      
        plots += [(sPlot, "strain.png"),
                  (epsPlot, "flaws.png")]

    elif DamageModelConstructor in (JohnsonCookDamageWeibull, JohnsonCookDamageGaussian):
        eps = damageModel.failureStrain
        epsl = ScalarFieldList()
        epsl.appendField(eps)
        epsPlot = plotFieldList(epsl, winTitle="JC failure strains",
                                plotStyle="o-")
        D1 = damageModel.D1
        D1l = ScalarFieldList()
        D1l.appendField(D1)
        D1Plot = plotFieldList(D1l, winTitle="JC D1",
                                plotStyle="o-")
        D2 = damageModel.D2
        D2l = ScalarFieldList()
        D2l.appendField(D2)
        D2Plot = plotFieldList(D2l, winTitle="JC D2",
                                plotStyle="o-")
        plots += [(epsPlot, "JC_flaw.png"),
                  (D1Plot, "D1.png"),
                  (D2Plot, "D2.png")]

    elif DamageModelConstructor is ProbabilisticDamageModel:
        ts = damageModel.strain
        s = ScalarField("strain", nodes)
        for i in range(nodes.numInternalNodes):
            s[i] = ts[i].xx
        sl = ScalarFieldList()
        sl.appendField(s)
        sPlot = plotFieldList(sl, winTitle="strain @ %g %i" % (control.time(), mpi.procs),
                              plotStyle="o-")
        epsMin = damageModel.minFlaw
        nflaws = damageModel.numFlaws
        epsl = ScalarFieldList()
        epsl.appendField(epsMin)
        epsPlot = plotFieldList(epsl, winTitle="Min flaw activation strains",
                                plotStyle="o-")
        nflawsl = IntFieldList()
        nflawsl.appendField(nflaws)
        nflawsPlot = plotFieldList(nflawsl, winTitle="Number of flaws",
                                   plotStyle="o-")
        plots += [(sPlot, "strain.png"),
                  (epsPlot, "flaws.png"),
                  (nflawsPlot, "numflaws.png")]

    fragPlot = plotFieldList(state.intFields(SolidFieldNames.fragmentIDs),
                             plotStyle = "o-",
                             winTitle = "Fragments @  %g %i" % (control.time(), mpi.procs))
    plots.append((fragPlot, "fragIDs.png"))

    # Save the figures.
    for p, fname in plots:
        savefig(p, os.path.join(dataDir, fname))

#-------------------------------------------------------------------------------
# If requested, write out the state in a global ordering to a file.
#-------------------------------------------------------------------------------
if outputFile:
    from SpheralTestUtilities import multiSort
    referenceFile = os.path.join(refDir, outputFile)
    outputFile = os.path.join(dataDir, outputFile)
    state = State(db, integrator.physicsPackages())
    pos = state.vectorFields(HydroFieldNames.position)
    rho = state.scalarFields(HydroFieldNames.massDensity)
    P = state.scalarFields(HydroFieldNames.pressure)
    vel = state.vectorFields(HydroFieldNames.velocity)
    eps = state.scalarFields(HydroFieldNames.specificThermalEnergy)
    Hfield = state.symTensorFields(HydroFieldNames.H)
    S = state.symTensorFields(SolidFieldNames.deviatoricStress)
    D = state.symTensorFields(SolidFieldNames.tensorDamage)
    xprof = mpi.reduce([x.x for x in internalValues(pos)], mpi.SUM)
    rhoprof = mpi.reduce(internalValues(rho), mpi.SUM)
    Pprof = mpi.reduce(internalValues(P), mpi.SUM)
    vprof = mpi.reduce([v.x for v in internalValues(vel)], mpi.SUM)
    epsprof = mpi.reduce(internalValues(eps), mpi.SUM)
    hprof = mpi.reduce([1.0/sqrt(H.Determinant()) for H in internalValues(Hfield)], mpi.SUM)
    sprof = mpi.reduce([x.xx for x in internalValues(S)], mpi.SUM)
    dprof = mpi.reduce([x.xx for x in internalValues(D)], mpi.SUM)
    mof = mortonOrderIndices(db)
    mo = mpi.reduce(internalValues(mof), mpi.SUM)
    if mpi.rank == 0:
        multiSort(mo, xprof, rhoprof, Pprof, vprof, epsprof, hprof, sprof, dprof)
        f = open(outputFile, "w")
        f.write(("#" + 8*" %16s" + "\n") % ("x", "rho", "P", "v", "eps", "h", "S", "D"))
        for (xi, rhoi, Pi, vi, epsi, hi, si, di) in zip(xprof, rhoprof, Pprof, vprof, epsprof, hprof, sprof, dprof):
            f.write((8*"%16.12e " + "\n") %
                    (xi, rhoi, Pi, vi, epsi, hi, si, di))
        f.close()

        if checkRef:   #if BuildData.cxx_compiler_id == "GNU":
            #---------------------------------------------------------------------------
            # Check the floating values for the state against reference data.
            #---------------------------------------------------------------------------
            import filearraycmp as fcomp
            assert fcomp.filearraycmp(outputFile, referenceFile, testtol, testtol)
            print("Floating point comparison test passed.")

            #---------------------------------------------------------------------------
            # Also we can optionally compare the current results with another file for
            # bit level consistency.
            #---------------------------------------------------------------------------
            if comparisonFile:
                comparisonFile = os.path.join(dataDir, comparisonFile)
                import filecmp
                assert filecmp.cmp(outputFile, comparisonFile)
