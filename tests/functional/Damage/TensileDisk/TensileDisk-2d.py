#-------------------------------------------------------------------------------
# A disk of stainless steel undergoing tensile strain.  This is intended as a
# test of the cracking/failure models.
#
# This test forces a constant radial expansion on the outermost nodes of the
# disk.
#-------------------------------------------------------------------------------
from SolidSpheral2d import *
from SpheralTestUtilities import *
from identifyFragments import identifyFragments, fragmentProperties
from GenerateNodeDistribution2d import GenerateNodeDistribution2d as GenerateNodeDistribution
from VoronoiDistributeNodes import distributeNodes2d as distributor
from math import *
import os, shutil
import mpi

#-------------------------------------------------------------------------------
# Identify ourselves!
#-------------------------------------------------------------------------------
title("2-D Tensile disk strength/damage model test")

#-------------------------------------------------------------------------------
# Build our units.
#-------------------------------------------------------------------------------
units = CGuS()  # (cm, gm, micorsec units)

#-------------------------------------------------------------------------------
# Generic problem parameters
# All (cm, gm, usec) units.
#-------------------------------------------------------------------------------
commandLine(

    # Geometry and initial conditions
    R0 = 1.0,
    nr = 100,
    vr0 = 5e-3,
    thetaFactor = 2.0,        # one of (0.5, 1.0, 2.0) -- how much of disk geometry to generate
    constantBoundary = True,  # force constant expansion on outer boundary nodes

    # Parameters for the time dependent strain and cracking.
    DamageModelConstructor = GradyKippTensorDamageOwen,
    useDamage = True,
    damageCoupling = ThreePointDamage,

    # Hydro
    nPerh = 4.01,
    crksph = False,     # Use CRK hydro?
    asph = False,        # Just the H tensor evolution -- applies to all hydros
    hminratio = 0.05,
    XSPH = False,
    densityUpdate = IntegrateDensity,
    compatibleEnergy = True,
    evolveTotalEnergy = False,
    correctionOrder = LinearOrder,
    volumeType = RKVoronoiVolume,

    # Time integration
    goalTime = 200.0,
    dt = 1e-10,
    dtMin = 1e-6,
    dtMax = 10.0,
    dtGrowth = 2.0,
    steps = None,
    maxSteps = None,
    domainIndependent = False,
    dtverbose = False,

    # Output
    dtSample = 1.0,
    statsStep = 10,
    redistributeStep = None,
    vizCycle = None,
    vizTime = 0.1,
    restartStep = 500,
    plotFlaws = False,
    clearDirectories = False,
    dataDirBase = "dumps-TensileDisk-2d",
    outputFile = None,

    # Should we restart (-1 => find most advanced available restart)
    restoreCycle = -1,
)

# Check input
assert thetaFactor in (0.5, 1.0, 2.0), "ERROR: thetaFactor must be one of (0.5, 1.0, 2.0)"
assert not (evolveTotalEnergy and compatibleEnergy), "ERROR, can only select one of (compatibleEnergy, evolveTotalEnergy)"

if crksph:
    hydroname = os.path.join("CRKSPH", str(volumeType))
else:
    hydroname = "SPH"
if asph:
    hydroname = "A" + hydroname

dataDir = os.path.join(dataDirBase,
                       "vr0=%g" % vr0,
                       "thetaFactor=%g" % thetaFactor,
                       hydroname,
                       "useDamage=%s" % useDamage,
                       DamageModelConstructor.__name__,
                       "damageCoupling=%s" % damageCoupling,
                       "nr=%i" % nr)
restartDir = os.path.join(dataDir, "restarts")
restartBaseName = os.path.join(restartDir, "TensileDisk-%i" % nr)
vizDir = os.path.join(dataDir, "visit")
vizBaseName = "TensileDisk-%i" % nr

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
rho0 = 7.9
etamin = 0.5
etamax = 1.5
eos = GruneisenEquationOfState("steel",
                               etamin = etamin,
                               etamax = etamax,
                               units = units)
rho0 = eos.referenceDensity
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
WT = TableKernel(WendlandC4Kernel(), 1000)
output("WT")

#-------------------------------------------------------------------------------
# Create the NodeLists.
#-------------------------------------------------------------------------------
hmin = 0.1*nPerh*R0/nr
hmax = 2.0*nPerh*R0/nr
nodes = makeSolidNodeList("Stainless steel", eos, strengthModel,
                          nPerh = nPerh,
                          hmin = hmin,
                          hmax = hmax,
                          hminratio = hminratio,
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
generator = GenerateNodeDistribution(nRadial = nr,
                                     nTheta = 1,     # Not used
                                     rho = rho0,
                                     distributionType = "constantDTheta",
                                     rmin = 0.0,
                                     rmax = R0,
                                     theta = thetaFactor * pi,
                                     nNodePerh = nPerh,
                                     SPH = not asph)
distributor((nodes, generator))
output("mpi.reduce(nodes.numInternalNodes, mpi.MIN)")
output("mpi.reduce(nodes.numInternalNodes, mpi.MAX)")
output("mpi.reduce(nodes.numInternalNodes, mpi.SUM)")

# # Set node specific thermal energies
# eps0 = eos.specificThermalEnergy(rho0, 300.0)
# nodes.specificThermalEnergy(ScalarField("tmp", nodes, eps0))

# Set node velocites to an initial linear (in radius) gradient
pos = nodes.positions()
vel = nodes.velocity()
for i in range(nodes.numInternalNodes):
    runit = pos[i].unitVector()
    vi = pos[i].magnitude()/R0 * vr0
    vel[i] = vi*runit

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
# Do we need any reflecting boundaries?
#-------------------------------------------------------------------------------
bcs = []
if thetaFactor in (0.5, 1.0):
    xbc = ReflectingBoundary(Plane(Vector(0, 0), Vector(1, 0)))
    bcs.append(xbc)
if thetaFactor == 0.5:
    ybc = ReflectingBoundary(Plane(Vector(0, 0), Vector(0, 1)))
    bcs.append(ybc)

#-------------------------------------------------------------------------------
# Construct constant velocity boundary conditions to be applied to the outer
# rind of the disk
#-------------------------------------------------------------------------------
if constantBoundary:
    dr = R0/nr
    constantNodes = vector_of_int([i for i in range(nodes.numInternalNodes)
                                   if pos[i].magnitude() > R0 - WT.kernelExtent*nPerh*dr])
    print("Selected %i constant velocity nodes." % mpi.allreduce(len(constantNodes)))

    # Set the nodes we're going to control to one single radial velocity 
    for i in constantNodes:
        vel[i] = vr0 * pos[i].unitVector()

    rbc = ConstantVelocityBoundary(nodes, constantNodes)
    bcs.append(rbc)

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
if crksph:
    hydro = CRKSPH(dataBase = db,
                   order = correctionOrder,
                   compatibleEnergyEvolution = compatibleEnergy,
                   evolveTotalEnergy = evolveTotalEnergy,
                   XSPH = XSPH,
                   densityUpdate = densityUpdate,
                   damageRelieveRubble = True,
                   ASPH = asph)
else:
    hydro = SPH(dataBase = db,
                W = WT,
                compatibleEnergyEvolution = compatibleEnergy,
                evolveTotalEnergy = evolveTotalEnergy,
                densityUpdate = densityUpdate,
                XSPH = XSPH,
                ASPH = asph)
output("hydro")
output("hydro.cfl")
output("hydro.useVelocityMagnitudeForDt")
output("hydro.HEvolution")
output("hydro.densityUpdate")
output("hydro.compatibleEnergyEvolution")

packages = [hydro]

#-------------------------------------------------------------------------------
# Construct a damage model.
#-------------------------------------------------------------------------------
# Parameters for Grady-Kipp models
#kWeibull = 8.8e4
#kWeibull = 6.52e3
kWeibull = 6.52e5
mWeibull = 2.63
numFlawsPerNode = 1
kWeibullFactor = 1.0
mWeibullFactor = 1.0
randomSeed = 548928513
strainType = PseudoPlasticStrain # BenzAsphaugStrain #
cullToWeakestFlaws = False
damageInCompression = False
negativePressureInDamage = False
volumeMultiplier = 1.0

# Johnson-Cook choices
D1 = 0.0
D2 = 2.0
D3 = -1.5
D4 = 0.0
D5 = 0.0
epsilondot0 = 0.01
Tcrit = -2.0
sigmamax = -3.0
efailmin = 0.01
sigmaD1 = 0.0  # Gaussian
aD1 = 0.0      # Weibull
bD1 = 0.0      # Weibull
eps0D1 = 0.0   # Weibull
sigmaD2 = 0.05 # Gaussian
aD2 = 0.065    # Weibull
bD2 = 2.0      # Weibull
eps0D2 = 0.165 # Weibull

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
                                         strainAlgorithm = strainType,
                                         minFlawsPerNode = numFlawsPerNode,
                                         damageCouplingAlgorithm = damageCoupling,
                                         damageInCompression = damageInCompression)

elif DamageModelConstructor is JohnsonCookDamageWeibull:
    damageModel = DamageModelConstructor(nodes,
                                         D1 = D1,
                                         D2 = D2,
                                         D3 = D3,
                                         D4 = D4,
                                         D5 = D5,
                                         aD1 = aD1,
                                         bD1 = bD1,
                                         eps0D1 = eps0D1,
                                         aD2 = aD2,
                                         bD2 = bD2,
                                         eps0D2 = eps0D2,
                                         epsilondot0 = epsilondot0,
                                         Tcrit = Tcrit,
                                         sigmamax = sigmamax,
                                         efailmin = efailmin,
                                         seed = randomSeed,
                                         domainIndependent = domainIndependent)

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

vizFields = []
if isinstance(damageModel, JohnsonCookDamage):
    vizFields += [damageModel.D1(), damageModel.D2()]

output("damageModel")

if cullToWeakestFlaws:
    damageModel.cullToWeakestFlaws()

# damageModel.excludeNodes = controlNodes
output("damageModel")

if useDamage:
    packages.append(damageModel)

#-------------------------------------------------------------------------------
# Construct a time integrator.
#-------------------------------------------------------------------------------
integrator = VerletIntegrator(db)
for package in packages:
    integrator.appendPhysicsPackage(package)
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
                            SPH = not ASPH,
                            vizFields = vizFields)
output("control")

#-------------------------------------------------------------------------------
# Monitor the evolution of the mass averaged strain.
#-------------------------------------------------------------------------------
if DamageModelConstructor in (GradyKippTensorDamageBenzAsphaug, GradyKippTensorDamageOwen):
    strainHistory = AverageStrain(damageModel,
                                      dataDir + "/strainhistory.txt")
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
# Plot the final flaw distribution, if requested.
#-------------------------------------------------------------------------------
if plotFlaws and DamageModelConstructor is StochasticWeibullTensorDamageModel2d:
    from SpheralGnuPlotUtilities import generateNewGnuPlot
    import Gnuplot
    flaws = [x for x in damageModel.sortedFlaws()]
    f = [i/volume for i in range(len(flaws))]
    fans = [kWeibull * x**mWeibull for x in flaws]
    d = Gnuplot.Data(flaws, f,
                     title = "Simulation",
                     inline = True)
    dans = Gnuplot.Data(flaws, fans,
                        with_ = "lines",
                        title = "Analytic",
                        inline = True)
    p = generateNewGnuPlot()
    p("set logscale xy")
    p("set key top left")
    p.plot(d)
    p.replot(dans)

#-------------------------------------------------------------------------------
# Compute the fragment properties and spit them out to a file.
#-------------------------------------------------------------------------------
## import pickle
## if DamageModelConstructor == GradyKippScalarDamage2d:
##     ufragIndex, dfragIndex = identifyFragments(nodes,
##                                                1.0,
##                                                damageModel,
##                                                nodesDamaged)
##     fragProps = fragmentProperties(nodes, ufragIndex,
##                                    nodesDamaged, dfragIndex)
## else:
##     ufragIndex = identifyFragments(nodes, 1.0)
##     fragProps = fragmentProperties(nodes, ufragIndex)

## if procID == 0:
##     f = open(dataDir + "/fragment-properties-t=%3.1fmicrosec" % (1.0e6*control.time()), "w")
##     pickle.dump(fragProps, f)
##     f.close()
