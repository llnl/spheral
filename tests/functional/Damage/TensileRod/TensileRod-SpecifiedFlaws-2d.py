#-------------------------------------------------------------------------------
# A rod of stainless steel undergoing tensile strain.  This is intended as a
# test of the cracking/failure models.
#
# This version imposes a constant velocity gradient across the rod, and uses
# a specified flaw distribution to create a known set of weakened nodes.
#
# See Benz & Asphaug (1994), Icarus, 107, 98
#-------------------------------------------------------------------------------
from Numeric import *
from Spheral import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from SpheralController import *
from findLastRestart import *
from SpheralVisitDump import dumpPhysicsState
from math import *

from Strength import *

# Load the mpi module if we're parallel.
mpi, procID, numProcs = loadmpi()

#-------------------------------------------------------------------------------
# Generic problem parameters
# All CGS units.
#-------------------------------------------------------------------------------
seed = 'lattice'

xlength, ylength = 20.0, 5.0
nx, ny = 96, 24
nPerh = 2.01

xmin = (-0.5*xlength, -0.5*ylength)
xmax = ( 0.5*xlength,  0.5*ylength)

rho0 = 7.9
m0 = (xlength*ylength)*rho0/(nx*ny)
dx = xlength/nx
dy = ylength/ny
hx = nPerh*dx
hy = nPerh*dy
H0 = SymTensor2d(1.0/hx, 0.0,
                 0.0, 1.0/hx)

# Parameters for the time dependent strain and cracking.
tfail = 20.0e-6
v0 = 1.0e3
xminFlaws = -1.0*dx
xmaxFlaws =  1.0*dx
yminFlaws = -1.0*dy
ymaxFlaws =  1.0*dy

Qconstructor = MonaghanGingoldViscosity2d
Cl, Cq = 1.0, 1.0
Qlimiter = False
balsaraCorrection = False
epsilon2 = 1e-4
negligibleSoundSpeed = 1e-5
csMultiplier = 1e-4
HsmoothMin, HsmoothMax = 1e-5, 0.5
cfl = 0.5
useVelocityMagnitudeForDt = False
XSPH = False
epsilonTensile = 0.0
nTensile = 4
hybridMassDensityThreshold = 0.01

goalTime = 500.0e-6
dtSample = 5.0e-6
dt = 1e-10
dtMin, dtMax = 1e-12, 1e-5
dtGrowth = 2.0
maxSteps = None
statsStep = 10
smoothIters = 0
HEvolution = Hydro2d.HEvolutionType.IdealH
sumForMassDensity = Hydro2d.MassDensityType.IntegrateDensity # HybridDensity # CorrectedSumDensity

restartStep = 1000
restartBaseName = "dumps-specifiedFlaws-2d/TensileRod-%i-%i" % (nx, ny)
restoreCycle = findLastRestart(restartBaseName)

#-------------------------------------------------------------------------------
# Identify ourselves!
#-------------------------------------------------------------------------------
title("2-D Tensile rod strength/damage model test")

#-------------------------------------------------------------------------------
# Stainless steel material properties.
#-------------------------------------------------------------------------------
eos = GruneisenEquationOfStateCGS2d(7.90,    # reference density  
                                    1e-4,    # etamin             
                                    6.0,     # etamax             
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
strengthModel = SteinbergGuinanStrengthCGS2d(eos,
                                             7.700000e11,        # G0
                                             2.2600e-12,         # A
                                             4.5500-0o4,          # B
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
                                           1e-4,    # etamin             
                                           6.0,     # etamax             
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
# Create the NodeLists.
#-------------------------------------------------------------------------------
nodes = SphSolidNodeList2d("Stainless steel", eos, strengthModel)
nodesDamaged = SphNodeList2d("Damaged stainless steel", eosDamaged)
for n in [nodes, nodesDamaged]:
    n.nodesPerSmoothingScale = nPerh
    n.epsilonTensile = epsilonTensile
    n.nTensile = nTensile
    n.XSPH = XSPH
    output('n.name()')
    output('n.nodesPerSmoothingScale')
    output('n.epsilonTensile')
    output('n.nTensile')
    output('n.XSPH')
del n

#-------------------------------------------------------------------------------
# Create our interpolation kernels -- one for normal hydro interactions, and
# one for use with the artificial viscosity
#-------------------------------------------------------------------------------
WT = TableKernel2d(BSplineKernel2d(), 1000)
WTPi = TableKernel2d(BSplineKernel2d(), 1000)
output('WT')
output('WTPi')
kernelExtent = WT.kernelExtent()

#-------------------------------------------------------------------------------
# Construct the neighbor objects and associate them with the node lists.
#-------------------------------------------------------------------------------
neighborTimer = SpheralTimer('Neighbor initialization.')
neighborTimer.start()
cache = []
for n in [nodes, nodesDamaged]:
    neighbor = TreeNeighbor2d(n,
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
    from GenerateNodeDistribution2d import *
    from DistributeNodes import distributeNodes2d
    generator = GenerateNodeDistribution2d(nx,
                                           ny,
                                           rho0,
                                           seed,
                                           xmin = xmin,
                                           xmax = xmax,
                                           nNodePerh = nPerh)
    n = generator.globalNumNodes()
    nodeInfo = distributeNodes2d([(nodes, n, generator)])
    output('mpi.reduce(nodes.numInternalNodes, mpi.MIN)')
    output('mpi.reduce(nodes.numInternalNodes, mpi.MAX)')
    output('mpi.reduce(nodes.numInternalNodes, mpi.SUM)')

    # Set the node masses.
    nodes.setMass(ScalarField2d("tmp", nodes, m0))

    # Set the smoothing scales.
    nodes.setHfield(SymTensorField2d("tmp", nodes, H0))

    # Set the node mass densities.
    nodes.setMassDensity(ScalarField2d("tmp", nodes, rho0))
    nodes.updateWeight()

    # Set node specific thermal energies
    u0 = eos.specificThermalEnergy(rho0, 300.0)
    nodes.setSpecificThermalEnergy(ScalarField2d("tmp", nodes, u0))

    # Set node velocites.
    for i in range(nodes.numInternalNodes):
        nodes.velocity()[i].x = nodes.positions()[i].x/(0.5*xlength)*v0

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase2d()
db.appendNodeList(nodes)
db.appendNodeList(nodesDamaged)
output('db')
output('db.numNodeLists')
output('db.numFluidNodeLists')

#-------------------------------------------------------------------------------
# Construct the artificial viscosities for the problem.
#-------------------------------------------------------------------------------
q = Qconstructor(Cl, Cq)
q.limiter = Qlimiter
q.balsaraShearCorrection = balsaraCorrection
q.epsilon2 = epsilon2
q.negligibleSoundSpeed = negligibleSoundSpeed
q.csMultiplier = csMultiplier
output('q')
output('q.Cl')
output('q.Cq')
output('q.limiter')
output('q.epsilon2')
output('q.negligibleSoundSpeed')
output('q.csMultiplier')
output('q.balsaraShearCorrection')

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
hydro = Hydro2d(WT, WTPi, q)
hydro.cfl = cfl
hydro.useVelocityMagnitudeForDt = True
hydro.HEvolution = HEvolution
hydro.sumForMassDensity = sumForMassDensity
hydro.HsmoothMin = HsmoothMin
hydro.HsmoothMax = HsmoothMax
hydro.hybridMassDensityThreshold = hybridMassDensityThreshold
output('hydro')
output('hydro.cfl')
output('hydro.useVelocityMagnitudeForDt')
output('hydro.HEvolution')
output('hydro.sumForMassDensity')
output('hydro.HsmoothMin')
output('hydro.HsmoothMax')
output('hydro.kernel()')
output('hydro.PiKernel()')
output('hydro.valid()')
output('hydro.hybridMassDensityThreshold')

#-------------------------------------------------------------------------------
# Construct a strength physics object.
#-------------------------------------------------------------------------------
strength = Strength2d()
output("strength")

#-------------------------------------------------------------------------------
# Construct a damage model.
#-------------------------------------------------------------------------------
# Compute the critical failure threshold for this strain rate.
K0 = eos.bulkModulus(rho0, 0.0)
mu0 = strengthModel.shearModulus(rho0, 0.0, 0.0)
assert distinctlyGreaterThan(3.0*K0 + mu0, 0.0)
E0 = 9.0*K0*mu0/(3.0*K0 + mu0)
assert distinctlyGreaterThan(E0, 0.0)
c0 = eos.C0
failureEnergy = (4.0*mu0*v0/(3.0*xlength)*tfail - rho0*c0*c0*(xlength/(xlength + v0*tfail) - 1.0))/E0

flawActivationStrain = failureEnergy
backgroundActivationStrain = 1.0e3*failureEnergy

initialFlaws = vector_of_vector_of_double()
numFlawedNodes = 0
for i in range(nodes.numInternalNodes):
    x = nodes.positions()[i].x
    y = nodes.positions()[i].y
    v = vector_of_double()
##    if (x > xminFlaws and x < xmaxFlaws and
##        (y > 0.5*ylength - 2.5*dy or
##         y < -0.5*ylength + 2.5*dy)):
##        v.append(flawActivationStrain)
##        numFlawedNodes += 1
##    else:
##        v.append(backgroundActivationStrain)
    if (x > xminFlaws and x < xmaxFlaws and
        y > yminFlaws and y < ymaxFlaws):
        v.append(flawActivationStrain)
        numFlawedNodes += 1
    else:
        v.append(backgroundActivationStrain)
    initialFlaws.append(v)
assert len(initialFlaws) == nodes.numInternalNodes
damageModel = SpecifiedFlawsDamage2d(nodes,
                                     nodesDamaged,
                                     kernelExtent,
                                     0.4,
                                     initialFlaws)
print("Selected %i flawed nodes" % mpi.allreduce(numFlawedNodes, mpi.SUM))
print(("Assigned (%g, %g) flaw activation strains to (flawed, regular) nodes." %
       (flawActivationStrain, backgroundActivationStrain)))
output("damageModel")

#-------------------------------------------------------------------------------
# Construct constant velocity boundary conditions to be applied to the rod ends.
#-------------------------------------------------------------------------------
xNodes = vector_of_ULL()
yNodes = vector_of_ULL()
for i in range(nodes.numInternalNodes):
    if (nodes.positions()[i].x < -0.5*xlength + 4*dx or
        nodes.positions()[i].x >  0.5*xlength - 4*dx):
        xNodes.append(i)
    if (nodes.positions()[i].y < -0.5*ylength + 2*dy or
        nodes.positions()[i].y >  0.5*ylength - 2*dy):
        yNodes.append(i)
print(("Selected (%i, %i) (x, y) constant velocity nodes." %
       (mpi.allreduce(len(xNodes), mpi.SUM),
        mpi.allreduce(len(yNodes), mpi.SUM))))
xbc = ConstantVelocityBoundary2d(nodes, xNodes)
ybc = ConstantYVelocityBoundary2d(nodes, yNodes)

#-------------------------------------------------------------------------------
# Construct a predictor corrector integrator.
#-------------------------------------------------------------------------------
integrator = PredictorCorrectorIntegrator2d(db)
integrator.appendPhysicsPackage(hydro)
integrator.appendPhysicsPackage(strength)
integrator.appendPhysicsPackage(damageModel)
integrator.lastDt = dt
if dtMin:
    integrator.dtMin = dtMin
if dtMax:
    integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
output('integrator')
output('integrator.havePhysicsPackage(hydro)')
output('integrator.havePhysicsPackage(strength)')
output('integrator.havePhysicsPackage(damageModel)')
output('integrator.valid()')
output('integrator.lastDt')
output('integrator.dtMin')
output('integrator.dtMax')
output('integrator.dtGrowth')

#-------------------------------------------------------------------------------
# Build the controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            boundaryConditions = [xbc],
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            initializeMassDensity = False)
output('control')

#-------------------------------------------------------------------------------
# Smooth the initial conditions/restore state.
#-------------------------------------------------------------------------------
if restoreCycle is not None:
    control.loadRestartFile(restoreCycle)
else:
    control.smoothState(smoothIters)

    # Viz the initial conditions.
    vx = ScalarField2d("x velocity", nodes)
    vy = ScalarField2d("y velocity", nodes)
    for i in range(nodes.numInternalNodes):
        vx[i] = nodes.velocity()[i].x
        vy[i] = nodes.velocity()[i].y
    dumpPhysicsState(integrator,
                     "TensileRod-2d-visit",
                     "dumps-specifiedFlaws-2d",
                     fields = [vx, vy,
                               damageModel.sumActivationEnergiesPerNode(),
                               damageModel.numFlawsPerNode()]
                     )

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
while control.time() < goalTime:
    nextGoalTime = min(control.time() + dtSample, goalTime)
    control.advance(nextGoalTime, maxSteps)
    control.dropRestartFile()

    # Viz the current state.
    vx = ScalarField2d("x velocity", nodes)
    vy = ScalarField2d("y velocity", nodes)
    for i in range(nodes.numInternalNodes):
        vx[i] = nodes.velocity()[i].x
        vy[i] = nodes.velocity()[i].y
    dumpPhysicsState(integrator,
                     "TensileRod-2d-visit",
                     "dumps-specifiedFlaws-2d",
                     fields = [vx, vy,
                               damageModel.sumActivationEnergiesPerNode(),
                               damageModel.numFlawsPerNode()]
                     )
