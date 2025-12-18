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

xlength = 20.0
nx = 96
nPerh = 2.01

xmin = -0.5*xlength
xmax =  0.5*xlength

rho0 = 7.9
m0 = xlength*rho0/nx
dx = xlength/nx
hx = nPerh*dx
H0 = SymTensor1d(1.0/hx)

# Parameters for the time dependent strain and cracking.
tfail = 20.0e-6
v0 = 1.0e3
xminFlaws = -1.0*dx
xmaxFlaws =  1.0*dx

Qconstructor = MonaghanGingoldViscosity1d
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

goalTime = 10.0e-6
dtSample = 5.0e-6
dt = 1e-10
dtMin, dtMax = 1e-12, 1e-5
dtGrowth = 2.0
maxSteps = 50
statsStep = 10
smoothIters = 0
HEvolution = Hydro1d.HEvolutionType.IdealH
sumForMassDensity = Hydro1d.MassDensityType.IntegrateDensity # HybridDensity # CorrectedSumDensity

restartStep = 1000
restartBaseName = "dumps-specifiedFlaws-1d/TensileRod-%i" % (nx)
restoreCycle = findLastRestart(restartBaseName)

#-------------------------------------------------------------------------------
# Identify ourselves!
#-------------------------------------------------------------------------------
title("1-D Tensile rod strength/damage model test")

#-------------------------------------------------------------------------------
# Stainless steel material properties.
#-------------------------------------------------------------------------------
eos = GruneisenEquationOfStateCGS1d(7.90,    # reference density  
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
strengthModel = SteinbergGuinanStrengthCGS1d(eos,
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
eosDamaged = GruneisenEquationOfStateCGS1d(7.90,    # reference density  
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
nodes = SphSolidNodeList1d("Stainless steel", eos, strengthModel)
nodesDamaged = SphNodeList1d("Damaged stainless steel", eosDamaged)
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
WT = TableKernel1d(BSplineKernel1d(), 1000)
WTPi = TableKernel1d(BSplineKernel1d(), 1000)
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
    neighbor = TreeNeighbor1d(n,
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
    from DistributeNodes import distributeNodes1d
    distributeNodes1d([(nodes, nx, (xmin, xmax))])
    output('mpi.reduce(nodes.numInternalNodes, mpi.MIN)')
    output('mpi.reduce(nodes.numInternalNodes, mpi.MAX)')
    output('mpi.reduce(nodes.numInternalNodes, mpi.SUM)')

    # Set the node masses.
    nodes.setMass(ScalarField1d("tmp", nodes, m0))

    # Set the smoothing scales.
    nodes.setHfield(SymTensorField1d("tmp", nodes, H0))

    # Set the node mass densities.
    nodes.setMassDensity(ScalarField1d("tmp", nodes, rho0))
    nodes.updateWeight()

    # Set node specific thermal energies
    u0 = eos.specificThermalEnergy(rho0, 300.0)
    nodes.setSpecificThermalEnergy(ScalarField1d("tmp", nodes, u0))

    # Set node velocites.
    for i in range(nodes.numInternalNodes):
        nodes.velocity()[i].x = nodes.positions()[i].x/(0.5*xlength)*v0

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase1d()
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
hydro = Hydro1d(WT, WTPi, q)
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
strength = Strength1d()
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
backgroundActivationStrain = 1.0e5*failureEnergy

initialFlaws = vector_of_vector_of_double()
numFlawedNodes = 0
for i in range(nodes.numInternalNodes):
    x = nodes.positions()[i].x
    v = vector_of_double()
    if x > xminFlaws and x < xmaxFlaws:
        v.append(flawActivationStrain)
        numFlawedNodes += 1
    else:
        v.append(backgroundActivationStrain)
    initialFlaws.append(v)
assert len(initialFlaws) == nodes.numInternalNodes
damageModel = SpecifiedFlawsDamage1d(nodes,
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
for i in range(nodes.numInternalNodes):
    if (nodes.positions()[i].x < -0.5*xlength + 4*dx or
        nodes.positions()[i].x >  0.5*xlength - 4*dx):
        xNodes.append(i)
print(("Selected %i constant velocity nodes." %
       (mpi.allreduce(len(xNodes), mpi.SUM))))
xbc = ConstantVelocityBoundary1d(nodes, xNodes)

# We need to tell the damage model not to damage the constant velocity nodes,
# or we'll lose our boundary condition!
damageModel.excludeNodes = xNodes

#-------------------------------------------------------------------------------
# Construct a predictor corrector integrator.
#-------------------------------------------------------------------------------
integrator = PredictorCorrectorIntegrator1d(db)
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

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
while control.time() < goalTime:
    nextGoalTime = min(control.time() + dtSample, goalTime)
    control.advance(nextGoalTime, maxSteps)
    control.dropRestartFile()

# Plot the state.
rhoPlot = plotFieldList(db.fluidMassDensity,
                        plotStyle="points",
                        winTitle="rho @ %g %i" % (control.time(), mpi.procs))
velPlot = plotFieldList(db.fluidVelocity,
                        yFunction = "%s.x",
                        plotStyle="points",
                        winTitle="vel @ %g %i" % (control.time(), mpi.procs))
mPlot = plotFieldList(db.fluidMass,
                      plotStyle="points",
                      winTitle="mass @ %g %i" % (control.time(), mpi.procs))
PPlot = plotFieldList(db.fluidPressure,
                      plotStyle="points",
                      winTitle="pressure @ %g %i" % (control.time(), mpi.procs))
HPlot = plotFieldList(db.fluidHfield,
                      yFunction = "%s.xx",
                      plotStyle="points",
                      winTitle="H @ %g %i" % (control.time(), mpi.procs))

s = damageModel.strain()
sl = ScalarFieldList1d()
sl.appendField(s)
sPlot = plotFieldList(sl, winTitle="strain @ %g %i" % (control.time(), mpi.procs))

d = damageModel.damage()
dl = ScalarFieldList1d()
dl.appendField(d)
dPlot = plotFieldList(dl, winTitle="damage @ %g %i" % (control.time(), mpi.procs))
