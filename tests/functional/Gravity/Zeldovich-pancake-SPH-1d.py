#-------------------------------------------------------------------------------
# The evolution of a Z'eldovich pancake in 1D.
#-------------------------------------------------------------------------------
from math import *
from Spheral import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from SpheralVisitDump import dumpPhysicsState
import mpi

title("1-D Zeldovich pancake collapse test")

# Plotter thingy.
class Plotter(object):
    def __init__(self, gravity, nodes):
        self.gravity = gravity
        self.nodes = nodes
    def __call__(self, cycle, time, dt):
        from pylab import plot, semilogy, ion, clf, savefig, \
             title, xlabel, ylabel
        ion()
        rho = self.nodes.massDensity().internalValues()
        x = [u.x for u in self.nodes.positions().internalValues()]
        semilogy(x, rho, 'k.')
        title('Mass density')
        xlabel('x')
        ylabel('log density')
        savefig('pancake-rho-%3d.png'%cycle)
        clf()
        v = [u.x for u in self.nodes.velocity().internalValues()]
        plot(x, v, 'k.')
        title('Velocity')
        savefig('pancake-v-%3d.png'%cycle)
        clf()
        psi = self.gravity.potential()[0].internalValues()
        acc = [u.x for u in self.nodes.DvelocityDt().internalValues()]
        plot(x, psi, 'k.', x, acc, 'r.')
        title('Gravitational potential, acceleration')
        xlabel('x')
        ylabel('potential')
        savefig('pancake-psi-%3d.png'%cycle)
        clf()


#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(NodeListConstructor = SphNodeList1d,

            nx = 101,
            
            rho0 = 1.0, # Initial background density
            u0 = 0.0, # Specific thermal energy
            vr1 = -1.0,
            nPerh = 1.25,

            gamma = 5.0/3.0,
            mu = 1.0,
            z1 = 50., # Starting redshift.
            z2 = 5., # Ending redshift.
            zc = 5., # Caustic redshift.
            L = 1.0, # Length scale of simulation 
            #Qconstructor = MonaghanGingoldViscosity1d,
            Qconstructor = TensorMonaghanGingoldViscosity1d,
            Cl = 1.0,
            Cq = 0.75,
            Qlimiter = True,
            balsaraCorrection = False,
            epsilon2 = 1e-2,
            negligibleSoundSpeed = 1e-5,
            csMultiplier = 1e-4,
            hmin = 1e-5,
            hmax = 1.0,
            hminratio = 0.05,
            HsmoothFraction = 0.0,
            cfl = 0.5,
            XSPH = False,
            epsilonTensile = 0.0,
            nTensile = 8,

            HEvolution = Hydro1d.HEvolutionType.IdealH,
            limitIdealH = False,

            dt = 0.0001,
            dtMin = 1.0e-5,
            dtMax = None,
            dtGrowth = 2.0,
            dtSample = 0.1,
            maxSteps = None,
            statsStep = 10,
            smoothIters = 0,
            sumForMassDensity = Hydro1d.MassDensityType.RigorousSumDensity,

            restoreCycle = None,
            plotEvery = None,

            L1v0 =   0.0889732,
            L1rho0 = 5.51975,
            L1u0 = 0.04701,
            L1P0 =   1.66301,
            L1A0 =   0.00344783,

            graphics = False,
            )

# nx must be odd!
if nx/2 != (nx-1)/2:
   raise ValueError('nx must be odd!')

# Boundaries of the simulation
xmin = -0.5*L
xmax =  0.5*L

import os.path
dumpName = "Zeldovich-pancake-SPH-%i"%nx
visitDir = os.path.join(dumpName, "visit")

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
import os, sys
if mpi.rank == 0:
    if not os.path.exists(visitDir):
        os.makedirs(visitDir)
mpi.barrier()

#-------------------------------------------------------------------------------
# Material properties.
#-------------------------------------------------------------------------------
#eos = GammaLawGasMKS1d(gamma, mu)
eos = PolytropicEquationOfStateMKS1d(0.0, 1.0, mu)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
WT = TableKernel1d(BSplineKernel1d(), 1000)
WTPi = TableKernel1d(BSplineKernel1d(), 1000)
output("WT")
output("WTPi")
kernelExtent = WT.kernelExtent()

#-------------------------------------------------------------------------------
# Make the NodeList.
#-------------------------------------------------------------------------------
nodes = NodeListConstructor("nodes", eos, WT, WTPi)
output("nodes")
nodes.HsmoothFraction = HsmoothFraction
nodes.XSPH = XSPH
nodes.nodesPerSmoothingScale = nPerh
nodes.epsilonTensile = epsilonTensile
nodes.nTensile = nTensile
nodes.hmin = hmin
nodes.hmax = hmax
nodes.hminratio = hminratio
output("nodes.HsmoothFraction")
output("nodes.nodesPerSmoothingScale")
output("nodes.epsilonTensile")
output("nodes.nTensile")
output("nodes.XSPH")
output("nodes.hmin")
output("nodes.hmax")
output("nodes.hminratio")

#-------------------------------------------------------------------------------
# Construct the neighbor object.
#-------------------------------------------------------------------------------
neighbor1 = TreeNeighbor1d(nodes,
                           kernelExtent = kernelExtent)
nodes.registerNeighbor(neighbor1)

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
from DistributeNodes import distributeNodesInRange1d
distributeNodesInRange1d([(nodes, nx, rho0, (xmin, xmax))])
nNodesThisDomain1 = nodes.numInternalNodes
output("nodes.numNodes")

# Back out the start and end times and expansion factors from the start 
# and end redshifts z1 and z2.  In a Friedmann cosmology the expansion 
# factor a(t) is proportional to t**2/3 in the matter-dominated phase.
# Of course, this applies directly only in co-moving coordianates.  We
# use Shandarin's coordinate transformation to solve the equations in 
# a stationary coordinate system.
a1 = 1.0/(1+z1)
t1 = -1.0/sqrt(a1)
a1dot = -2.0*t1**-3
a2 = 1.0/(1+z2)
t2 = -1.0/sqrt(a2)

# Set node specific thermal energies.
nodes.specificThermalEnergy(ScalarField1d("tmp", nodes, u0))

# Set node positions and velocities.
k = 2*pi/(xmax-xmin)
A = -5*(1+zc)/(2*k)
for nodeID in range(nodes.numNodes):
    q = nodes.positions()[nodeID].x
    nodes.positions()[nodeID] += Vector1d(0.4*a1*A*sin(k*q))
    nodes.velocity()[nodeID] += Vector1d(0.4*a1dot*A*sin(k*q))

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase1d()
output("db")
output("db.appendNodeList(nodes)")
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Construct the artificial viscosities for the problem.
#-------------------------------------------------------------------------------
#q = Qconstructor(Cl, Cq)
q = Qconstructor(0, 0)
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

##-------------------------------------------------------------------------------
## Construct the hydro physics object.
##-------------------------------------------------------------------------------
hydro = Hydro1d(WT, WTPi, q)
hydro.cfl = cfl
hydro.HEvolution = HEvolution
hydro.sumForMassDensity = sumForMassDensity
hydro.HsmoothMin = hmin
hydro.HsmoothMax = hmax
#output("hydro")
#output("hydro.cfl")
#output("hydro.HEvolution")
#output("hydro.sumForMassDensity")
#output("hydro.HsmoothMin")
#output("hydro.HsmoothMax")
#output("hydro.kernel()")
#output("hydro.PiKernel()")
#output("hydro.valid()")

#-------------------------------------------------------------------------------
# Construct the gravity physics object.
#-------------------------------------------------------------------------------
#from Spasmos import pause
#pause('Simulation paused...')
#G = 1.0
G = 4.3e-3  # G in units of megasolar masses and Megaparsecs, velocities in km/s.
gravity = SPHGravity1d(WT, G)

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
xPlane1 = Plane1d(Vector1d(xmin), Vector1d( 1.0))
xPlane2 = Plane1d(Vector1d(xmax), Vector1d(-1.0))
xbc = PeriodicBoundary1d(xPlane1, xPlane2)
hydro.appendBoundary(xbc)
gravity.appendBoundary(xbc)

#-------------------------------------------------------------------------------
# Construct a time integrator.
#-------------------------------------------------------------------------------
#integrator = ForwardEulerIntegrator1d(db)
integrator = SynchronousRK2Integrator1d(db)
integrator.appendPhysicsPackage(gravity)
integrator.appendPhysicsPackage(hydro)
integrator.lastDt = dt
if dtMin:
    integrator.dtMin = dtMin
if dtMax:
    integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
output("integrator")
output("integrator.havePhysicsPackage(gravity)")
output("integrator.havePhysicsPackage(hydro)")
output("integrator.valid()")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.dtGrowth")

#-------------------------------------------------------------------------------
# Build the controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            initializeMassDensity = True)
integrator.setCurrentTime(t1)
integrator.setVerbose(1)
plotter = Plotter(gravity, nodes)
if plotEvery is not None:
   control.appendPeriodicWork(plotter, plotEvery)
output("control")

# Initial dump.
if os.path.exists(visitDir):
   import shutil
   shutil.rmtree(visitDir)
dumpPhysicsState(integrator, dumpName, visitDir)

#-------------------------------------------------------------------------------
# Advance.
#-------------------------------------------------------------------------------
print('Running from t = %g to %g'%(t1, t2))
control.advance(t2, maxSteps)
dumpPhysicsState(integrator, dumpName, visitDir)
psi = gravity.potential()[0]
acc = nodes.DvelocityDt()
