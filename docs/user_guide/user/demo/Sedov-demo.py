#-------------------------------------------------------------------------------
# The cylindrical Sedov test case (2-D).
#-------------------------------------------------------------------------------
import os, shutil
import mpi
from Spheral2d import *
from SpheralTestUtilities import *
from GenerateNodeDistribution2d import *
from PeanoHilbertDistributeNodes import distributeNodes2d
from SpheralMatplotlib import *

title("2-D SPH hydrodynamics demonstration of the cylindrical Sedov problem")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(
    # Geometry of the problem
    rmin = 0.0,
    rmax = 1.0,

    # Resolution and how to lay down the points
    nRadial = 50,
    nPerh = 4.01,

    # Initial conditions
    rho0 = 1.0,
    eps0 = 0.0,
    Espike = 0.25*1.0,

    # Material equation of state options (gamma-law gas in this case)
    gamma = 5.0/3.0,
    mu = 1.0,

    # Time stepping
    goalTime = 0.5,
    dt = 1.0e-8,
    dtGrowth = 2.0,

    # IO
    vizCycle = None,
    vizTime = 0.1,
    restoreCycle = -1,
    restartStep = 1000,
    clearDirectories = False,
    dataDirBase = "dumps-cylindrical-Sedov",
    graphics = True,
)

#-------------------------------------------------------------------------------
# Path names.
#-------------------------------------------------------------------------------
dataDir = os.path.join(dataDirBase,
                       "nr={}".format(nRadial))
restartDir = os.path.join(dataDir, "restarts")
vizDir = os.path.join(dataDir, "visit")
restartBaseName = os.path.join(restartDir, "Sedov-cylindrical-2d-%i" % nRadial)

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
units = MKS()
eos = GammaLawGas(gamma = gamma,
                  mu = mu,
                  constants = units)

#-------------------------------------------------------------------------------
# Create our interpolation kernel
#-------------------------------------------------------------------------------
WT = TableKernel(WendlandC4Kernel())
output("WT")

#-------------------------------------------------------------------------------
# Create a NodeList and associated Neighbor object.
#-------------------------------------------------------------------------------
nodes = makeSolidNodeList(name = "gamma-law gas points",
                          eos = eos, 
                          kernelExtent = WT.kernelExtent,
                          nPerh = nPerh)

#-------------------------------------------------------------------------------
# Generate the initial node geometry (lay down the points).
#-------------------------------------------------------------------------------
generator = GenerateNodeDistribution2d(nRadial = nRadial,
                                       nTheta = 1,
                                       rho = rho0,
                                       distributionType = "constantDTheta",
                                       rmin = rmin,
                                       rmax = rmax,
                                       xmin = (0.0, 0.0),
                                       xmax = (rmax, rmax),
                                       theta = 0.5*pi,
                                       nNodePerh = nPerh)
distributeNodes2d((nodes, generator))
print("Point distribution across MPI ranks:")
output("  mpi.reduce(nodes.numInternalNodes, mpi.MIN)")
output("  mpi.reduce(nodes.numInternalNodes, mpi.MAX)")
output("  mpi.reduce(nodes.numInternalNodes, mpi.SUM)")

#-------------------------------------------------------------------------------
# Set the node properties.
# The above geometry initialization set the node positions, masses, mass
# density, and H tensors.  We still need to initialize the energy for the Sedov
# blast energy.
# In this case we're going to simply pick the innermost ring of points closest
# to the origin and dump the initial spike evenly across them.
#-------------------------------------------------------------------------------
pos = nodes.positions()
mass = nodes.mass()
eps = nodes.specificThermalEnergy()

# Find the points closest to the origin.
dr = rmax/nRadial
spikePoints = [i for i in range(nodes.numInternalNodes) if pos[i].magnitude() < 0.75*dr]
print("Selected a total of {} points for energy deposition.".format(mpi.allreduce(len(spikePoints))))

# Distribute the energy evenly in energy/mass.
Mspike = mpi.allreduce(sum([mass[i] for i in spikePoints]), mpi.SUM)
eps0 = Espike/Mspike
for i in spikePoints:
    eps[i] = eps0

# Verify we actually initialized the expected energy.
# We are using multiplication of Spheral Fields here, and the fact that Field.sumElements
# sums across all MPI domains by default
Esum = (mass*eps).sumElements()
print("Initialized a total energy of {}".format(Esum))
assert fuzzyEqual(Esum, Espike)

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase()
db.appendNodeList(nodes)
output("db")
output("  db.numNodeLists")
output("  db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
hydro = FSISPH(dataBase = db,
               W = WT)
packages = [hydro]

output("hydro")
output("  hydro.cfl")
output("  hydro.compatibleEnergyEvolution")
output("  hydro.densityUpdate")
output("  hydro.HEvolution")

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
xPlane0 = Plane(Vector(0.0, 0.0), Vector(1.0, 0.0))
yPlane0 = Plane(Vector(0.0, 0.0), Vector(0.0, 1.0))
xbc0 = ReflectingBoundary(xPlane0)
ybc0 = ReflectingBoundary(yPlane0)
for p in packages:
    p.appendBoundary(xbc0)
    p.appendBoundary(ybc0)

#-------------------------------------------------------------------------------
# Construct a time integrator, and add the one physics package.
#-------------------------------------------------------------------------------
integrator = CheapSynchronousRK2Integrator(db, packages)
integrator.lastDt = dt
integrator.allowDtCheck = True
output("integrator")
output("  integrator.havePhysicsPackage(hydro)")
output("  integrator.dtGrowth")
output("  integrator.lastDt")
output("  integrator.dtMin")
output("  integrator.dtMax")
output("  integrator.allowDtCheck")

#-------------------------------------------------------------------------------
# Build the controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator = integrator,
                            kernel = WT,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle,
                            vizBaseName = "Sedov-cylindrical-{}".format(nRadial),
                            vizDir = vizDir,
                            vizStep = vizCycle,
                            vizTime = vizTime,
                            SPH = True)
output("control")

#-------------------------------------------------------------------------------
# Finally run the problem and plot the results.
#-------------------------------------------------------------------------------
control.advance(goalTime)

#-------------------------------------------------------------------------------
# Plot the final state if desired.
#-------------------------------------------------------------------------------
if graphics:
    rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotRadialState(db)
