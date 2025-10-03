from Spheral import *
from SpheralTestUtilities import *
import os, shutil

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(geometry = "1d",    # (1d, 2d, 3d, RZ)
            
            nx1 = 100,
            nx2 = 100,

            rho1 = 1.0,
            rho2 = 1.0,
            vx0 = 1.0,

            x0 = -1.0,
            x1 = 0.0,
            x2 = 1.0,
            y0 = 0.0,
            y1 = 0.1,
            z0 = 0.0,
            z1 = 0.1,

            nPerh = 4.0,
            hmin = 1e-10,
            hmax = 1.0,

            dtMin = 1e-5,
            dtMax = 100.0,
            steps = None,
            goalTime = 1.0,

            graphics = True,
            vizStep = 1,
            clearDirectories = True,
            )

assert geometry in ("1d", "2d", "3d", "RZ")
exec("from Spheral%s import *" % geometry, globals())

vizDir = "dumps-InflowOutflow-%s" % geometry

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
import os, sys
if mpi.rank == 0:
    if clearDirectories and os.path.exists(vizDir):
        shutil.rmtree(vizDir)
    if not os.path.exists(vizDir):
        os.makedirs(vizDir)
mpi.barrier()

#-------------------------------------------------------------------------------
# Material properties.
#-------------------------------------------------------------------------------
units = CGS()
gamma = 5.0/3.0
mu = 1.0
eos = GammaLawGas(gamma, mu, units)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
WT = TableKernel(WendlandC4Kernel(), 1000)
kernelExtent = WT.kernelExtent
output("WT")

#-------------------------------------------------------------------------------
# Make the NodeList.
#-------------------------------------------------------------------------------
nodes1 = makeFluidNodeList("nodes1", eos, 
                           hmin = hmin,
                           hmax = hmax,
                           nPerh = nPerh,
                           kernelExtent = kernelExtent)
nodes2 = makeFluidNodeList("nodes2", eos, 
                           hmin = hmin,
                           hmax = hmax,
                           nPerh = nPerh,
                           kernelExtent = kernelExtent)
nodeSet = [nodes1, nodes2]
for nodes in nodeSet:
    output("nodes.name")
    output("    nodes.hmin")
    output("    nodes.hmax")
    output("    nodes.nodesPerSmoothingScale")

#-------------------------------------------------------------------------------
# Initial conditions
#-------------------------------------------------------------------------------
if geometry == "1d":
    from DistributeNodes import distributeNodesInRange1d
    distributeNodesInRange1d([(nodes1, nx1, rho1, (x0, x1)),
                              (nodes2, nx2, rho2, (x1, x2))],
                             nPerh = nPerh)
elif geometry == "2d":
    from GenerateNodeDistribution2d import GenerateNodeDistribution2d
    from PeanoHilbertDistributeNodes import distributeNodes2d
    ny1 = int((y1 - y0)/(x1 - x0)*nx1 + 0.5)
    ny2 = int((y1 - y0)/(x2 - x1)*nx2 + 0.5)
    gen1 = GenerateNodeDistribution2d(nx1, ny1, rho1, "lattice",
                                      xmin = (x0, y0),
                                      xmax = (x1, y1),
                                      nNodePerh = nPerh,
                                      SPH = True)
    gen2 = GenerateNodeDistribution2d(nx2, ny2, rho2, "lattice",
                                      xmin = (x1, y0),
                                      xmax = (x2, y1),
                                      nNodePerh = nPerh,
                                      SPH = True)
    distributeNodes2d((nodes1, gen1),
                      (nodes2, gen2))

output("nodes1.numNodes")
output("nodes2.numNodes")

# Set node velocities
for nodes in nodeSet:
    vel = nodes.velocity()
    for i in range(nodes.numInternalNodes):
        vel[i].x = vx0

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase()
output("db")
for nodes in nodeSet:
    output("db.appendNodeList(nodes)")
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
hydro = SPH(dataBase = db,
            W = WT)
output("hydro")
output("hydro.kernel")
output("hydro.PiKernel")
output("hydro.cfl")
output("hydro.compatibleEnergyEvolution")
output("hydro.densityUpdate")
output("hydro.XSPH")

packages = [hydro]

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
inplane = Plane(Vector(x0), Vector(1))
outplane = Plane(Vector(x2), Vector(-1))
inflow = InflowOutflowBoundary(db, inplane)
outflow = InflowOutflowBoundary(db, outplane)
bcs = [inflow, outflow]
packages += bcs

if geometry == "2d":
    yplane0 = Plane(Vector(x0, y0), Vector(0,  1))
    yplane1 = Plane(Vector(x0, y1), Vector(0, -1))
    ybc = PeriodicBoundary(yplane0, yplane1)
    bcs += [ybc]

for p in packages:
    for bc in bcs:
        p.appendBoundary(bc)

#-------------------------------------------------------------------------------
# Construct an integrator.
#-------------------------------------------------------------------------------
integrator = CheapSynchronousRK2Integrator(db)
for p in packages:
    integrator.appendPhysicsPackage(p)
del p
integrator.lastDt = dtMin
integrator.dtMin = dtMin
integrator.verbose = True
output("integrator")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.dtGrowth")
output("integrator.rigorousBoundaries")
output("integrator.updateBoundaryFrequency")
output("integrator.domainDecompositionIndependent")
output("integrator.cullGhostNodes")
output("integrator.verbose")

#-------------------------------------------------------------------------------
# Make the problem controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            vizBaseName = "inflowOutflow-%s" % geometry,
                            redistributeStep = 10,
                            vizDir = vizDir,
                            vizTime = 1e8,
                            vizStep = vizStep)
output("control")

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
if not steps is None:
    control.step(steps)
else:
    control.advance(goalTime)

#-------------------------------------------------------------------------------
# Plotting
#-------------------------------------------------------------------------------
if graphics:
    from SpheralMatplotlib import *
    rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotState(db)

    if geometry == "2d":
        posplot = plotNodePositions2d(db, plotGhosts=True)
