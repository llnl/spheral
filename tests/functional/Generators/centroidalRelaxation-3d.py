import os, shutil, mpi
from Spheral3d import *
from SpheralTestUtilities import *
from centroidalRelaxNodes import *
from GenerateNodeDistribution3d import *
from siloPointmeshDump import *
from fieldStatistics import *

title("3-D test of centroidal relaxation.")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(KernelConstructor = NBSplineKernel,
            order = 7,
            nPerh = 1.01,
            hmin = 1e-5,
            hmax = 1.0,
            hminratio = 1.0,

            # The initial density function coeffients:  rho(x) = a + bx*x + by*y + bz*z
            a = 1.0,
            bx = 0.0,
            by = 0.0,
            bz = 0.0,

            # Initial geometry
            nx = 20,
            ny = 20,
            nz = 20,
            x0 = 0.0,
            x1 = 1.0,
            y0 = 0.0,
            y1 = 1.0,
            z0 = 0.0,
            z1 = 1.0,
            ranfrac = 0.25,
            seed = 14892042,

            # Material properties
            gamma = 5.0/3.0,
            mu = 1.0,

            # Simulation control
            iterations = 100,
            tol = 1.0e-3,

            graphics = True,
            baseName = "centroidal_relaxation_3d",
            )

#-------------------------------------------------------------------------------
# Our density and gradient methods.
#-------------------------------------------------------------------------------
def rhofunc(posi):
   return a + bx*posi.x + by*posi.y + bz*posi.z

def gradrhofunc(posi):
   return Vector(bx, by, bz)

#-------------------------------------------------------------------------------
# Create a random number generator.
#-------------------------------------------------------------------------------
import random
random.seed(seed)

#-------------------------------------------------------------------------------
# Material properties.
#-------------------------------------------------------------------------------
eos = GammaLawGasMKS(gamma, mu)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
if KernelConstructor==NBSplineKernel:
    Wbase = NBSplineKernel(order)
else:
    Wbase = KernelConstructor()
WT = TableKernel(Wbase, 1000)
output("WT")

#-------------------------------------------------------------------------------
# Make the NodeList.
#-------------------------------------------------------------------------------
nodes = makeFluidNodeList("nodes", eos, 
                          hmin = hmin,
                          hmax = hmax,
                          hminratio = hminratio,
                          nPerh = nPerh,
                          kernelExtent = WT.kernelExtent)
    
output("nodes")
output("nodes.hmin")
output("nodes.hmax")
output("nodes.nodesPerSmoothingScale")

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
generator = GenerateNodeDistribution3d(nx, ny, nz, rhofunc, "lattice",
                                       xmin = (x0, y0, z0),
                                       xmax = (x1, y1, z1),
                                       nNodePerh = nPerh,
                                       SPH = True)

if mpi.procs > 1:
    from VoronoiDistributeNodes import distributeNodes3d
else:
    from DistributeNodes import distributeNodes3d

distributeNodes3d((nodes, generator))
output("mpi.reduce(nodes.numInternalNodes, mpi.MIN)")
output("mpi.reduce(nodes.numInternalNodes, mpi.MAX)")
output("mpi.reduce(nodes.numInternalNodes, mpi.SUM)")

# Randomly jitter the node positions.
dx = (x1 - x0)/nx
dy = (y1 - y0)/ny
dz = (z1 - z0)/nz
pos = nodes.positions()
for i in range(nodes.numInternalNodes):
   pos[i].x += ranfrac * dx * random.uniform(-1.0, 1.0)
   pos[i].y += ranfrac * dy * random.uniform(-1.0, 1.0)
   pos[i].z += ranfrac * dz * random.uniform(-1.0, 1.0)

# Initialize the mass and densities.
m = nodes.mass()
rho = nodes.massDensity()
for i in range(nodes.numNodes):
   rho[i] = rhofunc(pos[i])
   m[i] = rho[i]*dx*dy*dz

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase()
output("db")
output("db.appendNodeList(nodes)")
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
xPlane0 = Plane(Vector(x0, y0, z0), Vector( 1,  0,  0))
yPlane0 = Plane(Vector(x0, y0, z0), Vector( 0,  1,  0))
zPlane0 = Plane(Vector(x0, y0, z0), Vector( 0,  0,  1))

xPlane1 = Plane(Vector(x1, y1, z1), Vector(-1,  0,  0))
yPlane1 = Plane(Vector(x1, y1, z1), Vector( 0, -1,  0))
zPlane1 = Plane(Vector(x1, y1, z1), Vector( 0,  0, -1))

xbc0 = ReflectingBoundary(xPlane0)
ybc0 = ReflectingBoundary(yPlane0)
zbc0 = ReflectingBoundary(zPlane0)

xbc1 = ReflectingBoundary(xPlane1)
ybc1 = ReflectingBoundary(yPlane1)
zbc1 = ReflectingBoundary(zPlane1)

boundaries = []

#-------------------------------------------------------------------------------
# Call the centroidal relaxer.
#-------------------------------------------------------------------------------
# Report the initial mass matching.
print("Initial mass (min, max, avg, std dev) : ", fieldStatistics(m))

bcpoints = vector_of_Vector()
for p in [(x0, y0, z0), (x1, y0, z0), (x1, y1, z0), (x0, y1, z0),
          (x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1)]:
   bcpoints.append(Vector(*p))
boundary = Polyhedron(bcpoints)
vol, surfacePoint = centroidalRelaxNodes(nodeListsAndBounds = [(nodes, boundary)],
                                         W = WT,
                                         rho = rhofunc,
                                         gradrho = gradrhofunc,
                                         maxIterations = iterations,
                                         boundaries = boundaries,
                                         fracTol = tol,
                                         tessellationFileName = baseName)

# Report the final mass matching.
print("Final mass (min, max, avg, std dev) : ", fieldStatistics(m))

#-------------------------------------------------------------------------------
# Plot the final state.
#-------------------------------------------------------------------------------
if graphics:
    from SpheralGnuPlotUtilities import *

    # rho
    rhoPlot = plotFieldList(db.fluidMassDensity,
                            winTitle = "Density",
                            plotStyle = "points",
                            colorNodeLists = False,
                            plotGhosts = False)

    # mass
    massPlot = plotFieldList(db.fluidMass,
                             winTitle = "Mass",
                             plotStyle = "points",
                             colorNodeLists = False,
                             plotGhosts = False)

