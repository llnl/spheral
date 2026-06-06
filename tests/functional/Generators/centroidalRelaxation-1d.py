import os, shutil
from Spheral1d import *
from SpheralTestUtilities import *
from centroidalRelaxNodes import *
from fieldStatistics import *

title("1-D test of centroidal relaxation.")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(KernelConstructor = NBSplineKernel,
            order = 7,
            nPerh = 1.01,
            hmin = 1e-5,
            hmax = 1.0,

            # The initial density function coeffients:  rho(x) = a + b*x + c*x^2
            a = 1.0,
            b = 0.0,
            c = 0.0,

            # Initial geometry
            nx = 100,
            x0 = 0.0,
            x1 = 1.0,
            ranfrac = 0.25,
            seed = 14892042,

            # Material properties
            gamma = 5.0/3.0,
            mu = 1.0,

            # Simulation control
            iterations = 100,
            tol = 1.0e-3,

            graphics = True,
            )

#-------------------------------------------------------------------------------
# Our density and gradient methods.
#-------------------------------------------------------------------------------
def rhofunc(posi):
   if isinstance(posi, Vector):
      xi = posi.x
   else:
      xi = posi
   return a + b*xi + c*xi*xi

def gradrhofunc(posi):
   xi = posi.x
   return Vector(b + 2.0*c*xi)

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
                           nPerh = nPerh,
                           kernelExtent = WT.kernelExtent)
    
output("nodes")
output("nodes.hmin")
output("nodes.hmax")
output("nodes.nodesPerSmoothingScale")

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
from DistributeNodes import distributeNodesInRange1d
distributeNodesInRange1d([(nodes, nx, a, (x0, x1))],
                         nPerh = nPerh)
output("nodes.numNodes")

# Randomly jitter the node positions.
dx = (x1 - x0)/nx
pos = nodes.positions()
for i in range(nodes.numInternalNodes):
   pos[i].x += ranfrac * dx * random.uniform(-1.0, 1.0)

# Initialize the mass and densities.
m = nodes.mass()
rho = nodes.massDensity()
for i in range(nodes.numNodes):
   rho[i] = rhofunc(pos[i])
   m[i] = rho[i]*dx

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
xPlane0 = Plane(Vector(x0), Vector(1.0))
xPlane1 = Plane(Vector(x1), Vector(-1.0))
xbc0 = ReflectingBoundary(xPlane0)
xbc1 = ReflectingBoundary(xPlane1)

boundaries = [] # [xbc0, xbc1]

#-------------------------------------------------------------------------------
# Call the centroidal relaxer.
#-------------------------------------------------------------------------------
# Report the initial mass matching.
print("Initial mass (min, max, avg, std dev) : ", fieldStatistics(m))

vol, surfacePoint = centroidalRelaxNodes(nodeListsAndBounds = [(nodes, Box1d(Vector(0.5*(x0 + x1)), 0.5*(x1 - x0)))],
                                         W = WT,
                                         rho = rhofunc,
                                         gradrho = gradrhofunc,
                                         maxIterations = iterations,
                                         boundaries = boundaries,
                                         fracTol = tol)

# Report the final mass matching.
print("Final mass (min, max, avg, std dev) : ", fieldStatistics(m))

#-------------------------------------------------------------------------------
# Plot the final state.
#-------------------------------------------------------------------------------
if graphics:
    from SpheralGnuPlotUtilities import *

    xprof = mpi.allreduce([x.x for x in pos.internalValues()], mpi.SUM)
    xprof.sort()

    # rho
    rhoPlot = plotFieldList(db.fluidMassDensity,
                            winTitle = "Density",
                            plotStyle = "points",
                            colorNodeLists = False,
                            plotGhosts = False)
    rhoAns = Gnuplot.Data(xprof, [rhofunc(xi) for xi in xprof],
                            with_ = "lines",
                            title = "Solution",
                            inline = True)
    rhoPlot.replot(rhoAns)

    # mass
    massPlot = plotFieldList(db.fluidMass,
                             winTitle = "Mass",
                             plotStyle = "points",
                             colorNodeLists = False,
                             plotGhosts = False)

    # delta points
    deltaData = Gnuplot.Data([0.5*(xprof[i] + xprof[i+1]) for i in range(len(xprof)-1)],
                             [xprof[i+1] - xprof[i] for i in range(len(xprof)-1)],
                             with_ = "linespoints",
                             title = "Delta spacing",
                             inline = True)
    deltaPlot = generateNewGnuPlot()
    deltaPlot.replot(deltaData)
