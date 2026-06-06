#ATS:test(SELF, "--graphics False", label="Voronoi hourglass filter tests -- 1-D (serial)")
#-------------------------------------------------------------------------------
# Test out the Voronoi hourglass filtering algorithm.
#-------------------------------------------------------------------------------
from Spheral2d import *
from SpheralTestUtilities import *
from generateMesh import generateLineMesh
from SpheralVoronoiSiloDump import dumpPhysicsState

import random, numpy, Gnuplot

#-------------------------------------------------------------------------------
# Command line parameters.
#-------------------------------------------------------------------------------
commandLine(n = 100,
            rho0 = 2.0,
            rhoSlope = -1.0,      # drho/dr

            x0 = 0.0,
            y0 = 0.0,
            x1 = 1.0,
            y1 = 1.0,

            nPerh = 1.01,

            gammaGas = 5.0/3.0,
            mu = 1.0,
            hmin = 1e-10,
            hmax = 1.0,
            hourglassOrder = 1,
            hourglassLimiter = 1,

            iterations = 100,
            graphics = True,
            silodir = "VoronoiHourglassControTest",
            silofile = "visit",
            )

#-------------------------------------------------------------------------------
# Material.
#-------------------------------------------------------------------------------
eos = GammaLawGasMKS(gammaGas, mu)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
WT = TableKernel(BSplineKernel(), 1000)

#-------------------------------------------------------------------------------
# Make the NodeLists.
#-------------------------------------------------------------------------------
nodes = makeFluidNodeList("nodes1", eos, 
                          hmin = hmin,
                          hmax = hmax,
                          nPerh = nPerh)

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
M = (x1 - x0)*(y1 - y0)*rho0
nodes.numInternalNodes = n
mass = nodes.mass()
pos = nodes.positions()
H = nodes.Hfield()

nxcells = 2**16
ncells = nxcells**2
dx = (x1 - x0)/nxcells
dy = (y1 - y0)/nxcells
H0 = SymTensor(1.0/1.0, 0.0,
               0.0, 1.0/1.0)
cells = random.sample(range(ncells), n)
i = 0
for cell in cells:
    ix = cell % nxcells
    iy = cell // nxcells
    mass[i] = M/n
    pos[i] = Vector((ix + 0.5)*dx, (iy + 0.5)*dy)
    H[i] = H0
    i += 1
print(i)
assert i == n

## nx = int(sqrt(float(n)))
## assert nx*nx == n
## dx = (x1 - x0)/nx
## dy = (y1 - y0)/nx
## for iy in xrange(nx):
##     for ix in xrange(nx):
##         i = ix + nx*iy
##         mass[i] = M/n
##         pos[i] = Vector((ix + 0.5)*dx, (iy + 0.5)*dy)
##         H[i] = H0

nodes.neighbor().updateNodes()

def setRho():
    pos = nodes.positions()
    rho = nodes.massDensity()
    for i in range(nodes.numInternalNodes):
        rho[i] = rho0 + rhoSlope*pos[i].magnitude()
    return

setRho()

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase()
output("db")
output("db.appendNodeList(nodes)")
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Construct an hourglass control object.
#-------------------------------------------------------------------------------
hg = VoronoiHourglassControl(WT, hourglassOrder, hourglassLimiter)
output("hg")
output("hg.order")
output("hg.limiter")
packages = vector_of_Physics()
packages.append(hg)

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
xPlane0 = Plane(Vector(x0, 0.0), Vector( 1.0,  0.0))
xPlane1 = Plane(Vector(x1, 0.0), Vector(-1.0,  0.0))
yPlane0 = Plane(Vector(0.0, y0), Vector( 0.0,  1.0))
yPlane1 = Plane(Vector(0.0, y1), Vector( 0.0, -1.0))
xbc0 = ReflectingBoundary(xPlane0)
xbc1 = ReflectingBoundary(xPlane1)
ybc0 = ReflectingBoundary(yPlane0)
ybc1 = ReflectingBoundary(yPlane1)
bcs = vector_of_Boundary()
for bc in [xbc0, ybc0]:
    bcs.append(bc)
    hg.appendBoundary(bc)

#-------------------------------------------------------------------------------
# A helpful method for enforcing boundaries.
#-------------------------------------------------------------------------------
def enforceBoundaries(state, derivs):
    nodes.numGhostNodes = 0
    nodes.neighbor().updateNodes()
    for bc in bcs:
        bc.setGhostNodes(nodes)
        bc.finalizeGhostBoundary()
        nodes.neighbor().updateNodes()
    for f in [nodes.mass(), nodes.massDensity()]:
        for bc in bcs:
            bc.applyGhostBoundary(f)
    hg.applyGhostBoundaries(state, derivs)
    return

#-------------------------------------------------------------------------------
# Provide a method of convering the H's.
#-------------------------------------------------------------------------------
def convergeH():
    W = TableKernel(BSplineKernel(), 1000)
    meth = SPHSmoothingScale()
    iterateIdealH(db, bcs, W, meth,
                  tolerance = 1.0e-4)
    return
convergeH()

#-------------------------------------------------------------------------------
# Iteratively let the hourglass control adjust the point positions, and lets
# see where we converge to.
#-------------------------------------------------------------------------------
mass = nodes.mass()
pos = nodes.positions()
H = nodes.Hfield()
rho = nodes.massDensity()

packages = vector_of_Physics()
packages.append(hg)

for iter in range(iterations):
    setRho()
    state = State()
    for f in (pos, mass, rho, H):
        state.enroll(f)
    derivs = StateDerivatives(db, packages)
    hg.registerState(db, state)
    hg.registerDerivatives(db, derivs)
    enforceBoundaries(state, derivs)
    state.update(derivs, 1.0, 0.0, 1.0)
    hg.finalize(0.0, 0.0, db, state, derivs)
    convergeH()
    if graphics:
        for f in [nodes.mass(), nodes.massDensity(), nodes.Hfield()]:
            state.enroll(f)
        for fl in [hg.mask(), hg.gradRho(), hg.A(), hg.B(), hg.gradA(), hg.gradB(),
                   hg.weight()]:
            state.enrollFieldList(fl)
        dumpPhysicsState(state, silofile, silodir, 
                         currentTime = iter,
                         currentCycle = iter)
