from Spheral import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *

################################################################################
nx, ny = 50, 50
n1 = nx*ny
rho1 = 1.0
m1 = 1.0*rho1/n1
eps0 = 0.0
epsMultiplier = 1.0
v1 = Vector2d(0.0)
x0, x1 = 0.0, 1.0
y0, y1 = 0.0, 1.0

gamma = 2.0
mu = 1.0

neighborSearchType = Neighbor2d.NeighborSearchType.GatherScatter
numGridLevels = 10
topGridCellSize = 2.0
origin = Vector2d(0.0, 0.0)

################################################################################
def epsFunction(x):
    return eps0 + epsMultiplier*x**-1

################################################################################
title('2-D FieldList Test')

eos = GammaLawGasMKS2d(gamma, mu)

nodes1 = SphNodeList2d(eos, n1)
output('nodes1.numNodes')

W = BSplineKernel2d()
#W = W4SplineKernel2d()
#W = GaussianKernel2d()
#W = SuperGaussianKernel2d()
#W = PiGaussianKernel2d(1.0)
output('W')
kernelExtent = W.kernelExtent()

# Set the table kernel for the FieldList divergence.
WT = TableKernel2d()
WT.setTableData(W, 100)
output('WT')

import random
generator = random.Random()
ranMag = 0.1

dx = (x1 - x0)/nx
dy = (y1 - y0)/ny
for iy in xrange(ny):
    for ix in xrange(ny):
        nodes1.positions()[ix + iy*nx] = Vector2d((ix + 0.5)*dx, (iy + 0.5)*dy)
        #nodes1.positions[i] = (i + 0.5 + ranMag*generator.uniform(-1.0,1.0))*dx1

nodes1.setMass(ScalarField2d(nodes1, m1))

for nodeID in xrange(nodes1.numNodes):
    nodes1.specificThermalEnergy()[nodeID] = epsFunction(nodes1.positions()[nodeID].magnitude())
nodes1.setVelocity(VectorField2d(nodes1, v1))

h1 = 1.0/(2.01*dx)
H0 = SymTensor2d(h1, 0.0,
                 0.0, h1)
nodes1.setHfield(SymTensorField2d(nodes1, H0))

neighbor1 = NestedGridNeighbor2d(nodes1,
                                 neighborSearchType,
                                 numGridLevels,
                                 topGridCellSize,
                                 origin,
                                 kernelExtent)
nodes1.registerNeighbor(neighbor1)

nodes1.setMassDensity(ScalarField2d(nodes1, rho1))

nodes1.updateWeight()

db = DataBase2d()
output('db')
output('db.appendNodeList(nodes1)')
output('db.numNodeLists')
output('db.numFluidNodeLists')

pressure = db.fluidPressure

fluidPosition = db.fluidPosition
fluidWeight = db.fluidWeight
fluidMass = db.fluidMass
fluidRho = db.fluidMassDensity
fluidHfield = db.fluidHfield

# Enforce boundary conditions.
xPlane = Plane2d(Vector2d(0.0, 0.0), Vector2d(1.0, 0.0))
yPlane = Plane2d(Vector2d(0.0, 0.0), Vector2d(0.0, 1.0))
xbc = ReflectingBoundary2d(xPlane)
ybc = ReflectingBoundary2d(yPlane)
output('xbc.setGhostNodes(db)')
output('ybc.setGhostNodes(db)')
output('nodes1.neighbor().updateNodes()')

output('xbc.applyScalarFieldListGhostBoundary(pressure)')
output('xbc.applyScalarFieldListGhostBoundary(fluidMass)')
output('xbc.applyScalarFieldListGhostBoundary(fluidRho)')
output('xbc.applyScalarFieldListGhostBoundary(fluidWeight)')
output('xbc.applySymTensorFieldListGhostBoundary(fluidHfield)')
output('xbc.applyScalarFieldListGhostBoundary(fluidWeight)')

output('ybc.applyScalarFieldListGhostBoundary(pressure)')
output('ybc.applyScalarFieldListGhostBoundary(fluidMass)')
output('ybc.applyScalarFieldListGhostBoundary(fluidRho)')
output('ybc.applyScalarFieldListGhostBoundary(fluidWeight)')
output('ybc.applySymTensorFieldListGhostBoundary(fluidHfield)')
output('ybc.applyScalarFieldListGhostBoundary(fluidWeight)')

################################################################################
# Generate the analytic answer
import Gnuplot
xans = array([0.0]*n1)
i = 0
for nodeList in db.nodeLists:
    for r in nodeList.positions()[:nodeList.numInternalNodes]:
        xans[i] = r.magnitude()
        i = i + 1
yans = -1.0/xans**2
multiSort(xans, yans)
ansData = Gnuplot.Data(xans, yans, with='lines', title='Analytic solution')

################################################################################
# Test the gradient function of the FieldList
title('2-D FieldList.gradient() test')

gradPressure = gradientScalar2d(pressure,
                               fluidPosition,
                               fluidWeight,
                               fluidMass,
                               fluidRho,
                               fluidHfield,
                               WT)

# Plot the results of the grad operation.
gradPressurePlot = plotFieldList(gradPressure,
                                 xFunction = '%s.magnitude()',
                                 yFunction = '-%s.magnitude()',
                                 plotStyle = 'points',
                                 winTitle='Standard SPH gradient')
gradPressurePlot.replot(ansData)

gradPerssurePlot2 = plotVectorField2d(db, gradPressure,
                                      vectorMultiplier = 1e-3,
                                      title = "Standard SPH gradient")

################################################################################
# Test the pair wise form of the gradient.
title('2-D FieldList gradientPairWise test')

gradPressurePW = gradientPairWiseScalar2d(pressure,
                                          fluidPosition, fluidMass,
                                          fluidMass, fluidRho, fluidHfield, WT)

# Plot the results of the grad operation.
gradPressurePlotPW = plotFieldList(gradPressurePW,
                                   xFunction = '%s.magnitude()',
                                   yFunction = '-%s.magnitude()',
                                   plotStyle = 'points',
                                   winTitle='Pairwise gradient')
gradPressurePlotPW.replot(ansData)

gradPerssurePlotPW2 = plotVectorField2d(db, gradPressurePW,
                                        vectorMultiplier = 1e-3,
                                        title = "Pairwise gradient")

