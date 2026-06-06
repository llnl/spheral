from Spheral import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *

################################################################################
n1, n2 = 50, 50
rho1, rho2 = 1.0, 1.0
m1, m2 = 0.5*rho1/n1, 0.5*rho2/n2
eps0 = 0.0
epsMultiplier = 1.0
v1, v2 = Vector1d(1.0), Vector1d(-1.0)
x0, x1 = 0.0, 1.0

gamma = 2.0
mu = 1.0

neighborSearchType = Neighbor2d.NeighborSearchType.GatherScatter
numGridLevels = 10
topGridCellSize = 2.0
origin = Vector1d(0.0)

################################################################################
def epsFunction(x):
    return eps0 + epsMultiplier*x**3

################################################################################
title('1-D FieldList Test')

eos = GammaLawGasMKS1d(gamma, mu)

nodes1 = SphNodeList1d(eos, n1)
nodes2 = SphNodeList1d(eos, n2)
output('nodes1.numNodes')
output('nodes2.numNodes')

W = BSplineKernel1d()
#W = W4SplineKernel1d()
#W = GaussianKernel1d()
#W = SuperGaussianKernel1d()
#W = PiGaussianKernel1d(1.0)
output('W')
kernelExtent = W.kernelExtent()

# Set the table kernel for the FieldList divergence.
WT = TableKernel1d()
WT.setTableData(W, 100)
output('WT')

import random
generator = random.Random()
ranMag = 0.1

dx1 = 0.5*(x1 - x0)/n1
dx2 = 0.5*(x1 - x0)/n2
for i in xrange(n1):
    nodes1.positions()[i] = Vector1d((i + 0.5)*dx1)
    #nodes1.positions()[i] = Vector1d((i + 0.5 + ranMag*generator.uniform(-1.0,1.0))*dx1)
for i in xrange(n2):
    nodes2.positions()[i] = Vector1d(0.5 + (i + 0.5)*dx2)
    #nodes2.positions()[i] = Vector1d(0.5 + (i + 0.5 + ranMag*generator.uniform(-1.0,1.0))*dx2)

nodes1.setMass(ScalarField1d(nodes1, m1))
nodes2.setMass(ScalarField1d(nodes2, m2))

for nodeID in xrange(nodes1.numNodes):
    nodes1.specificThermalEnergy()[nodeID] = epsFunction(nodes1.positions()[nodeID].x)
for nodeID in xrange(nodes2.numNodes):
    nodes2.specificThermalEnergy()[nodeID] = epsFunction(nodes2.positions()[nodeID].x)
nodes1.setVelocity(VectorField1d(nodes1, v1))
nodes2.setVelocity(VectorField1d(nodes2, v2))

H1 = SymTensor1d(1.0/(2.01*dx1))
H2 = SymTensor1d(1.0/(2.01*dx2))
nodes1.setHfield(SymTensorField1d(nodes1, H1))
nodes2.setHfield(SymTensorField1d(nodes2, H2))

neighbor1 = NestedGridNeighbor1d(nodes1,
                                 neighborSearchType,
                                 numGridLevels,
                                 topGridCellSize,
                                 origin,
                                 kernelExtent)
nodes1.registerNeighbor(neighbor1)

neighbor2 = NestedGridNeighbor1d(nodes2,
                                 neighborSearchType,
                                 numGridLevels,
                                 topGridCellSize,
                                 origin,
                                 kernelExtent)
nodes2.registerNeighbor(neighbor2)

nodes1.setMassDensity(ScalarField1d(nodes1, rho1))
nodes2.setMassDensity(ScalarField1d(nodes2, rho2))

nodes1.updateWeight()
nodes2.updateWeight()

db = DataBase1d()
output('db')
output('db.appendNodeList(nodes1)')
output('db.appendNodeList(nodes2)')
output('db.numNodeLists')
output('db.numFluidNodeLists')

pressure = db.fluidPressure

fluidPosition = db.fluidPosition
fluidWeight = db.fluidWeight
fluidMass = db.fluidMass
fluidRho = db.fluidMassDensity
fluidHfield = db.fluidHfield

################################################################################
# Generate the analytic answer, grad(P) = 2.0*x
import Gnuplot
xans = array([0.0]*(n1 + n2))
i = 0
for nodeList in db.nodeLists:
    for r in nodeList.positions()[:nodeList.numInternalNodes]:
        xans[i] = r.x
        i = i + 1
yans = 3.0*xans**2 #2.0*xans
ansData = Gnuplot.Data(xans, yans, with='lines', title='Analytic solution')

################################################################################
# Test the gradient function of the FieldList
title('1-D FieldList.gradient() test')

gradPressure = gradientScalar1d(pressure,
                                fluidPosition,
                                fluidWeight,
                                fluidMass,
                                fluidRho,
                                fluidHfield,
                                WT)

# Plot the results of the grad operation.
gradPressurePlot = plotFieldList(gradPressure,
                                 yFunction = '%s.x',
                                 plotStyle = 'points',
                                 winTitle='1-D FieldList gradient test -- no Boundary conditions')
gradPressurePlot.replot(ansData)

################################################################################
# Try applying boundary conditions to the pressure, and redo the gradient.
title('1-D FieldList.gradient() test with reflecting boundaries.')

xPlane0, xPlane1 = Plane1d(Vector1d(0.0), Vector1d(1.0)), Plane1d(Vector1d(1.0), Vector1d(-1.0))
xbc0, xbc1 = ReflectingBoundary1d(xPlane0), ReflectingBoundary1d(xPlane1)
output('xbc0, xbc1')
output('xbc0.setGhostNodes(db)')
output('xbc1.setGhostNodes(db)')

output('xbc0.applyScalarFieldListGhostBoundary(pressure)')
output('xbc0.applyScalarFieldListGhostBoundary(fluidMass)')
output('xbc0.applyScalarFieldListGhostBoundary(fluidRho)')
output('xbc0.applySymTensorFieldListGhostBoundary(fluidHfield)')
output('xbc0.applyScalarFieldListGhostBoundary(fluidWeight)')

output('xbc1.applyScalarFieldListGhostBoundary(pressure)')
output('xbc1.applyScalarFieldListGhostBoundary(fluidMass)')
output('xbc1.applyScalarFieldListGhostBoundary(fluidRho)')
output('xbc1.applySymTensorFieldListGhostBoundary(fluidHfield)')
output('xbc1.applyScalarFieldListGhostBoundary(fluidWeight)')

output('nodes1.neighbor().updateNodes()')
output('nodes2.neighbor().updateNodes()')

gradPressure = gradientScalar1d(pressure,
                                fluidPosition,
                                fluidWeight,
                                fluidMass,
                                fluidRho,
                                fluidHfield,
                                WT)

# Plot the results of the grad operation.
gradPressurePlot2 = plotFieldList(gradPressure,
                                  yFunction = '%s.x',
                                  plotStyle = 'points',
                                  winTitle='Reflecting Boundary conditions')
gradPressurePlot2.replot(ansData)

################################################################################
# Now extend the right answer over the ghosts, and replot the gradient
title('1-D FieldList.gradient() test with reflecting boundaries, and forced ghost values.')
for nodeID in xrange(nodes1.firstGhostNode, nodes1.numNodes):
    pressure[0][nodeID] = (gamma - 1.0)*rho1*epsFunction(nodes1.positions()[nodeID].x)
for nodeID in xrange(nodes2.firstGhostNode, nodes2.numNodes):
    pressure[1][nodeID] = (gamma - 1.0)*rho2*epsFunction(nodes2.positions()[nodeID].x)

gradPressure2 = gradientScalar1d(pressure,
                                 fluidPosition,
                                 fluidWeight,
                                 fluidMass,
                                 fluidRho,
                                 fluidHfield,
                                 WT)

# Plot the results of the grad operation.
gradPressurePlot3 = plotFieldList(gradPressure2,
                                  yFunction = '%s.x',
                                  plotStyle = 'points',
                                  winTitle='Analytic Boundary conditions')
gradPressurePlot3.replot(ansData)

################################################################################
# Test the pair wise form of the gradient.
title('1-D FieldList gradientPairWise test with reflecting boundaries, and forced ghost values.')
gradPressurePW = gradientPairWiseScalar1d(pressure,
                                          fluidPosition, fluidMass,
                                          fluidMass, fluidRho, fluidHfield, WT)

# Plot the results of the grad operation.
gradPressurePlotPW = plotFieldList(gradPressurePW,
                                   yFunction = '%s.x',
                                   plotStyle = 'points',
                                   winTitle='Pairwise gradient')
gradPressurePlotPW.replot(ansData)

################################################################################
# Test the MASH form of the gradient.
title('1-D FieldList gradientMash test with reflecting boundaries, and forced ghost values.')
gradPressureMash = gradientMashScalar1d(pressure,
                                        fluidPosition,
                                        fluidWeight,
                                        fluidHfield,
                                        WT)

# Plot the results of the grad operation.
gradPressurePlotMash = plotFieldList(gradPressureMash,
                                     yFunction = '%s.x',
                                     plotStyle = 'points',
                                     winTitle='MASH gradient')
gradPressurePlotMash.replot(ansData)
