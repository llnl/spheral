from Spheral import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *

################################################################################
n1, n2 = 50, 50
rho1, rho2 = 1.0, 1.0
m1, m2 = 0.5*rho1/n1, 0.5*rho2/n2
eps0 = 0.0

v0 = 0.0
vMultiplier = 1.0

x0, x1 = 0.0, 1.0

gamma = 2.0
mu = 1.0

neighborSearchType = 3 # GatherScatter
numGridLevels = 10
topGridCellSize = 2.0
origin = Vector1d(0.0)

################################################################################
def vxFunction(x):
#    return v0 + vMultiplier*x**3
#    return v0 + vMultiplier*x**4
    return v0 + vMultiplier*x**5

################################################################################
title('1-D FieldList Test')

eos = GammaLawGasMKS1d(gamma, mu)

nodes1 = SphNodeList1d(eos, n1)
nodes2 = SphNodeList1d(eos, n2)
output('nodes1.numNodes')
output('nodes2.numNodes')

#W = NBSplineKernel1d(5)
W = BSplineKernel1d()
#W = W4SplineKernel1d()
#W = GaussianKernel1d()
#W = SuperGaussianKernel1d()
#W = PiGaussianKernel1d(1.0)
#W = NSincPolynomialKernel1d(5)
#W = QuarticSplineKernel1d()
output('W')
kernelExtent = W.kernelExtent

# Set the table kernel for the FieldList divergence.
WT = TableKernel1d()
WT.setTableData(W, 1000)
output('WT')

import random
generator = random.Random()
ranMag = 0.1

dx1 = 0.5*(x1 - x0)/n1
dx2 = 0.5*(x1 - x0)/n2
for i in xrange(n1):
    nodes1.positions[i] = (i + 0.5)*dx1
    #nodes1.positions[i] = (i + 0.5 + ranMag*generator.uniform(-1.0,1.0))*dx1
for i in xrange(n2):
    nodes2.positions[i] = 0.5 + (i + 0.5)*dx2
    #nodes2.positions[i] = 0.5 + (i + 0.5 + ranMag*generator.uniform(-1.0,1.0))*dx2

output('nodes1.positions[:]')
output('nodes2.positions[:]')

nodes1.mass[:] = [m1]*nodes1.numNodes
nodes2.mass[:] = [m2]*nodes2.numNodes

for nodeID in xrange(nodes1.numNodes):
    nodes1.velocity[nodeID].x = vxFunction(nodes1.positions[nodeID].x)
for nodeID in xrange(nodes2.numNodes):
    nodes2.velocity[nodeID].x = vxFunction(nodes2.positions[nodeID].x)

h1 = 1.0/(2.01*dx1)
h2 = 1.0/(2.01*dx2)
for H in nodes1.Hfield:
    H.xx = h1
for H in nodes2.Hfield:
    H.xx = h2
output('nodes1.Hfield[:]')
output('nodes2.Hfield[:]')

neighbor1 = NestedGridNeighbor1d(nodes1,
                                 neighborSearchType,
                                 numGridLevels,
                                 topGridCellSize,
                                 origin,
                                 kernelExtent)
nodes1.neighbor = neighbor1

neighbor2 = NestedGridNeighbor1d(nodes2,
                                 neighborSearchType,
                                 numGridLevels,
                                 topGridCellSize,
                                 origin,
                                 kernelExtent)
nodes2.neighbor = neighbor2

nodes1.massDensity[:] = [rho1]*nodes1.numNodes
nodes2.massDensity[:] = [rho2]*nodes2.numNodes
output('nodes1.massDensity[:]')
output('nodes2.massDensity[:]')

#output('nodes.updateMassDensity(W)')
output('nodes1.updateWeight()')
output('nodes2.updateWeight()')

db = DataBase1d()
output('db')
output('db.appendNodeList(nodes1)')
output('db.appendNodeList(nodes2)')
output('db.numNodeLists')
output('db.numFluidNodeLists')

output('db.globalMass[:]')
output('db.fluidMass[:]')
output('db.globalPosition[:]')
output('db.fluidMassDensity[:]')
output('db.fluidSpecificThermalEnergy[:]')
output('db.fluidVelocity[:]')
output('db.fluidWeight[:]')
output('db.fluidHfield[:]')

velocity = db.fluidVelocity
##output('velocity[:]')
##output('velocity[0][:]')
##output('velocity[1][:]')

fluidPosition = db.fluidPosition
fluidWeight = db.fluidWeight
fluidMass = db.fluidMass
fluidRho = db.fluidMassDensity
fluidHfield = db.fluidHfield

# Create boundary conditions.  We need at least this much to create the initial
# mass density field.
xPlane0 = Plane1d((x0), ( 1.0))
xPlane1 = Plane1d((x1), (-1.0))
xbc0 = ReflectingBoundary1d(xPlane0)
xbc1 = ReflectingBoundary1d(xPlane1)
boundaryConditions = [xbc0, xbc1]

# Enforce the boundary conditions.
for i in xrange(10):
    for bc in boundaryConditions:
        bc.setGhostNodes(db)
        bc.applyFieldListGhostBoundary(fluidWeight)
        bc.applyFieldListGhostBoundary(fluidMass)
        bc.applyFieldListGhostBoundary(fluidRho)
        bc.applyFieldListGhostBoundary(fluidHfield)
    db.updateFluidMassDensity()
    for nodes in [nodes1, nodes2]:
        nodes.numGhostNodes = 0
        nodes.neighbor.updateNodes()
        nodes.updateWeight()
for bc in boundaryConditions:
    bc.setGhostNodes(db)
    bc.applyFieldListGhostBoundary(fluidWeight)
    bc.applyFieldListGhostBoundary(fluidMass)
    bc.applyFieldListGhostBoundary(fluidRho)
    bc.applyFieldListGhostBoundary(fluidHfield)
for nodes in [nodes1, nodes2]:
    nodes.neighbor.updateNodes()
    for i in xrange(nodes.firstGhostNode, nodes.numNodes):
        nodes.velocity[i].x = vxFunction(nodes.positions[i].x)

################################################################################
# Generate the analytic answer, grad(P) = 2.0*x
import Gnuplot
xans = array([0.0]*(n1 + n2))
i = 0
for nodeList in db.nodeLists:
    for r in nodeList.positions[:nodeList.numInternalNodes]:
        xans[i] = r.x
        i = i + 1
#yans = 6.0*vMultiplier*xans
#yans = 12.0*vMultiplier*xans**2
yans = 20.0*vMultiplier*xans**3
ansData = Gnuplot.Data(xans, yans, with='lines', title='Analytic answer')

################################################################################
# Plot the direct, successive first derivatives against the known solution.
dummy = FieldFunctions()

# directGradDivVel = dummy.gradDivVectorFieldList1d(velocity,
#                                                   fluidPosition,
#                                                   fluidWeight,
#                                                   fluidMass,
#                                                   fluidRho,
#                                                   fluidHfield,
#                                                   WT,
#                                                   boundaryConditions)

# plotDirect = plotFieldList(directGradDivVel, yFunction='%s.x',
#                            plotStyle='points',
#                            winTitle = 'Direct second derivative of velocity.')
# plotDirect.replot(ansData)

# ################################################################################
# # Plot the simple second derivative method against the known solution.
# simpleGradDivVel = dummy.gradDivVectorFieldListSimple1d(velocity,
#                                                         fluidPosition,
#                                                         fluidWeight,
#                                                         fluidMass,
#                                                         fluidRho,
#                                                         fluidHfield,
#                                                         WT)

# plotSimple = plotFieldList(simpleGradDivVel, yFunction='%s.x',
#                            plotStyle='points',
#                            winTitle = 'Simple second derivative of velocity.')
# plotSimple.replot(ansData)

################################################################################
# Plot the golden second derivative method against the known solution.
goldenGradDivVel = dummy.gradDivVectorFieldListGolden1d(velocity,
                                                        fluidPosition,
                                                        fluidMass,
                                                        fluidMass,
                                                        fluidRho,
                                                        fluidHfield,
                                                        WT)

plotGolden = plotFieldList(goldenGradDivVel, yFunction='%s.x',
                           plotStyle='points',
                           winTitle = 'Golden second derivative of velocity.')
plotGolden.replot(ansData)

# ################################################################################
# # Plot the golden2 second derivative method against the known solution.
# golden2GradDivVel = dummy.gradDivVectorFieldListGolden21d(velocity,
#                                                         fluidPosition,
#                                                         fluidWeight,
#                                                         fluidMass,
#                                                         fluidRho,
#                                                         fluidHfield,
#                                                         WT)

# plotGolden2 = plotFieldList(golden2GradDivVel, yFunction='%s.x',
#                            plotStyle='points',
#                            winTitle = 'Golden2 second derivative of velocity.')
# plotGolden2.replot(ansData)

# ################################################################################
# # Plot the mash second derivative method against the known solution.
# mashGradDivVel = dummy.gradDivVectorFieldListMash1d(velocity,
#                                                     fluidPosition,
#                                                     fluidWeight,
#                                                     fluidMass,
#                                                     fluidRho,
#                                                     fluidHfield,
#                                                     WT)

# plotMash = plotFieldList(mashGradDivVel, yFunction='%s.x',
#                          plotStyle='points',
#                          winTitle = 'Mash second derivative of velocity.')
# plotMash.replot(ansData)

################################################################################
# Plot the pair wise direct second derivative method against the known solution.
pwGradDivVel = dummy.gradDivVectorFieldListPairWise1d(velocity,
                                                      fluidPosition,
                                                      fluidWeight,
                                                      fluidMass,
                                                      fluidRho,
                                                      fluidHfield,
                                                      WT)

plotPW = plotFieldList(pwGradDivVel, yFunction='%s.x',
                       plotStyle='points',
                       winTitle = 'Pair wise direct second derivative of velocity.')
plotPW.replot(ansData)
plotPW.refresh()

