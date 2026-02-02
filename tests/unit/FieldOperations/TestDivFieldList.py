from Spheral import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *

################################################################################
n1, n2 = 50, 50
rho1, rho2 = 1.0, 1.0
m1, m2 = 0.5*rho1/n1, 0.5*rho2/n2
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
def vFunction(x):
    return v0 + vMultiplier*x**3

################################################################################
title('1-D FieldList Divergence Test')

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
kernelExtent = W.kernelExtent

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
    #nodes1.positions[i] = (i + 0.5)*dx1
    nodes1.positions[i] = (i + 0.5 + ranMag*generator.uniform(-1.0,1.0))*dx1
for i in xrange(n2):
    #nodes2.positions[i] = 0.5 + (i + 0.5)*dx2
    nodes2.positions[i] = 0.5 + (i + 0.5 + ranMag*generator.uniform(-1.0,1.0))*dx2

nodes1.mass[:] = [m1]*nodes1.numNodes
nodes2.mass[:] = [m2]*nodes2.numNodes

for nodeID in xrange(nodes1.numNodes):
    nodes1.velocity[nodeID].x = vFunction(nodes1.positions[nodeID].x)
for nodeID in xrange(nodes2.numNodes):
    nodes2.velocity[nodeID].x = vFunction(nodes2.positions[nodeID].x)

h1 = 1.0/(2.01*dx1)
h2 = 1.0/(2.01*dx2)
for H in nodes1.Hfield:
    H.xx = h1
for H in nodes2.Hfield:
    H.xx = h2

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

nodes1.updateWeight()
nodes2.updateWeight()

db = DataBase1d()
output('db')
output('db.appendNodeList(nodes1)')
output('db.appendNodeList(nodes2)')
output('db.numNodeLists')
output('db.numFluidNodeLists')

velocity = db.fluidVelocity

fluidPosition = db.fluidPosition
fluidWeight = db.fluidWeight
fluidMass = db.fluidMass
fluidRho = db.fluidMassDensity
fluidHfield = db.fluidHfield

################################################################################
# Generate the analytic answer
import Gnuplot
xans = array([0.0]*(n1 + n2))
i = 0
for nodeList in db.nodeLists:
    for r in nodeList.positions:
        xans[i] = r.x
        i = i + 1
yans = 3.0*xans**2 #2.0*xans
ansData = Gnuplot.Data(xans, yans, with='lines', title='Analytic solution')

################################################################################
# Test the divergence function of the FieldList
title('1-D FieldList.divergence() test')
dummy = FieldFunctions()

divVelocity = dummy.divergence1d(velocity, fluidPosition, fluidWeight,
                                 fluidMass, fluidRho, fluidHfield, WT)

# Plot the results of the divergence operation.
divVelocityPlot = plotFieldList(divVelocity,
                                plotStyle = 'points',
                                winTitle='1-D FieldList divergence test -- no Boundary conditions')
divVelocityPlot.replot(ansData)

################################################################################
# Test the pair wise form of the divergence.
title('1-D FieldList divergencePairWise test.')
divVelocityPW = dummy.divergencePairWise1d(velocity,
                                           fluidPosition, fluidMass,
                                           fluidMass, fluidRho, fluidHfield, WT)

# Plot the results of the divergence operation.
divVelocityPlotPW = plotFieldList(divVelocityPW,
                                  plotStyle = 'points',
                                  winTitle='Pairwise divergence')
divVelocityPlotPW.replot(ansData)

##################################################################################
### Test the MASH form of the divergence.
##title('1-D FieldList divergence test.')
##divVelocityMash = dummy.divergenceMash1d(velocity,
##                                         fluidPosition,
##                                         fluidWeight,
##                                         fluidHfield,
##                                         WT)

### Plot the results of the divergence operation.
##divVelocityPlotMash = plotFieldList(divVelocityMash,
##                                    plotStyle = 'points',
##                                    winTitle='MASH divergence')
##divVelocityPlotMash.replot(ansData)
