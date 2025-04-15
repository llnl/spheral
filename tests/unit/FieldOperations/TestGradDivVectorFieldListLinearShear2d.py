from Spheral import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from GenerateNodeDistribution2d import *

################################################################################
nx1, ny1 = 40, 20
nx2, ny2 = 40, 20
nPerh = 4.01
seed = 'lattice'
rho0 = 1.0
eps0 = 1.0

offset = 0.2

vx0, vxslope = 1.0, 1.0
vy0 = 0.0
dvy = 0.0 #0.005
perturbLength = 0.1

xmax = nx1/float(ny1 + ny2) * 1.0
rmin1, rmax1 = (0.0, 0.0), (xmax, 0.5)
rmin2, rmax2 = (0.0, 0.5), (xmax, 1.0)

gamma = 2.0
mu = 1.0

neighborSearchType = 3 # GatherScatter
#neighborSearchType = 2 # Scatter
#neighborSearchType = 1 # Gather
numGridLevels = 10
topGridCellSize = 2.0
origin = Vector2d(0.0, 0.0)

################################################################################
title('2-D FieldList grad div test')

eos = GammaLawGasMKS2d(gamma, mu)

# Set node positions, masses, and H's for this domain.
nodes1 = SphNodeList2d(0, eos)
nodes2 = SphNodeList2d(0, eos)
from DistributeNodes import distributeNodes2d
generator1 = GenerateNodeDistribution2d(nx1, ny1, rho0, seed,
                                        rmin = rmin1,
                                        rmax = rmax1,
                                        nNodePerh = nPerh)
generator2 = GenerateNodeDistribution2d(nx2, ny2, rho0, seed,
                                        rmin = rmin2,
                                        rmax = rmax2,
                                        nNodePerh = nPerh)
n1 = generator1.globalNumNodes()
n2 = generator2.globalNumNodes()
nodeTuples = [(nodes1, n1, generator1),
              (nodes2, n2, generator2)]
nodeInfo = distributeNodes2d(nodeTuples)
output('nodes1.numInternalNodes')
output('nodes2.numInternalNodes')
assert len(nodeInfo[nodes1]['globalNodeListID']) == nodes1.numInternalNodes
assert len(nodeInfo[nodes2]['globalNodeListID']) == nodes2.numInternalNodes

for tup in nodeTuples:
    nodes = tup[0]
    gen = tup[2]
    for i in xrange(nodes.numInternalNodes):
        globalID = nodeInfo[nodes]['globalNodeListID'][i]
        nodes.mass[i] = gen.mass(globalID)
        Hi = gen.Htensor(globalID)
        h = sqrt(Hi.Determinant())
        Hi = SymTensor2d(h, 0, 0, h)
        nodes.Hfield[i] = Hi

# Set node specific thermal energies
nodes1.specificThermalEnergy[:] = [eps0]*nodes1.numNodes
nodes2.specificThermalEnergy[:] = [eps0]*nodes2.numNodes

# Set node velocities
for nodes in [nodes1, nodes2]:
    for nodeID in xrange(nodes.numNodes):
        nodes.velocity[nodeID] = (vx0 + vxslope*nodes.positions[nodeID].y, vy0)

# If requested, perturb the velocity across the interface.
xlength = rmax1[0] - rmin1[0]
ylength = rmax2[1] - rmin1[1]
assert xlength > 0.0
assert ylength > 0.0
for nodeID in xrange(nodes1.numNodes):
    r = nodes1.positions[nodeID]
    dvyi = dvy*(sin(2.0*pi*(r.x - rmin1[0])/xlength)*
                exp(-((r.y - rmax1[1])/perturbLength)**2))
#                cos(pi*(r.y - rmin1[1])/ylength - 0.5))
    nodes1.velocity[nodeID].y += dvyi
for nodeID in xrange(nodes2.numNodes):
    r = nodes2.positions[nodeID]
    dvyi = dvy*(sin(2.0*pi*(r.x - rmin2[0])/xlength)*
                exp(-((r.y - rmax1[1])/perturbLength)**2))
#                cos(pi*(r.y - rmin2[1])/ylength - 0.5))
    nodes2.velocity[nodeID].y += dvyi

# Offset the nodes as prescribed.
for nodes in [nodes1, nodes2]:
    for i in xrange(nodes.numNodes):
        nodes.positions[i].x += nodes.velocity[i].x*offset

# If we're integrating for mass density, 'ya have to initialize the density.
nodes1.massDensity[:] = [rho0]*nodes1.numNodes
nodes2.massDensity[:] = [rho0]*nodes2.numNodes

# Build a kernel.
#W = BSplineKernel2d()
#W = W4SplineKernel2d()
W = GaussianKernel2d()
#W = SuperGaussianKernel2d()
#W = PiGaussianKernel2d(1.0)
output('W')
kernelExtent = W.kernelExtent

# Set the table kernel for the FieldList divergence.
WT = TableKernel2d()
WT.setTableData(W, 1000)
output('WT')

neighbor1 = NestedGridNeighbor2d(nodes1,
                                 neighborSearchType,
                                 numGridLevels,
                                 topGridCellSize,
                                 origin,
                                 kernelExtent)
nodes1.neighbor = neighbor1

neighbor2 = NestedGridNeighbor2d(nodes2,
                                 neighborSearchType,
                                 numGridLevels,
                                 topGridCellSize,
                                 origin,
                                 kernelExtent)
nodes2.neighbor = neighbor2

output('nodes1.updateWeight()')
output('nodes2.updateWeight()')

db = DataBase2d()
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

fluidPosition = db.fluidPosition
fluidWeight = db.fluidWeight
fluidMass = db.fluidMass
fluidRho = db.fluidMassDensity
fluidHfield = db.fluidHfield

# Create boundary conditions.  We will enforce reflecting boundaries at
# y=0 and y=1, and a periodic boundary wrapping around x=0 and x=1.
xPlane0 = Plane2d(rmin1, (1.0, 0.0))
xPlane1 = Plane2d(rmax2, (-1.0, 0.0))
yPlane0 = Plane2d(rmin1, (0.0, 1.0))
yPlane1 = Plane2d(rmax2, (0.0, -1.0))
xbc0 = PeriodicBoundary2d(xPlane0, xPlane1)
ybc0 = ReflectingBoundary2d(yPlane0)
ybc1 = ReflectingBoundary2d(yPlane1)
boundaryConditions = [xbc0, ybc0, ybc1]

# Enforce the boundary conditions.
for bc in boundaryConditions:
    bc.setViolationNodes(db)
for bc in boundaryConditions:
    bc.setGhostNodes(db)
    bc.applyFieldListGhostBoundary(fluidWeight)
    bc.applyFieldListGhostBoundary(fluidMass)
    bc.applyFieldListGhostBoundary(fluidRho)
    bc.applyFieldListGhostBoundary(velocity)
    bc.applyFieldListGhostBoundary(fluidHfield)
for nodes in [nodes1, nodes2]:
    nodes.neighbor.updateNodes()
    for i in xrange(nodes.firstGhostNode, nodes.numNodes):
        nodes.velocity[i] = (vx0 + vxslope*nodes.positions[i].y, vy0)

################################################################################
# Plot the node positions.
rPlot = plotNodePositions2d(db)

################################################################################
# Generate the analytic answer.
import Gnuplot
n = n1 + n2
xans = array([0.0]*n)
yans = array([0.0]*n)
i = 0
for nodeList in db.nodeLists:
    for r in nodeList.positions[:nodeList.numInternalNodes]:
        xans[i] = r.y
        yans[i] = 0.0
        i = i + 1
ansData = Gnuplot.Data(xans, yans, with='lines', title='Analytic answer')

################################################################################
# Plot the direct, successive first derivatives against the known solution.
dummy = FieldFunctions()

directGradDivVel = dummy.gradDivVectorFieldList2d(velocity,
                                                  fluidPosition,
                                                  fluidWeight,
                                                  fluidMass,
                                                  fluidRho,
                                                  fluidHfield,
                                                  WT,
                                                  boundaryConditions)
x = max([qq.magnitude() for qq in directGradDivVel[0][:] + directGradDivVel[0][:]])
print mpi.reduce(x, mpi.MAX)

plotDirect = plotVectorField2d(db, directGradDivVel, vectorMultiplier=10.0,
                               title = 'Direct second derivative of velocity.')
plotDirect2 = plotFieldList(directGradDivVel,
                            xFunction = '%s.y',
                            yFunction = '%s.y',
                            plotStyle = 'points',
                            winTitle = 'Direct second derivative of velocity.')
plotDirect2.replot(ansData)

# ##################################################################################
# # Plot the simple second derivative method against the known solution.
# simpleGradDivVel = dummy.gradDivVectorFieldListSimple2d(velocity,
#                                                         fluidPosition,
#                                                         fluidWeight,
#                                                         fluidMass,
#                                                         fluidRho,
#                                                         fluidHfield,
#                                                         WT)
# x = max([qq.magnitude() for qq in simpleGradDivVel[0][:] + simpleGradDivVel[0][:]])
# print mpi.reduce(x, mpi.MAX)

# plotSimple = plotVectorField2d(db, simpleGradDivVel, vectorMultiplier=10.0,
#                                title = 'Simple second derivative of velocity.')
# plotSimple2 = plotFieldList(simpleGradDivVel,
#                             xFunction = '%s.y',
#                             yFunction = '%s.y',
#                             plotStyle = 'points',
#                             winTitle = 'Simple second derivative of velocity.')
# plotSimple2.replot(ansData)

##################################################################################
# Plot the golden second derivative method against the known solution.
goldenGradDivVel = dummy.gradDivVectorFieldListGolden2d(velocity,
                                                        fluidPosition,
                                                        fluidWeight,
                                                        fluidMass,
                                                        fluidRho,
                                                        fluidHfield,
                                                        WT)
x = max([qq.magnitude() for qq in goldenGradDivVel[0][:] + goldenGradDivVel[0][:]])
print mpi.reduce(x, mpi.MAX)

plotGolden = plotVectorField2d(db, goldenGradDivVel, vectorMultiplier=10.0,
                               title = 'Golden second derivative of velocity.')
plotGolden2 = plotFieldList(goldenGradDivVel,
                            xFunction = '%s.y',
                            yFunction = '%s.y',
                            plotStyle = 'points',
                            winTitle = 'Golden second derivative of velocity.')
plotGolden2.replot(ansData)

# ##################################################################################
# # Plot the mash second derivative method against the known solution.
# mashGradDivVel = dummy.gradDivVectorFieldListMash2d(velocity,
#                                                         fluidPosition,
#                                                         fluidWeight,
#                                                         fluidMass,
#                                                         fluidRho,
#                                                         fluidHfield,
#                                                         WT)
# x = max([qq.magnitude() for qq in mashGradDivVel[0][:] + mashGradDivVel[0][:]])
# print mpi.reduce(x, mpi.MAX)

# plotMash = plotVectorField2d(db, mashGradDivVel, vectorMultiplier=10.0,
#                                title = 'Mash second derivative of velocity.')
# plotMash2 = plotFieldList(mashGradDivVel,
#                             xFunction = '%s.y',
#                             yFunction = '%s.y',
#                             plotStyle = 'points',
#                             winTitle = 'Mash second derivative of velocity.')
# plotMash2.replot(ansData)
