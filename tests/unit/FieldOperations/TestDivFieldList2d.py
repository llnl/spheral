from Spheral import *
from SpheralTestUtilities import *

################################################################################
nx1, nx2 = 25, 25
nx = nx1 + nx2
ny = 50
n1, n2 = nx1*ny, nx2*ny
rho1, rho2 = 1.0, 1.0
m1, m2 = rho1/(n1 + n2), rho2/(n1 + n2)
eps0 = 0.0
epsMultiplier = 1.0
vr = -1.0
x0, x1 = 0.0, 1.0
y0, y1 = 0.0, 1.0

gamma = 2.0
mu = 1.0

neighborSearchType = 3 # GatherScatter
numGridLevels = 10
topGridCellSize = 0.5
origin = Vector2d(0.0)

################################################################################
title('2-D FieldList Test')

eos = GammaLawGasMKS2d(gamma, mu)

nodes1 = SphNodeList2d(eos, n1)
nodes2 = SphNodeList2d(eos, n2)
output('nodes1.numNodes')
output('nodes2.numNodes')

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

# Set node positions
dx, dy = (x1 - x0)/nx, (y1 - y0)/ny
for iy in xrange(ny):
    for ix in xrange(nx1):
        nodeID = ix + nx1*iy
        nodes1.positions()[nodeID] = Vector2d((ix + 0.5)*dx, (iy + 0.5)*dy)
for iy in xrange(ny):
    for ix in xrange(nx2):
        nodeID = ix + nx2*iy
        nodes2.positions()[nodeID] = Vector2d((ix + 0.5)*dx + 0.5, (iy + 0.5)*dy)

# Set node masses
nodes1.setMass(ScalarField2d(nodes1, m1))
nodes2.setMass(ScalarField2d(nodes2, m2))

# Set node specific thermal energies
for nodeID in xrange(nodes1.numNodes):
    nodes1.specificThermalEnergy()[nodeID] = eps0 + epsMultiplier*(nodes1.positions()[nodeID].magnitude2())
for nodeID in xrange(nodes2.numNodes):
    nodes2.specificThermalEnergy()[nodeID] = eps0 + epsMultiplier*(nodes2.positions()[nodeID].magnitude2())

# Set node velocities
for nodeID in xrange(nodes1.numNodes):
    unit = nodes1.positions()[nodeID].unitVector()
    nodes1.velocity()[nodeID] = unit*vr
for nodeID in xrange(nodes2.numNodes):
    unit = nodes2.positions()[nodeID].unitVector()
    nodes2.velocity()[nodeID] = unit*vr

# Set the smoothing scales.
h = 1.0/(2.01*dx)
H = SymTensor2d(h, 0.0, 0.0, h)
nodes1.setHfield(SymTensorField2d(nodes1, H))
nodes2.setHfield(SymTensorField2d(nodes2, H))

# Set the mass densities.
nodes1.setMassDensity(ScalarField2d(nodes1, rho1))
nodes2.setMassDensity(ScalarField2d(nodes2, rho2))

nodes1.updateWeight()
nodes2.updateWeight()

# Update the nodal weights
updateTimer = SpheralTimer('Update node weights')
updateTimer.start()
nodes1.updateWeight()
nodes2.updateWeight()
updateTimer.stop()
updateTimer.printStatus()

# Construct neighbor objects and associate them with the node lists.
neighborTimer = SpheralTimer('Neighbor initialization.')
neighborTimer.start()
neighbor1 = NestedGridNeighbor2d(nodes1,
                                 neighborSearchType,
                                 numGridLevels,
                                 topGridCellSize,
                                 origin,
                                 kernelExtent)
nodes1.registerNeighbor(neighbor1)

neighbor2 = NestedGridNeighbor2d(nodes2,
                                 neighborSearchType,
                                 numGridLevels,
                                 topGridCellSize,
                                 origin,
                                 kernelExtent)
nodes2.registerNeighbor(neighbor2)
neighborTimer.stop()
neighborTimer.printStatus()

# Construct a DataBase to hold our two node lists
db = DataBase2d()
output('db')
output('db.appendNodeList(nodes1)')
output('db.appendNodeList(nodes2)')
output('db.numNodeLists')
output('db.numFluidNodeLists')

# Plot the velocity field.
from SpheralGnuPlotUtilities import *
p = plotVelocityField2d(db, title )

################################################################################
# Test the divergence function of the FieldList
title('2-D FieldList.divergence() test')

velocity = db.fluidVelocity
fluidPosition = db.fluidPosition
fluidWeight = db.fluidWeight
fluidMass = db.fluidMass
fluidDensity = db.fluidMassDensity
fluidHfield = db.fluidHfield
divTimer = SpheralTimer('Calculate divergence')
divTimer.start()
##divVelocity = divergenceVector2d(velocity, fluidPosition, fluidWeight, fluidMass,
##                                   fluidDensity, fluidHfield, WT)
divVelocity = divergencePairWise2d(velocity, fluidPosition, fluidWeight, fluidMass,
                                   fluidDensity, fluidHfield, WT)
divTimer.stop()
divTimer.printStatus()

# Plot the results of the div operation.
p1 = plotFieldList(divVelocity,
                   xFunction = "%s.magnitude()",
                   plotStyle = "points",
                   winTitle = "Radial divergence without boundaries")

# Plot the analytic solution, which in this case is 1/r.
xans = []
yans = []
for nodeList in db.nodeLists:
    for r in nodeList.positions()[:nodeList.numInternalNodes]:
        xans.append(r.magnitude())
        yans.append(-1.0/xans[-1])
multiSort(xans, yans)
ansData = Gnuplot.Data(xans, yans,
                       with = "lines",
                       title = "Analytic solution")
p1.replot(ansData)

##################################################################################
### Try applying boundary conditions to the pressure, and redo the divergence.
##title('2-D FieldList.divergence() test with reflecting boundaries.')

##xPlane0 = Plane2d((0.0, 0.0), (1.0, 0.0))
##xPlane1 = Plane2d((1.0, 0.0), (-1.0, 0.0))
##yPlane0 = Plane2d((0.0, 0.0), (0.0, 1.0))
##yPlane1 = Plane2d((0.0, 1.0), (0.0, -1.0))
##xbc0 = ReflectingBoundary2d(xPlane0)
##xbc1 = ReflectingBoundary2d(xPlane1)
##ybc0 = ReflectingBoundary2d(yPlane0)
##ybc1 = ReflectingBoundary2d(yPlane1)
##output('xbc0, xbc1')
##output('ybc0, ybc1')

##bcTimer = SpheralTimer('Apply boundary conditions')
##bcTimer.start()
##for bc in [xbc0, ybc0]:
##    output('bc.setGhostNodes(nodes1)')
##    output('bc.applyVectorGhostBoundary(velocity[0])')
##    output('bc.applySymTensorGhostBoundary(fluidHfield[0])')
##    output('bc.applyScalarGhostBoundary(fluidWeight[0])')

##    output('bc.setGhostNodes(nodes2)')
##    output('bc.applyVectorGhostBoundary(velocity[1])')
##    output('bc.applySymTensorGhostBoundary(fluidHfield[1])')
##    output('bc.applyScalarGhostBoundary(fluidWeight[1])')
##bcTimer.stop()
##bcTimer.printStatus()

##neighborTimer.start()
##output('nodes1.neighbor.updateNodes()')
##output('nodes2.neighbor.updateNodes()')
##neighborTimer.stop()
##neighborTimer.printStatus()

##divTimer.start()
##divVelocity = dummy.divergencePairWise2d(velocity, fluidPosition, fluidWeight, fluidMass,
##                                         fluidDensity, fluidHfield, WT)
##divTimer.stop()
##output('divVelocity[:]')
##divTimer.printStatus()

### Plot the results of the divergence operation.
##p2 = plotFieldList(divVelocity,
##                   xFunction = "%s.magnitude()",
##                   plotStyle = "points",
##                   winTitle = "Radial divergence with reflecting boundaries")

### Plot the analytic solution, which in this case is 1/r.
##p2.replot(ansData)
