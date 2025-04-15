from Spheral import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from GenerateNodeDistribution3d import *

################################################################################
nx, ny, nz = 10, 10, 10
rmin, rmax = 0.0, 1.0
xmin = Vector3d(0, 0, 0)
xmax = Vector3d(1, 1, 1)
rho0 = 1.0
eps0 = 0.0
r0, r1, r2 = 0.0, 0.0, 1.0
nPerh = 1.51

v0 = 0.0
v1 = -1.0

gamma = 2.0
mu = 1.0

neighborSearchType = 3 # GatherScatter
numGridLevels = 20
topGridCellSize = 4.0
origin = Vector3d(0, 0, 0)

################################################################################
def vFunction(position):
    r = position.magnitude()
    rhat = position.unitVector()
    if r < r1:
        return Vector3d(v0*rhat.x, v0*rhat.y, v0*rhat.z)
    else:
        return Vector3d(v1*rhat.x, v1*rhat.y, v1*rhat.z)

################################################################################
def gradDivVelFunction(r):
    if r < r1:
        return 2.0*abs(v0)/(r*r)
    else:
        return 2.0*abs(v1)/(r*r)

################################################################################
title('3-D FieldList grad div test')

eos = GammaLawGasMKS3d(gamma, mu)

# Create the cylindrical node distribution.
nodes1 = SphNodeList3d(eos)
generator = GenerateNodeDistribution3d(nx, ny, nz, rho0, 'lattice',
                                       xmin = xmin,
                                       xmax = xmax,
                                       rmin = rmin,
                                       rmax = rmax,
                                       nNodePerh = nPerh)
n1 = generator.globalNumNodes()
nodes1.numInternalNodes = n1
for i in xrange(n1):
    nodes1.positions[i] = generator.position(i)
    h = (generator.Htensor(i).Determinant())**(1.0/3.0)
    nodes1.Hfield[i] = SymTensor3d(h, 0, 0,
                                   0, h, 0,
                                   0, 0, h)
    nodes1.mass[i] = generator.m[i]
    nodes1.massDensity[i] = rho0
output('nodes1.numNodes')

# Build a kernel.
#W = BSplineKernel3d()
#W = W4SplineKernel3d()
#W = GaussianKernel3d()
#W = SuperGaussianKernel3d()
#W = PiGaussianKernel3d(1.0)
#W = QuarticSplineKernel3d()
W = NBSplineKernel3d(5)
output('W')
kernelExtent = W.kernelExtent

# Set the table kernel for the FieldList divergence.
WT = TableKernel3d()
WT.setTableData(W, 1000)
output('WT')

for nodeID in xrange(nodes1.numNodes):
    nodes1.velocity[nodeID] = vFunction(nodes1.positions[nodeID])

neighbor1 = NestedGridNeighbor3d(nodes1,
                                 neighborSearchType,
                                 numGridLevels,
                                 topGridCellSize,
                                 origin,
                                 kernelExtent)
nodes1.neighbor = neighbor1

#output('nodes.updateMassDensity(W)')
output('nodes1.updateWeight()')

db = DataBase3d()
output('db')
output('db.appendNodeList(nodes1)')
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

# Create boundary conditions.  We need at least this much to create the initial
# mass density field.
xPlane0 = Plane3d((0.0, 0.0, 0.0), (1.0, 0.0, 0.0))
yPlane0 = Plane3d((0.0, 0.0, 0.0), (0.0, 1.0, 0.0))
zPlane0 = Plane3d((0.0, 0.0, 0.0), (0.0, 0.0, 1.0))
xbc0 = ReflectingBoundary3d(xPlane0)
ybc0 = ReflectingBoundary3d(yPlane0)
zbc0 = ReflectingBoundary3d(zPlane0)
boundaryConditions = [xbc0, ybc0, zbc0]

# If this is a parallel run, then build communication boundary condition.
try:
    import mpi
    if mpi.procs > 1:
        domainbc = NestedGridDistributedBoundary3d()
        boundaryConditions.insert(0, domainbc)
except:
    print "Serial test."
    pass

# Enforce the boundary conditions.
for bc in boundaryConditions:
    bc.setGhostNodes(db)
    nodes1.neighbor.updateNodes()
    bc.applyFieldListGhostBoundary(fluidWeight)
    bc.applyFieldListGhostBoundary(fluidMass)
    bc.applyFieldListGhostBoundary(fluidRho)
    bc.applyFieldListGhostBoundary(fluidHfield)
db.updateFluidMassDensity()
nodes1.numGhostNodes = 0
nodes1.neighbor.updateNodes()
nodes1.updateWeight()
for bc in boundaryConditions:
    bc.setGhostNodes(db)
    nodes1.neighbor.updateNodes()
    bc.applyFieldListGhostBoundary(fluidWeight)
    bc.applyFieldListGhostBoundary(fluidMass)
    bc.applyFieldListGhostBoundary(fluidRho)
    bc.applyFieldListGhostBoundary(fluidHfield)
for i in xrange(nodes1.firstGhostNode, nodes1.numNodes):
    nodes1.velocity[i] = vFunction(nodes1.positions[i])

################################################################################
# Generate the analytic answer.
import Gnuplot
xl = []
for nodeList in db.nodeLists:
    for r in nodeList.positions[:nodeList.numInternalNodes]:
        xl.append(r.magnitude())
if mpi:
    xll = mpi.allreduce(xl, mpi.SUM)
    xl = xll
xl.sort()
xans = array(xl)
yans = gradDivVelFunction(xans)
ansData = Gnuplot.Data(xans, yans, with='lines', title='Analytic answer')

################################################################################
# Plot the direct, successive first derivatives against the known solution.
dummy = FieldFunctions()

directGradDivVel = dummy.gradDivVectorFieldList3d(velocity,
                                                  fluidPosition,
                                                  fluidWeight,
                                                  fluidMass,
                                                  fluidRho,
                                                  fluidHfield,
                                                  WT,
                                                  boundaryConditions)

plotDirect1 = plotVectorField3d(db, directGradDivVel, vectorMultiplier=1e-3,
                               title = 'Direct second derivative of velocity.')
plotDirect = plotFieldList(directGradDivVel,
                           xFunction = '%s.magnitude()',
                           yFunction = '%s.magnitude()',
                           plotStyle = 'points',
                           winTitle = 'Direct second derivative of velocity.')
plotDirect.replot(ansData)

# ##################################################################################
# # Plot the simple second derivative method against the known solution.
# simpleGradDivVel = dummy.gradDivVectorFieldListSimple3d(velocity,
#                                                         fluidPosition,
#                                                         fluidWeight,
#                                                         fluidMass,
#                                                         fluidRho,
#                                                         fluidHfield,
#                                                         WT)

# ##plotSimple = plotVectorField3d(db, simpleGradDivVel, vectorMultiplier=1e-4,
# ##                               title = 'Simple second derivative of velocity.')
# plotSimple1 = plotVectorField3d(db, simpleGradDivVel, vectorMultiplier=1e-4,
#                                title = 'Simple second derivative of velocity.')
# plotSimple = plotFieldList(simpleGradDivVel,
#                            xFunction = '%s.magnitude()',
#                            yFunction = '%s.magnitude()',
#                            plotStyle = 'points',
#                            winTitle = 'Simple second derivative of velocity.')
# plotSimple.replot(ansData)

##################################################################################
# Plot the golden second derivative method against the known solution.
goldenGradDivVel = dummy.gradDivVectorFieldListGolden3d(velocity,
                                                        fluidPosition,
                                                        fluidWeight,
                                                        fluidMass,
                                                        fluidRho,
                                                        fluidHfield,
                                                        WT)

integrator = PredictorCorrectorIntegrator3d(db)
integrator.timerSummary()

plotGolden = plotFieldList(goldenGradDivVel,
                           xFunction = '%s.magnitude()',
                           yFunction = '%s.magnitude()',
                           plotStyle = 'points',
                           winTitle = 'Golden second derivative of velocity.')
plotGolden.replot(ansData)

####################################################################################
### Plot the mash second derivative method against the known solution.
##mashGradDivVel = dummy.gradDivVectorFieldListMash3d(velocity,
##                                                        fluidPosition,
##                                                        fluidWeight,
##                                                        fluidMass,
##                                                        fluidRho,
##                                                        fluidHfield,
##                                                        WT)

##plotMash1 = plotVectorField3d(db, mashGradDivVel, vectorMultiplier=1e-4,
##                               title = 'Mash second derivative of velocity.')
##plotMash = plotFieldList(mashGradDivVel,
##                           xFunction = '%s.magnitude()',
##                           yFunction = '%s.magnitude()',
##                           plotStyle = 'points',
##                           winTitle = 'Mash second derivative of velocity.')
##plotMash.replot(ansData)
