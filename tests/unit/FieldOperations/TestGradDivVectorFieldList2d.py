from Spheral import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from GenerateNodeDistribution2d import *

################################################################################
nRadial, nTheta = 40, 40
theta = 0.5*pi
rho0 = 1.0
eps0 = 0.0
r0, r1, r2 = 0.0, 0.0, 1.0
nPerh = 2.01

v0 = 0.0
v1 = -1.0

gamma = 2.0
mu = 1.0

neighborSearchType = 3 # GatherScatter
numGridLevels = 20
topGridCellSize = 4.0
origin = Vector2d(0.0, 0.0)

################################################################################
def vFunction(position):
    r = position.magnitude()
    rhat = position.unitVector()
    if r < r1:
        return Vector2d(v0*rhat.x, v0*rhat.y)
    else:
        return Vector2d(v1*rhat.x, v1*rhat.y)

################################################################################
def gradDivVelFunction(position):
    r2 = position.magnitude2()
    if r < r1:
        return abs(v0)/r2
    else:
        return abs(v1)/r2

################################################################################
title('2-D FieldList grad div test')

eos = GammaLawGasMKS2d(gamma, mu)

# Create the cylindrical node distribution.
nodes1 = SphNodeList2d(eos)
output('nodes1.numNodes')
generator = GenerateNodeDistribution2d(nRadial, nTheta, rho0, 'optimal',
                                       rmin = 0.0,
                                       rmax = 1.0,
                                       theta = theta,
                                       nNodePerh = nPerh)
n1 = generator.globalNumNodes()
nodes1.numInternalNodes = n1
for i in xrange(n1):
    nodes1.positions[i].x = generator.x[i]
    nodes1.positions[i].y = generator.y[i]
    h = sqrt(generator.H[i].Determinant())
    nodes1.Hfield[i] = SymTensor2d(h, 0, 0, h)
    nodes1.mass[i] = generator.m[i]
    nodes1.massDensity[i] = rho0

# Build a kernel.
W = BSplineKernel2d()
#W = W4SplineKernel2d()
#W = GaussianKernel2d()
#W = SuperGaussianKernel2d()
#W = PiGaussianKernel2d(1.0)
#W = QuarticSplineKernel2d()
#W = NBSplineKernel2d(5)
output('W')
kernelExtent = W.kernelExtent

# Set the table kernel for the FieldList divergence.
WT = TableKernel2d()
WT.setTableData(W, 1000)
output('WT')

for nodeID in xrange(nodes1.numNodes):
    nodes1.velocity[nodeID] = vFunction(nodes1.positions[nodeID])

neighbor1 = NestedGridNeighbor2d(nodes1,
                                 neighborSearchType,
                                 numGridLevels,
                                 topGridCellSize,
                                 origin,
                                 kernelExtent)
nodes1.neighbor = neighbor1

#output('nodes.updateMassDensity(W)')
output('nodes1.updateWeight()')

db = DataBase2d()
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
xPlane0 = Plane2d((0.0, 0.0), (1.0, 0.0))
yPlane0 = Plane2d((0.0, 0.0), (0.0, 1.0))
xbc0 = ReflectingBoundary2d(xPlane0)
ybc0 = ReflectingBoundary2d(yPlane0)
boundaryConditions = [xbc0, ybc0]

# Enforce the boundary conditions.
for bc in boundaryConditions:
    bc.setGhostNodes(db)
    nodes1.neighbor.updateNodes()
    bc.applyFieldListGhostBoundary(fluidWeight)
    bc.applyFieldListGhostBoundary(fluidMass)
    bc.applyFieldListGhostBoundary(fluidRho)
    bc.applyFieldListGhostBoundary(velocity)
    bc.applyFieldListGhostBoundary(fluidHfield)

################################################################################
# Generate the analytic answer.
import Gnuplot
xans = array([0.0]*n1)
yans = array([0.0]*n1)
i = 0
for nodeList in db.nodeLists:
    for r in nodeList.positions[:nodeList.numInternalNodes]:
        xans[i] = r.magnitude()
        yans[i] = gradDivVelFunction(r)
        i = i + 1
ansData = Gnuplot.Data(xans, yans, with='lines', title='Analytic answer')

################################################################################
# Plot the direct, successive first derivatives against the known solution.
dummy = FieldFunctions()

# directGradDivVel = dummy.gradDivVectorFieldList2d(velocity,
#                                                   fluidPosition,
#                                                   fluidWeight,
#                                                   fluidMass,
#                                                   fluidRho,
#                                                   fluidHfield,
#                                                   WT,
#                                                   boundaryConditions)

# plotDirect1 = plotVectorField2d(db, directGradDivVel, vectorMultiplier=1e-3,
#                                title = 'Direct second derivative of velocity.')
# plotDirect = plotFieldList(directGradDivVel,
#                            xFunction = '%s.magnitude()',
#                            yFunction = '%s.magnitude()',
#                            plotStyle = 'points',
#                            winTitle = 'Direct second derivative of velocity.')
# plotDirect.replot(ansData)

# ##################################################################################
# # Plot the simple second derivative method against the known solution.
# simpleGradDivVel = dummy.gradDivVectorFieldListSimple2d(velocity,
#                                                         fluidPosition,
#                                                         fluidWeight,
#                                                         fluidMass,
#                                                         fluidRho,
#                                                         fluidHfield,
#                                                         WT)

# plotSimple1 = plotVectorField2d(db, simpleGradDivVel, vectorMultiplier=1e-4,
#                                title = 'Simple second derivative of velocity.')
# plotSimple = plotFieldList(simpleGradDivVel,
#                            xFunction = '%s.magnitude()',
#                            yFunction = '%s.magnitude()',
#                            plotStyle = 'points',
#                            winTitle = 'Simple second derivative of velocity.')
# plotSimple.replot(ansData)

##################################################################################
# Plot the golden second derivative method against the known solution.
goldenGradDivVel = dummy.gradDivVectorFieldListGolden2d(velocity,
                                                        fluidPosition,
                                                        fluidWeight,
                                                        fluidMass,
                                                        fluidRho,
                                                        fluidHfield,
                                                        WT)
integrator = PredictorCorrectorIntegrator2d(db)
integrator.timerSummary()

plotGolden1 = plotVectorField2d(db, goldenGradDivVel, vectorMultiplier=1e-4,
                               title = 'Golden second derivative of velocity.')

plotGolden = plotFieldList(goldenGradDivVel,
                           xFunction = '%s.magnitude()',
                           yFunction = '%s.magnitude()',
                           plotStyle = 'points',
                           winTitle = 'Golden second derivative of velocity.')
plotGolden.replot(ansData)
plotGolden("set yrange [1:]; set logscale y; replot")

# ##################################################################################
# # Plot the golden2 second derivative method against the known solution.
# golden2GradDivVel = dummy.gradDivVectorFieldListGolden22d(velocity,
#                                                         fluidPosition,
#                                                         fluidWeight,
#                                                         fluidMass,
#                                                         fluidRho,
#                                                         fluidHfield,
#                                                         WT)

# plotGolden21 = plotVectorField2d(db, golden2GradDivVel, vectorMultiplier=1e-4,
#                                title = 'Golden2 second derivative of velocity.')
# plotGolden2 = plotFieldList(golden2GradDivVel,
#                            xFunction = '%s.magnitude()',
#                            yFunction = '%s.magnitude()',
#                            plotStyle = 'points',
#                            winTitle = 'Golden2 second derivative of velocity.')
# plotGolden2.replot(ansData)

# ##################################################################################
# # Plot the mash second derivative method against the known solution.
# mashGradDivVel = dummy.gradDivVectorFieldListMash2d(velocity,
#                                                         fluidPosition,
#                                                         fluidWeight,
#                                                         fluidMass,
#                                                         fluidRho,
#                                                         fluidHfield,
#                                                         WT)

# plotMash1 = plotVectorField2d(db, mashGradDivVel, vectorMultiplier=1e-4,
#                                title = 'Mash second derivative of velocity.')
# plotMash = plotFieldList(mashGradDivVel,
#                            xFunction = '%s.magnitude()',
#                            yFunction = '%s.magnitude()',
#                            plotStyle = 'points',
#                            winTitle = 'Mash second derivative of velocity.')
# plotMash.replot(ansData)

##################################################################################
# Plot the pairwise second derivative method against the known solution.
pwGradDivVel = dummy.gradDivVectorFieldListPairWise2d(velocity,
                                                        fluidPosition,
                                                        fluidMass,
                                                        fluidMass,
                                                        fluidRho,
                                                        fluidHfield,
                                                        WT)

plotPW1 = plotVectorField2d(db, pwGradDivVel, vectorMultiplier=1e-3,
                               title = 'Pairwise second derivative of velocity.')

plotPW = plotFieldList(pwGradDivVel,
                           xFunction = '%s.magnitude()',
                           yFunction = '%s.magnitude()',
                           plotStyle = 'points',
                           winTitle = 'Pairwise second derivative of velocity.')
plotPW.replot(ansData)
plotPW("set yrange [1:]; set logscale y; replot")
