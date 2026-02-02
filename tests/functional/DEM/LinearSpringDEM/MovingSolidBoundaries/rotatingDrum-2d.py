import os, shutil, mpi
from math import *
from Spheral2d import *
from SpheralTestUtilities import *
from findLastRestart import *
from GenerateNodeDistribution2d import *
from GenerateDEMfromSPHGenerator import GenerateDEMfromSPHGenerator2d

if mpi.procs > 1:
    from PeanoHilbertDistributeNodes import distributeNodes2d
else:
    from DistributeNodes import distributeNodes2d

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(omegaDrum = 0.25,                         # angular velocity of drum
            radiusDrum = 10.0,                        # radius of the drum
            yThresh = 0.2,                            # level to initial fill to
            yClip = 0.0,                              # level to clip at after settling
            g0 = 5.0,                                 # grav acceleration

            radius = 0.5,                            # particle radius
            normalSpringConstant=10000.0,             # spring constant for LDS model
            normalRestitutionCoefficient=0.55,        # restitution coefficient to get damping const
            tangentialSpringConstant=2857.0,          # spring constant for LDS model
            tangentialRestitutionCoefficient=0.55,    # restitution coefficient to get damping const
            dynamicFriction = 1.0,                    # static friction coefficient sliding
            staticFriction = 1.0,                     # dynamic friction coefficient sliding
            rollingFriction = 1.05,                   # static friction coefficient for rolling
            torsionalFriction = 1.3,                  # static friction coefficient for torsion
            cohesiveTensileStrength = 0.0,            # units of pressure
            shapeFactor = 0.1,                        # in [0,1] shape factor from Zhang 2018, 0 - no torsion or rolling

            neighborSearchBuffer = 0.1,             # multiplicative buffer to radius for neighbor search algo

            # integration
            IntegratorConstructor = VerletIntegrator, # Verlet one integrator to garenteee conservation
            stepsPerCollision = 50,                   # replaces CFL for DEM
            updateBoundaryFrequency = 10,             # CAREFUL: make sure fast time stepping is off for DEM
            settleTime = 5.0,                         # time simulated before we start spinning
            goalTime = 15.0,                          # duration of spin cycle
            dt = 1e-8,
            dtMin = 1.0e-8, 
            dtMax = 0.1,
            dtGrowth = 2.0,
            steps = None,
            maxSteps = None,
            statsStep = 10,
            domainIndependent = False,
            rigorousBoundaries = False,
            dtverbose = False,
            
            # output control
            vizCycle = None,
            vizTime = 0.1,
            clearDirectories = False,
            restoreCycle = None,
            restartStep = 1000,
            redistributeStep = 500,
            dataDir = "dumps-DEM-rotating-drum-2d", 

             # ats parameters
            boolCheckRestitutionCoefficient=False, # turn on error checking for restitution coefficient
            boolCheckSlidingFrictionX=False,       # checks sliding friction reduces relative rotation
            boolCheckSlidingFrictionY=False,       # checks rolling friction reduces relative rotation 
            boolCheckTorsionalFriction=False,      # checks torsional friction reduces relative rotation
            restitutionErrorThreshold = 0.02,      # relative error actual restitution vs nominal
            omegaThreshold = 1e-14,                # theshold for perpendicular components that should stay zero
            )

#-------------------------------------------------------------------------------
# check for bad inputs
#-------------------------------------------------------------------------------
assert radius < radiusDrum/2.0
assert shapeFactor <= 1.0 and shapeFactor >= 0.0
assert dynamicFriction >= 0.0
assert staticFriction >= 0.0
assert torsionalFriction >= 0.0
assert rollingFriction >= 0.0
assert cohesiveTensileStrength >= 0.0
    
#-------------------------------------------------------------------------------
# file things
#-------------------------------------------------------------------------------
testName = "DEM-RotatingDrum-2d"

dataDir = os.path.join(dataDir,
                  "restitutionCoefficient=%s" % normalRestitutionCoefficient)

restartDir = os.path.join(dataDir, "restarts")
vizDir = os.path.join(dataDir, "visit")
restartBaseName = os.path.join(restartDir, testName)
vizBaseName = testName

if vizCycle is None and vizTime is None:
    vizBaseName=None

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
if mpi.rank == 0:
    if clearDirectories and os.path.exists(dataDir):
        shutil.rmtree(dataDir)
    if not os.path.exists(restartDir):
        os.makedirs(restartDir)
    if not os.path.exists(vizDir):
        os.makedirs(vizDir)
mpi.barrier()

#-------------------------------------------------------------------------------
# If we're restarting, find the set of most recent restart files.
#-------------------------------------------------------------------------------
if restoreCycle is None:
    restoreCycle = findLastRestart(restartBaseName)

#-------------------------------------------------------------------------------
# This doesn't really matter kernel filler for neighbor algo
#-------------------------------------------------------------------------------
WT = TableKernel(WendlandC2Kernel(), 1000)

#-------------------------------------------------------------------------------
# Make the NodeList.
#-------------------------------------------------------------------------------
units = CGuS()
nodes1 = makeDEMNodeList("nodeList1",
                          hmin = 1.0e-30,
                          hmax = 1.0e30,
                          hminratio = 100.0,
                          neighborSearchBuffer = neighborSearchBuffer,
                          kernelExtent = WT.kernelExtent)
nodeSet = [nodes1]
for nodes in nodeSet:
    output("nodes.name")
    output("nodes.hmin")
    output("nodes.hmax")
    output("nodes.hminratio")
    output("nodes.nodesPerSmoothingScale")

#-------------------------------------------------------------------------------
# Set the node properties. (gen 2 particles visit doesn't like just one)
#-------------------------------------------------------------------------------
def rejecter(x,y,m,H):
    xnew,ynew,mnew,Hnew = [],[],[],[]
    for i in range(len(x)):
        if y[i] < yThresh and sqrt(x[i]**2+y[i]**2) < radiusDrum - 1.10*radius:
            xnew.append(x[i])
            ynew.append(y[i])
            mnew.append(m[i])
            Hnew.append(H[i])

    return xnew,ynew,mnew,Hnew

nx = int(2*radiusDrum / (radius * 1.02))
ny = nx

generator0 = GenerateNodeDistribution2d(nx, ny,
                                        rho = 1.0,
                                        distributionType = "lattice",
                                        xmin = (-radiusDrum,  -radiusDrum),
                                        xmax = ( radiusDrum,   radiusDrum),
                                        rejecter=rejecter)

generator1 = GenerateDEMfromSPHGenerator2d(WT,
                                           generator0)
distributeNodes2d((nodes1, generator1))
 
#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase()
output("db")
for nodes in nodeSet:
    db.appendNodeList(nodes)
output("db.numNodeLists")
output("db.numDEMNodeLists")
output("db.numFluidNodeLists")


#-------------------------------------------------------------------------------
# Physics Package: DEM
#-------------------------------------------------------------------------------
dem = DEM(db,
          normalSpringConstant = normalSpringConstant,
          normalRestitutionCoefficient = normalRestitutionCoefficient,
          tangentialSpringConstant = tangentialSpringConstant,
          tangentialRestitutionCoefficient = tangentialRestitutionCoefficient,
          dynamicFrictionCoefficient = dynamicFriction,
          staticFrictionCoefficient = staticFriction,
          rollingFrictionCoefficient = rollingFriction,
          torsionalFrictionCoefficient = torsionalFriction,
          cohesiveTensileStrength =cohesiveTensileStrength,
          shapeFactor = shapeFactor,
          stepsPerCollision = stepsPerCollision,
          enableFastTimeStepping = False)

packages = [dem]

solidWall = SphereSolidBoundary(center = Vector(0.0, 0.0), 
                                 radius = radiusDrum,
                                 angularVelocity = 0.0)

dem.appendSolidBoundary(solidWall)

#-------------------------------------------------------------------------------
# Gravity: DEM
#-------------------------------------------------------------------------------
gravity = ConstantAcceleration(a0 = Vector(0.0,-g0),
                               nodeList = nodes1)
packages += [gravity]

#-------------------------------------------------------------------------------
# initial conditions
#-------------------------------------------------------------------------------
# velocity = nodes1.velocity()
# particleRadius = nodes1.particleRadius()

# velocity[0] = Vector(0.0,0.0,-vImpact)
# particleRadius[0] = radius
    
#-------------------------------------------------------------------------------
# Construct a time integrator, and add the physics packages.
#-------------------------------------------------------------------------------
integrator = IntegratorConstructor(db)
for p in packages:
    integrator.appendPhysicsPackage(p)
integrator.lastDt = dt
integrator.dtMin = dtMin
integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
integrator.domainDecompositionIndependent = domainIndependent
integrator.verbose = dtverbose
integrator.rigorousBoundaries = rigorousBoundaries
integrator.updateBoundaryFrequency = 10
integrator.cullGhostNodes = False

output("integrator")
output("integrator.havePhysicsPackage(dem)")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.dtGrowth")
output("integrator.domainDecompositionIndependent")
output("integrator.rigorousBoundaries")
output("integrator.verbose")

#-------------------------------------------------------------------------------
# Make the problem controller.
#-------------------------------------------------------------------------------
from SpheralPointmeshSiloDump import dumpPhysicsState
control = SpheralController(integrator, WT,
                            iterateInitialH = False,
                            initializeDerivatives = True,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle,
                            vizBaseName = vizBaseName,
                            vizMethod=dumpPhysicsState,
                            vizDir = vizDir,
                            vizStep = vizCycle,
                            vizTime = vizTime)
output("control")

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
if not steps is None:
    control.step(steps)
else:
    control.advance(settleTime, maxSteps)

solidWall.angularVelocity = omegaDrum

if not steps is None:
    control.step(steps)
else:
    control.advance(settleTime+goalTime, maxSteps)
