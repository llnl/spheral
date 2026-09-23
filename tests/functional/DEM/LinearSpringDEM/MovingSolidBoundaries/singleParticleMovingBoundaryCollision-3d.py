#ATS:DEM3dMSPBC1 = test( SELF, "--clearDirectories True  --boolCheckRestitutionCoefficient True  --normalRestitutionCoefficient 1.0 --steps 100", label="DEM perfectly elastic collision with solid boundary -- 3-D (serial)")
#ATS:DEM3dMSPBC2 = test( SELF, "--clearDirectories True  --boolCheckRestitutionCoefficient True  --normalRestitutionCoefficient 0.5 --steps 100", label="DEM inelastic collision with solid boundary -- 3-D (serial)")

import os, shutil, mpi
from math import *
from Spheral3d import *
from SpheralTestUtilities import *
from findLastRestart import *
from GenerateNodeDistribution3d import *
from GenerateDEMfromSPHGenerator import GenerateDEMfromSPHGenerator3d

if mpi.procs > 1:
    from PeanoHilbertDistributeNodes import distributeNodes3d
else:
    from DistributeNodes import distributeNodes3d

title("DEM Moving Boundary Restitution Coefficient Test")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(vImpact = 1.0,                            # impact velocity (negative z direction)
            vbc = 0.2,                                # velocity of bc  (positive z direction)

            h0 = 1.00,                                # initial height above the solid bc plane
            radius = 0.95,                            # particle radius
            normalSpringConstant=10000.0,             # spring constant for LDS model
            normalRestitutionCoefficient=1.00,        # restitution coefficient to get damping const
            tangentialSpringConstant=2857.0,          # spring constant for LDS model
            tangentialRestitutionCoefficient=0.55,    # restitution coefficient to get damping const
            dynamicFriction = 1.0,                    # static friction coefficient sliding
            staticFriction = 1.0,                     # dynamic friction coefficient sliding
            rollingFriction = 1.05,                   # static friction coefficient for rolling
            torsionalFriction = 1.3,                  # static friction coefficient for torsion
            cohesiveTensileStrength = 0.0,            # units of pressure
            shapeFactor = 0.5,                        # in [0,1] shape factor from Zhang 2018, 0 - no torsion or rolling

            neighborSearchBuffer = 0.1,             # multiplicative buffer to radius for neighbor search algo

            # integration
            IntegratorConstructor = VerletIntegrator, # Verlet one integrator to garenteee conservation
            stepsPerCollision = 50,                   # replaces CFL for DEM
            goalTime = 3.0,
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
            dataDir = "dumps-DEM-particle-boundary-3d", 

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
assert mpi.procs == 1 
assert h0 > radius
assert shapeFactor <= 1.0 and shapeFactor >= 0.0
assert dynamicFriction >= 0.0
assert staticFriction >= 0.0
assert torsionalFriction >= 0.0
assert rollingFriction >= 0.0
assert cohesiveTensileStrength >= 0.0
    
#-------------------------------------------------------------------------------
# file things
#-------------------------------------------------------------------------------
testName = "DEM-SingleParticleMovingBoundaryCollision-3d"

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
generator0 = GenerateNodeDistribution3d(1, 1, 1,
                                        rho = 1.0,
                                        distributionType = "lattice",
                                        xmin = (-1.0,  -1.0, -1+h0),
                                        xmax = (1.0,  1.0, 1+h0))

generator1 = GenerateDEMfromSPHGenerator3d(WT,
                                           generator0)
distributeNodes3d((nodes1, generator1))
 
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
# DEM
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
          stepsPerCollision = stepsPerCollision)

packages = [dem]

# set the solid bcs and add to dem package
solidWall = InfinitePlaneSolidBoundary(Vector(0.0, 0.0, 0.0), 
                                       Vector(0.0, 0.0, 1.0),
                                       Vector(0.0, 0.0, 0.0))
solidWall.velocity = Vector(0.0,0.0,vbc)

dem.appendSolidBoundary(solidWall)

#-------------------------------------------------------------------------------
# initial conditions
#-------------------------------------------------------------------------------
velocity = nodes1.velocity()
particleRadius = nodes1.particleRadius()

velocity[0] = Vector(0.0,0.0,-vImpact)
particleRadius[0] = radius
    
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
    control.advance(goalTime, maxSteps)

#-------------------------------------------------------------------------------
# Great success?
#-------------------------------------------------------------------------------
if boolCheckRestitutionCoefficient:
# check our restitution coefficient is correct
#-------------------------------------------------------------
    vijPostImpact = -velocity[0].z + vbc
    vijPreImpact = vImpact + vbc
    restitutionEff = vijPostImpact/vijPreImpact
    restitutionError = abs(restitutionEff + normalRestitutionCoefficient)/normalRestitutionCoefficient
    if  restitutionError > restitutionErrorThreshold:
        print("    final velocity = {0}".format(vijPostImpact))
        print("  initial velocity = {0}".format(vijPreImpact))
        raise ValueError("  relative restitution coefficient error, %g, exceeds bounds" % restitutionError)