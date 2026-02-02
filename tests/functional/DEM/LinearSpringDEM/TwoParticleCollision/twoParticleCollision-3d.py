#ATS:DEM3d0 = test(        SELF, "--clearDirectories True  --checkError True  --checkConservation True --normalRestitutionCoefficient 1.0 --steps 100", label="DEM perfectly elastic 2 particle collision -- 3-D (serial)")
#ATS:DEM3d1 = test(        SELF, "--clearDirectories True  --checkError True  --checkConservation True --normalRestitutionCoefficient 0.8 --steps 100", label="DEM perfectly inelastic 2 particle collision -- 3-D (serial)")
#ATS:DEM3d2 = test(        SELF, "--clearDirectories True  --checkError True --boolCheckSlidingFriction True --checkConservation True --normalRestitutionCoefficient 0.8 --steps 100", label="DEM inelastic 2 particle collision - sliding friction -- 3-D (serial)")
#ATS:DEM3d3 = test(        SELF, "--clearDirectories True  --checkError True --boolCheckRollingFriction True --checkConservation True --normalRestitutionCoefficient 0.8 --steps 100", label="DEM inelastic 2  particle collision - rolling friction -- 3-D (serial)")
#ATS:DEM3d4 = test(        SELF, "--clearDirectories True  --checkError True --boolCheckTorsionalFriction True --checkConservation True --normalRestitutionCoefficient 0.8 --steps 100", label="DEM inelastic 2  particle collision - torsional friction -- 3-D (serial)")
#ATS:DEM3d5 = test(        SELF, "--clearDirectories True  --checkError True --boolCheckTorsionalObjectivity True --checkConservation True --normalRestitutionCoefficient 0.8 --steps 100", label="DEM inelastic 2  particle collision - torsional objectivity -- 3-D (serial)")

import os, sys, shutil, mpi
from math import *
from Spheral3d import *
from SpheralTestUtilities import *
from findLastRestart import *
from GenerateNodeDistribution3d import *
from GenerateDEMfromSPHGenerator import GenerateDEMfromSPHGenerator3d

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from DEMConservationTracker import TrackConservation3d as TrackConservation

if mpi.procs > 1:
    from PeanoHilbertDistributeNodes import distributeNodes3d
else:
    from DistributeNodes import distributeNodes3d

title("DEM Restitution Coefficient Test")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(vImpact = 1.0,                            # impact velocity
            omega0 = 0.1,                             # initial angular velocity it we're doing that

            radius = 0.25,                            # particle radius
            normalSpringConstant=10000.0,             # spring constant for LDS model
            normalRestitutionCoefficient=0.55,        # restitution coefficient to get damping const
            tangentialSpringConstant=2857.0,          # spring constant for LDS model
            tangentialRestitutionCoefficient=0.55,    # restitution coefficient to get damping const
            dynamicFriction = 1.0,                    # static friction coefficient sliding
            staticFriction = 1.0,                     # dynamic friction coefficient sliding
            rollingFriction = 1.05,                   # static friction coefficient for rolling
            torsionalFriction = 1.3,                  # static friction coefficient for torsion
            cohesiveTensileStrength = 0.0,            # units of pressure
            shapeFactor = 0.5,                        # in [0,1] shape factor from Zhang 2018, 0 - no torsion or rolling

            neighborSearchBuffer = 0.1,               # multiplicative buffer to radius for neighbor search algo

            # integration
            IntegratorConstructor = VerletIntegrator, # Verlet one integrator to garenteee conservation
            stepsPerCollision = 50,                   # replaces CFL for DEM
            goalTime = None,
            dt = 1e-8,
            dtMin = 1.0e-8, 
            dtMax = 0.1,
            dtGrowth = 2.0,
            steps = 500,
            maxSteps = None,
            statsStep = 10,
            domainIndependent = False,
            rigorousBoundaries = False,
            dtverbose = False,
            
            # output control
            vizCycle = None,
            vizTime = None,
            clearDirectories = False,
            restoreCycle = None,
            restartStep = 1000,
            redistributeStep = 500,
            dataDir = "dumps-DEM-2particle-3d", 

             # ats parameters
            checkError = False,                    # turn on error checking for restitution coefficient
            boolCheckSlidingFriction=False,        # checks sliding friction reduces relative rotation
            boolCheckRollingFriction=False,        # checks rolling friction reduces relative rotation 
            boolCheckTorsionalFriction=False,      # checks torsional friction reduces relative rotation
            boolCheckTorsionalObjectivity=False,   # checks to make sure torsion is objective
            checkRestart = False,                  # turn on error checking for restartability
            checkConservation = False,             # turn on error checking for momentum conservation
            restitutionErrorThreshold = 0.02,      # relative error actual restitution vs nominal
            conservationErrorThreshold = 1e-15,    # relative error for momentum conservation
            omegaThreshold = 1e-14,                # tolerance for erroneous spin up in inactive directions
            torsionalObjectivityThreshold = 1e-10  # relative error bounds on torsion objectivity test
            )

#------------------------------------------------------------------------------
# check for bad inputs
#-------------------------------------------------------------------------------
assert mpi.procs == 1 
assert shapeFactor <= 1.0 and shapeFactor >= 0.0
assert dynamicFriction >= 0.0
assert staticFriction >= 0.0
assert torsionalFriction >= 0.0
assert rollingFriction >= 0.0
assert cohesiveTensileStrength >= 0.0
assert sum([boolCheckSlidingFriction,
            boolCheckRollingFriction,
            boolCheckTorsionalFriction,
            boolCheckTorsionalObjectivity]) <= 1

if boolCheckSlidingFriction:
    shapeFactor = 0.0
    
#-------------------------------------------------------------------------------
# file things
#-------------------------------------------------------------------------------
testName = "DEM-twoParticleCollision-3d"

dataDir = os.path.join(dataDir,
                  "restitutionCoefficient=%s" % normalRestitutionCoefficient,
                  "boolCheckSlidingFriction=%s" % boolCheckSlidingFriction,
                  "boolCheckRollingFriction=%s" % boolCheckRollingFriction,
                  "boolCheckTorsionalFriction=%s" % boolCheckTorsionalFriction,
                  "boolCheckTorsionalObjectivity=%s" % boolCheckTorsionalObjectivity)

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
# Set the node properties.
#-------------------------------------------------------------------------------
generator0 = GenerateNodeDistribution3d(2, 1, 1,
                                        rho = 1.0,
                                        distributionType = "lattice",
                                        xmin = (0.0,  0.0, 0.0),
                                        xmax = (1.0,  0.5, 0.5))

generator1 = GenerateDEMfromSPHGenerator3d(WT,
                                           generator0,
                                           particleRadius=radius)
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


#-------------------------------------------------------------------------------
# initial conditions
#-------------------------------------------------------------------------------

velocity = nodes1.velocity()
velocity[0] = Vector(vImpact,0.0,0.0)
velocity[1] = Vector(-vImpact,0.0,0.0)

particleRadius = nodes1.particleRadius()
particleRadius[0] = radius
particleRadius[1] = radius

bonusSpace = radius
position = nodes1.positions()
position[0].x -= bonusSpace
position[1].x += bonusSpace

omega = dem.omega
if boolCheckSlidingFriction:
    omega[0][0] = Vector(0.0,0.0,omega0)
    omega[0][1] = Vector(0.0,0.0,omega0)
elif boolCheckRollingFriction:  
    omega[0][0] = Vector(0.0,0.0, omega0)
    omega[0][1] = Vector(0.0,0.0,-omega0)
elif boolCheckTorsionalFriction:
    omega[0][0] = Vector( omega0,0.0,0.0)
    omega[0][1] = Vector(-omega0,0.0,0.0)
elif boolCheckTorsionalObjectivity:
    omega[0][0] = Vector( omega0,0.0,0.0)
    omega[0][1] = Vector( omega0,0.0,0.0)
    

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
# Periodic Work Function : track conservation
#-------------------------------------------------------------------------------
conservation = TrackConservation(db,
                                  dem,
                                  verbose=False)
                                  
periodicWork = [(conservation.periodicWorkFunction,1)]

#-------------------------------------------------------------------------------
# Make the problem controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            iterateInitialH = False,
                            initializeDerivatives = True,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle,
                            vizBaseName = vizBaseName,
                            vizDir = vizDir,
                            vizStep = vizCycle,
                            vizTime = vizTime,
                            periodicWork=periodicWork)
output("control")

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------

if not steps is None:
    if checkRestart:
        control.setRestartBaseName(restartBaseName + "_CHECK")
    control.step(steps)
else:
    control.advance(goalTime, maxSteps)

#-------------------------------------------------------------------------------
# Great success?
#-------------------------------------------------------------------------------
if checkRestart:
# check reproducibility when restarted
#-------------------------------------------------------------
    control.setRestartBaseName(restartBaseName)
    state0 = State(db, integrator.physicsPackages())
    state0.copyState()
    control.loadRestartFile(control.totalSteps)
    state1 = State(db, integrator.physicsPackages())
    if not state1 == state0:
        raise ValueError("The restarted state does not match!")
    else:
        print("Restart check PASSED.")

if checkError:
# check our restitution coefficient is correct
#-------------------------------------------------------------
    vijPostImpact = velocity[0].x - velocity[1].x
    vijPreImpact = 2.0*vImpact
    restitutionEff = vijPostImpact/vijPreImpact
    restitutionError = abs(restitutionEff + normalRestitutionCoefficient)/normalRestitutionCoefficient
    if  restitutionError > restitutionErrorThreshold:
        raise ValueError("relative restitution coefficient error, %g, exceeds bounds" % restitutionError)

if checkConservation:
# check momentum conservation
#-------------------------------------------------------------
    if  conservation.deltaLinearMomentumX() > conservationErrorThreshold:
        raise ValueError("linear momentum - x conservation error, %g, exceeds bounds" % conservation.deltaLinearMomentumX())
    if  conservation.deltaLinearMomentumY() > conservationErrorThreshold:
        raise ValueError("linear momentum - y conservation error, %g, exceeds bounds" % conservation.deltaLinearMomentumY())
    if  conservation.deltaLinearMomentumZ() > conservationErrorThreshold:
        raise ValueError("linear momentum - z conservation error, %g, exceeds bounds" % conservation.deltaLinearMomentumZ())
    if  conservation.deltaRotationalMomentumX() > conservationErrorThreshold:
        raise ValueError("rotational momentum - x conservation error, %g, exceeds bounds" % conservation.deltaRotationalMomentumX())
    if  conservation.deltaRotationalMomentumY() > conservationErrorThreshold:
        raise ValueError("rotational momentum - y conservation error, %g, exceeds bounds" % conservation.deltaRotationalMomentumY())
    if  conservation.deltaRotationalMomentumZ() > conservationErrorThreshold:
        raise ValueError("rotational momentum -z conservation error, %g, exceeds bounds" % conservation.deltaRotationalMomentumZ())

if boolCheckSlidingFriction or boolCheckRollingFriction:
# check for non-physical behavior
#-------------------------------------------------------------
    if omega[0][0].magnitude()+omega[0][1].magnitude() > 2*omega0:
        raise ValueError("particles are rotating faster post-collision")
    if abs(omega[0][0].x) > omegaThreshold or abs(omega[0][0].y) > omegaThreshold:
        raise ValueError("erroneous spin-up of particle 0 in perpendicular direction")
    if abs(omega[0][1].x) > omegaThreshold or abs(omega[0][1].y) > omegaThreshold:
        raise ValueError("erroneous spin-up of particle 1 in perpendicular direction")
    
if  boolCheckTorsionalFriction:
# check for non-physical behavior
#-------------------------------------------------------------
    if omega[0][0].magnitude()+omega[0][1].magnitude() > 2*omega0:
        raise ValueError("particles are rotating faster post-collision")
    if abs(omega[0][0].z) > omegaThreshold or abs(omega[0][0].y) > omegaThreshold:
        raise ValueError("erroneous spin-up of particle 0 in perpendicular direction")
    if abs(omega[0][1].z) > omegaThreshold or abs(omega[0][1].y) > omegaThreshold:
        raise ValueError("erroneous spin-up of particle 1 in perpendicular direction")
    
if boolCheckTorsionalObjectivity:
# to satify objectivity omega (along axis) should not change when equal
#-------------------------------------------------------------
    omegaError = (2*omega0 - omega[0][0][0] - omega[0][1][0]) / (2*omega0)
    if omegaError > torsionalObjectivityThreshold:
        raise ValueError("torsional objectivity failure with relative angular velocity error, %g, exceeds bounds" % omegaError)
    if abs(omega[0][0].z) > omegaThreshold or abs(omega[0][0].y) > omegaThreshold:
        raise ValueError("erroneous spin-up of particle 0 in perpendicular direction")
    if abs(omega[0][1].z) > omegaThreshold or abs(omega[0][1].y) > omegaThreshold:
        raise ValueError("erroneous spin-up of particle 1 in perpendicular direction")
