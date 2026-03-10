#ATS:DEM2d5part0 = test(SELF, "--clearDirectories True --checkContacts True --checkAngularMomentum True", label="DEM 5 particle collision tests for pairwise storage -- 2-D (serial)")
#ATS:DEM2d5part1 = test(SELF, "--clearDirectories True --checkContacts True --checkAngularMomentum True", label="DEM 5 particle collision tests for pairwise storage -- 2-D (parallel", np=5)

import os, sys, shutil, mpi
from math import *
from Spheral2d import *
from SpheralTestUtilities import *
from findLastRestart import *
from GenerateNodeDistribution2d import *
from GenerateDEMfromSPHGenerator import GenerateDEMfromSPHGenerator2d

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from DEMConservationTracker import TrackConservation2d as TrackConservation

if mpi.procs > 1:
    from PeanoHilbertDistributeNodes import distributeNodes2d
else:
    from DistributeNodes import distributeNodes2d

title("DEM Neighbor Tracking Check")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(vImpact = 1.0,                            # impact velocity

            radius = 0.25,                            # particle radius
            normalSpringConstant=10000.0,             # spring constant for LDS model
            normalRestitutionCoefficient=0.55,        # restitution coefficient to get damping const
            tangentialSpringConstant=2857.0,          # spring constant for LDS model
            tangentialRestitutionCoefficient=0.55,    # restitution coefficient to get damping const
            dynamicFriction = 1.0,                    # dynamic sliding friction coefficient
            staticFriction = 1.0,                     # static sliding friction coefficient
            rollingFriction = 1.05,                   # rolling friction coefficient
            torsionalFriction = 1.3,                  # torisional friction coefficient
            cohesiveTensileStrength =0.0,             # units of pressure
            shapeFactor = 0.5,                        # shape irregularity parameter 0-1 (1 most irregular)
            
            neighborSearchBuffer = 0.001,             # multiplicative buffer to radius for neighbor search algo

            # integration
            IntegratorConstructor = VerletIntegrator,    # Verlet only integrator that garentees conservation of Rot Mom w/ DEM
            stepsPerCollision = 50,                      # replaces CFL for DEM
            goalTime = None,
            dt = 1e-8,
            dtMin = 1.0e-8, 
            dtMax = 0.1,
            dtGrowth = 2.0,
            steps = 10,
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
            restartStep = 10,
            redistributeStep = 500,
            dataDir = "dumps-DEM-5particle-2d",

             # ats parameters
            checkAngularMomentum = True,       # soft check to make sure angular velocity of central particle drops
            checkContacts = True,              # turn on error checking for our contacts
            checkRestart = False,              # turn on error checking for restartability
            checkConservation = False,         # turn on error checking for momentum conservation
            conservationErrorThreshold = 1e-15 # relative error for momentum conservation
            )

#-------------------------------------------------------------------------------
# check for bad inputs
#-------------------------------------------------------------------------------
assert mpi.procs == 1 or mpi.procs==5
assert shapeFactor <= 1.0 and shapeFactor >= 0.0
assert dynamicFriction >= 0.0
assert staticFriction >= 0.0
assert torsionalFriction >= 0.0
assert rollingFriction >= 0.0
assert cohesiveTensileStrength >= 0.0

#-------------------------------------------------------------------------------
# file things
#-------------------------------------------------------------------------------
testName = "DEM-fiveParticleCollision-2d"
dataDir = os.path.join(dataDir,
                       "mpiprocs=%s" % mpi.procs)
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

# used to generate 5 points, locations and properties changed later
generator0 = GenerateNodeDistribution2d(5, 1,
                                        rho = 1.0,
                                        distributionType = "lattice",
                                        xmin = (-0.5,  0.0),
                                        xmax = (2.0,  0.5))

# make it a DEM generator
generator1 = GenerateDEMfromSPHGenerator2d(WT,
                                           generator0,
                                           particleRadius=radius)
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
# Initial Conditions
#-------------------------------------------------------------------------------

velocity = nodes1.velocity()
position = nodes1.positions()
particleRadius = nodes1.particleRadius()

uniqueIndices = db.DEMUniqueIndex
omega = dem.omega 
neighborIndices = dem.neighborIndices

for i in range(nodes.numInternalNodes):
    if position[i][0]<0.0:
        velocity[i] = Vector(0,0.0)
        position[i] = Vector( 0.25, 0.25)
        omega[0][i]=1
        uniqueIndices[0][i]=1
    elif position[i][0]<0.5:
        velocity[i] = Vector( vImpact,0.0)
        position[i] = Vector(-0.25, 0.25)
        uniqueIndices[0][i]=2
    elif position[i][0]<1.0:
        velocity[i] = Vector(-vImpact,0.0)
        position[i] = Vector( 0.75, 0.25)
        uniqueIndices[0][i]=3
    elif position[i][0]<1.5:
        velocity[i] = Vector(0.0, vImpact)
        position[i] = Vector( 0.25,-0.25)
        uniqueIndices[0][i]=4
    elif position[i][0]<2.0:
        velocity[i] = Vector(0.0,-vImpact)
        position[i] = Vector( 0.25, 0.75)
        uniqueIndices[0][i]=5



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
                            periodicWork = periodicWork)
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
    control.setRestartBaseName(restartBaseName)
    state0 = State(db, integrator.physicsPackages())
    state0.copyState()
    control.loadRestartFile(control.totalSteps)
    state1 = State(db, integrator.physicsPackages())
    if not state1 == state0:
        raise ValueError("The restarted state does not match!")
    else:
        print("Restart check PASSED.")

if checkContacts and mpi.procs==5:
    for i in range(nodes.numInternalNodes):
        if  uniqueIndices[0][i]==1:
            if set(neighborIndices[0][i]) != set([2,3,4,5]):
                print("-----------------------------")
                print("unique index   : %s" % uniqueIndices[0][i])
                print("mpi rank       : %s" % mpi.rank)
                print("neighborIndices: %s" % neighborIndices[0][i])
                print("-----------------------------")
                raise ValueError("Central node contact list is incorrect ")
        elif uniqueIndices[0][i]==2:
            if set(neighborIndices[0][i]) != set([1]):
                print("-----------------------------")
                print("unique index   : %s" % uniqueIndices[0][i])
                print("mpi rank       : %s" % mpi.rank)
                print("neighborIndices: %s" % neighborIndices[0][i])
                print("-----------------------------")
                raise ValueError("node 1 contact list is incorrect ")
        elif uniqueIndices[0][i]==3:
            if set(neighborIndices[0][i]) != set([1]):
                print("-----------------------------")
                print("unique index   : %s" % uniqueIndices[0][i])
                print("mpi rank       : %s" % mpi.rank)
                print("neighborIndices: %s" % neighborIndices[0][i])
                print("-----------------------------")
                raise ValueError("node 2 contact list is incorrect ")
        elif uniqueIndices[0][i]==4:
            if set(neighborIndices[0][i]) != set([1]):
                print("-----------------------------")
                print("unique index   : %s" % uniqueIndices[0][i])
                print("mpi rank       : %s" % mpi.rank)
                print("neighborIndices: %s" % neighborIndices[0][i])
                print("-----------------------------")
                raise ValueError("node 3 contact list is incorrect ")
        elif uniqueIndices[0][i]==5:
            if set(neighborIndices[0][i]) != set([1]):
                print("-----------------------------")
                print("unique index   : %s" % uniqueIndices[0][i])
                print("mpi rank       : %s" % mpi.rank)
                print("neighborIndices: %s" % neighborIndices[0][i])
                print("-----------------------------")
                raise ValueError("node 4 contact list is incorrect ")


if checkContacts and mpi.procs==1:
    for i in range(nodes.numInternalNodes):
        if  uniqueIndices[0][i]==1:
            if set(neighborIndices[0][i]) != set([2,3,4,5]):
                print("-----------------------------")
                print("unique index   : %s" % uniqueIndices[0][i])
                print("mpi rank       : %s" % mpi.rank)
                print("neighborIndices: %s" % neighborIndices[0][i])
                print("-----------------------------")
                raise ValueError("Central node contact list is incorrect ")
        elif uniqueIndices[0][i]==2:
            if set(neighborIndices[0][i]) != set([]):
                print("-----------------------------")
                print("unique index   : %s" % uniqueIndices[0][i])
                print("mpi rank       : %s" % mpi.rank)
                print("neighborIndices: %s" % neighborIndices[0][i])
                print("-----------------------------")
                raise ValueError("node 1 contact list is incorrect ")
        elif uniqueIndices[0][i]==3:
            if set(neighborIndices[0][i]) != set([]):
                print("-----------------------------")
                print("unique index   : %s" % uniqueIndices[0][i])
                print("mpi rank       : %s" % mpi.rank)
                print("neighborIndices: %s" % neighborIndices[0][i])
                print("-----------------------------")
                raise ValueError("node 2 contact list is incorrect ")
        elif uniqueIndices[0][i]==4:
            if set(neighborIndices[0][i]) != set([]):
                print("-----------------------------")
                print("unique index   : %s" % uniqueIndices[0][i])
                print("mpi rank       : %s" % mpi.rank)
                print("neighborIndices: %s" % neighborIndices[0][i])
                print("-----------------------------")
                raise ValueError("node 3 contact list is incorrect ")
        elif uniqueIndices[0][i]==5:
            if set(neighborIndices[0][i]) != set([]):
                print("-----------------------------")
                print("unique index   : %s" % uniqueIndices[0][i])
                print("mpi rank       : %s" % mpi.rank)
                print("neighborIndices: %s" % neighborIndices[0][i])
                print("-----------------------------")
                raise ValueError("node 4 contact list is incorrect ")

if checkAngularMomentum:
    if  uniqueIndices[0][i]==0:
            if omega[0][i] > 1:
                print("-----------------------------")
                print("Angular momentum : %s" % omega[0][i])
                print("mpi rank       : %s" % mpi.rank)
                print("-----------------------------")
                raise ValueError("Central particles spin increased!")

if checkConservation:
    if  conservation.deltaLinearMomentumX() > conservationErrorThreshold:
        raise ValueError("linear momentum - x conservation error, %g, exceeds bounds" % conservation.deltaLinearMomentumX())
    if  conservation.deltaLinearMomentumY() > conservationErrorThreshold:
        raise ValueError("linear momentum - y conservation error, %g, exceeds bounds" % conservation.deltaLinearMomentumY())
    if  conservation.deltaRotationalMomentumZ() > conservationErrorThreshold:
        raise ValueError("rotational momentum -z conservation error, %g, exceeds bounds" % conservation.deltaRotationalMomentumZ())
