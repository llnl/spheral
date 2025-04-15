#ATS:DEM2dImpact = test(          SELF, "--clearDirectories True  --checkConservation True --goalTime 1.0", label="DEM impacting squares -- 2-D (parallel)", np=8)

import os, sys, shutil, mpi

from math import *
from Spheral2d import *
from SpheralTestUtilities import *
from findLastRestart import *
from GenerateNodeDistribution2d import *
from GenerateDEMfromSPHGenerator import GenerateDEMfromSPHGenerator2d

sys.path.insert(0, '..')
from DEMConservationTracker import TrackConservation2d as TrackConservation

if mpi.procs > 1:
    from PeanoHilbertDistributeNodes import distributeNodes2d
else:
    from DistributeNodes import distributeNodes2d

title("DEM 2d Impacting Squares")
# This tests the conservation properties of the DEM package when
# distribution across multiple processors

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(numParticlePerLength = 10,                # number of particles on a side of the box
            radius = 0.25,                            # particle radius
            normalSpringConstant=1000.0,              # spring constant for LDS model
            normalRestitutionCoefficient=0.55,        # restitution coefficient to get damping const
            tangentialSpringConstant=285.70,          # spring constant for LDS model
            tangentialRestitutionCoefficient=0.55,    # restitution coefficient to get damping const
            dynamicFriction = 1.0,                    # static friction coefficient sliding
            staticFriction = 1.0,                     # dynamic friction coefficient sliding
            rollingFriction = 1.05,                   # static friction coefficient for rolling
            torsionalFriction = 1.3,                  # static friction coefficient for torsion
            cohesiveTensileStrength = 0.0,            # units of pressure
            shapeFactor = 0.1,                        # in [0,1] shape factor from Zhang 2018, 0 - no torsion or rolling

            neighborSearchBuffer = 0.1,             # multiplicative buffer to radius for neighbor search algo

            # integration
            IntegratorConstructor = VerletIntegrator,
            stepsPerCollision = 50,  # replaces CFL for DEM
            goalTime = 25.0,
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
            vizTime = None, 
            clearDirectories = False,
            restoreCycle = None,
            restartStep = 10000,
            redistributeStep = 1000000000,
            dataDir = "dumps-DEM-impactingSquares-2d",
            
            # ats
            checkRestart = False,
            checkConservation = False,             # turn on error checking for momentum conservation
            conservationErrorThreshold = 1e-14,    # relative error for momentum conservation
            
            )

#-------------------------------------------------------------------------------
# file things
#-------------------------------------------------------------------------------
testName = "DEM-ImpactingSquares-2d"
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
if restoreCycle is None:
    generator0 = GenerateNodeDistribution2d(2* numParticlePerLength, numParticlePerLength,
                                            rho = 1.0,
                                            distributionType = "lattice",
                                            xmin = (-1.0,  0.0),
                                            xmax = ( 1.0,  1.0))
    
    # really simple bar shaped particle
    def DEMParticleGenerator(xi,yi,Hi,mi,Ri):
        xout = [xi]
        yout = [yi]
        mout = [mi/1.1]
        Rout = [Ri/1.1]
        return xout,yout,mout,Rout

    generator1 = GenerateDEMfromSPHGenerator2d(WT,
                                               generator0,
                                               particleRadius= 0.5/float(numParticlePerLength+1),
                                               DEMParticleGenerator=DEMParticleGenerator)

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
# PhysicsPackage : DEM
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
          cohesiveTensileStrength = cohesiveTensileStrength,
          shapeFactor = shapeFactor,
          stepsPerCollision = stepsPerCollision)

packages = [dem]


#-------------------------------------------------------------------------------
# Initial Conditions
#-------------------------------------------------------------------------------

numNodeLists = db.numNodeLists
nodeLists = db.nodeLists
for i in range(db.numNodeLists):
    nodeListi = nodeLists[i]
    numNodes = nodeListi.numInternalNodes
    v = nodeListi.velocity()
    p = nodeListi.positions()
    for j in range(numNodes):
        if p[j][0] > 0.0:
            v[j][0]= -0.1
            p[j][1]+=0.25/numParticlePerLength
        else:
            v[j][0]=  0.1

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
# Periodic Work Function: Track conseravation
#-------------------------------------------------------------------------------

conservation = TrackConservation(db,
                                  dem,
                                  verbose=True)
                                  
periodicWork = [(conservation.periodicWorkFunction,1)]

#-------------------------------------------------------------------------------
# Make the problem controller.
#-------------------------------------------------------------------------------
from SpheralPointmeshSiloDump import dumpPhysicsState
control = SpheralController(integrator, WT,
                            iterateInitialH = False,
                            initializeDerivatives = True,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            redistributeStep=redistributeStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle,
                            vizBaseName = vizBaseName,
                            vizMethod = dumpPhysicsState,
                            vizGhosts=True,
                            vizDir = vizDir,
                            vizStep = vizCycle,
                            vizTime = vizTime,
                            SPH = SPH,
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

if checkConservation:
    if  conservation.deltaLinearMomentumX() > conservationErrorThreshold:
        raise ValueError("linear momentum - x conservation error, %g, exceeds bounds" % conservation.deltaLinearMomentumX())
    if  conservation.deltaLinearMomentumY() > conservationErrorThreshold:
        raise ValueError("linear momentum - y conservation error, %g, exceeds bounds" % conservation.deltaLinearMomentumY())
    if  conservation.deltaRotationalMomentumZ() > conservationErrorThreshold:
        raise ValueError("rotational momentum -z conservation error, %g, exceeds bounds" % conservation.deltaRotationalMomentumZ())
