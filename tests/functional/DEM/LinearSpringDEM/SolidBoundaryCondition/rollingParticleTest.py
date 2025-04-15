#ATS:DEM3dImpact = test(          SELF, "--clearDirectories True --checkConservation True   --goalTime 1.0", label="DEM impacting squares -- 3-D (parallel)", np=8)

import os, sys, shutil, mpi, random
from math import *
from Spheral3d import *
from SpheralTestUtilities import *
from findLastRestart import *
from GenerateNodeDistribution3d import *
from GenerateDEMfromSPHGenerator import GenerateDEMfromSPHGenerator3d

import numpy as np

sys.path.insert(0, '..')
from DEMConservationTracker import TrackConservation3d as TrackConservation

if mpi.procs > 1:
    from PeanoHilbertDistributeNodes import distributeNodes3d
else:
    from DistributeNodes import distributeNodes3d

title("DEM 3d Drop Test with Particle Generation")
# this tests the ability to generate particle on the fly
# during the course of a simulation using a periodic 
# work function. It also tests the solid boundary condition

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(numParticlePerLength = 3,                 # number of particles on a side of the box
            normalSpringConstant=1.0,          # spring constant for LDS model
            normalRestitutionCoefficient=0.55,        # restitution coefficient to get damping const
            tangentialSpringConstant=0.2857,      # spring constant for LDS model
            tangentialRestitutionCoefficient=0.55,    # restitution coefficient to get damping const
            dynamicFriction = 1.0,                    # static friction coefficient sliding
            staticFriction = 1.0,                     # dynamic friction coefficient sliding
            rollingFriction = 1.05,                   # static friction coefficient for rolling
            torsionalFriction = 1.3,                  # static friction coefficient for torsion
            cohesiveTensileStrength = 0.0,            # units of pressure
            shapeFactor = 0.1,                        # in [0,1] shape factor from Zhang 2018, 0 - no torsion or rolling
            
            particleRadius = 0.10,                    # particle radius
            particleDensity = 2.60,
            particleVelocity = 0.1,

            neighborSearchBuffer = 0.1,             # multiplicative buffer to radius for neighbor search algo
            nPerh = 1.01,

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
            vizTime = 1.0, 
            clearDirectories = False,
            restoreCycle = None,
            restartStep = 10000,
            redistributeStep = 100000000000000,
            dataDir = "dumps-DEM-impactingSquares-3d",

            # ats
            checkRestart = False,
            checkConservation = False,             # turn on error checking for momentum conservation
            conservationErrorThreshold = 2e-14,    # relative error for momentum conservation  
            )

#-------------------------------------------------------------------------------
# file things
#-------------------------------------------------------------------------------
testName = "DEM-ImpactingSquares-3d"
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
if restoreCycle is None:
    generator0 = GenerateNodeDistribution3d(2, 2, 1,
                                            rho = 1.0,
                                            distributionType = "lattice",
                                            xmin = (0, 0.0, 1.0),
                                            xmax = ( 1.0,1.0, 2.0),
                                            nNodePerh = nPerh)
    
    # # really simple bar shaped particle
    def DEMParticleGenerator(xi,yi,zi,Hi,mi,Ri):
        xout = [xi]
        yout = [yi]
        zout = [particleRadius]
        mout = [particleDensity * 4.0/3.0*np.pi*particleRadius**3]
        Rout = [particleRadius]
        return xout,yout,zout,mout,Rout

    generator1 = GenerateDEMfromSPHGenerator3d(WT,
                                               generator0,
                                               DEMParticleGenerator=DEMParticleGenerator)

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

solidWall = InfinitePlaneSolidBoundary(Vector(0.0, 0.0, 0.0), Vector( 0.0, 0.0, 1.0))
#solidWall2 = CylinderSolidBoundary(Vector(0.0, 0.0, -10.0),Vector(0.0, 0.0,1.0),5.0,15.0)

dem.appendSolidBoundary(solidWall)
#dem.appendSolidBoundary(solidWall2)
# #-------------------------------------------------------------------------------
# # PhysicsPackage : gravity
# #-------------------------------------------------------------------------------
gravity = ConstantAcceleration(a0 = Vector(0.0, 0.0, -1.00),
                               nodeList = nodes1)
packages += [gravity]

# #-------------------------------------------------------------------------------
# # Create boundary conditions.
# #-------------------------------------------------------------------------------
# plane1 = Plane(Vector(0.0, 0.0, 0.0), Vector(  0.0, 0.0, 1.0))
# bc1 = ReflectingBoundary(plane1)
# bcSet = [bc1]

# for p in packages:
#     for bc in bcSet:
#         p.appendBoundary(bc)

#-------------------------------------------------------------------------------
# Fields and Variables
#-------------------------------------------------------------------------------
numNodeLists = db.numNodeLists
nodeLists = db.nodeLists

position = db.DEMPosition
mass = db.DEMMass
velocity = db.DEMVelocity
H = db.DEMHfield
radius = db.DEMParticleRadius
compositeParticleIndex = db.DEMCompositeParticleIndex

uniqueIndex = db.DEMUniqueIndex
omega = dem.omega

# pure rolling
if mpi.rank == 0 :
    for i in range(nodes.numInternalNodes):
        if i > 0:
            velocity[0][i]=Vector(particleVelocity,0,0)
            omega[0][i]=Vector(0,particleVelocity/particleRadius,0)

#-------------------------------------------------------------------------------
# Initial Conditions
#-------------------------------------------------------------------------------
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
control = SpheralController(integrator,
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
                            periodicWork=periodicWork)
output("control")

#control.redistribute = PeanoHilbertOrderRedistributeNodes(db.maxKernelExtent,workBalance=False)
#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------

if not steps is None:
    control.step(steps)
else:
    control.advance(goalTime, maxSteps)

