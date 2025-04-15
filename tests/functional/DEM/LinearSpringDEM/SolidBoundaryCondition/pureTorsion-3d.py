#ATS:DEM3dPT = test( SELF, "--clearDirectories True --checkRestitutionCoefficient True --checkNaturalFrequency True", label="DEM pure torsion test -- 3-D (serial)")

import os, sys, shutil, mpi
from math import *
from Spheral3d import *
from SpheralTestUtilities import *
from findLastRestart import *
from GenerateNodeDistribution3d import *
from GenerateDEMfromSPHGenerator import GenerateDEMfromSPHGenerator3d

import numpy as np

if mpi.procs > 1:
    from PeanoHilbertDistributeNodes import distributeNodes3d
else:
    from DistributeNodes import distributeNodes3d

title("DEM 3d Pure Torsion Test")

# This tests the natural freq and restitution coefficient for Zhang et.al. formulation
#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(numParticlePerLength = 3,                 # number of particles on a side of the box
            normalSpringConstant=100.0,               # spring constant for LDS model
            normalRestitutionCoefficient=0.55,        # restitution coefficient to get damping const
            tangentialSpringConstant=20.857,          # spring constant for LDS model
            tangentialRestitutionCoefficient=0.55,    # restitution coefficient to get damping const
            dynamicFriction = 1.0,                    # static friction coefficient sliding
            staticFriction = 1.0,                     # dynamic friction coefficient sliding
            rollingFriction = 1.05,                   # static friction coefficient for rolling
            torsionalFriction = 1.3,                  # static friction coefficient for torsion
            cohesiveTensileStrength = 0.0,            # units of pressure
            shapeFactor = 0.1,                        # in [0,1] shape factor from Zhang 2018, 0 - no torsion or rolling
            
            particleRadius = 1.0,                     # particle radius
            particleDensity = 2.60,
            particleVelocity = 0.1,

            neighborSearchBuffer = 0.1,               # multiplicative buffer to radius for neighbor search algo
            nPerh = 1.01,


            useSolidBoundary = True,

            # integration
            IntegratorConstructor = VerletIntegrator,
            stepsPerCollision = 50,
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
            dataDir = "dumps-DEM-PureTorsionTest-3d",

            # ats
            checkRestitutionCoefficient = False,
            threshRestitutionCoefficient = 0.05,
            checkNaturalFrequency = False,
            threshNaturalFrequency = 0.1,
            )

#-------------------------------------------------------------------------------
# file things
#-------------------------------------------------------------------------------
testName = "DEM-PureTorsionTest-3d"
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
                                            xmin = (0, 0.0, 0),
                                            xmax = ( particleRadius*5,particleRadius*5, particleRadius*5),
                                            nNodePerh = nPerh)
    
    # transforms particle properties from generator
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


#-------------------------------------------------------------------------------
# PhysicsPackage : gravity
#-------------------------------------------------------------------------------
gravity = ConstantAcceleration(a0 = Vector(0.0, 0.0, -1.00),
                               nodeList = nodes1)
packages.append(gravity)
#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
# implement boundary condition using the DEM packages solid wall feature
if useSolidBoundary:
    
    solidWall = InfinitePlaneSolidBoundary(Vector(0.0, 0.0, 0.0), Vector( 0.0, 0.0, 1.0))
    dem.appendSolidBoundary(solidWall)

# implement boundary condition using Spheral's ghost particle reflection
else:
    bcs = [ReflectingBoundary(Plane(Vector(0.0, 0.0, 0.0), Vector( 0.0, 0.0, 1.0)))]
    for package in packages:
        for bc in bcs:
            package.appendBoundary(bc)

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
# Periodic Work Function: track resitution coefficient
#-------------------------------------------------------------------------------
class OmegaTracker:
    def __init__(self):
        self.maxOmega = 0.0
        self.period = 0.0
        self.omegan = 0.0
    def periodicWorkFunction(self,cycle,time,dt):
        omegai = omega[0][1][2]
        if omegai < self.maxOmega:
            self.maxOmega = omegai
            self.period = time - 15
            self.omegan =  pi / self.period
omegaTracker = OmegaTracker()
periodicWork=[(omegaTracker.periodicWorkFunction,1)]

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

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------

if not steps is None:
    control.step(steps)
else:

    # settle dem particles on solid bc in grav field
    control.advance(15.0, maxSteps)

    # give 3 particles torsional rotation
    for i in range(nodes.numInternalNodes):
        if i > 0:
            velocity[0][i]=Vector(0,0,0)
            omega[0][i]=Vector(0,0,particleVelocity/particleRadius)

    # run to the goal time
    control.advance(15.0+goalTime, maxSteps)

if checkRestitutionCoefficient or checkNaturalFrequency:

    # get the sliding damping constant
    beta = pi/log(min(max(tangentialRestitutionCoefficient,1e-4),0.9999))
    nu = 2*tangentialSpringConstant/(1+beta*beta)
    mi = particleDensity * 4.0/3.0*np.pi*particleRadius**3 
    C = sqrt(2*mi*nu)

    # get the natural frequency and effective resitution coefficient of torsion
    omegan = sqrt(5*tangentialSpringConstant*shapeFactor*shapeFactor/mi)
    alphan = 5 * C*shapeFactor*shapeFactor/ (2 * mi)
    analyticRestitution = exp(-alphan*pi/sqrt(omegan**2-alphan**2))

    print("")
    print("==============================================================")

    if checkNaturalFrequency:
        print("")
        print(" Checking Torsional Natural Frequency ")
        print("")
        print("    analytic  natural freq : %g" % omegan)
        print("    numerical natural freq : %g" % omegaTracker.omegan)

        relativeErrorNaturalFrequency = (abs(omegaTracker.omegan-omegan)/omegan)

        print("    relative error   : %g" % relativeErrorNaturalFrequency)

        if relativeErrorNaturalFrequency > threshNaturalFrequency:
            raise ValueError(" natural frequency is not within error bounds ")
    
    if checkRestitutionCoefficient:

        numericalRestitutionCoefficient = (-omegaTracker.maxOmega/particleVelocity*particleRadius)
        
        print("")
        print(" Checking Torsional Restitution Coefficient ")
        print("")
        print("    analytic  restitution coefficient : %g" % analyticRestitution)
        print("    numerical restitution coefficient : %g" % numericalRestitutionCoefficient)

        relativeErrorRestitution = (abs(numericalRestitutionCoefficient-analyticRestitution)/analyticRestitution)

        print("    relative error                    : %g" % relativeErrorRestitution)

        if relativeErrorRestitution > threshRestitutionCoefficient:
            raise ValueError(" restitution coefficient is not within error bounds ")
    
    print("==============================================================")
