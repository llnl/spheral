#-------------------------------------------------------------------------------
# Diametral Compression Test 2D with inflow BC
#-------------------------------------------------------------------------------
#
# Solid FSISPH
#
#ATS:t100 = test(        SELF, "--clearDirectories True --checkError True --fsisph True --goalTime 2.0 --nrSpecimen 15 ", label="Diametral Compression Test FSISPH -- 2-D", np=1, fsisph=True)

from Spheral2d import *
import sys, os
import shutil
from math import *
import mpi
import numpy as np
from SpheralTestUtilities import *
from findLastRestart import *

from CompositeNodeDistribution import *
from GenerateNodeDistribution2d import *
from LatticeSampler import *
from HerzianSolution import HerzianSolution

if mpi.procs > 1:
    from PeanoHilbertDistributeNodes import distributeNodes2d as distributeNodes
else:
    from DistributeNodes import distributeNodes2d as distributeNodes

# Note -- all units in (cm, gm, microsec)

commandLine(

    compressionSpeed = 1.0e-6,   # downward velocity of top bc (cm/usec)

    # ats settings
    checkError=False,             # check error rel to analytic
    tol = 5.0,                   # error toleration (%) average tensile stress
    leaveNoTrace=False,           # delete output dirs

    # Specimen
    rho0 = 3.24,    # g/cc

    # Clamps
    compressionWithInflowBC = True,      # inflow vs periodic work function
    cylindricalClamps = True,           # cylinder vs lattice
    incompressibleClamps = True,         # makes them out of Cu (10x bulk modulus)
    rClamp = 0.165,                      # radius of our clamp

    # Specimen geometry.
    SpecimenDistribution = "conformal",  # node distribution for the rock
    rSpecimen =  0.165,                  # 3.3 mm DIA (Zaytzev 2021)
    nrSpecimen = 20,     	             # number of radial bins
    rotation = 0.0,                      # rotates Specimen distribution (radians)

    # Material parameters 
    eosChoice = "tillotson",       # (gruneisen, tillotson)
    isPressureEquilibrium = False, # initialize eps to give pressure equilibrium
    PoissonsRatio = 0.25,          # possion's ratio
    Y0 = 1000.0 / 1e5,             # yield strength  (Zaytsev 2021)
    Yultimate = 0.372 ,            # ultimate strength  (Zaytsev 2021)
    etaMin = 0.1,                  # rho/rho0 limit
    etaMax = 6.0,                  # rho/rho0 limit

    # Damage Model
    strengthModel = "collins",      # ("collins", "constant", "null")
    fullyDamagedSpecimen=False,     # treat Specimen as full damaged to start?
    useDamage = True,               # turns damage model on
    strengthInDamage = False,       # this isn't ready right now
    strainType = BenzAsphaugStrain, # (BenzAsphaugStrain, PseudoPlasticStrain, MeloshRyanAsphaugStrain, PlasticStrain, PseudoPlasticStrain)

    # hydro type
    crksph = False,
    fsisph = False,

    # Hydro parameters  
    correctVelocityGradient=True,        # linear order velocity gradient for PSPH SPH or FSISPH
    asph = False,                        # turns on elliptic kernels
    xsph = False,                        # updates position based on averaged velocity
    epsilonTensile = 0.00,               # term to fight the tensile instability 
    nTensile = 4,                        # term to fight the tensile instability 
    densityUpdate = IntegrateDensity,    # (IntegrateDensity, RigorousSumDensity) update method for the density
    compatibleEnergyEvolution = True,    # activate secondary "accounting" loop for rigorous energy conservation
    totalEnergyEvolution = False,        # evolve the total energy instead of specific thermal energy
    
    # FSISPH parameters
    fsiSurfaceCoefficient = 0.00,          # adds additional repulsive force to material interfaces)
    fsiRhoStabilizeCoeff = 0.1,            # coefficient that smooths the density field
    fsiEpsDiffuseCoeff = 0.1,              # explicit diiffusion of the thermal energy
    fsiXSPHCoeff = 0.00,                   # fsi uses multiplier for XSPH instead of binary switch
    fsiInterfaceMethod = ModulusInterface, # (HLLCInterface, ModulusInterface)
    
    # CRKSPH parameters
    correctionOrder = LinearOrder,   # for CRKSPH higher order field approximations
    
    # Artificial Viscosity
    Cl = 1.0,                      # Cl linear coefficient for av -- None allows spheral to decide
    Cq = 2.0,                      # Cq quadratic coefficient for av -- None allows spheral to decide
    epsilon2 = 1e-30,              # denominator term in balsara correct to prevent singularity
    balsaraCorrection=False,       # do we want the balsara shear correction
    linearInExpansion = None,      # can run the linear portion when nodes are expanding
    quadraticInExpansion = None,   # same for quadratic term
    etaCritFrac = None,            

    # smoothing scale parameters
    kernelOrder = None,       # for b spline kernels, None defaults to WendlandC2 
    nPerh = 3.01,             # number of neighbors / smoothing length     
    HEvolution =  IdealH, # (IdealH , IntegrateH) update method for smoothing kernel
    iterateInitialH = False,  # to calc initial ideal H in controller constructor
     
    # Times, and simulation control.
    cfl = 0.35,               #
    steps = None,             #
    goalTime = 10.0,          # usec -- final time
    dt = 1e-4,                # usec -- initial time step
    dtMin = 1e-4,             # usec -- minimum time step size
    dtMax = 1000.0,           # usec -- maximum time step size
    dtGrowth = 2.0,           # multiplication factor for timestep growth
    dtcheck = True,           # Should the time integrator check the timestep throughout a cycle?
    verbosedt = False,        # will print stats about what controls time-step
    maxSteps = None,          # do we want to limit the total number of cycles
    statsStep = 10,           # 
    redistributeStep = 100,   # cycle for redistribution amongst processors to get more even loading
    restartStep = 5000,       # cycle to write a restart file
    restoreCycle = -1,        # set a cycle to read in an old restart file
    inflowvelfreq = 50,       # How often we measure the Specimen velocity and slow inflow
    historyfreq = 1000,       # how often we write to the history file

    # Output
    graphics = False,
    vizDerivs = True,                    # output derivatives in viz dump
    vizTime = 10.0,                      # usec -- visit output 
    vizCycle = None,                     # can also output visit silo files based on cycles
    vizName = "DiametralTest",           # name for silo files
    baseDir = "dumps-DiametralTest",     # name for output directory this should be /p/lustre/username/baseName
    clearDirectories = False,            # do we want to delete all old files in baseDir upon restart
)


title("2d diametral compression test.")

#--------------------------------------------------------------------------------
# limits on inputs
#--------------------------------------------------------------------------------

eosChoice = eosChoice.lower()
strengthModel = strengthModel.lower()
SpecimenDistribution = SpecimenDistribution.lower()

# valid options for things
assert not (fsisph and crksph)
assert eosChoice in (["gruneisen","tillotson"])
assert strengthModel in (["constant", "collins"]) 
assert SpecimenDistribution in (["lattice","conformal"])

assert not (useDamage and strengthModel=="null") # if were using damage strength must be used

assert fsiXSPHCoeff >= 0.0 
assert fsiEpsDiffuseCoeff >= 0.0 
assert fsiRhoStabilizeCoeff >= 0.0 
assert fsiSurfaceCoefficient >= 0.0 

assert useDamage

#--------------------------------------------------------------------------------
# Path and File Name
#--------------------------------------------------------------------------------
hydroname = "SPH"
if crksph:
    hydroname = "CRK"+hydroname
elif fsisph:
    hydroname = "FSI"+hydroname
if asph:
    hydroname = "A" + hydroname


# Restart and output files.
dataDir = os.path.join(baseDir,
                       "fullyDamage=%s" % fullyDamagedSpecimen,
                       "nrSpecimen%s" % nrSpecimen,
                       "rot%s" % rotation,
                       hydroname)
restartDir = os.path.join(dataDir, "restarts", "proc-%04i" % mpi.rank)
vizDir = os.path.join(dataDir, "viz")
restartBaseName = os.path.join(restartDir, "Airburst-TestProblem02")

if mpi.rank == 0:
    if clearDirectories and os.path.exists(dataDir):
        shutil.rmtree(dataDir)
    if not os.path.exists(dataDir):
        os.makedirs(dataDir)
    if not os.path.exists(vizDir):
        os.makedirs(vizDir)
    if not os.path.exists(restartDir):
        os.makedirs(restartDir)
mpi.barrier()
if not os.path.exists(restartDir):
    os.makedirs(restartDir)
mpi.barrier()

#-------------------------------------------------------------------------------
# Our Kernel 
#-------------------------------------------------------------------------------
if kernelOrder is not None:
    WT = TableKernel(NBSplineKernel(kernelOrder), 1000)
else:
    WT = TableKernel(WendlandC2Kernel(), 1000)
output("WT")


#-------------------------------------------------------------------------------
# SiO2 material properties.
#-------------------------------------------------------------------------------
units = CGuS()

if eosChoice == "gruneisen":
    rho0Granite = 2.68               # g/cc
    eosSolidGranite = GruneisenEquationOfState(rho0,                      # ref density (g/cc)
                                               etaMin,      
                                               etaMax,      
                                               C0 = 0.36839,              # cm/usec      
                                               S1 = 1.8954,               # dimensionless
                                               S2 = 0.0,                  # dimensionless
                                               S3 = 0.0,                  # dimensionless
                                               gamma0 = 0.9,              # dimensionless
                                               b = 1.0,                   # dimensionless
                                               atomicWeight = 60.0843,    # atomic weight
                                               constants = units,
                                               minimumPressure = minimumPressure,
                                               )
    T0Granite, eps0Granite = 0.0, 0.0
elif eosChoice == "tillotson":
    eosSolidGranite = TillotsonEquationOfState("granite",
                                               etaMin,
                                               etaMax,
                                               units)
    eosSolidGranite.referenceDensity = rho0
    T0Granite, eps0Granite = 0.0, 0.0 
    rho0Granite=rho0

cs0 = eosSolidGranite.soundSpeed(rho0Granite, eps0Granite)
Ks0 = rho0Granite*cs0*cs0
Gs = 3.0 * Ks0 * (1.0-2.0*PoissonsRatio)/(2.0*(1.0+PoissonsRatio))
Ys = 9.0*Ks0*Gs/(3.0*Ks0+Gs)

print("Reference (rho0, T0, eps0, c0) = (%g, %g, %g, %g)" % (rho0Granite, T0Granite, eps0Granite,cs0))
print("Reference (K, Y, G, nu) = (%g, %g, %g, %g)" % (Ks0, Ys, Gs, PoissonsRatio))

if strengthModel == 'constant':
    strengthSolidGranite = ConstantStrength(Gs,  # mu = 3.2 GPa
                                            Y0)  # Y0 = 15 MPa
elif strengthModel == 'collins':
    shearModulusModel = ConstantStrength(Gs,     # mu = 3.2 GPa
                                         Y0)     # Y0 = 15 MPa
    strengthSolidGranite=CollinsStrength(shearModulusModel, 
                                         2.0,            # intact coeff of friction
                                         0.7,            # damaged coeff of friction
                                         Y0,             # strength at P=0    (15 MPa)
                                         Yultimate);     # strength as P->oo  (1.5GPa)
else:
    strengthSolidGranite = NullStrength()


eosSpecimen = eosSolidGranite
strengthSpecimen = strengthSolidGranite

rho0Specimen = rho0Granite
c0Specimen = cs0

# Clamps default to granite
eosClamps = eosSolidGranite
strengthClamps = strengthSolidGranite
rho0Clamps = rho0Granite
cS0Clamps = cs0
PoissonsRatioClamps = PoissonsRatio

# set up incompressible Clampss as copper with ridic poisson
if incompressibleClamps:
    eosSolidCu = TillotsonEquationOfState("copper",
                                       etaMin,
                                       etaMax,
                                       units)
    PoissonsRatioCu=0.33
    rho0Cu = eosSolidCu.referenceDensity 
    cs0Cu = eosSolidCu.soundSpeed(rho0Cu, 0.0)
    Ks0Cu = rho0Cu*cs0Cu*cs0Cu
    GsCu = Ks0Cu / ((2.0*(1.0+PoissonsRatioCu))/(3.0*(1.0-2.0*PoissonsRatioCu)))
    print(((2.0*(1.0+PoissonsRatioCu))/(3.0*(1.0-2.0*PoissonsRatioCu))))
    strengthCu = ConstantStrength(GsCu, 10.0*Y0) 

    print("Reference K   (granite, copper) = (%g, %g)" % (Ks0, Ks0Cu))
    print("Reference G   (granite, copper) = (%g, %g)" % (Gs, GsCu))
    print("Reference rho (granite, copper) = (%g, %g)" % (rho0Specimen, rho0Cu))
    print("Reference c   (granite, copper) = (%g, %g)" % (c0Specimen, cs0Cu))

    eosClamps = eosSolidCu
    strengthClamps = strengthCu
    rho0Clamps = rho0Cu
    cS0Clamps = cs0Cu
    PoissonsRatioClamps = PoissonsRatioCu
#-------------------------------------------------------------------------------
# Create the NodeLists. 
#-------------------------------------------------------------------------------
boxsize = 2.0*rSpecimen
hmaxSpecimen = boxsize
nodesSpecimen = makeSolidNodeList("Test Specimen",
                          eosSpecimen,
                          strengthSpecimen,
                          xmin = -boxsize * Vector.one,
                          xmax =  boxsize * Vector.one,
                          hmax = hmaxSpecimen,
                          nPerh = nPerh,
                          rhoMin = etaMin*rho0Granite,
                          rhoMax = etaMax*rho0Granite,
                          kernelExtent = WT.kernelExtent)

nodesDriver = makeSolidNodeList("Test Driver",
                                eosClamps,
                                strengthClamps,
                                xmin = -boxsize * Vector.one,
                                xmax =  boxsize * Vector.one,
                                hmax = hmaxSpecimen,
                                nPerh = nPerh,
                                rhoMin = etaMin*rho0Granite,
                                rhoMax = etaMax*rho0Granite,
                                kernelExtent = WT.kernelExtent)

nodesBase = makeSolidNodeList("Test Base",
                              eosClamps,
                              strengthClamps,
                              xmin = -boxsize * Vector.one,
                              xmax =  boxsize * Vector.one,
                              hmax = hmaxSpecimen,
                              nPerh = nPerh,
                              rhoMin = etaMin*rho0Granite,
                              rhoMax = etaMax*rho0Granite,
                              kernelExtent = WT.kernelExtent)


nodeListSet = [nodesSpecimen,nodesDriver,nodesBase]

#-------------------------------------------------------------------------------
# Generate and distribute the nodes. 
#-------------------------------------------------------------------------------
if compressionWithInflowBC:
    theta = pi
    ymaxLattice = 3*rSpecimen
    nx = 4*nrSpecimen
    ny = 3*nrSpecimen
else:
    theta = 2*pi
    ymaxLattice = 3*rSpecimen
    nx = 4*nrSpecimen
    ny = 3*nrSpecimen

# Driver
#-----------------------------
if cylindricalClamps:
    radRatio = rClamp/rSpecimen

    # we need to know where inflow is if were doing that
    yLowerBC = -rSpecimen*(1.0+radRatio)
    yUpperBC =  rSpecimen*(1.0+radRatio)

    generatorDriver = GenerateNodeDistribution2d(nRadial = int(nrSpecimen*radRatio), 
                                                 nTheta = int(nrSpecimen*radRatio),
                                                        rho = rho0Clamps,
                                                        distributionType = "constantDTheta",
                                                        theta = theta,
                                                        rotation = pi,
                                                        rmax = rClamp,
                                                        rmin = 0.0,
                                                        offset=[0.0,(1.0+radRatio)*rSpecimen],
                                                        nNodePerh = nPerh)

    generatorBase = GenerateNodeDistribution2d(nRadial = int(nrSpecimen*radRatio), 
                                               nTheta = int(nrSpecimen*radRatio),
                                                        rho = rho0Clamps,
                                                        distributionType = "constantDTheta",
                                                        theta = theta,
                                                        rotation = 0.0,
                                                        rmax = rClamp,
                                                        rmin = 0.0,
                                                        offset=[0.0,-(1.0+radRatio)*rSpecimen],
                                                        nNodePerh = nPerh)
else:
    yLowerBC = -rSpecimen-ymaxLattice
    yUpperBC =  rSpecimen+ymaxLattice
    radRatio = 1.0
    generatorDriver = GenerateNodeDistribution2d(nx, ny,
                                                        rho = rho0Clamps,
                                                        distributionType = "lattice",
                                                        xmax = [2*rSpecimen,ymaxLattice],
                                                        xmin = [-2*rSpecimen,0.0],
                                                        offset=[0.0,rSpecimen],
                                                        nNodePerh = nPerh)

    generatorBase = GenerateNodeDistribution2d(nx, ny,
                                                        rho = rho0Clamps,
                                                        distributionType = "lattice",
                                                        xmax = [2*rSpecimen,0.0],
                                                        xmin = [-2*rSpecimen,-ymaxLattice],
                                                        offset=[0.0,-rSpecimen],
                                                        nNodePerh = nPerh)
# Specimen
#-----------------------------
generatorSpecimen = GenerateNodeDistribution2d(nRadial = nrSpecimen, nTheta = nrSpecimen,
                                                        rho = rho0Granite,
                                                        distributionType = "constantDTheta",
                                                        theta = 2.0*pi,
                                                        rotation = rotation,
                                                        rmax = rSpecimen,
                                                        rmin = 0.0,
                                                        nNodePerh = nPerh)

distributeNodes((nodesSpecimen,generatorSpecimen),(nodesDriver,generatorDriver),(nodesBase,generatorBase))


#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node lists.
#-------------------------------------------------------------------------------
db = DataBase()
for n in nodeListSet:
    db.appendNodeList(n)
    print(n.name)
del n
nodeLists = db.nodeLists

output("db")
output("db.nodeLists")
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Physics Package : construct the hydro physics object.
#-------------------------------------------------------------------------------
packages = []
periodicWork = []

if crksph:
    hydro = CRKSPH(dataBase = db,
                   cfl = cfl,
                   useVelocityMagnitudeForDt = False,
                   densityUpdate = densityUpdate,
                   order = correctionOrder,
                   compatibleEnergyEvolution = compatibleEnergyEvolution,
                   evolveTotalEnergy = totalEnergyEvolution,
                   XSPH = xsph,
                   ASPH = asph)

elif fsisph: 
    q = LimitedMonaghanGingoldViscosity(Cl,Cq)   
    hydro = FSISPH(dataBase = db,
                   Q=q,
                   W = WT,
                   cfl = cfl,
                   surfaceForceCoefficient = fsiSurfaceCoefficient,                       
                   densityStabilizationCoefficient = fsiRhoStabilizeCoeff,                 
                   specificThermalEnergyDiffusionCoefficient = fsiEpsDiffuseCoeff,  
                   xsphCoefficient = fsiXSPHCoeff,
                   interfaceMethod = fsiInterfaceMethod,      
                   linearCorrectGradients = correctVelocityGradient,
                   compatibleEnergyEvolution = compatibleEnergyEvolution,  
                   evolveTotalEnergy = totalEnergyEvolution,         
                   ASPH = asph,
                   HUpdate=HEvolution,
                   epsTensile = epsilonTensile,
                   nTensile = nTensile)

else:
    hydro = SPH(dataBase = db,
                W = WT,
                cfl = cfl,
                useVelocityMagnitudeForDt = False,
                densityUpdate = densityUpdate,
                compatibleEnergyEvolution = compatibleEnergyEvolution,
                evolveTotalEnergy = totalEnergyEvolution,
                XSPH = xsph,
                ASPH = asph,
                HUpdate=HEvolution,
                epsTensile = epsilonTensile,
                nTensile = nTensile)

if not epsilon2 is None:
    hydro.Q.epsilon2=epsilon2              
if not Cl is None:
    hydro.Q.Cl = Cl
if not Cq is None:
    hydro.Q.Cq = Cq
if balsaraCorrection:
    hydro.Q.balsaraShearCorrection=True
if not linearInExpansion is None:
    hydro.Q.linearInExpansion = linearInExpansion
if not quadraticInExpansion is None:
    hydro.Q.quadraticInExpansion = quadraticInExpansion
if hasattr(hydro.Q, "etaCritFrac") and not etaCritFrac is None:
    hydro.Q.etaCritFrac = etaCritFrac

output("hydro")
output("hydro.cfl")
output("hydro.useVelocityMagnitudeForDt")
output("hydro.densityUpdate")
output("hydro.HEvolution")
if hasattr(hydro, "correctionOrder"):
    output("hydro.correctionOrder")
if hasattr(hydro, "volumeType"):
    output("hydro.volumeType")
output("hydro.Q")
output("hydro.Q.Cl")
output("hydro.Q.epsilon2")
output("hydro.Q.Cq")
output("hydro.Q.balsaraShearCorrection")
output("hydro.Q.linearInExpansion")
output("hydro.Q.quadraticInExpansion")
if hasattr(hydro.Q, "etaCritFrac"):
    output("hydro.Q.etaCritFrac")

packages += [hydro]


#-------------------------------------------------------------------------------
# Physics Package :  construct the damage model.
#-------------------------------------------------------------------------------
damageModelSpecimen = GradyKippTensorDamage("granite",
                                            nodeList = nodesSpecimen,
                                            strainAlgorithm = strainType,
                                            kernel = WT,
                                            units = units)
            
packages.append(damageModelSpecimen)

#-------------------------------------------------------------------------------
# Create initial conditions.
#-------------------------------------------------------------------------------
nodesDriver.velocity(VectorField("initial velocity driver", nodesDriver,Vector(0.0,-compressionSpeed)))
nodesBase.velocity(VectorField("initial velocity base", nodesBase,Vector(0.0,compressionSpeed)))


#-------------------------------------------------------------------------------
# How are we compressing the specimen? (BC vs periodic work)
#-------------------------------------------------------------------------------
# the bc method uses hemispherical clamps with an inflow conditions on the
# top chuck to drive the system. The lower chuck has a reflection condition. 
# This method results in a predictable total displacement from the compression
# and is prefered when the total displacement is small. For large compression
# turn this feature off. 
#-------------------------------------------------------------------------------
if compressionWithInflowBC:
    yp1 = Plane(Vector(0.0,yUpperBC), Vector( 0.0, -1.0))
    yp2 = Plane(Vector(0.0,yLowerBC), Vector( 0.0, 1.0))
    bcUpper = InflowOutflowBoundary(db,yp1)
    bcLower = InflowOutflowBoundary(db,yp2)

    packages.extend([bcUpper,bcLower])

    bcSet = [bcUpper,bcLower]
    for p in packages:
        for bc in bcSet:
            p.appendBoundary(bc)
else:
    class Squeezer:
        def __init__(self,nodes,vNodes,criterion):
            self.nodes = nodes
            self.vNodes = vNodes
            self.criterion = criterion

        def squeeze(self,cycle,time,dt):
            vel = self.nodes.velocity()
            pos = self.nodes.positions()
            for i in range(self.nodes.numInternalNodes):
                if self.criterion(pos[i].y):
                    vel[i].y = self.vNodes

    def critUpper(yPos):
        crit = False
        if yPos>rSpecimen*(1.0+radRatio):
            crit = True
        return crit

    def critLower(yPos):
        crit = False
        if yPos<-rSpecimen*(1.0+radRatio):
            crit = True
        return crit

    res = rSpecimen/nrSpecimen
    unstoppableForce = Squeezer(nodesDriver,-compressionSpeed,critUpper)
    immovableObject = Squeezer(nodesBase,compressionSpeed,critLower)

    periodicWork += [(unstoppableForce.squeeze,1)]
    periodicWork += [(immovableObject.squeeze,1)]


#-------------------------------------------------------------------------------
# Construct a time integrator.
#-------------------------------------------------------------------------------
integrator = CheapSynchronousRK2Integrator(db)
for p in packages:
    integrator.appendPhysicsPackage(p)
integrator.lastDt = dt
integrator.verbose = verbosedt
integrator.allowDtCheck = dtcheck
if dtMin:
    integrator.dtMin = dtMin
if dtMax:
    integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
integrator.cullGhostNodes = False  #  <--- Fix me!
output("integrator")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.dtGrowth")
output("integrator.rigorousBoundaries")
output("integrator.verbose")
output("integrator.allowDtCheck")

#-------------------------------------------------------------------------------
# PeriodicWork: zap the net x velocity.
#-------------------------------------------------------------------------------
class removeTransverseVelocity:
    def __init__(self,nodes):
        self.nodes = nodes
    def apply(self,cycle,time,dt):
        vel = self.nodes.velocity()
        sumVelx = 0.0
        for i in range(self.nodes.numInternalNodes):
            sumVelx += vel[i].x

        totalSumVelx =  mpi.allreduce(sumVelx,mpi.SUM)
        totalNodes = mpi.allreduce(self.nodes.numInternalNodes,mpi.SUM)
        if totalNodes>0:
            averageVelx = totalSumVelx/totalNodes
            vel -= Vector(averageVelx,0.0)

killVelx = removeTransverseVelocity(nodesSpecimen)
periodicWork += [(killVelx.apply,1)]


#-------------------------------------------------------------------------------
# PeriodicWork: writeOut our strain and oxx
#-------------------------------------------------------------------------------
class loadCurveStorage:

    def __init__(self,
                 vel,
                 dia,
                 db,
                 nodes,
                 pressure):

        self.hist = []
        self.strainHist2 = [0.0]
        self.strainHist = [0.0]
        self.stressHist = [0.0]
        self.timeHist = [0.0]

        self.vel = vel
        self.diameter = dia
        self.delta = 0.0
        self.timeLast = 0.0
        self.sigma = nodes.deviatoricStress()
        self.position = nodes.positions()
        self.pressure = pressure
        self.nodes = nodes
        self.db = db

        self.maxy0 = self.maxY()
        self.miny0 = self.minY()

    def maxY(self):
        maxy = 0.0
        nodeLists = self.db.nodeLists
        for nodeList in nodeLists:
            positions = nodeList.positions()
            for position in positions:
                maxy = max(maxy,position.y)
        maxy = mpi.allreduce(maxy,mpi.MAX)
        return maxy

    def minY(self):
        miny = 0.0
        nodeLists = self.db.nodeLists
        for nodeList in nodeLists:
            positions = nodeList.positions()
            for position in positions:
                miny = min(miny,position.y)
        miny = mpi.allreduce(miny,mpi.MIN)
        return miny

    def readSigmaxx(self):
        Syy = []
        for i in range(self.nodes.numInternalNodes):
            if abs(self.position[i].x)<2.0*rSpecimen/nrSpecimen and abs(self.position[i].y)<0.5*rSpecimen:
                Syy.append(float(self.sigma[i].xx-self.pressure[i]))

        Syyreduced = mpi.allreduce(Syy,mpi.SUM)
        sigmaxx = sum(Syyreduced)/len(Syyreduced)
        return sigmaxx

    def storeLoadCurve(self,cycle,time,dt):

        dtNew = time-self.timeLast
        self.delta = 2.0*time*self.vel
        straini = self.delta/self.diameter
        sigmaxx = self.readSigmaxx()

        # other form of strain calc
        maxy = self.maxY()
        miny = self.minY()
        delta2 = (self.maxy0 - self.miny0) - (maxy-miny)
        straini2 = delta2/self.diameter

        self.timeLast = time
        self.timeHist.append(time)
        self.strainHist2.append(100.0*straini2)
        self.strainHist.append(100.0*straini)
        self.stressHist.append(sigmaxx*1e5)
        self.hist.append([time,straini,sigmaxx])


LCS = loadCurveStorage(compressionSpeed, 
                       2.0*rSpecimen, 
                       db,
                       nodesSpecimen,
                       hydro.pressure[2])

periodicWork += [(LCS.storeLoadCurve ,1)]


header_label='"x(km)" "y(km)" "z(km)"  "P (psi)" "Sxx" "Syy" "Sxy"'
P = hydro.pressure
Sigma =db.solidDeviatoricStress
def eulerianSampleVars(j,i):
    
    return (P(j,i),
           Sigma(j,i).xx,
           Sigma(j,i).yy,
           Sigma(j,i).xy)

res_clip = 2.0*rSpecimen/float(nrSpecimen)

#-------------------------------------------------------------------------------
# PeriodicWork: sampling function
#-------------------------------------------------------------------------------
samplers = [] # list of sample types (currently only line segments are supported)
LSxx = LineSegment(p0=Vector(0.0,-rSpecimen), p1=Vector(0.0,rSpecimen),
           resolution=rSpecimen/float(nrSpecimen),
           N_Fields = 4,
           histDir=dataDir+'/vertical_',
           header_label=header_label)
samplers.append(LSxx)

LSyy = LineSegment(p0=Vector(-rSpecimen,0.0), p1=Vector(rSpecimen,0.0),
           resolution=rSpecimen/float(nrSpecimen),
           N_Fields = 4,
           histDir=dataDir+'/horizontal_',
           header_label=header_label)
samplers.append(LSyy)

Sample = Sample(db, eulerianSampleVars, samplers )
periodicWork += [(Sample.sample,100)]


#-------------------------------------------------------------------------------
# Build the controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            iterateInitialH=iterateInitialH,
                            statsStep = statsStep,
                            restoreCycle = restoreCycle,
                            restartStep = restartStep,
                            redistributeStep = redistributeStep,
                            restartBaseName = restartBaseName,
                            vizBaseName = vizName,
                            vizDir = vizDir,
                            vizDerivs = vizDerivs,
                            vizStep = vizCycle,
                            vizTime = vizTime,
                            periodicWork=periodicWork)

output("control")

#-------------------------------------------------------------------------------
# Advance to completion.
#-------------------------------------------------------------------------------
if not steps is None:
    control.step(steps)
else:
    control.advance(goalTime, maxSteps)


#-------------------------------------------------------------------------------
# Comparison to analytic solution.
#-------------------------------------------------------------------------------

# get our force from sampling
#----------------------------
Sigmaxx=LSyy.DATA_Tot[:,5]-LSyy.DATA_Tot[:,3]
x_sample = LSyy.DATA_Tot[:,0]
y_sample = LSyy.DATA_Tot[:,1]
z_sample = LSyy.DATA_Tot[:,2]
xdiff = x_sample[:-1]-x_sample[1:]
forceSampled = sum((0.5*Sigmaxx[1:]+0.5*Sigmaxx[:-1])*xdiff)


# estimate plane-strain force from displacement
#----------------------------------------------
Ks = rho0Specimen*c0Specimen**2
Kc = rho0Clamps*cS0Clamps**2
cSs = sqrt((Gs)/rho0Specimen)
cSc = sqrt((GsCu)/rho0Clamps)
Es = 3.0*Ks*(1.0-2.0*PoissonsRatio)
Ec = 3.0*Kc*(1.0-2.0*PoissonsRatioClamps)
pi=3.1415

oneOverEstar = (1.0-PoissonsRatio**2)/Es + (1.0-PoissonsRatioClamps**2)/Ec
Estar = 1.0/oneOverEstar

delta = compressionSpeed*goalTime
forceEstimated = 1.0*pi/4.0*Estar*delta

print(" ")
print("force from displacement : %g" % forceEstimated)
print("Force from sampling     : %g" % forceSampled)
print(" ")

# analytic soln
#--------------------------------
analyticSolution = HerzianSolution(forceSampled,rSpecimen)

#-----------------------
# Plot along the x-axis  
#-----------------------
P = []
Sxx = []
Syy = []
Sxy = []
Y = []
X = []
pressure = hydro.pressure
sigma = db.solidDeviatoricStress
pos = db.globalPosition

for i in range(nodesSpecimen.numInternalNodes):
    if abs(pos(2,i).y)<2.0*rSpecimen/nrSpecimen and abs(pos(2,i).x)<0.95*rSpecimen:
        P.append(float(pressure(2,i)))
        Sxx.append(float(sigma(2,i).xx-P[-1]))
        Syy.append(float(sigma(2,i).yy-P[-1]))
        Sxy.append(float(sigma(2,i).xy))
        X.append(float(pos(2,i).x))
        Y.append(float(pos(2,i).y))
    

Preduced = mpi.allreduce(P,mpi.SUM)
Sxxreduced = mpi.allreduce(Sxx,mpi.SUM)
Sxyreduced = mpi.allreduce(Sxy,mpi.SUM)
Syyreduced = mpi.allreduce(Syy,mpi.SUM)
Yreduced = mpi.allreduce(Y,mpi.SUM)
Xreduced = mpi.allreduce(X,mpi.SUM)

FS=14
Nsteps = 100
SxxY0 = Sxxreduced
SyyY0 = Syyreduced
xY0 = [(-rSpecimen + i/float(Nsteps)*(2.0*rSpecimen)) for i in range(Nsteps) ]
yY0 = [0.0 for i in range(Nsteps)]
SxxAnalyticY0 = analyticSolution.sigmaxx(xY0,yY0)
SyyAnalyticY0 = analyticSolution.sigmayy(xY0,yY0)


if mpi.rank==0 and graphics:
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6, 6), dpi=200)

    ax.plot(Xreduced,Sxxreduced,'c.',label='$\sigma_{xx}$ -- FSISPH')
    ax.plot(xY0,SxxAnalyticY0,'b--',label='$\sigma_{xx}$ -- Analytic')
    ax.plot(Xreduced,Syyreduced,'y.',label='$\sigma_{yy}$ -- FSISPH')
    ax.plot(xY0,SyyAnalyticY0,'k-',label='$\sigma_{yy}$ -- Analytic')
    plt.xlabel(r'x-$cm$',fontsize=FS)
    plt.ylabel(r'Stress-$CGuS$',fontsize=FS)
    plt.savefig('DiametralCompression-xaxis.png')
    plt.show()

#-----------------------
# Plot along the y-axis  
#-----------------------
P = []
Sxx = []
Syy = []
Sxy = []
Y = []
X = []
pressure = hydro.pressure
sigma = db.solidDeviatoricStress
pos = db.globalPosition

for i in range(nodesSpecimen.numInternalNodes):
    if abs(pos(2,i).x)<2.0*rSpecimen/nrSpecimen and abs(pos(2,i).y)<0.9*rSpecimen:
        P.append(float(pressure(2,i)))
        Sxx.append(float(sigma(2,i).xx-P[-1]))
        Syy.append(float(sigma(2,i).yy-P[-1]))
        Sxy.append(float(sigma(2,i).xy))
        X.append(float(pos(2,i).x))
        Y.append(float(pos(2,i).y))
    

Preduced = mpi.allreduce(P,mpi.SUM)
Sxxreduced = mpi.allreduce(Sxx,mpi.SUM)
Sxyreduced = mpi.allreduce(Sxy,mpi.SUM)
Syyreduced = mpi.allreduce(Syy,mpi.SUM)
Yreduced = mpi.allreduce(Y,mpi.SUM)
Xreduced = mpi.allreduce(X,mpi.SUM)
Sxxanalytic = analyticSolution.sigmaxx(Xreduced,Yreduced)
Syyanalytic = analyticSolution.sigmayy(Xreduced,Yreduced)

Nsteps = 100
SxxY0 = Sxxreduced
SyyY0 = Syyreduced
yY0 = [(-rSpecimen*0.9 + i/float(Nsteps)*(2.0*rSpecimen*0.9)) for i in range(Nsteps) ]
xY0 = [0.0 for i in range(Nsteps)]
SxxAnalyticY0 = analyticSolution.sigmaxx(xY0,yY0)
SyyAnalyticY0 = analyticSolution.sigmayy(xY0,yY0)


if mpi.rank==0:
    if graphics:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(6, 6), dpi=200)

        ax.plot(Yreduced,Sxxreduced,'c.',label='$\sigma_{xx}$ -- FSISPH')
        ax.plot(yY0,SxxAnalyticY0,'b--',label='$\sigma_{xx}$ -- Analytic')
        #ax.plot(Yreduced,Syyreduced,'y.',label='$\sigma_{yy}$ -- FSISPH')
        #ax.plot(yY0,SyyAnalyticY0,'k-',label='$\sigma_{yy}$ -- Analytic')
        ax.legend(fontsize=FS)
        plt.xlabel(r'y-$cm$',fontsize=FS)
        plt.ylabel(r'Stress-$CGuS$',fontsize=FS)
        plt.savefig('DiametralCompression-yaxis.png')
        plt.show()

    if checkError:

        avgSxxanalytic = sum(Sxxanalytic)/len(Sxxanalytic)
        avgSxxreduced = sum(Sxxreduced)/len(Sxxreduced)
        print(avgSxxanalytic)
        print(avgSxxreduced)
        error = 100.0*abs(avgSxxanalytic - avgSxxreduced)/max(abs(avgSxxanalytic),1e-30)
        
        if error > tol:
            raise ValueError("tensile stress error bounds violated (error, error tolerance) = (%g,%g)." % (error,tol))

    if leaveNoTrace:
        os.system("rm -rf "+baseDir)



