#-------------------------------------------------------------------------------
# Diametral Compression Test 2D
#-------------------------------------------------------------------------------
#
# Solid FSISPH
#
#ATS:t100 = test(        SELF, "--clearDirectories True --checkError True --goalTime 5.0 --fsisph True --nrSpecimen 15 ", label="Diametral Compression Test FSISPH -- 2-D", np=8, fsisph=True)

from Spheral2d import *

import sys, os, shutil, mpi
from math import *
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
    checkError=True,             # check error rel to analytic
    tol = 10.0,                  # error toleration (%) average tensile stress
    leaveNoTrace=True,           # delete output dirs

    # Specimen geometry.
    SpecimenDistribution = "conformal",  # node distribution for the rock
    rSpecimen =  0.165,                  # 3.3 mm DIA (Zaytzev 2021)
    nrSpecimen = 20,     	             # number of radial bins
    rotation = 0.0,                      # rotates Specimen distribution (radians)

    # Specimen material 
    materialSpecimen = "granite",
    isPressureEquilibrium = False, # initialize eps to give pressure equilibrium
    PoissonsRatio = 0.25,          # possion's ratio
    Y0 = 1000.0 / 1e5,             # yield strength 
    etaMin = 0.1,                  # rho/rho0 limit
    etaMax = 6.0,                  # rho/rho0 limit

    # Clamps geometry
    rClamp = 0.0825,

    # Clamps material
    materialClamps = "granite",
    PoissonsRatioClamps = 0.25,    # possion's ratio

    # hydro selection
    crksph = False,
    fsisph = False,

    # Hydro parameters
    cfl = 0.25,   
    correctVelocityGradient=True,        # linear order velocity gradient for PSPH SPH or FSISPH
    asph = False,                        # turns on elliptic kernels
    xsph = False,                        # updates position based on averaged velocity
    epsilonTensile = 0.00,               # term to fight the tensile instability 
    nTensile = 4,                        # term to fight the tensile instability 
    densityUpdate = IntegrateDensity,    # (IntegrateDensity, RigorousSumDensity) update method for the density
    compatibleEnergyEvolution = True,    # activate secondary "accounting" loop for rigorous energy conservation
    totalEnergyEvolution = False,        # evolve the total energy instead of specific thermal energy
    
    # FSISPH parameters
    fsiSurfaceCoefficient = 0.00,            # adds additional repulsive force to material interfaces)
    fsiRhoStabilizeCoeff = 0.0,              # coefficient that smooths the density field
    fsiEpsDiffuseCoeff = 0.0,                # explicit diiffusion of the thermal energy
    fsiXSPHCoeff = 0.00,                     # fsi uses multiplier for XSPH instead of binary switch
    fsiInterfaceMethod = ModulusInterface,   # (HLLCInterface, ModulusInterface)

    # CRKSPH parameters
    correctionOrder = LinearOrder,   # for CRKSPH higher order field approximations
    
    # Artificial Viscosity
    Cl = 0.5,                     # Cl linear coefficient for av -- None allows spheral to decide
    Cq = 0.5,                     # Cq quadratic coefficient for av -- None allows spheral to decide
    epsilon2 = 1e-30,              # denominator term in balsara correct to prevent singularity
    balsaraCorrection=False,       # do we want the balsara shear correction
    linearInExpansion = None,      # can run the linear portion when nodes are expanding
    quadraticInExpansion = None,   # same for quadratic term
    etaCritFrac = None,            

    # smoothing scale parameters
    KernelConstructor = WendlandC2Kernel, # WendlandC2Kernel, WendlandC4Kernel, WendlandC6Kernel, NBSplineKernel,
    kernelOrder = None,                   # 3,5,7 for NBSpline Kernels
    nPerh = 3.01,                         # number of neighbors / smoothing length     
    HEvolution =  IntegrateH,             # (IdealH , IntegrateH) update method for smoothing kernel
    iterateInitialH = False,              # to calc initial ideal H in controller constructor
     
    # Times, and simulation control.
    steps = None,
    goalTime = 100.0,         # usec -- final time
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
    vizDerivs = False,                   # output derivatives in viz dump
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

SpecimenDistribution = SpecimenDistribution.lower()

# valid options for things
assert not (fsisph and crksph)
assert SpecimenDistribution in (["lattice","conformal"])
assert fsiXSPHCoeff >= 0.0 
assert fsiEpsDiffuseCoeff >= 0.0 
assert fsiRhoStabilizeCoeff >= 0.0 
assert fsiSurfaceCoefficient >= 0.0 


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
# Interpolation kernels.
#-------------------------------------------------------------------------------
if KernelConstructor==NBSplineKernel:
    assert order in (3,5,7)
    WT = TableKernel(KernelConstructor(order), 1000)
else:
    WT = TableKernel(KernelConstructor(), 1000) 
kernelExtent = WT.kernelExtent
output("WT")

#-------------------------------------------------------------------------------
# SiO2 material properties.
#-------------------------------------------------------------------------------
units = CGuS()

# specimen
eos = TillotsonEquationOfState(materialSpecimen,
                                       etaMin,
                                       etaMax,
                                       units)

rho0 = eos.referenceDensity
eps0 = 0.0
c0 = eos.soundSpeed(rho0, eps0)
K0 = rho0*c0*c0
G = 3.0 * K0 * (1.0-2.0*PoissonsRatio)/(2.0*(1.0+PoissonsRatio))

strength = ConstantStrength(G,  Y0)  

# clamps
eosClamps= TillotsonEquationOfState(materialClamps,
                                     etaMin,
                                     etaMax,
                                     units)

rho0Clamps = eosClamps.referenceDensity 
c0Clamps = eosClamps.soundSpeed(rho0Clamps, 0.0)
K0Clamps = rho0Clamps*c0Clamps*c0Clamps
GClamps = K0Clamps / ((2.0*(1.0+PoissonsRatioClamps))/(3.0*(1.0-2.0*PoissonsRatioClamps)))
strengthClamps = ConstantStrength(GClamps, 10.0*Y0) 

print("Reference K   (specimen, clamps) = (%g, %g)" % (K0, K0Clamps))
print("Reference G   (specimen, clamps) = (%g, %g)" % (G, GClamps))
print("Reference rho (specimen, clamps) = (%g, %g)" % (rho0, rho0Clamps))
print("Reference c   (specimen, clamps) = (%g, %g)" % (c0, c0Clamps))

#-------------------------------------------------------------------------------
# Create the NodeLists. 
#-------------------------------------------------------------------------------
boxsize = 2.0*rSpecimen
hmaxSpecimen = boxsize

nodesSpecimen = makeSolidNodeList("Test Specimen",
                          eos,
                          strength,
                          xmin = -boxsize * Vector.one,
                          xmax =  boxsize * Vector.one,
                          hmax = hmaxSpecimen,
                          nPerh = nPerh,
                          rhoMin = etaMin*rho0,
                          rhoMax = etaMax*rho0,
                          kernelExtent = WT.kernelExtent)

nodesDriver = makeSolidNodeList("Test Driver",
                                eosClamps,
                                strengthClamps,
                                xmin = -boxsize * Vector.one,
                                xmax =  boxsize * Vector.one,
                                hmax = hmaxSpecimen,
                                nPerh = nPerh,
                                rhoMin = etaMin*rho0Clamps,
                                rhoMax = etaMax*rho0Clamps,
                                kernelExtent = WT.kernelExtent)

nodesBase = makeSolidNodeList("Test Base",
                              eosClamps,
                              strengthClamps,
                              xmin = -boxsize * Vector.one,
                              xmax =  boxsize * Vector.one,
                              hmax = hmaxSpecimen,
                              nPerh = nPerh,
                              rhoMin = etaMin*rho0Clamps,
                              rhoMax = etaMax*rho0Clamps,
                              kernelExtent = WT.kernelExtent)


nodeListSet = [nodesSpecimen,nodesDriver,nodesBase]

#-------------------------------------------------------------------------------
# Generate and distribute the nodes. 
#-------------------------------------------------------------------------------

# Driver
#-----------------------------
radRatio = rClamp/rSpecimen

generatorDriver = GenerateNodeDistribution2d(nRadial = int(nrSpecimen*radRatio), nTheta = int(nrSpecimen*radRatio),
                                                        rho = rho0Clamps,
                                                        distributionType = "constantDTheta",
                                                        theta = 2.0*pi,
                                                        rotation = pi,
                                                        rmax = rClamp,
                                                        rmin = 0.0,
                                                        offset=[0.0,(1.0+radRatio)*rSpecimen],
                                                        nNodePerh = nPerh)

generatorBase = GenerateNodeDistribution2d(nRadial = int(nrSpecimen*radRatio), nTheta = int(nrSpecimen*radRatio),
                                                        rho = rho0Clamps,
                                                        distributionType = "constantDTheta",
                                                        theta = 2.0*pi,
                                                        rotation = 0.0,
                                                        rmax = rClamp,
                                                        rmin = 0.0,
                                                        offset=[0.0,-(1.0+radRatio)*rSpecimen],
                                                        nNodePerh = nPerh)

# Specimen
#-----------------------------
generatorSpecimen = GenerateNodeDistribution2d(nRadial = nrSpecimen, nTheta = nrSpecimen,
                                                        rho = rho0,
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
packages=[]

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
    hydro = FSISPH(dataBase = db,
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
                   HUpdate=HEvolution)

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
output("hydro._smoothingScaleMethod.HEvolution")
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
# Create initial conditions.
#-------------------------------------------------------------------------------
nodesDriver.velocity(VectorField("initial velocity driver", nodesDriver,Vector(0.0,-compressionSpeed)))
nodesBase.velocity(VectorField("initial velocity base", nodesBase,Vector(0.0,compressionSpeed)))

#-------------------------------------------------------------------------------
# Construct a time integrator.
#-------------------------------------------------------------------------------
integrator = VerletIntegrator(db)
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
# PeriodicWork:  zap the net x velocity.
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
periodicWork = [(killVelx.apply,1)]

#-------------------------------------------------------------------------------
# PeriodicWork: squeeze the specimen
#-------------------------------------------------------------------------------
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
        positions = self.nodes.positions()
        for position in positions:
            maxy = max(maxy,position.y)
        maxy = mpi.allreduce(maxy,mpi.MAX)
        return maxy

    def minY(self):
        miny = 0.0
        positions = self.nodes.positions()
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


#-------------------------------------------------------------------------------
# PeriodicWork: Sampling function
#-------------------------------------------------------------------------------
header_label='"x(km)" "y(km)" "z(km)"  "P (psi)" "Sxx" "Syy" "Sxy"'
P = hydro.pressure
Sigma =db.solidDeviatoricStress
def eulerianSampleVars(j,i):
    
    return (P(j,i),
           Sigma(j,i).xx,
           Sigma(j,i).yy,
           Sigma(j,i).xy)

res_clip = 2.0*rSpecimen/float(nrSpecimen)

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

# sampled force 
#-------------------
Sigmaxx= LSyy.DATA_Tot[:,5]-LSyy.DATA_Tot[:,3]
x_sample = LSyy.DATA_Tot[:,0]
y_sample = LSyy.DATA_Tot[:,1]
z_sample = LSyy.DATA_Tot[:,2]
xdiff = x_sample[:-1]-x_sample[1:]
forceSampled = sum((0.5*Sigmaxx[1:]+0.5*Sigmaxx[:-1])*xdiff)

# plane strain estimate of force
#--------------------------------
Es = 3.0*K0*(1-2.0*PoissonsRatio)
Ec = 3.0*K0Clamps*(1-2.0*PoissonsRatioClamps)
oneOverEstar = (1.0-PoissonsRatio**2)/Es + (1.0-PoissonsRatioClamps**2)/Ec
Estar = 2.0/oneOverEstar
# penetration depth
delta = compressionSpeed*goalTime
# herzian penetration depth
forceEstimated = pi/4.0*Es*delta
forceEstimated2 = pi/4.0*Es*LCS.strainHist2[-1]*rSpecimen/100.0
forceEstimated1 = pi/4.0*Es*LCS.strainHist[-1]*rSpecimen/100.0


print("   ")
print("---------------------------------------------------------------------------")
print("Specimen  Bulk Modulus: %g" % K0)
print("Clamps    Bulk Modulus: %g" % K0Clamps)
print("Specimen  Poisson's Ratio: %g" % PoissonsRatio)
print("Clamps    Poisson's Ratio: %g" % PoissonsRatioClamps)
print("Specimen  Elastic Modulus: %g" % (Es/(1.0-PoissonsRatio**2)))
print("Clamps    Elastic Modulus: %g" % (Ec/(1.0-PoissonsRatioClamps**2)))
print("Effective Elastic Modulus: %g" % Estar)

print("---------------------------------------------------------------------------")
print("Displacement from nominal compression rate: %g" % delta)
print("Displacement from nominal compression rate: %g" % (LCS.strainHist[-1]*rSpecimen/100.0))
print("measured displacement of specimen         : %g" % (LCS.strainHist2[-1]*rSpecimen/100.0))


print("---------------------------------------------------------------------------")
print("Force from Displacement based on nominal compression rate: %g" % forceEstimated)
print("Force from Displacement based on nominal compression rate: %g" % forceEstimated1)
print("Force from Displacement based on the measured compression: %g" % forceEstimated2)
print("Force from integrating the stress : %g" % forceSampled)
print("---------------------------------------------------------------------------")
print("   ")

# analytic herzian soln for stress field
#----------------------------------------
analyticSolution = HerzianSolution(forceSampled,rSpecimen)

# numerical Solution
#--------------------------------
P = []
Sxx = []
Syy = []
Sxy = []
Y = []
X = []

pressure = hydro.pressure.fieldForNodeList(nodesSpecimen)

for i in range(nodesSpecimen.numInternalNodes):
    pos = nodesSpecimen.positions()
    S = nodesSpecimen.deviatoricStress()
    if abs(pos[i].x)<2.0*rSpecimen/float(nrSpecimen) and abs(pos[i].y)<0.75*rSpecimen:
        P.append(float(pressure[i]))
        Sxx.append(float(S[i].xx-P[-1]))
        Syy.append(float(S[i].yy-P[-1]))
        Sxy.append(float(S[i].xy))
        X.append(float(pos[i].x))
        Y.append(float(pos[i].y))
    

Preduced = mpi.allreduce(P,mpi.SUM)
Sxxreduced = mpi.allreduce(Sxx,mpi.SUM)
Sxyreduced = mpi.allreduce(Sxy,mpi.SUM)
Syyreduced = mpi.allreduce(Syy,mpi.SUM)
Yreduced = mpi.allreduce(Y,mpi.SUM)
Xreduced = mpi.allreduce(X,mpi.SUM)
Sxxanalytic = analyticSolution.sigmaxx(Xreduced,Yreduced)
Syyanalytic = analyticSolution.sigmayy(Xreduced,Yreduced)

if mpi.rank==0:
    if graphics:
        import matplotlib.pyplot as plt
        
        fig, ax = plt.subplots(figsize=(6, 6), dpi=200)
        FS = 14
        factor = 1e5 # cgus -> MPa
        ax.plot(Yreduced,[entry*factor for entry in Sxxreduced],'c.',label='$\sigma_{xx}$ -- FSISPH')
        ax.plot(Yreduced,[entry*factor for entry in Sxxanalytic],'b.',label='$\sigma_{xx}$ -- Analytic')
        ax.plot(Yreduced,[entry*factor for entry in Syyreduced],'y.',label='$\sigma_{yy}$ -- FSISPH')
        ax.plot(Yreduced,[entry*factor for entry in Syyanalytic],'k.',label='$\sigma_{yy}$ -- Analytic')
        ax.legend(fontsize=FS)
        plt.xlabel(r'y-$cm$',fontsize=FS)
        plt.ylabel(r'Stress-$MPa$',fontsize=FS)
        plt.savefig('DiametralCompression-yaxis.png')
        plt.show()

    if checkError:

        avgSxxanalytic = sum(Sxxanalytic)/len(Sxxanalytic)
        avgSxxreduced = sum(Sxxreduced)/len(Sxxreduced)
        error = 100.0*abs(avgSxxanalytic - avgSxxreduced)/max(abs(avgSxxanalytic),1e-30)
        
        print("average analytic Sxx  = %s" % avgSxxanalytic)
        print("average simulated Sxx = %s" % avgSxxreduced)
        if error > tol:
            raise ValueError("tensile stress error bounds violated (error, error tolerance) = (%g,%g)." % (error,tol))

    if leaveNoTrace:
        os.system("rm -rf "+baseDir)
