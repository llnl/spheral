#-------------------------------------------------------------------------------
# Model of the shock water column experiment of Sembian 2016 which was also
# modelled by Xiang and Wang 2017
#-------------------------------------------------------------------------------

import sys, os, shutil
import mpi
import pickle
from math import *
import numpy as np

from Spheral2d import *
from SpheralTestUtilities import *
from findLastRestart import *
from CompositeNodeDistribution import *
from GenerateNodeDistribution2d import *

from LatticeSampler import * 
from NormalShockRelations import *

if mpi.procs > 1:
    from PeanoHilbertDistributeNodes import distributeNodes2d as distributeNodes
else:
    from DistributeNodes import distributeNodes2d as distributeNodes


# Note -- all units in (m, kg, sec)
commandLine(

    # Air Parameters
    rhoAir = 1.17,            # unshocked air density  (Xiang 2017 Tab. 1)
    PAir = 101000.0,          # unshocked air pressure (Xiang 2017 Tab. 1)
    gammaAir = 1.4,           # specific heat ratio air
    M = 2.40,                 # mach number

    # Water Parameters
    rhoWater = 1000.0,        # water density
    PWater = 100000.0,        # water initial pressure
    gammaWater = 6.2,         # water gamma stiffed gas   (Xiang 2017)
    P0Water = 3.43e8,         # water stiffened parameter (Xiang 2017)
    CvWater = 4184,           # water specific heat
    etaMax = 1.1,             # max density ratio
    etaMin = 0.9,             # min density ratio

    # Geometry
    rColumn = 0.011,          # (m) outer radius of bolide  (Sembian 2016 22mm column dia)
    nrColumn = 50,     	      # number of radial bins
    boxwidthmult = 27.0,      # multiplier of rColumn for simulation box width (Sembian 2016 74mm test section)
    boxheightmult = 6.7,      # multiplier of rColumn for simulation box height
    halfDomain = False,       # for 2d are we using symmetry?

    # Kernel 
    KernelConstructor = WendlandC2Kernel, # WendlandC2Kernel, WendlandC4Kernel, WendlandC6Kernel, NBSplineKernel,
    nPerh = 2.75,                         # neighbors per smoothing length       
    kernelOrder = None,                   # 3,5,7 order for NBSpline
    HEvolution =  IdealH,                 # (IdealH, IntegrateH)

    # Hydro type
    crksph = False,
    psph = False,
    gsph = False,
    fsisph = False,

    # general hydro parameters
    asph = False,                            # turn on ellipitic kernels                         
    xsph = False,                            # smoothed velocity for position integration                       
    epsilonTensile = 0.00,                   # coeff for tensile instability term
    nTensile = 4,                            # exponent of tensile instability term
    densityUpdate = IntegrateDensity,        # (IntegrateDensity,RigorousSumDensity),
    useVelocityMagnitudeForDt=False,         # 
    compatibleEnergyEvolution = True,        # turn on secondary loop to rigorously conserve energy 
    totalEnergyEvolution = False,            # evolve total instead of specific thermal energy
    correctVelocityGradient=True,            # linear corrected velocity gradient in SPH (corrected kernels for GSPH and FSISPH)
    
    # fsisph-specific parameters
    fsiSurfaceCoefficient = 0.00,                # magnificaiton factor for interface face
    fsiRhoStabilizeCoeff = 0.1,                  # coeff for HLLC-pressure term style diffusion
    fsiXSPHCoeff=0.00,                           # xsph coeff
    fsiEpsDiffuseCoeff=0.1,                      # coeff second order artificial conduction
    fsiInterfaceMethod = HLLCInterface,          # how do we handle multimat? (HLLCInterface,ModulusInterface,NoInterface)
    fsiKernelMethod = NeverAverageKernels,       # should we avg the kernels? (AlwaysAverageKernels,NeverAverageKernels,AverageInterfaceKernels)
    fsiSumDensity = True,                        # apply sum density to air nodes
    fsiSlides = True,                             # apply slide condition to interface

    # gsph-specific parameters
    gsphEpsDiffuseCoeff = 0.0,                # 1st order artificial conduction coeff
    gsphLinearCorrect = True,                 # True - second order False - first order HLLC
    gsphGradientMethod = MixedMethodGradient, # gradient def (RiemannGradient,HydroAccelerationGradient,SPHGradient,MixedMethodGradient)
    
    # crksph-specific parameters
    volumeType = RKVolumeType.RKVoronoiVolume,   
    correctionOrder = LinearOrder,
    
    # Artificial Viscosity
    Cl = 1.0,                           # linear coeff
    Cq = 2.0,                           # quadratic coeff
    epsilon2 = 1e-30,                   # factor to keep denomenator non-zero
    balsaraCorrection=False,            # balsara on/off
    linearInExpansion = None,           # turn on linear coeff when particles are separating
    quadraticInExpansion = None,        # turn on quadratic coeff when particles are separating
    etaCritFrac = None,
    
    boolReduceViscosity =True,
    nh = 1.00,
    aMin = 0.05,
    aMax = 2.0,

    # Times, and simulation control.
    cfl = 0.2,
    goalTime = 13.6,          # tstar units
    dt = 1e-9,                # sec - mixed units i know lock me up
    dtMin = 1e-12,            # sec
    dtMax = 1000.0,           # sec
    dtGrowth = 2.0,           # multiplication factor
    dtcheck = True,           # Should the time integrator check the timestep throughout a cycle?
    verbosedt = False,        # output time step info
    steps = None,             # should we do a fixed num of steps?
    maxSteps = None,          # stop at this num steps
    statsStep = 10,           # cycle freq to dump stats
    redistributeStep = 100,   # load balance step
    restartStep = 2000,       # save state
    restoreCycle = -1,        # cycle to load state from
    historyfreq = 1000,       # how often we write to the history file

     # Output
    outputTime = 0.1,                     # tstar units dump pkl files
    vizTime = 0.1,                        # tstar units dump silo
    vizCycle = None,                      # cycle freq to dump silos
    vizDerivs = False,                    # dump our derivative fields
    vizName = "Sembian2016_WaterColumn",  # name for our silo files
    restartName = "Sembian2016",          # name for our rst files
    baseDir = "dumps-Sembian2016",        # name for dump dir
    clearDirectories = False,             # wipe out existing files in dump dir on start
)

# grams - centimeters - micro seconds
units = CGuS()

#--------------------------------------------------------------------------------
# File stuff
#--------------------------------------------------------------------------------
title("Sembian2016__ShockWaterColumn")

if crksph:
    hydroname = os.path.join("CRKSPH", str(volumeType), str(correctionOrder))
elif psph:
    hydroname = "PSPH"
elif fsisph:
    hydroname = "FSISPH"
elif gsph:
    hydroname = "GSPH"
else:
    hydroname = "SPH"
if asph:
    hydroname = "A" + hydroname

# Restart and output files.
dataDir = os.path.join(baseDir,
                       "nrColumn%s" % nrColumn,
                       hydroname)
restartDir = os.path.join(dataDir, "restarts", "proc-%04i" % mpi.rank)
vizDir = os.path.join(dataDir, "viz")
restartBaseName = os.path.join(restartDir, restartName)

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
if KernelConstructor==NBSplineKernel:
    assert kernelOrder in (3,5,7)
    WT = TableKernel(KernelConstructor(kernelOrder), 1000)
else:
    WT = TableKernel(KernelConstructor(), 1000) 
kernelExtent = WT.kernelExtent
output("WT")

#--------------------------------------------------------------------------------
# EOS s
#--------------------------------------------------------------------------------
eosAir = GammaLawGas(gamma=gammaAir, 
                     mu=2.0, 
                     constants=units,
                     minimumPressure = 0.0)

eosWater = StiffenedGas(gamma=gammaWater, 
                        P0=P0Water, 
                        Cv=CvWater,
                        constants=units)

#--------------------------------------------------------------------------------
# Simulation Domain -- were going to set up our initial bounds. These will
#                      act as the max bounds within which we will vary things
#                      dynamically
#--------------------------------------------------------------------------------
boxwidth = boxwidthmult * rColumn
boxheight = boxheightmult * rColumn
boxsize = 2.0*max(boxwidth, boxheight)

# use analytic normal shock relations to get post shock state
NSR = NormalShockRelations(M,gammaAir)
rhoRatio = NSR.densityRatio()
pRatio = NSR.pressureRatio()
vRatio = NSR.velocityRatio()

# 1D compression
lengthRatio = pow(rhoRatio,0.5)

# post shock P and rho
PAirShocked = PAir*pRatio
rhoAirShocked = rhoAir*rhoRatio

#calc thermal energies for pressure equilibrium
epsUnshockedAir = PAir/((gammaAir-1.0)*rhoAir)
epsShockedAir = PAirShocked/((gammaAir-1.0)*rhoAirShocked)
epsWater0 = (PWater + gammaWater * P0Water)/(rhoWater*(gammaWater-1.0))

# calc our post shock velocity
c0S = eosAir.soundSpeed(rhoAir,epsUnshockedAir) # un shocked air sound speed
v0S = -c0S*M                                    # velocity of shock
v1S = v0S*vRatio                                # velocity of air behind shock in shock frame
velShockedAir = v1S-v0S                         # air velocity in lab-frame

print("----------------------------------")
print('Shock Properties')
print("(M...............%s" % M)
print("(P/P0............%s" % pRatio)
print("(rho/rho0........%s" % rhoRatio)
print(("inflow velocity..%s" % velShockedAir))

xmax = [  0.4*boxwidth,  0.5*boxheight]
xmin = [ -0.6*boxwidth, -0.5*boxheight]

nx = int((xmax[0]-xmin[0])/rColumn*nrColumn/lengthRatio)
ny = int((xmax[1]-xmin[1])/rColumn*nrColumn/lengthRatio)

ShockX = -1.5*rColumn
xmaxS = [  ShockX,  0.5*boxheight]
xminS = [ -0.6*boxwidth, -0.5*boxheight]
xmaxU = [  0.4*boxwidth,  0.5*boxheight]
xminU = [ ShockX, -0.5*boxheight]

# using a little bit of a fudge factor instead of calc the proper reflect shock state
nxS = int((xmaxS[0]-xminS[0])/rColumn*nrColumn/lengthRatio*1.2)
nyS = int((xmaxS[1]-xminS[1])/rColumn*nrColumn/lengthRatio*1.2)
nxU = int((xmaxU[0]-xminU[0])/rColumn*nrColumn/lengthRatio/lengthRatio*1.2)
nyU = int((xmaxU[1]-xminU[1])/rColumn*nrColumn/lengthRatio/lengthRatio*1.2)

print("----------------------------------")
print(("shock air lattice  : %s, %s" % (nxS, nyS)))
print(("unshock air lattice: %s, %s" % (nxU, nyU)))

if halfDomain:
    xmin[1]=0.0
    ny1=ny1//2

dx1 = (xmax[0]-xmin[0])/nx
dy1 = (xmax[1]-xmin[1])/ny

hmaxColumn = rColumn
hmaxAir = boxheight

#-------------------------------------------------------------------------------
# set our "real" time based on specified t-star
#-------------------------------------------------------------------------------
tstar = 2*rColumn/(velShockedAir)
print("----------------------------------")
print(("t* = %s" % tstar))

print(("goalTime = %s - tstar" % goalTime))
goalTime *= tstar
vizTime *= tstar
outputTime *= tstar
print(("goalTime = %s - t" % goalTime))
print(("vizTime = %s - t" % vizTime))
print("----------------------------------")

#-------------------------------------------------------------------------------
# Create the NodeLists.
#-------------------------------------------------------------------------------
nodesWater = makeSolidNodeList("Water",
                                eosWater,
                                xmin = -boxsize * Vector.one,
                                xmax =  boxsize * Vector.one,
                                hmax = hmaxColumn,
                                nPerh = nPerh,
                                rhoMin = etaMin*rhoWater,
                                rhoMax = etaMax*rhoWater,
                                kernelExtent = WT.kernelExtent)

nodesAir = makeSolidNodeList("Air", eosAir,
                             xmin = -boxsize * Vector.one,
                             xmax =  boxsize * Vector.one,
                             hmax = hmaxAir,
                             nPerh = nPerh,
                             rhoMin = 0.0,
                             rhoMax = 1000.0,
                             kernelExtent = WT.kernelExtent)

nodeListSet = [nodesWater, nodesAir]

#-------------------------------------------------------------------------------
# Create distributions and add them to the nodelists
#-------------------------------------------------------------------------------
class CircularRejecter:
    def __init__(self, origin, radius):
        self.origin = origin
        self.radius = radius
        return
    def __call__(self, x, y, m, H):
        n = len(x)
        assert (len(y) == n and len(m) == n and len(H) == n)
        xnew, ynew, mnew, Hnew = [], [], [], []
        R2 = self.radius**2
        for i in range(n):
            if ((x[i] - self.origin[0])**2 +
                (y[i] - self.origin[1])**2 > R2):
                xnew.append(x[i])
                ynew.append(y[i])
                mnew.append(m[i])
                Hnew.append(H[i])
        return xnew, ynew, mnew, Hnew


# shocked air ahead of water column
generatorAirS = GenerateNodeDistribution2d(nxS,2*(nyS//2),rhoAirShocked,
                                        distributionType='lattice',
                                        xmin = xminS,
                                        xmax = xmaxS,
                                        nNodePerh = nPerh,
                                        rejecter = CircularRejecter(origin = (0.0,0.0),
                                                                radius = (1.01*rColumn)))

# unshocked air surrounding column
generatorAirU = GenerateNodeDistribution2d(nxU,2*(nyU//2),rhoAir,
                                        distributionType='lattice',
                                        xmin = xminU,
                                        xmax = xmaxU,
                                        nNodePerh = nPerh,
                                        rejecter = CircularRejecter(origin = (0.0,0.0),
                                                                radius = (1.01*rColumn)))

# stitch'em
generatorAir = CompositeNodeDistribution(generatorAirU,generatorAirS)

# create our column
generatorWater= GenerateNodeDistribution2d(nRadial = nrColumn, nTheta = nrColumn,
                                            rho = rhoWater,
                                            distributionType = "constantDTheta",
                                            theta = 2.0*pi,
                                            rotation = 0.0,
                                            rmax = rColumn,
                                            rmin = 0.0,
                                            nNodePerh = nPerh)

generators=[(nodesWater,generatorWater),(nodesAir,generatorAir)]
distributeNodes(*tuple(generators))

#-------------------------------------------------------------------------------
# set some initial conditions
#-------------------------------------------------------------------------------

# fields for air
posFieldAir = nodesAir.positions()
epsFieldAir = nodesAir.specificThermalEnergy()
velFieldAir = nodesAir.velocity()
n = nodesAir.numInternalNodes

# air ICs
for i in range((n)):
    if posFieldAir[i].x<ShockX:
        velFieldAir[i].x = velShockedAir
        epsFieldAir[i] = epsShockedAir
    else:
        epsFieldAir[i] = epsUnshockedAir

# water IC
nodesWater.specificThermalEnergy(ScalarField("initial thermal energy", nodesWater, epsWater0))

# jostle to prevent columnation
f = 5e-2*dy1
posFieldAir = [posi + f * Vector(0.0,np.random.rand()-0.5,np.random.rand()-0.5) for posi in posFieldAir]

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node lists.
#-------------------------------------------------------------------------------
db = DataBase()
for n in nodeListSet:
    print((n.name))
    db.appendNodeList(n)
del n
output("db")
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Physics Package: build our hydro
#-------------------------------------------------------------------------------
packages =[]

if crksph:
    hydro = CRKSPH(dataBase = db,
                   cfl = cfl,
                   densityUpdate = RigorousSumDensity,
                   order = correctionOrder,
                   compatibleEnergyEvolution = compatibleEnergyEvolution,
                   evolveTotalEnergy = totalEnergyEvolution,
                   XSPH = xsph,
                   ASPH = asph)
elif fsisph: 
    sumDensityNodeLists = []
    slidePairs = []
    if fsiSumDensity:
        sumDensityNodeLists = [nodesAir]
    if fsiSlides:
        slidePairs = [(nodesAir,nodesWater)]

    slides = makeSlideSurfaces(db,slidePairs)

    hydro = FSISPH(dataBase = db,
                slides=slides,
                W = WT,
                cfl = cfl,
                surfaceForceCoefficient = fsiSurfaceCoefficient,                   
                densityStabilizationCoefficient = fsiRhoStabilizeCoeff,         
                specificThermalEnergyDiffusionCoefficient = fsiEpsDiffuseCoeff,  
                xsphCoefficient = fsiXSPHCoeff,
                interfaceMethod = fsiInterfaceMethod,  
                kernelAveragingMethod = fsiKernelMethod,
                sumDensityNodeLists = sumDensityNodeLists,
                linearCorrectGradients = correctVelocityGradient,
                compatibleEnergyEvolution = compatibleEnergyEvolution,
                evolveTotalEnergy = totalEnergyEvolution,
                HUpdate = HEvolution,
                ASPH = asph,
                epsTensile = epsilonTensile,
                nTensile = nTensile)

elif gsph:
    limiter = VanLeerLimiter()
    waveSpeed = DavisWaveSpeed()
    solver = HLLC(limiter,waveSpeed,gsphLinearCorrect)
    hydro = GSPH(dataBase = db,
                riemannSolver = solver,
                W = WT,
                cfl=cfl,
                specificThermalEnergyDiffusionCoefficient = gsphEpsDiffuseCoeff,
                compatibleEnergyEvolution = compatibleEnergyEvolution,
                correctVelocityGradient= correctVelocityGradient,
                evolveTotalEnergy = evolveTotalEnergy,
                XSPH = xsph,
                ASPH = asph,
                epsTensile = epsilonTensile,
                nTensile = nTensile)         
else: 
    hydro = SPH(dataBase = db,
                W = WT,
                cfl = cfl,
                densityUpdate = IntegrateDensity,
                compatibleEnergyEvolution = compatibleEnergyEvolution,
                evolveTotalEnergy = totalEnergyEvolution,
                XSPH = xsph,
                ASPH = asph,
                epsTensile = epsilonTensile,
                nTensile = nTensile)

if not gsph:
    if not epsilon2 is None:
        hydro.Q.epsilon2=epsilon2              
    if not Cl is None:
        hydro.Q.Cl = Cl
    if not Cq is None:
        hydro.Q.Cq = Cq
    if balsaraCorrection:
        hydro.Q.balsaraShearCorrection=balsaraCorrection
    if not linearInExpansion is None:
        hydro.Q.linearInExpansion = linearInExpansion
    if not quadraticInExpansion is None:
        hydro.Q.quadraticInExpansion = quadraticInExpansion
    if hasattr(hydro.Q, "etaCritFrac") and not etaCritFrac is None:
        hydro.Q.etaCritFrac = etaCritFrac

    output("hydro.Q")
    output("hydro.Q.Cl")
    output("hydro.Q.epsilon2")
    output("hydro.Q.Cq")
    output("hydro.Q.linearInExpansion")
    output("hydro.Q.quadraticInExpansion")
    output("hydro.densityUpdate")

output("hydro")
output("hydro.cfl")
output("hydro.useVelocityMagnitudeForDt")
output("hydro.HEvolution")

packages += [hydro]

#-------------------------------------------------------------------------------
# Physics Package : AV limiters
#-------------------------------------------------------------------------------
if boolReduceViscosity:
    evolveReducingViscosityMultiplier = MorrisMonaghanReducingViscosity(nh,nh,aMin,aMax)
    packages.append(evolveReducingViscosityMultiplier)

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
# We're going to set up reflection conditions on the top and bottom,
# an inlet conditions to inflow nodes and we're going to do a hack
# extrapolation condtion on the outflow BC. The extrapolation BC
# is handled through a periodicWorkFunction
#-------------------------------------------------------------------------------
xp1 = Plane(Vector(*xmin), Vector(  1.0,  0.0))
yp1 = Plane(Vector(*xmin), Vector(  0.0,  1.0))
xp2 = Plane(Vector(*xmax), Vector( -1.0,  0.0))
yp2 = Plane(Vector(*xmax), Vector(  0.0, -1.0))

bcInflow = InflowOutflowBoundary(db, xp1)

packages += [bcInflow]

bcSet = [bcInflow,
         ReflectingBoundary(yp1),
         ReflectingBoundary(yp2)]

for p in packages:
    for bc in bcSet:
        p.appendBoundary(bc)

# bootleg extrapolation condition
class clipOutflow:
    def __init__(self,
                 database,
                 plane):
        self.database = database
        self.plane = plane
    def periodicWorkFunction(self,cycle,time,dt):
        nodeLists = self.database.nodeLists
        for nodeList in nodeLists:
            cullList=[]
            pos = nodeList.positions()
            ni = nodeList.numInternalNodes
            for i in range(ni):
                
                if self.plane.signedDistance(pos[i])<0.0:
                    cullList.append(i)

            nodeList.deleteNodes(vector_of_int(cullList))

extrapolationBoundary=clipOutflow(db,xp2)

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
# Periodic Work Function : track the heigh and width of the droplet
#-------------------------------------------------------------------------------
# track the height and width ignoring stripped mist
class TrackDropletShape:
    def __init__(self,
                 interfaceFraction,
                 dropletNodeList,
                 rDroplet,
                 nrDroplet,
                 tstar,
                 dataDir):

        self.restart = RestartableObject(self)
        self.dataDir = dataDir
        self.sampleWidth = 5.0 * rDroplet/float(nrDroplet) 
        self.tstar = tstar

        self.interfaceFlag = interfaceFraction
        self.dropletNodeList = dropletNodeList
        self.rDroplet = rDroplet

        self.hist_time=[]
        self.hist_tstar=[]
        self.hist_xCOM=[]
        self.hist_width=[]
        self.hist_height=[]
        self.hist_ratio=[]

    def calcDropletWidth(self):

        ni = self.dropletNodeList.numInternalNodes
        pos = self.dropletNodeList.positions()

        xmax = -1.0e9
        xmin =  1.0e9
        for i in range(ni):
            fraci = self.frac[i]
            posi = pos[i]
            if (posi.y < self.sampleWidth and posi.y > -self.sampleWidth) and fraci < self.fracThresh:
                xmax = max(xmax,posi.x)
                xmin = min(xmin,posi.x)
        xmax = mpi.allreduce(xmax,mpi.MAX)
        xmin = mpi.allreduce(xmin,mpi.MIN)
        return xmax-xmin

    def calcDropletCOM(self):

        ni = self.dropletNodeList.numInternalNodes
        pos = self.dropletNodeList.positions()
        mass = self.dropletNodeList.mass()

        msum = 0.0
        xmsum = 0.0
        for i in range(ni):
            msum += mass[i]
            xmsum += pos[i].x*mass[i]

        xmsum = mpi.allreduce(xmsum,mpi.SUM)
        msum = mpi.allreduce(msum,mpi.SUM)
        return xmsum/msum

    def calcDropletHeight(self):

        pos = self.dropletNodeList.positions()

        ymax = -1.0e9
        ymin =  1.0e9
        for i in range(pos.numInternalElements):
            posi = pos[i]
            if self.frac[i] < self.fracThresh:
                ymax = max(ymax,posi.y)
                ymin = min(ymin,posi.y)
        ymax = mpi.allreduce(ymax,mpi.MAX)
        ymin = mpi.allreduce(ymin,mpi.MIN)
        return ymax-ymin

    def setInterfaceFrac(self):
        self.frac = self.interfaceFraction.fieldForNodeList(self.dropletNodeList)

    def storeHistory(self,cycle,time,dt):
        self.setInterfaceFrac()
        xCOM = self.calcDropletCOM()
        width = self.calcDropletWidth()
        height = self.calcDropletHeight()

        self.hist_time.append(time)
        self.hist_tstar.append(time/self.tstar)
        self.hist_xCOM.append(xCOM)
        self.hist_width.append(width)
        self.hist_height.append(height)
        self.hist_ratio.append(height/width)

    def write(self):
        if mpi.rank==0:
            historyFile = os.path.join(self.dataDir,'dropletShapeHistory.pkl')
            f = open(historyFile,'w')
            pickle.dump([self.hist_time,
                         self.hist_tstar, 
                         self.hist_xCOM,
                         self.hist_width,
                         self.hist_height,
                         self.hist_ratio],f)
            f.close()

    def label(self):
        return "DropShapeHistory"

    def dumpState(self, file, path):
        file.writeObject(self.hist_time, path + "/hist_time")
        file.writeObject(self.hist_tstar, path + "/hist_tstar")
        file.writeObject(self.hist_xCOM, path + "/hist_xCOM")
        file.writeObject(self.hist_width, path + "/hist_width")
        file.writeObject(self.hist_height, path + "/hist_height")
        file.writeObject(self.hist_ratio, path + "/hist_ratio")

    def restoreState(self, file, path):
        self.hist_time = file.readObject(path + "/hist_time")
        self.hist_tstar = file.readObject(path + "/hist_tstar")
        self.hist_xCOM = file.readObject(path + "/hist_xCOM")
        self.hist_width = file.readObject(path + "/hist_width")
        self.hist_height = file.readObject(path + "/hist_height")
        self.hist_ratio = file.readObject(path + "/hist_ratio")


#-------------------------------------------------------------------------------
# Periodic Work Function : write out the interface nodes
#-------------------------------------------------------------------------------
class trackInterface:
    def __init__(self,
                 interfaceFlag,
                 dropletNodeList,
                 tstar,
                 dataDir):
        self.tstar = tstar
        self.interfaceFlag = interfaceFlag
        self.dropletNodeList = dropletNodeList
        self.dataDir = os.path.join(dataDir,"interface")
        if not os.path.exists(self.dataDir) and mpi.rank==0:
            os.makedirs(self.dataDir)

    def setInterfaceFlag(self):
        self.flag = self.interfaceFlag.fieldForNodeList(self.dropletNodeList)

    def write(self,cycle,time,dt):
 
        self.setInterfaceFlag()

        ni = self.dropletNodeList.numInternalNodes
        pos = self.dropletNodeList.positions()

        flag = []
        x = []
        y = []
        z = []

        for i in range(ni):
            flagi = self.flag[i]
            if flagi==4:
                flag.append(flagi)
                x.append(pos[i].x)
                y.append(pos[i].y)
                z.append(pos[i].z)

        # share info
        flag = mpi.allreduce(flag,mpi.SUM)
        x = mpi.allreduce(x,mpi.SUM)
        y = mpi.allreduce(y,mpi.SUM)
        z = mpi.allreduce(z,mpi.SUM)

        # write it
        if mpi.rank == 0:
            filename = os.path.join(self.dataDir,'interface_'+str(cycle))
            f=open(filename+'.pkl','w')
            pickle.dump([cycle, time, time/self.tstar, x, y, z, flag], f)    
            f.close()

#-------------------------------------------------------------------------------
# Periodic Work Function : write the pressure along the centerline 
#-------------------------------------------------------------------------------
# were going to write out the centerline pressure
class trackCenterlinePressure:
    def __init__(self,
                 pressure,
                 database,
                 rDroplet,
                 nrDroplet,
                 tstar,
                 dataDir):

        self.sampleWidth = 5.0 * rDroplet/float(nrDroplet) 
        self.tstar = tstar

        self.pressure = pressure
        self.database = database
        self.position = database.fluidPosition
        self.massDensity = database.fluidMassDensity
        self.dataDir = os.path.join(dataDir,"centerlinePressure")
        if not os.path.exists(self.dataDir) and mpi.rank==0:
            os.makedirs(self.dataDir)
    def write(self,cycle,time,dt):
        
        numNodeLists = self.database.numNodeLists
        nodeLists = self.database.nodeLists

        P   = []
        rho = []
        x   = []
        y   = []
        z   = []

        # collect our interface nodes
        for i in range(numNodeLists):

            nodeListi = nodeLists[i]
            ni = nodeListi.numInternalNodes
            
            for j in range(ni):

                yi = self.position(i,j).y                
                if yi < self.sampleWidth and yi > -self.sampleWidth:
                    P.append(self.pressure(i,j))
                    rho.append(self.massDensity(i,j))
                    x.append(self.position(i,j).x)
                    y.append(yi)
                    z.append(self.position(i,j).z)

        # share info
        rho = mpi.allreduce(rho,mpi.SUM)
        P = mpi.allreduce(P,mpi.SUM)
        x = mpi.allreduce(x,mpi.SUM)
        y = mpi.allreduce(y,mpi.SUM)
        z = mpi.allreduce(z,mpi.SUM)

        # write it
        if mpi.rank == 0:
            filename = os.path.join(self.dataDir,'centerlinePressure_'+str(cycle))
            f=open(filename+'.pkl','w')
            pickle.dump([cycle, time, time/self.tstar, x, y, z, P, rho], f)
            f.close()


# create the periodic work functions. These are time-based work functions so we append them
# after the control is constructe (these require fields from fsi's slide package)
if False: #fsisph:
    interfaceFlags=hydro.interfaceFlags
    interfaceWriter = trackInterface(interfaceFlags,nodesWater,tstar,dataDir)
    #shapeTracker = TrackDropletShape(interfaceFlag,nodesWater,rColumn,nrColumn,tstar,dataDir)   
    centerlineWriter = trackCenterlinePressure(hydro.pressure,db,rColumn,nrColumn,tstar,dataDir)

#-------------------------------------------------------------------------------
# Build the controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
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
                            periodicWork = [(extrapolationBoundary.periodicWorkFunction,1)])

control.redistribute = PeanoHilbertOrderRedistributeNodes(db.maxKernelExtent,workBalance=False)

if False: #fsisph:
    control.appendPeriodicTimeWork(interfaceWriter.write, outputTime)
    control.appendPeriodicTimeWork(centerlineWriter.write, outputTime)
    #control.appendPeriodicTimeWork(shapeTracker.storeHistory, outputTime)

output("control")

#-------------------------------------------------------------------------------
# Advance to completion.
#-------------------------------------------------------------------------------
if not steps is None:
    control.step(steps)
else:
    control.advance(goalTime, maxSteps)

# shape tracker dumps width/height at end
#if fsisph:
#    shapeTracker.write()
