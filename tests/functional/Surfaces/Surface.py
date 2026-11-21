#ATS:t1 = test(      SELF, "--CRKSPH True --cfl 0.25 --clearDirectories True  --steps=2 --nx=50 --ny=50 --checkAnswer=True --detectSurfaces=True --detectThreshold=0.99 --sweepAngle=pi/4.0 --detectRange=2.0", label="Surface Detection Test -- 2-D (serial)")

from math import *
import mpi
import os, sys, shutil
from Spheral2d import *
from SpheralTestUtilities import *
from findLastRestart import *
import SpheralPointmeshSiloDump
from GenerateNodeDistribution2d import *

title("Surface Detection Test")

class Rejecter(object):
    def __init__(self,radius):
        self.radius = radius
    def __call__(self,x,y,m,H):
        nX = []
        nY = []
        nM = []
        nH = []
        for i in range(len(x)):
            ri = sqrt(x[i]*x[i]+y[i]*y[i])
            if (ri > self.radius):
                nX.append(x[i])
                nY.append(y[i])
                nM.append(m[i])
                nH.append(H[i])
        return nX,nY,nM,nH

class dSurface(object):
    def __init__(self,nodes,db,Kern,Bf,Sf,hydro,file):
        self.nodes  = nodes
        self.db     = db
        self.Kern   = Kern
        self.Bf     = Bf
        self.Sf     = Sf
        self.hydro  = hydro
        self.file   = file
    def __call__(self,cycle,time,dt):
        #self.renormMat()
        self.momentNorm()
    def renormMat(self):
        f = open(self.file, 'w')
        f.write("i\tSi\txi\n")
        self.db.updateNeighbors()
        self.db.updateConnectivityMap(True)
        cm = self.db.connectivityMap()
        for i in range(self.nodes.numInternalNodes):
            xi = self.nodes.positions()[i]
            Hi = self.nodes.Hfield()[i]
            neighbors = cm.connectivityForNode(self.nodes, i)
            Bi = Tensor.zero
            Vi = self.hydro.volume()[0][i]
            for j in neighbors[0]:
                xj  = self.nodes.positions()[j]
                xij = xj-xi
                Hj  = self.nodes.Hfield()[j]
                Vj  = self.hydro.volume()[0][j] # this could be done better
                gWj = Hj*xij.unitVector()*self.Kern.gradValue((Hj*xij).magnitude(),Hj.Determinant())
                Bij = gWj.dyad(xij)*Vj
                Bi += Bij
            Bi = Bi.Inverse()
            Ei = Bi.eigenValues()
            Si = min(abs(Ei[0]),abs(Ei[1]))
            f.write("%d\t%f\t%f\n" % (i,Si,xi.magnitude()))
    def momentNorm(self):
        f = open(self.file, 'w')
        f.write("i\tSi\txi\tSSi\n")
        for i in range(self.nodes.numInternalNodes):
            xi = self.nodes.positions()[i]
            m0i = self.hydro.m0()[0][i]
            m1i = self.hydro.m1()[0][i]
            f.write("%d\t%f\t%f\t%f\n" %(i,m0i,xi.magnitude(),m1i.magnitude()))
                                     
#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(lattice = True,
            nx      = 50,
            ny      = 50,
            rmin    = 0.0,
            rmax    = 1.0,
            nPerh   = 1.01,
            
            rho0    = 1.0,
            eps0    = 0.0,
            gamma   = 5.0/3.0,
            mu      = 1.0,
            rhomin  = 1.0e-8,
            holeRadius = 0.5,
            
            ASPH    = False,
            CRKSPH  = True,
            SPH     = True,
            XSPH    = False,
            filter  = 0,
            
            KernelConstructor = NBSplineKernel,
            order   = 7,
            
            # Hydro
            Qconstructor        = MonaghanGingoldViscosity2d,
            correctionOrder     = LinearOrder,
            Cl                  = 1.0,
            Cq                  = 2.0,
            Qlimiter            = False,
            balsaraCorrection   = False,
            epsilon2            = 1e-4,
            negligibleSoundSpeed = 1e-5,
            csMultiplier        = 0.1,
            hmin                = 0.004,
            hmax                = 10.0,
            hminratio           = 0.1,
            compatibleEnergy    = False,
            gradhCorrection     = False,
            HEvolution          = IdealH,
            sumForMassDensity   = RigorousSumDensity,
            densityUpdate       = RigorousSumDensity,
            HUpdate             = IdealH,
            linearInExpansion   = False,
            volumeType          = RKVoronoiVolume,
            
            # Timestep constraints
            cfl                 = 0.5,
            deltaPhi            = 0.01,
            domainIndependent   = False,
            
            # Integrator
            IntegratorConstructor = CheapSynchronousRK2Integrator,
            goalTime            = 1.0,
            dt                  = 0.0001,
            dtMin               = 1.0e-5,
            dtMax               = 1.0e5,
            dtGrowth            = 2.0,
            maxSteps            = None,
            steps               = None,
            statsStep           = 10,
            redistributeStep    = 500,
            restartStep         = 500,
            restoreCycle        = None,
            smoothIters         = 0,
            rigorousBoundaries  = True,
            dtverbose           = False,
            
            vizCycle            = 1,
            vizTime             = 1.0e5,
            vizMethod           = SpheralPointmeshSiloDump.dumpPhysicsState,
            
            clearDirectories    = False,
            renormFile          = "renorm.txt",
            
            detectSurfaces      = False,
            detectRange         = 2.0,
            sweepAngle          = pi/4.0,
            detectThreshold     = 0.99,
            checkAnswer         = False,
            )

if CRKSPH:
    Qconstructor = LimitedMonaghanGingoldViscosity2d
    if ASPH:
        HydroConstructor = ACRKSPHHydro
    else:
        HydroConstructor = CRKSPHHydro
else:
    if ASPH:
        HydroConstructor = ASPHHydro
    else:
        HydroConstructor = SPHHydro

dataDir = "surface-%i-%i" % (nx,ny)
dataDir = os.path.join(dataDir, "CRK=%s-nPerh=%f" % (CRKSPH,nPerh))
dataDir = os.path.join(dataDir, "Cl=%f-Cq=%f" % (Cl,Cq))
restartBaseName = "%s/SurfaceTest-%i-%i" % (dataDir,nx,ny)

vizDir = os.path.join(dataDir, "visit")
vizBaseName = "SurfaceTest"

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
import os, sys
if mpi.rank == 0:
    if clearDirectories and os.path.exists(dataDir):
        shutil.rmtree(dataDir)
    if not os.path.exists(dataDir):
        os.makedirs(dataDir)
    if not os.path.exists(vizDir):
        os.makedirs(vizDir)
mpi.barrier()

#-------------------------------------------------------------------------------
# If we're restarting, find the set of most recent restart files.
#-------------------------------------------------------------------------------
if restoreCycle is None:
    restoreCycle = findLastRestart(restartBaseName)

#-------------------------------------------------------------------------------
# Material properties.
#-------------------------------------------------------------------------------
eos = GammaLawGasMKS(gamma, mu)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
if KernelConstructor==NBSplineKernel:
    Wbase = NBSplineKernel(order)
else:
    Wbase = KernelConstructor()
WT = TableKernel(KernelConstructor(order), 1000)
WTPi = TableKernel(KernelConstructor(order), 1000)
output('WT')
output('WTPi')
kernelExtent = WT.kernelExtent
output("WT")

nodes1 = makeFluidNodeList("nodes1", eos,
                           hmin = hmin,
                           hmax = hmax,
                           nPerh = nPerh,
                           kernelExtent = kernelExtent,
                           rhoMin = rhomin)

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
pos     = nodes1.positions()
vel     = nodes1.velocity()
mass    = nodes1.mass()
eps     = nodes1.specificThermalEnergy()
H       = nodes1.Hfield()

if restoreCycle is None:
    if lattice == True:
        xmin = (-1.0, -1.0)
        xmax = (1.0, 1.0)
        myRejecter = Rejecter(holeRadius)
        generator = GenerateNodeDistribution2d(nx,ny,rho0,"lattice",
                                               rmin = rmin,
                                               rmax = rmax,
                                               xmin = xmin,
                                               xmax = xmax,
                                               theta = 2*pi,
                                               nNodePerh = nPerh,
                                               SPH = (not ASPH),
                                               rejecter = myRejecter)
    if mpi.procs > 1:
        from VoronoiDistributeNodes import distribueNodes2d
    else:
        from DistributeNodes import distributeNodes2d

    distributeNodes2d((nodes1,generator))
    output("mpi.reduce(nodes1.numInternalNodes, mpi.MIN)")
    output("mpi.reduce(nodes1.numInternalNodes, mpi.MAX)")
    output("mpi.reduce(nodes1.numInternalNodes, mpi.SUM)")

    for nodeID in range(nodes1.numInternalNodes):
        eps[nodeID] = eps0

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase()
output("db")
output("db.appendNodeList(nodes1)")
output("db.numNodeLists")
output("db.numFluidNodeLists")

Bf = db.newFluidTensorFieldList(Tensor.zero, "Normalization")
Sf = db.newFluidScalarFieldList(0.0, "Surface")

#-------------------------------------------------------------------------------
# Construct the artificial viscosity.
#-------------------------------------------------------------------------------
q = Qconstructor(Cl, Cq, linearInExpansion)
q.epsilon2 = epsilon2
q.limiter = Qlimiter
q.balsaraShearCorrection = balsaraCorrection
output("q")
output("q.Cl")
output("q.Cq")
output("q.epsilon2")
output("q.limiter")
output("q.balsaraShearCorrection")
output("q.linearInExpansion")
output("q.quadraticInExpansion")

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
if CRKSPH:
    hydro = HydroConstructor(W = WT,
                             WPi = WTPi,
                             Q = q,
                             filter = filter,
                             cfl = cfl,
                             compatibleEnergyEvolution = compatibleEnergy,
                             XSPH = XSPH,
                             densityUpdate = densityUpdate,
                             correctionOrder = correctionOrder,
                             volumeType = volumeType,
                             HUpdate = HEvolution,
                             detectSurfaces = detectSurfaces,
                             detectThreshold = detectThreshold,
                             sweepAngle = sweepAngle,
                             detectRange = detectRange)
else:
    hydro = HydroConstructor(W = WT,
                             Q = q,
                             cfl = cfl,
                             compatibleEnergyEvolution = compatibleEnergy,
                             evolveTotalEnergy = evolveTotalEnergy,
                             gradhCorrection = gradhCorrection,
                             correctVelocityGradient = correctVelocityGradient,
                             densityUpdate = densityUpdate,
                             XSPH = XSPH,
                             HUpdate = HEvolution)
output("hydro")
output("hydro.kernel()")
output("hydro.PiKernel()")
output("hydro.cfl")
output("hydro.compatibleEnergyEvolution")
output("hydro.XSPH")
output("hydro.densityUpdate")
output("hydro.HEvolution")

packages = [hydro]

#-------------------------------------------------------------------------------
# Construct the surface detection periodic work function
#-------------------------------------------------------------------------------

#ds      = detectSurface(nodes1,db,WT,Bf,Sf,hydro,renormFile)
#dsFreq  = 1

#-------------------------------------------------------------------------------
# Construct a time integrator, and add the one physics package.
#-------------------------------------------------------------------------------
integrator = IntegratorConstructor(db)
for p in packages:
    integrator.appendPhysicsPackage(p)
integrator.lastDt = dt
if dtMin:
    integrator.dtMin = dtMin
if dtMax:
    integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
output("integrator")
output("integrator.havePhysicsPackage(hydro)")
output("integrator.dtGrowth")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")

#-------------------------------------------------------------------------------
# Build the controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle,
                            vizMethod = vizMethod,
                            vizBaseName = "surface-test-%ix%i" % (nx, ny),
                            vizDir = vizDir,
                            vizStep = vizCycle,
                            vizTime = vizTime,
                            SPH = (not ASPH))
output("control")

#control.appendPeriodicWork(ds,dsFreq)

#-------------------------------------------------------------------------------
# Finally run the problem and plot the results.
#-------------------------------------------------------------------------------
if not steps is None:
    control.step(steps)
else:
    control.advance(goalTime,maxSteps)

if checkAnswer:
    sp = hydro.surfacePoint()
    count = 0
    for i in range(nodes1.numInternalNodes):
        if sp[0][i] == 1:
            count += 1
    if not count == 212:
        raise ValueError("The surface detection algorithm failed!")
    else:
        print("Surface Detection PASSED.")
