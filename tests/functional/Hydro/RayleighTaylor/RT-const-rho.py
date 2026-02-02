#ATS:test(SELF, "--CRKSPH=True --nx1=128 --nx2=128 --ny1=256 --ny2=512 --cfl=0.25 --Cl=1.0 --Cq=1.0 --clearDirectories=False --filter=0 --nPerh=1.51 --serialDump=True", label="RT CRK, nPerh=1.5", np=16)
#ATS:test(SELF, "--CRKSPH=True --nx1=128 --nx2=128 --ny1=256 --ny2=512 --cfl=0.25 --Cl=1.0 --Cq=1.0 --clearDirectories=False --filter=0 --nPerh=2.01 --serialDump=True", label="RT CRK, nPerh=2.0", np=16)
#ATS:test(SELF, "--CRKSPH=False --nx1=128 --nx2=128 --ny1=256 --ny2=512 --cfl=0.25 --Cl=1.0 --Cq=1.0 --clearDirectories=False --filter=0 --nPerh=1.51 --serialDump=True", label="RT Spheral, nPerh=1.5", np=16)
#ATS:test(SELF, "--CRKSPH=False --nx1=128 --nx2=128 --ny1=256 --ny2=512 --cfl=0.25 --Cl=0.0 --Cq=0.0 --clearDirectories=False --filter=0 --nPerh=1.51 --serialDump=True", label="RT Spheral-NoQ, nPerh=1.5", np=16)
#ATS:test(SELF, "--CRKSPH=False --nx1=128 --nx2=128 --ny1=256 --ny2=512 --cfl=0.25 --Cl=0.0 --Cq=0.0 --clearDirectories=False --filter=0 --nPerh=1.51  --serialDump=True --compatibleEnergy=False", label="RT TSPH-NoQ, nPerh=1.5", np=16)

#-------------------------------------------------------------------------------
# This is the basic Rayleigh-Taylor Problem
#-------------------------------------------------------------------------------
import shutil
from math import *
from Spheral2d import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from findLastRestart import *
from GenerateNodeDistribution2d import *
from CompositeNodeDistribution import *
from CentroidalVoronoiRelaxation import *

import mpi
import DistributeNodes

title("Rayleigh-Taylor test problem in 2D")

class ExponentialDensity:
    def __init__(self,
                 y1,
                 rho0,
                 alpha):
        self.y1 = y1
        self.rho0 = rho0
        self.alpha = alpha
        return
    def __call__(self, r):
        return self.rho0*exp(self.alpha*(r.y - self.y1))

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(nx1     = 50,
            ny1     = 100,
            nx2     = 50,
            ny2     = 200,
            reso    = 1,    # optional scale modifier for the resolution in all directions
            rho0    = 1.0,
            eps0    = 1.0,
            x0      = 0.0,
            x1      = 1.0,
            y0      = 0.0,
            y1      = 2.0,  # position of the interface
            y2      = 6.0,
            P0      = 1.0,  # pressure at top of simulation (y2)
            freq    = 1.0,
            alpha   = 0.01, # amplitude of displacement
            beta    = 5.0,  # speed at which displacement decays away from midline
            S       = 2.0,  # density jump at surface
            g0      = -2.0, # gravitational acceleration
            
            gamma   = 5.0/3.0,
            mu      = 1.0,
            
            nPerh   = 1.51,
            
            SVPH    = False,
            CRKSPH  = False,
            ASPH    = False,
            SPH     = True,   # This just chooses the H algorithm -- you can use this with CRKSPH for instance.
            filter  = 0.0,   # CRKSPH filtering
            Qconstructor = MonaghanGingoldViscosity,
            #Qconstructor = TensorMonaghanGingoldViscosity,
            linearConsistent = False,
            fcentroidal = 0.0,
            fcellPressure = 0.0,
            boolReduceViscosity = False,
            nh      = 5.0,
            aMin    = 0.1,
            aMax    = 2.0,
            Qhmult  = 1.0,
            Cl      = 1.0,
            Cq      = 1.0,
            linearInExpansion = False,
            Qlimiter = False,
            balsaraCorrection = False,
            epsilon2 = 1e-2,
            hmin    = 0.0001,
            hmax    = 0.5,
            hminratio = 0.1,
            cfl     = 0.5,
            useVelocityMagnitudeForDt = False,
            XSPH = False,
            epsilonTensile = 0.0,
            nTensile = 8,
            
            IntegratorConstructor = CheapSynchronousRK2Integrator,
            goalTime = 5.0,
            steps = None,
            vizCycle = None,
            vizTime = 0.01,
            dt = 0.0001,
            dtMin = 1.0e-8,
            dtMax = 0.1,
            dtGrowth = 2.0,
            maxSteps = None,
            statsStep = 10,
            smoothIters = 0,
            HUpdate = IdealH,
            domainIndependent = False,
            rigorousBoundaries = False,
            dtverbose = False,
            
            densityUpdate = RigorousSumDensity, # VolumeScaledDensity,
            compatibleEnergy = True,            # <--- Important!  rigorousBoundaries does not work with the compatibleEnergy algorithm currently.
            gradhCorrection = False,
            
            useVoronoiOutput = False,
            clearDirectories = False,
            restoreCycle = None,
            restartStep = 100,
            redistributeStep = 500,
            checkRestart = False,
            dataDir = "dumps-Rayleigh-Taylor-2d-constRho",
            outputFile = None,
            comparisonFile = None,
            
            serialDump = False, #whether to dump a serial ascii file at the end for viz
            
            bArtificialConduction = False,
            arCondAlpha = 0.5,
            )

nx1 = nx1*reso
nx2 = nx2*reso
ny1 = ny1*reso
ny2 = ny2*reso

#-------------------------------------------------------------------------------
# Computing and printing the growth rate
#-------------------------------------------------------------------------------
atwood  = (S-1.0)/(S+1.0)
zdot    = sqrt(freq*atwood*abs(g0))

print("\n\n\nzdot = exp({0:3.3e}*t)  <-<-<-<-<-<-<-<-<-<------\n\n\n".format(zdot))



# Decide on our hydro algorithm.
if SVPH:
    if ASPH:
        HydroConstructor = ASVPHFacetedHydro
    else:
        HydroConstructor = SVPHFacetedHydro
elif CRKSPH:
    if ASPH:
        HydroConstructor = ACRKSPHHydro
    else:
        HydroConstructor = CRKSPHHydro
else:
    if ASPH:
        HydroConstructor = ASPHHydro
    else:
        HydroConstructor = SPHHydro

dataDir = os.path.join(dataDir,
                       "S=%g" % (S),
                       "CRKSPH=%s" % CRKSPH,
                       str(HydroConstructor).split("'")[1].split(".")[-1],
                       "densityUpdate=%s" % (densityUpdate),
                       "XSPH=%s" % XSPH,
                       "filter=%s" % filter,
                       "%s-Cl=%g-Cq=%g" % (str(Qconstructor).split("'")[1].split(".")[-1], Cl, Cq),
                       "%ix%i" % (nx1, ny1 + ny2),
                       "nPerh=%g-Qhmult=%g" % (nPerh, Qhmult))
restartDir = os.path.join(dataDir, "restarts")
vizDir = os.path.join(dataDir, "visit")
restartBaseName = os.path.join(restartDir, "Rayleigh-Taylor-2d")
vizBaseName = "Rayleigh-Taylor-2d"

#-------------------------------------------------------------------------------
# CRKSPH Switches to ensure consistency
#-------------------------------------------------------------------------------
if CRKSPH:
    Qconstructor = LimitedMonaghanGingoldViscosity

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
import os, sys
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
# Material properties.
#-------------------------------------------------------------------------------
eos = GammaLawGasMKS(gamma, mu)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
WT = TableKernel(BSplineKernel(), 1000)
WTPi = TableKernel(BSplineKernel(), 1000, Qhmult)
output("WT")
output("WTPi")
kernelExtent = WT.kernelExtent

#-------------------------------------------------------------------------------
# Make the NodeList.
#-------------------------------------------------------------------------------
nodes1 = makeFluidNodeList("High density gas", eos,
                           hmin = hmin,
                           hmax = hmax,
                           hminratio = hminratio,
                           nPerh = nPerh)
nodes2 = makeFluidNodeList("Low density gas", eos,
                           hmin = hmin,
                           hmax = hmax,
                           hminratio = hminratio,
                           nPerh = nPerh)
nodeSet = [nodes1, nodes2]
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
    generator1 = GenerateNodeDistribution2d(nx1, ny1,
                                            rho = rho0/S,
                                            distributionType = "lattice",
                                            xmin = (x0,y0),
                                            xmax = (x1,y1),

                                            nNodePerh = nPerh,
                                            SPH = SPH)
    generator2 = GenerateNodeDistribution2d(nx2, ny2,
                                            rho = rho0,
                                            distributionType = "lattice",
                                            xmin = (x0,y1),
                                            xmax = (x1,y2),
                                            nNodePerh = nPerh,
                                            SPH = SPH)

    if mpi.procs > 1:
        from VoronoiDistributeNodes import distributeNodes2d
    else:
        from DistributeNodes import distributeNodes2d

    distributeNodes2d((nodes1, generator1),
                      (nodes2, generator2))

    # A helpful method for setting y displacement.
    def dy(ri):
        thpt = alpha*cos(2.0*pi*ri.x*freq)
        return thpt*exp(-beta*abs(ri.y-y1))

    # Finish initial conditions.
    eps1 = nodes1.specificThermalEnergy()
    eps2 = nodes2.specificThermalEnergy()
    pos1 = nodes1.positions()
    pos2 = nodes2.positions()
    
    rho1 = rho0/S
    rho2 = rho0
    P01  = P0 + g0*(y1-y2)*(rho2-rho1)
    P02  = P0

    
    for i in range(nodes1.numInternalNodes):
        y = pos1[i].y
        eps1[i] = (P01+g0*rho1*(y-y2))/((gamma-1.0)*rho1)
    for i in range(nodes2.numInternalNodes):
        y = pos2[i].y
        eps2[i] = (P02+g0*rho2*(y-y2))/((gamma-1.0)*rho2)
    #nodes1.specificThermalEnergy(ScalarField("tmp", nodes1, eps0))
    #nodes2.specificThermalEnergy(ScalarField("tmp", nodes2, eps0/S))
    for nodes in (nodes1,nodes2):
        pos = nodes.positions()
        vel = nodes.velocity()
        for i in range(nodes.numInternalNodes):
            pos[i].y += dy(pos[i])

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase()
output("db")
for nodes in nodeSet:
    db.appendNodeList(nodes)
output("db.numNodeLists")
output("db.numFluidNodeLists")

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
if SVPH:
    hydro = HydroConstructor(W = WT, 
                             Q = q,
                             cfl = cfl,
                             useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                             compatibleEnergyEvolution = compatibleEnergy,
                             densityUpdate = densityUpdate,
                             XSVPH = XSPH,
                             linearConsistent = linearConsistent,
                             generateVoid = False,
                             HUpdate = HUpdate,
                             fcentroidal = fcentroidal,
                             fcellPressure = fcellPressure,
                             xmin = Vector(-2.0, -2.0),
                             xmax = Vector(3.0, 3.0))
# xmin = Vector(x0 - 0.5*(x2 - x0), y0 - 0.5*(y2 - y0)),
# xmax = Vector(x2 + 0.5*(x2 - x0), y2 + 0.5*(y2 - y0)))
elif CRKSPH:
    hydro = HydroConstructor(W = WT,
                             WPi = WTPi,
                             Q = q,
                             filter = filter,
                             cfl = cfl,
                             useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                             compatibleEnergyEvolution = compatibleEnergy,
                             XSPH = XSPH,
                             densityUpdate = densityUpdate,
                             HUpdate = HUpdate)
else:
    hydro = HydroConstructor(W = WT,
                             WPi = WTPi,
                             Q = q,
                             cfl = cfl,
                             useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                             compatibleEnergyEvolution = compatibleEnergy,
                             gradhCorrection = gradhCorrection,
                             XSPH = XSPH,
                             densityUpdate = densityUpdate,
                             HUpdate = HUpdate,
                             epsTensile = epsilonTensile,
                             nTensile = nTensile)
output("hydro")
output("hydro.kernel()")
output("hydro.PiKernel()")
output("hydro.cfl")
output("hydro.compatibleEnergyEvolution")
output("hydro.densityUpdate")
output("hydro.HEvolution")

packages = [hydro]

#-------------------------------------------------------------------------------
# Construct the MMRV physics object.
#-------------------------------------------------------------------------------

if boolReduceViscosity:
    evolveReducingViscosityMultiplier = MorrisMonaghanReducingViscosity(nh,aMin,aMax)
    
    packages.append(evolveReducingViscosityMultiplier)

#-------------------------------------------------------------------------------
# Construct the Artificial Conduction physics object.
#-------------------------------------------------------------------------------

if bArtificialConduction:
    #q.reducingViscosityCorrection = True
    ArtyCond = ArtificialConduction(WT,arCondAlpha)
    
    packages.append(ArtyCond)

#-------------------------------------------------------------------------------
# Construct the gravitational acceleration object.
#-------------------------------------------------------------------------------

gravity1 = ConstantAcceleration2d(Vector2d(0.0, g0),
                                  nodes1)
gravity2 = ConstantAcceleration2d(Vector2d(0.0, g0),
                                  nodes2)

packages.append(gravity1)
packages.append(gravity2)

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
xp1 = Plane(Vector(x0, y0), Vector( 1.0, 0.0))
xp2 = Plane(Vector(x1, y0), Vector(-1.0, 0.0))
yp1 = Plane(Vector(x0, y0), Vector(0.0,  1.0))
yp2 = Plane(Vector(x0, y2), Vector(0.0, -1.0))
xbc = PeriodicBoundary(xp1, xp2)
#ybc = PeriodicBoundary(yp1, yp2)
ybc1 = ReflectingBoundary(yp1)
ybc2 = ReflectingBoundary(yp2)
bcSet = [xbc, ybc1, ybc2]
#bcSet = [xbc,ybc1]

for bc in bcSet:
    for p in packages:
        p.appendBoundary(bc)

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

output("integrator")
output("integrator.havePhysicsPackage(hydro)")
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
if useVoronoiOutput:
    import SpheralVoronoiSiloDump
    vizMethod = SpheralVoronoiSiloDump.dumpPhysicsState
else:
    import SpheralPointmeshSiloDump
    vizMethod = SpheralPointmeshSiloDump.dumpPhysicsState
control = SpheralController(integrator, WT,
                            initializeDerivatives = True,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle,
                            redistributeStep = redistributeStep,
                            vizMethod = vizMethod,
                            vizBaseName = vizBaseName,
                            vizDir = vizDir,
                            vizStep = vizCycle,
                            vizTime = vizTime,
                            SPH = SPH)
output("control")

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
if not steps is None:
    control.step(steps)

else:
    control.advance(goalTime, maxSteps)
    control.updateViz(control.totalSteps, integrator.currentTime, 0.0)
    control.dropRestartFile()

if serialDump:
    procs = mpi.procs
    rank = mpi.rank
    serialData = []
    i,j = 0,0
    for i in range(procs):
        for nodeL in nodeSet:
            if rank == i:
                for j in range(nodeL.numInternalNodes):
                    serialData.append([nodeL.positions()[j],3.0/(nodeL.Hfield()[j].Trace()),nodeL.mass()[j],nodeL.massDensity()[j],nodeL.specificThermalEnergy()[j]])
    serialData = mpi.reduce(serialData,mpi.SUM)
    if rank == 0:
        f = open(dataDir + "/serialDump.ascii",'w')
        for i in range(len(serialData)):
            f.write("{0} {1} {2} {3} {4} {5} {6} {7}\n".format(i,serialData[i][0][0],serialData[i][0][1],0.0,serialData[i][1],serialData[i][2],serialData[i][3],serialData[i][4]))
        f.close()
