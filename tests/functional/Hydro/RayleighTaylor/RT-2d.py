#-------------------------------------------------------------------------------
# This is the basic Rayleigh-Taylor Problem
#-------------------------------------------------------------------------------
import shutil, os, sys, mpi
from math import *
from Spheral2d import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from findLastRestart import *
from GenerateNodeDistribution2d import *
from CompositeNodeDistribution import *
from CentroidalVoronoiRelaxation import *
from HydrostaticReflectingBoundary import HydrostaticReflectingBoundary2d as HydrostaticReflectingBoundary

if mpi.procs > 1:
    from PeanoHilbertDistributeNodes import distributeNodes2d
else:
    from DistributeNodes import distributeNodes2d

title("Rayleigh-Taylor test problem in 2D")

class ExponentialProfile:
    def __init__(self,
                 y1,
                 rho0,
                 alpha):
        self.y1 = y1
        self.rho0 = rho0
        self.alpha = alpha
        return
    def __call__(self, r):
        #if r.y > 1.0:
        #    print self.rho0*exp(self.alpha*(r.y - self.y1))
        return self.rho0*exp(self.alpha*(r.y - self.y1))

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(nx1 = 100,
            ny1 = 100,
            nx2 = 100,
            ny2 = 100,
            refineFactor = 1,
            nybound = 8,       # number of layers in const node bc

            rho0 = 1.0,
            eps0 = 1.0,
            x0 = 0.0,
            x1 = 1.0,
            y0 = 0.0,
            y1 = 1.0,# position of the interface
            y2 = 2.0,
            P1 = 2.5,
            P2 = 2.5,
            vx1 = 0.0,
            vx2 = 0.0,
            freq = 1.0,
            alpha = 0.0025,   # amplitude of displacement
            beta = 5.0,     # speed at which displacement decays away from midline
            S = 10.0,        # density jump at surface
            g0 = -0.5,
            w0 = 0.005,
            sigma = 0.05/sqrt(2.0),
            
            gamma = 5.0/3.0,
            mu = 1.0,
            
            useHydrostaticBoundary = True,

            # kernel options
            HUpdate = IdealH,
            KernelConstructor = WendlandC2Kernel,
            order = 3,
            nPerh = 3.0,
            hmin = 0.0001,
            hmax = 0.5,
            hminratio = 0.1,

            # hydro type
            svph = False,
            crksph = False,
            psph = False,
            fsisph = False,
            gsph = False,
            mfm = False,

            # hydro options
            asph = False,
            xsph = False,
            solid = False,
            filter = 0.0,
            densityUpdate = IntegrateDensity,
            compatibleEnergy = True,
            evolveTotalEnergy = False,
            useVelocityMagnitudeForDt = False,
            correctVelocityGradient = False,
            epsilonTensile = 0.0,
            nTensile = 8,

            # SPH/PSPH options
            gradhCorrection = False,
            
            # svph options
            fcentroidal = 0.0,
            fcellPressure = 0.0,
            linearConsistent = False,

            # FSISPH parameters
            fsiSurfaceCoefficient = 0.00,           # adds additional repulsive force to material interfaces)
            fsiRhoStabilizeCoeff = 0.1,             # coefficient that smooths the density field
            fsiEpsDiffuseCoeff = 0.1,               # explicit diiffusion of the thermal energy
            fsiXSPHCoeff = 0.00,                    # fsi uses multiplier for XSPH instead of binary switch
            fsiInterfaceMethod = ModulusInterface,     # (HLLCInterface, ModulusInterface)
            fsiKernelMethod  = NeverAverageKernels, # (NeverAverageKernels, AlwaysAverageKernels, AverageInterfaceKernels)
    
            # GSPH/MFM parameters
            gsphEpsDiffuseCoeff = 0.0,
            gsphLinearCorrect = True,
            LimiterConstructor = VanLeerLimiter,
            WaveSpeedConstructor = DavisWaveSpeed,
            
            # artificial viscosity options
            Qconstructor = LimitedMonaghanGingoldViscosity,
            Cl = 1.0,
            Cq = 1.0,
            linearInExpansion = False,
            Qlimiter = False,
            balsaraCorrection = False,
            epsilon2 = 1e-2,

            boolReduceViscosity = False,
            nh = 5.0,
            aMin = 0.1,
            aMax = 2.0,
            Qhmult = 1.0,
            
            # artificial conduction options
            bArtificialConduction = False,
            arCondAlpha = 0.5,

            # integrator & options
            IntegratorConstructor = CheapSynchronousRK2Integrator,
            cfl = 0.25,
            goalTime = 2.0,
            steps = None,
            dt = 0.0001,
            dtMin = 1.0e-8,
            dtMax = 0.1,
            dtGrowth = 2.0,
            maxSteps = None,
            statsStep = 10,
            smoothIters = 0,
            domainIndependent = False,
            rigorousBoundaries = False,
            dtverbose = False,
            
            # output options
            vizCycle = None,
            vizTime = 0.1,
            useVoronoiOutput = False,
            clearDirectories = False,
            restoreCycle = None,
            restartStep = 100,
            redistributeStep = 50000,
            checkRestart = False,
            dataDir = "dumps-Rayleigh-Taylor-2d",
            outputFile = None,
            comparisonFile = None,
            
            serialDump = False, #whether to dump a serial ascii file at the end for viz
            )


assert not (compatibleEnergy and evolveTotalEnergy)
assert sum([fsisph,psph,gsph,crksph,svph,mfm])<=1
assert not (fsisph and not solid)
assert not ((mfm or gsph) and (boolReduceViscosity))

# Decide on our hydro algorithm.
hydroname = 'SPH'
useArtificialViscosity=True

if svph:
    hydroname = "SVPH"
elif crksph:
    Qconstructor = LimitedMonaghanGingoldViscosity
    hydroname = "CRK"+hydroname
elif psph:
    hydroname = "P"+hydroname
elif fsisph:
    hydroname = "FSI"+hydroname
elif gsph:
    hydroname = "G"+hydroname
    useArtificialViscosity=False
elif mfm:
    hydroname = "MFM"
    useArtificialViscosity=False
if asph: 
    hydorname = "A"+hydroname
if solid: 
    hydroname = "solid"+hydroname

dataDir = os.path.join(dataDir,
                       "S=%g" % (S),
                       "vx1=%g-vx2=%g" % (abs(vx1), abs(vx2)),
                       hydroname,
                       "densityUpdate=%s" % (densityUpdate),
                       "XSPH=%s" % xsph,
                       "filter=%s" % filter,
                       "%s-Cl=%g-Cq=%g" % (str(Qconstructor).split("'")[1].split(".")[-1], Cl, Cq),
                       "%ix%i" % (nx1, ny1 + ny2),
                       "nPerh=%g-Qhmult=%g" % (nPerh, Qhmult))
restartDir = os.path.join(dataDir, "restarts")
vizDir = os.path.join(dataDir, "visit")
restartBaseName = os.path.join(restartDir, "Rayleigh-Taylor-2d")
vizBaseName = "Rayleigh-Taylor-2d"


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
# Material properties.
#-------------------------------------------------------------------------------
eos = GammaLawGasMKS(gamma, mu)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
if KernelConstructor==NBSplineKernel:
  WT = TableKernel(NBSplineKernel(order), 1000)
else:
  WT = TableKernel(KernelConstructor(), 1000)
output("WT")
kernelExtent = WT.kernelExtent

#-------------------------------------------------------------------------------
# Make the NodeList.
#-------------------------------------------------------------------------------
if solid:
    nodeListConstructor=makeSolidNodeList
else:
    nodeListConstructor=makeFluidNodeList

nodes1 = nodeListConstructor("Low density gas", eos,
                           hmin = hmin,
                           hmax = hmax,
                           hminratio = hminratio,
                           nPerh = nPerh)
nodes2 = nodeListConstructor("High density gas", eos,
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
# functions for ICs.
#-------------------------------------------------------------------------------
eps1=eps0
eps2=eps1/S

lowerDensity = ExponentialProfile(y1,
                                rho0/S,
                                 g0/((gamma - 1.0)*eps1))
upperDensity = ExponentialProfile(y1,
                                rho0,
                                g0/((gamma - 1.0)*eps2))


#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------

if restoreCycle is None:
    nx1 *= refineFactor
    ny1 *= refineFactor
    nx2 *= refineFactor
    ny2 *= refineFactor
    dy = (y1 - y0)/ny1
    
    generator1 = GenerateNodeDistribution2d(nx1, ny1+nybound,
                                            rho = lowerDensity,
                                            distributionType = "xstaggeredLattice",
                                            xmin = (x0,y0-nybound*dy),
                                            xmax = (x1,y1),
                                            nNodePerh = nPerh,
                                            SPH = SPH)
    generator2 = GenerateNodeDistribution2d(nx2, ny2+nybound,
                                            rho = upperDensity,
                                            distributionType = "xstaggeredLattice",
                                            xmin = (x0,y1),
                                            xmax = (x1,y2+nybound*dy),
                                            nNodePerh = nPerh,
                                            SPH = SPH)

    distributeNodes2d((nodes1, generator1),
                      (nodes2, generator2))

    # A helpful method for setting y displacement.
    def dy(ri):
        thpt = alpha*cos(2.0*pi*ri.x*freq)
        return thpt*exp(-beta*abs(ri.y-y1))

    # Finish initial conditions.
    nodes1.specificThermalEnergy(ScalarField("tmp", nodes1, eps0))
    nodes2.specificThermalEnergy(ScalarField("tmp", nodes2, eps0/S))
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
if useArtificialViscosity:
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
if crksph:
    hydro = CRKSPH(dataBase = db,
                   cfl = cfl,
                   filter = filter,
                   epsTensile = epsilonTensile,
                   nTensile = nTensile,
                   compatibleEnergyEvolution = compatibleEnergy,
                   evolveTotalEnergy = evolveTotalEnergy,
                   XSPH = xsph,
                   densityUpdate = densityUpdate,
                   HUpdate = HUpdate)
elif psph:
    hydro = PSPH(dataBase = db,
                 cfl = cfl,
                 W = WT,
                 Q = q,
                 filter = filter,
                 compatibleEnergyEvolution = compatibleEnergy,
                 evolveTotalEnergy = evolveTotalEnergy,
                 densityUpdate = densityUpdate,
                 HUpdate = HUpdate,
                 XSPH = xsph,
                 ASPH = asph)
if fsisph:
    sumDensityNodeListSwitch =[nodes1,nodes2]  
    hydro = FSISPH(dataBase = db,
                   Q=q, 
                   W = WT,
                   cfl = cfl,
                   surfaceForceCoefficient = fsiSurfaceCoefficient,              
                   densityStabilizationCoefficient = fsiRhoStabilizeCoeff,         
                   specificThermalEnergyDiffusionCoefficient = fsiEpsDiffuseCoeff,     
                   xsphCoefficient = fsiXSPHCoeff,
                   interfaceMethod = fsiInterfaceMethod,
                   kernelAveragingMethod = fsiKernelMethod,
                   sumDensityNodeLists = sumDensityNodeListSwitch,
                   correctVelocityGradient = correctVelocityGradient,
                   compatibleEnergyEvolution = compatibleEnergy,  
                   evolveTotalEnergy = evolveTotalEnergy,         
                   ASPH = asph,
                   epsTensile = epsilonTensile)
elif gsph:
    limiter = LimiterConstructor()
    waveSpeed = WaveSpeedConstructor()
    solver = HLLC(limiter,waveSpeed,gsphLinearCorrect)
    hydro = GSPH(dataBase = db,
                riemannSolver = solver,
                W = WT,
                cfl=cfl,
                compatibleEnergyEvolution = compatibleEnergy,
                correctVelocityGradient= correctVelocityGradient,
                evolveTotalEnergy = evolveTotalEnergy,
                densityUpdate=densityUpdate,
                XSPH = xsph,
                ASPH = asph,
                epsTensile = epsilonTensile,
                nTensile = nTensile)
elif mfm:
    limiter = LimiterConstructor()
    waveSpeed = WaveSpeedConstructor()
    solver = HLLC(limiter,waveSpeed,gsphLinearCorrect)
    hydro = MFM(dataBase = db,
                riemannSolver = solver,
                W = WT,
                cfl=cfl,
                compatibleEnergyEvolution = compatibleEnergy,
                correctVelocityGradient= correctVelocityGradient,
                evolveTotalEnergy = evolveTotalEnergy,
                densityUpdate=densityUpdate,
                XSPH = xsph,
                ASPH = asph,
                epsTensile = epsilonTensile,
                nTensile = nTensile)
else:
    hydro = SPH(dataBase = db,
                cfl = cfl,
                W = WT,
                Q = q,
                compatibleEnergyEvolution = compatibleEnergy,
                evolveTotalEnergy = evolveTotalEnergy,
                gradhCorrection = gradhCorrection,
                correctVelocityGradient = correctVelocityGradient,
                XSPH = xsph,
                ASPH = asph,
                densityUpdate = densityUpdate,
                HUpdate = HUpdate,
                epsTensile = epsilonTensile,
                nTensile = nTensile)

output("hydro")
output("hydro.kernel()")
output("hydro.cfl")
output("hydro.compatibleEnergyEvolution")
output("hydro.densityUpdate")
output("hydro.HEvolution")

packages = [hydro]

#-------------------------------------------------------------------------------
# Construct the MMRV physics object.
#-------------------------------------------------------------------------------

if boolReduceViscosity and useArtificialViscosity:
    evolveReducingViscosityMultiplier = MorrisMonaghanReducingViscosity(nh,aMin,aMax)
    packages.append(evolveReducingViscosityMultiplier)

#-------------------------------------------------------------------------------
# Construct the Artificial Conduction physics object.
#-------------------------------------------------------------------------------

if bArtificialConduction:
    ArtyCond = ArtificialConduction(WT,arCondAlpha) 
    packages.append(ArtyCond)

#-------------------------------------------------------------------------------
# Construct the gravitational acceleration object.
#-------------------------------------------------------------------------------
nodeIndicies1 = vector_of_int()
nodeIndicies2 = vector_of_int()

for i in range(nodes1.numInternalNodes):
    nodeIndicies1.append(i)
for i in range(nodes2.numInternalNodes):
    nodeIndicies2.append(i)

gravity1 = ConstantAcceleration2d(Vector2d(0.0, g0),
                                  nodes1,
                                  nodeIndicies1)
gravity2 = ConstantAcceleration2d(Vector2d(0.0, g0),
                                  nodes2,
                                  nodeIndicies2)

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

pos = nodes1.positions()
ylow, yhigh = vector_of_int(), vector_of_int()
for i in xrange(nodes1.numInternalNodes):
    if pos[i].y < y0:
        ylow.append(i)

pos = nodes2.positions()
for i in xrange(nodes2.numInternalNodes):
    if pos[i].y > y2:
        yhigh.append(i)

print(yhigh)
ybc1 = ConstantBoundary(db, nodes1, ylow, yp1)
ybc2 = ConstantBoundary(db, nodes2, yhigh, yp2)

bcSet = [ybc1, ybc2, xbc]

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
                            vizGhosts=True,
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
