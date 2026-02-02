#ATS:test(SELF, "--CRKSPH=True --nx1=256 --nx2=256 --ny1=128 --ny2=128 --cfl=0.25 --Cl=1.0 --Cq=1.0 --clearDirectories=True --filter=0 --nPerh=1.51", label="KH CRK, nPerh=1.5", np=16)
#ATS:test(SELF, "--CRKSPH=True --nx1=256 --nx2=256 --ny1=128 --ny2=128 --cfl=0.25 --Cl=1.0 --Cq=1.0 --clearDirectories=True --filter=0 --nPerh=2.01", label="KH CRK, nPerh=2.0", np=16)
#-------------------------------------------------------------------------------
# This is the basic Kelvin-Helmholtz problem as discussed in
# Springel 2010, MNRAS, 401, 791-851.
#-------------------------------------------------------------------------------
import shutil, os, sys
from math import *
from Spheral2d import *
from SpheralTestUtilities import *
#from SpheralGnuPlotUtilities import *
from findLastRestart import *
from GenerateNodeDistribution2d import *
from CompositeNodeDistribution import *
from CentroidalVoronoiRelaxation import *

import mpi
import DistributeNodes

if mpi.procs > 1:
    from PeanoHilbertDistributeNodes import distributeNodes2d
else:
    from DistributeNodes import distributeNodes2d

title("Kelvin-Helmholtz test problem in 2D")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(nx1 = 100,
            ny1 =  50,
            nx2 = 100,
            ny2 =  50,
            
            rho1 = 2.0,
            rho2 = 1.0,
            P1 = 2.5,
            P2 = 2.5,
            vx1 = 0.5,
            vx2 = -0.5,
            vxboost = 0.0,
            vyboost = 0.0,
            freq = 4.0,
            w0 = 0.1,
            sigma = 0.05/sqrt(2.0),

            numNodeLists = 2,  # If 2, makes this a two material problem.

            gamma = 5.0/3.0,
            mu = 1.0,

            # kernel
            HUpdate = IdealH,
            nPerh = 3.0,
            KernelConstructor = WendlandC2Kernel,
            order = 3,
            hmin = 0.0001, 
            hmax = 0.5,
            hminratio = 0.1,

            # hydro 
            svph = False,
            psph = False,
            crksph = False,
            fsisph = False,
            gsph = False,
            mfm = False,
            mfv = False,

            # hydro options
            solid = False,                        # fluid limit of the solid hydro
            asph = False,                         # elliptic kernels
            xsph = False,                         # smoothed position update
            filter = 0.0,                         #
            useVelocityMagnitudeForDt = False,    #
            epsilonTensile = 0.0,                 # tensile correction
            nTensile = 8,                         # tensile correction
            densityUpdate = RigorousSumDensity,   # (RigorousSumDensity, IntegrateDensity)
            compatibleEnergy = True,              # 2nd loop to correct spec. therm. energy for conservation
            evolveTotalEnergy = False,            # integrate total instead of specific energy
            gradhCorrection = True,               # correct for temporal variation in h
            correctVelocityGradient = True,       # linear exact velocity gradient (M correction) (corrected kernesl for GSPH and FSISPH)
            
            # SVPH parameters
            fcentroidal = 0.0,
            fcellPressure = 0.0,

            # FSISPH parameters
            fsiSurfaceCoefficient = 0.00,           # adds additional repulsive force to material interfaces)
            fsiRhoStabilizeCoeff = 0.1,             # coefficient that smooths the density field
            fsiEpsDiffuseCoeff = 0.1,               # explicit diiffusion of the thermal energy
            fsiXSPHCoeff = 0.00,                    # fsi uses multiplier for XSPH instead of binary switch
            fsiInterfaceMethod = HLLCInterface,     # (HLLCInterface, ModulusInterface)
            fsiKernelMethod  = NeverAverageKernels, # (NeverAverageKernels, AlwaysAverageKernels, AverageInterfaceKernels)
    
            # GSPH/MFM/MFV parameters 
            gsphEpsDiffuseCoeff = 0.0,
            gsphLinearCorrect = True,
            LimiterConstructor = VanLeerLimiter,
            WaveSpeedConstructor = DavisWaveSpeed,
            nodeMotionCoefficient = 0.2,

            # artificial viscosity
            Qconstructor = LimitedMonaghanGingoldViscosity,
            Cl = 1.0, 
            Cq = 1.0,
            linearConsistent = False,
            boolReduceViscosity = False,
            nh = 5.0,
            aMin = 0.1,
            aMax = 2.0,
            Qhmult = 1.0,
            linearInExpansion = False,
            Qlimiter = False,
            balsaraCorrection = False,
            epsilon2 = 1e-2,

            # artificial conduction
            bArtificialConduction = False,
            arCondAlpha = 0.5,

            # integrator
            cfl = 0.25,
            IntegratorConstructor = CheapSynchronousRK2Integrator,
            goalTime = 2.0,
            steps = None,
            vizCycle = None,
            vizTime = 0.1,
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

            # outputs
            vizDerivs = False,
            useVoronoiOutput = False,
            clearDirectories = False,
            restoreCycle = None,
            restartStep = 100,
            redistributeStep = 500,
            checkRestart = False,
            dataDir = "dumps-KelvinHelmholtz-2d",

            serialDump = False, #whether to dump a serial ascii file at the end for viz
            
            )

assert numNodeLists in (1, 2)

assert not svph 
assert not (compatibleEnergy and evolveTotalEnergy)
assert sum([fsisph,psph,gsph,crksph,svph,mfm])<=1
assert not (fsisph and not solid)
assert not ((mfm or gsph or mfv) and ( boolReduceViscosity))

# Decide on our hydro algorithm.
hydroname = 'SPH'
useArtificialViscosity=True

if svph:
    hydroname = "SVPH"
elif crksph:
    hydroname = "CRK"+hydroname
    Qconstructor = LimitedMonaghanGingoldViscosity
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
elif mfv:
    hydroname = "MFV"
    useArtificialViscosity=False
if asph: 
    hydorname = "A"+hydroname
if solid: 
    hydroname = "solid"+hydroname

dataDir = os.path.join(dataDir,
                       "rho1=%g-rho2=%g" % (rho1, rho2),
                       "vx1=%g-vx2=%g" % (abs(vx1), abs(vx2)),
                       "vxboost=%g-vyboost=%g" % (vxboost, vyboost),
                       hydroname,
                       "densityUpdate=%s" % (densityUpdate),
                       "XSPH=%s" % xsph,
                       "filter=%s" % filter,
                       "%s-Cl=%g-Cq=%g" % (str(Qconstructor).split("'")[1].split(".")[-1], Cl, Cq),
                       "%ix%i" % (nx1, ny1 + ny2),
                       "nPerh=%g-Qhmult=%g" % (nPerh, Qhmult))
restartDir = os.path.join(dataDir, "restarts")
vizDir = os.path.join(dataDir, "visit")
restartBaseName = os.path.join(restartDir, "KelvinHelmholtz-2d")
vizBaseName = "KelvinHelmholtz-2d"

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
if KernelConstructor == NBSplineKernel:
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

nodes1 = nodeListConstructor("High density gas", eos,
                           hmin = hmin,
                           hmax = hmax,
                           hminratio = hminratio,
                           nPerh = nPerh)
nodes2 = nodeListConstructor("Low density gas", eos,
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
                                            rho = rho1,
                                            distributionType = "lattice",
                                            xmin = (0.0,  0.25),
                                            xmax = (1.0,  0.75),
                                            nNodePerh = nPerh,
                                            SPH = SPH)
    generator21 = GenerateNodeDistribution2d(nx2, int(0.5*ny2 + 0.5),
                                             rho = rho2,
                                             distributionType = "lattice",
                                             xmin = (0.0, 0.0),
                                             xmax = (1.0, 0.25),
                                             nNodePerh = nPerh,
                                             SPH = SPH)
    generator22 = GenerateNodeDistribution2d(nx2, int(0.5*ny2 + 0.5),
                                             rho = rho2,
                                             distributionType = "lattice",
                                             xmin = (0.0, 0.75),
                                             xmax = (1.0, 1.0),
                                             nNodePerh = nPerh,
                                             SPH = SPH)
    generator2 = CompositeNodeDistribution(generator21, generator22)

    if numNodeLists == 2:
        distributeNodes2d((nodes1, generator1),
                          (nodes2, generator2))
    else:
        gen = CompositeNodeDistribution(generator1, generator2)
        distributeNodes2d((nodes1, gen))

    # A helpful method for setting y velocities.
    def vy(ri):
        thpt = 1.0/(2.0*sigma*sigma)
        return (w0*sin(freq*pi*ri.x) *
                (exp(-((ri.y - 0.25)**2 * thpt)) +
                 exp(-((ri.y - 0.75)**2 * thpt))))*abs(0.5 - ri.y)

    # Finish initial conditions.
    eps1 = P1/((gamma - 1.0)*rho1)
    eps2 = P2/((gamma - 1.0)*rho2)
    if numNodeLists == 2:
        nodes1.specificThermalEnergy(ScalarField("tmp", nodes1, eps1))
        nodes2.specificThermalEnergy(ScalarField("tmp", nodes2, eps2))
        for (nodes, vx) in ((nodes1, vx1),
                            (nodes2, vx2)):
            pos = nodes.positions()
            vel = nodes.velocity()
            for i in range(nodes.numInternalNodes):
                vel[i] = Vector(vx + vxboost, vy(pos[i]) + vyboost)
    else:
        pos = nodes1.positions()
        vel = nodes1.velocity()
        eps = nodes1.specificThermalEnergy()
        for i in range(nodes1.numInternalNodes):
            if pos[i].y > 0.25 and pos[i].y < 0.75:
                eps[i] = eps1
                vel[i] = Vector(vx1 + vxboost, vy(pos[i]) + vyboost)
            else:
                eps[i] = eps2
                vel[i] = Vector(vx2 + vxboost, vy(pos[i]) + vyboost)

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
                   Q=q,
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
                   linearCorrectGradients = correctVelocityGradient,
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
                gradientType = SPHSameTimeGradient,
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
                gradientType = SPHSameTimeGradient,
                densityUpdate=densityUpdate,
                XSPH = xsph,
                ASPH = asph,
                epsTensile = epsilonTensile,
                nTensile = nTensile)
elif mfv:
    limiter = LimiterConstructor()
    waveSpeed = WaveSpeedConstructor()
    solver = HLLC(limiter,waveSpeed,gsphLinearCorrect)
    hydro = MFV(dataBase = db,
                riemannSolver = solver,
                W = WT,
                cfl=cfl,
                compatibleEnergyEvolution = compatibleEnergy,
                correctVelocityGradient= correctVelocityGradient,
                nodeMotionCoefficient = nodeMotionCoefficient,
                nodeMotionType = NodeMotionType.Lagrangian,
                gradientType = SPHSameTimeGradient,
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
# Create boundary conditions.
#-------------------------------------------------------------------------------
xp1 = Plane(Vector(0.0, 0.0), Vector( 1.0, 0.0))
xp2 = Plane(Vector(1.0, 0.0), Vector(-1.0, 0.0))
yp1 = Plane(Vector(0.0, 0.0), Vector(0.0,  1.0))
yp2 = Plane(Vector(0.0, 1.0), Vector(0.0, -1.0))
xbc = PeriodicBoundary(xp1, xp2)
ybc = PeriodicBoundary(yp1, yp2)
bcSet = [xbc, ybc]

for p in packages:
    for bc in bcSet:
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

# Blago!  Currently a problem with periodic boundaries.
# integrator.cullGhostNodes = False

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
#if not useVoronoiOutput:
#    import SpheralPointmeshSiloDump
#    vizMethod = SpheralPointmeshSiloDump.dumpPhysicsState
control = SpheralController(integrator, WT,
                            initializeDerivatives=True,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle,
                            redistributeStep = redistributeStep,
                            vizBaseName = vizBaseName,
                            vizDerivs = vizDerivs,
                            vizDir = vizDir,
                            vizStep = vizCycle,
                            vizTime = vizTime)
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

