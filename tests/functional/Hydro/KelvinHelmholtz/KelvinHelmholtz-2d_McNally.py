#ATS:test(SELF, "--crksph=True --nx1=256 --nx2=256 --ny1=128 --ny2=128 --cfl=0.25 --Cl=1.0 --Cq=1.0 --clearDirectories=True --KernelConstructor NBSplineKernel --filter=0.0 --nPerh=1.51 --graphMixing True --mixFile KH_CRK_256x256.gnu --serialDump=False", label="KH CRK 256^2, nPerh=1.5", np=10)
#ATS:test(SELF, "--crksph=True --nx1=512 --nx2=512 --ny1=256 --ny2=256 --cfl=0.25 --Cl=1.0 --Cq=1.0 --clearDirectories=True --KernelConstructor NBSplineKernel --filter=0.0 --nPerh=1.51 --graphMixing True --mixFile KH_CRK_512x512.gnu --serialDump=False", label="KH CRK 512^2, nPerh=1.5", np=70)

#-------------------------------------------------------------------------------
# This is the basic Kelvin-Helmholtz problem as discussed in
# Springel 2010, MNRAS, 401, 791-851.
#-------------------------------------------------------------------------------
import os, sys, shutil
import mpi
from math import *

from Spheral2d import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from GenerateNodeDistribution2d import *
from CompositeNodeDistribution import *
from CentroidalVoronoiRelaxation import *

import DistributeNodes
if mpi.procs > 1:
    from PeanoHilbertDistributeNodes import distributeNodes2d
else:
    from DistributeNodes import distributeNodes2d

title("Kelvin-Helmholtz test problem in 2D")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(nx1 = 256,
            ny1 = 128,
            nx2 = 256,
            ny2 = 128,
            
            rho1 = 1.0,
            rho2 = 2.0,
            P1 = 2.5,
            vx1 = 0.5,
            vx2 = -0.5,
            vxboost = 0.0,
            vyboost = 0.0,
            smooth = 0.025,
            delta = 0.01,
            freq = 4.0,
            w0 = 0.1,
            sigma = 0.05/sqrt(2.0),

            numNodeLists = 2,  # If 2, makes this a two material problem.

            gamma = 5.0/3.0,
            mu = 1.0,

            # kernel
            KernelConstructor = WendlandC4Kernel,
            nbSplineOrder = 7,
            nPerh = 4.01,
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
            solid = False,
            asph = False,   # Just for choosing the H algorithm
            XSPH = False,
            useVelocityMagnitudeForDt = False,
            epsilonTensile = 0.0,
            nTensile = 8,
            filter = 0.0,
            densityUpdate = RigorousSumDensity,
            compatibleEnergy = True,         
            gradhCorrection = True,
            correctVelocityGradient = True,
            evolveTotalEnergy = False,     

            # artificial viscosity 
            Cl = None, 
            Cq = None,
            linearConsistent = False,
            boolReduceViscosity = False,
            nh = 5.0,
            aMin = 0.1,
            aMax = 2.0,
            Qhmult = 1.0,
            boolCullenViscosity = False,
            alphMax = 2.0,
            alphMin = 0.02,
            betaC = 0.7,
            betaD = 0.05,
            betaE = 1.0,
            fKern = 1.0/3.0,
            boolHopkinsCorrection = True,
            linearInExpansion = False,
            Qlimiter = None,
            balsaraCorrection = False,
            epsilon2 = None,
            
            # CRKSPH parameters
            correctionOrder = LinearOrder,
            volumeType = RKSumVolume,

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
    
            # GSPH parameters
            gsphEpsDiffuseCoeff = 0.0,
            gsphLinearCorrect = True,
            LimiterConstructor = VanLeerLimiter,
            WaveSpeedConstructor = DavisWaveSpeed,

            ## integrator
            cfl = 0.5,
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
            HUpdate = IdealH,
            domainIndependent = False,
            rigorousBoundaries = False,
            dtverbose = False,

            # outputs
            useVoronoiOutput = True,
            clearDirectories = False,
            restoreCycle = -1,
            restartStep = 100,
            redistributeStep = None,
            checkRestart = False,
            dataDir = "dumps-KelvinHelmholtz-2d_McNally",
            graphMixing = False,
            mixInterval = 0.02,
            mixFile = "MixingModeAmp.gnu",
            
            serialDump = False, #whether to dump a serial ascii file at the end for viz
            
            bArtificialConduction = False,
            arCondAlpha = 0.5,
            )

assert not(boolReduceViscosity and boolCullenViscosity)
assert numNodeLists in (1, 2)
assert not svph 
assert not (compatibleEnergy and evolveTotalEnergy)
assert sum([fsisph,psph,gsph,crksph,svph,mfm])<=1
assert not (fsisph and not solid)
assert not ((mfm or gsph) and (boolCullenViscosity or boolReduceViscosity))

# hydro algorithm label
useArtificialViscosity = True
if svph:
    hydroname = "SVPH"
elif crksph:
    hydroname = os.path.join("CRKSPH",
                             "correctionOrder=%s" % correctionOrder,
                             "volumeType=%s" % volumeType)
elif psph:
    hydroname = "PSPH"
elif fsisph:
    hydroname = "FSISPH"
elif gsph:
    hydroname = "GSPH"
    useArtificialViscosity=False
elif mfm:
    hydroname = "MFM"
    useArtificialViscosity=False
else:
    hydroname = "SPH"

if asph:
    hydroname = "A" + hydroname

if solid:
    hydroname = "solid" + hydroname

dataDir = os.path.join(dataDir,
                       "rho1=%g-rho2=%g" % (rho1, rho2),
                       "vx1=%g-vx2=%g" % (abs(vx1), abs(vx2)),
                       "vxboost=%g-vyboost=%g" % (vxboost, vyboost),
                       hydroname,
                       "densityUpdate=%s" % (densityUpdate),
                       "compatibleEnergy=%s" % (compatibleEnergy),
                       "Cullen=%s" % (boolCullenViscosity),
                       "filter=%g" % filter,
                       "Cl=%s-Cq=%s" % (Cl, Cq),
                       "%ix%i" % (nx1, ny1 + ny2),
                       "nPerh=%g-Qhmult=%g" % (nPerh, Qhmult))
restartDir = os.path.join(dataDir, "restarts")
vizDir = os.path.join(dataDir, "visit")
restartBaseName = os.path.join(restartDir, "KelvinHelmholtz-2d_McNally")
vizBaseName = "KelvinHelmholtz-2d_McNally"

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
# Material properties.
#-------------------------------------------------------------------------------
eos = GammaLawGasMKS(gamma, mu)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
if KernelConstructor is NBSplineKernel:
    WT = TableKernel(KernelConstructor(nbSplineOrder), 1000)
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
                           kernelExtent = kernelExtent,
                           nPerh = nPerh)
nodes2 = nodeListConstructor("Low density gas", eos,
                           hmin = hmin,
                           hmax = hmax,
                           hminratio = hminratio,
                           kernelExtent = kernelExtent,
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
generator1 = GenerateNodeDistribution2d(nx1, ny1,
                                        rho = rho2,
                                        distributionType = "lattice",
                                        xmin = (0.0,  0.25),
                                        xmax = (1.0,  0.75),
                                        nNodePerh = nPerh,
                                        SPH = (not ASPH))
generator21 = GenerateNodeDistribution2d(nx2, int(0.5*ny2 + 0.5),
                                         rho = rho1,
                                         distributionType = "lattice",
                                         xmin = (0.0, 0.0),
                                         xmax = (1.0, 0.25),
                                         nNodePerh = nPerh,
                                         SPH = (not ASPH))
generator22 = GenerateNodeDistribution2d(nx2, int(0.5*ny2 + 0.5),
                                         rho = rho1,
                                         distributionType = "lattice",
                                         xmin = (0.0, 0.75),
                                         xmax = (1.0, 1.0),
                                         nNodePerh = nPerh,
                                         SPH = (not ASPH))
generator2 = CompositeNodeDistribution(generator21, generator22)

if numNodeLists == 2:
    distributeNodes2d((nodes1, generator1),
                      (nodes2, generator2))
else:
    gen = CompositeNodeDistribution(generator1, generator2)
    distributeNodes2d((nodes1, gen))

# Finish initial conditions.
rhom = 0.5*(rho1-rho2)
vxm = 0.5*(vx1-vx2)
if numNodeLists == 2:
    for (nodes, vx) in ((nodes1, vx1),
                        (nodes2, vx2)):
        pos = nodes.positions()
        vel = nodes.velocity()
        rho = nodes.massDensity()
        eps = nodes.specificThermalEnergy()
        mass = nodes.mass()
        for i in range(nodes.numInternalNodes):
            yval = pos[i].y
            xval = pos[i].x
            velx = 0.0
            rho[i] = 0.0
            vely = delta*sin(4*pi*xval)
            if yval >= 0 and yval < 0.25:
               rho[i]=rho1 - rhom*exp((yval-0.25)/smooth)
               mass[i] *= rho[i]/rho1
               velx = vx1 - vxm*exp((yval-0.25)/smooth)
            elif yval >= 0.25 and yval < 0.5:
               rho[i]=rho2 + rhom*exp((0.25-yval)/smooth)
               mass[i] *= rho[i]/rho2
               velx = vx2 + vxm*exp((0.25-yval)/smooth)
            elif yval >= 0.5 and yval < 0.75:
               rho[i]=rho2 + rhom*exp((yval-0.75)/smooth)
               mass[i] *= rho[i]/rho2
               velx = vx2 + vxm*exp((yval-0.75)/smooth)
            else:
               rho[i]=rho1 - rhom*exp((0.75-yval)/smooth)
               mass[i] *= rho[i]/rho1
               velx = vx1 - vxm*exp((0.75-yval)/smooth)
            vel[i] = Vector(velx + vxboost, vely + vyboost)
            eps[i] = P1/((gamma - 1.0)*rho[i])
else:
    pos = nodes1.positions()
    vel = nodes1.velocity()
    rho = nodes1.massDensity()
    eps = nodes1.specificThermalEnergy()
    mass = nodes1.mass()
    for i in range(nodes1.numInternalNodes):
            yval = pos[i].y
            xval = pos[i].x
            velx = 0.0
            vely = delta*sin(4*pi*xval)
            if yval >= 0 and yval < 0.25:
               rho[i]=rho1 - rhom*exp((yval-0.25)/smooth)
               mass[i] *= rho[i]/rho1
               velx = vx1 - vxm*exp((yval-0.25)/smooth)
            elif yval >= 0.25 and yval < 0.5:
               rho[i]=rho2 + rhom*exp((0.25-yval)/smooth)
               mass[i] *= rho[i]/rho2
               velx = vx2 + vxm*exp((0.25-yval)/smooth)
            elif yval >= 0.5 and yval < 0.75:
               rho[i]=rho2 + rhom*exp((yval-0.75)/smooth)
               mass[i] *= rho[i]/rho2
               velx = vx2 + vxm*exp((yval-0.75)/smooth)
            else:
               rho[i]=rho1 - rhom*exp((0.75-yval)/smooth)
               mass[i] *= rho[i]/rho1
               velx = vx1 - vxm*exp((0.75-yval)/smooth)

            vel[i] = Vector(velx + vxboost, vely + vyboost)
            eps[i] = P1/((gamma - 1.0)*rho[i])

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
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
if svph:
    hydro = SVPH(dataBase = db,
                 W = WT, 
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
                 xmax = Vector(3.0, 3.0),
                 ASPH = asph)
    # xmin = Vector(x0 - 0.5*(x2 - x0), y0 - 0.5*(y2 - y0)),
    # xmax = Vector(x2 + 0.5*(x2 - x0), y2 + 0.5*(y2 - y0)))
elif crksph:
    hydro = CRKSPH(dataBase = db,
                   order = correctionOrder,
                   filter = filter,
                   cfl = cfl,
                   useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                   compatibleEnergyEvolution = compatibleEnergy,
                   evolveTotalEnergy = evolveTotalEnergy,
                   XSPH = XSPH,
                   densityUpdate = densityUpdate,
                   HUpdate = HUpdate,
                   ASPH = asph)
elif psph:
    hydro = PSPH(dataBase = db,
                 W = WT,
                 filter = filter,
                 cfl = cfl,
                 compatibleEnergyEvolution = compatibleEnergy,
                 evolveTotalEnergy = evolveTotalEnergy,
                 correctVelocityGradient = correctVelocityGradient,
                 densityUpdate = densityUpdate,
                 HUpdate = HUpdate,
                 XSPH = XSPH,
                 ASPH = asph)
if fsisph:
    sumDensityNodeListSwitch =[nodes1,nodes2]  
    hydro = FSISPH(dataBase = db,
                   W = WT,
                   cfl = cfl,
                   surfaceForceCoefficient = fsiSurfaceCoefficient,              
                   densityStabilizationCoefficient = fsiRhoStabilizeCoeff,         
                   specificThermalEnergyDiffusionCoefficient = fsiEpsDiffuseCoeff,     
                   xsphCoefficient = fsiXSPHCoeff,
                   interfaceMethod = fsiInterfaceMethod,
                   kernelAveragingMethod = fsiKernelMethod,
                   sumDensityNodeLists = sumDensityNodeListSwitch,
                   linearCorectGradients = correctVelocityGradient,
                   compatibleEnergyEvolution = compatibleEnergy,  
                   evolveTotalEnergy = evolveTotalEnergy,         
                   ASPH = asph,
                   epsTensile = epsilonTensile)
elif gsph:
    limiter = LimiterConstructor
    waveSpeed = WaveSpeedConstructor
    solver = HLLC(limiter,waveSpeed,gsphLinearCorrect)
    hydro = GSPH(dataBase = db,
                riemannSolver = solver,
                W = WT,
                cfl=cfl,
                compatibleEnergyEvolution = compatibleEnergy,
                correctVelocityGradient= correctVelocityGradient,
                evolveTotalEnergy = evolveTotalEnergy,
                densityUpdate=densityUpdate,
                XSPH = XSPH,
                ASPH = asph,
                epsTensile = epsilonTensile,
                nTensile = nTensile)
elif mfm:
    limiter = LimiterConstructor
    waveSpeed = WaveSpeedConstructor
    solver = HLLC(limiter,waveSpeed,gsphLinearCorrect)
    hydro = GSPH(dataBase = db,
                riemannSolver = solver,
                W = WT,
                cfl=cfl,
                compatibleEnergyEvolution = compatibleEnergy,
                correctVelocityGradient= correctVelocityGradient,
                evolveTotalEnergy = evolveTotalEnergy,
                densityUpdate=densityUpdate,
                XSPH = XSPH,
                ASPH = asph,
                epsTensile = epsilonTensile,
                nTensile = nTensile)
else:
    hydro = SPH(dataBase = db,
                W = WT,
                cfl = cfl,
                useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                compatibleEnergyEvolution = compatibleEnergy,
                evolveTotalEnergy = evolveTotalEnergy,
                gradhCorrection = gradhCorrection,
                correctVelocityGradient = correctVelocityGradient,
                XSPH = XSPH,
                densityUpdate = densityUpdate,
                HUpdate = HUpdate,
                epsTensile = epsilonTensile,
                nTensile = nTensile,
                ASPH = asph)
output("hydro")
output("hydro.cfl")
output("hydro.compatibleEnergyEvolution")
output("hydro.densityUpdate")
output("hydro.HEvolution")

packages = [hydro]

#-------------------------------------------------------------------------------
# Set the artificial viscosity parameters.
#-------------------------------------------------------------------------------
if useArtificialViscosity:
    q = hydro.Q
    if not Cl is None:
        q.Cl = Cl
    if not Cq is None:
        q.Cq = Cq
    if not epsilon2 is None:
        q.epsilon2 = epsilon2
    if not Qlimiter is None:
        q.limiter = Qlimiter
    if not balsaraCorrection is None:
        q.balsaraShearCorrection = balsaraCorrection
    output("q")
    output("q.Cl")
    output("q.Cq")  
    output("q.epsilon2")
    output("q.limiter")
    output("q.balsaraShearCorrection")
    try:
        output("q.linearInExpansion")
        output("q.quadraticInExpansion")
    except:
        pass

#-------------------------------------------------------------------------------
# Construct the MMRV physics object.
#-------------------------------------------------------------------------------
if boolReduceViscosity and useArtificialViscosity:
    evolveReducingViscosityMultiplier = MorrisMonaghanReducingViscosity(nh,aMin,aMax)
    packages.append(evolveReducingViscosityMultiplier)
elif boolCullenViscosity and useArtificialViscosity:
    evolveCullenViscosityMultiplier = CullenDehnenViscosity(WT,alphMax,alphMin,betaC,betaD,betaE,fKern,boolHopkinsCorrection)
    packages.append(evolveCullenViscosityMultiplier)

#-------------------------------------------------------------------------------
# Construct the Artificial Conduction physics object.
#-------------------------------------------------------------------------------
if bArtificialConduction:
    #q.reducingViscosityCorrection = True
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
#ybc1 = ReflectingBoundary(yp1)
#ybc2 = ReflectingBoundary(yp2)
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
control = SpheralController(integrator, WT,
                            volumeType = volumeType,
                            initializeDerivatives = True,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle,
                            redistributeStep = redistributeStep,
                            vizBaseName = vizBaseName,
                            vizDir = vizDir,
                            vizStep = vizCycle,
                            vizTime = vizTime,
                            SPH = (not asph))
output("control")

#-------------------------------------------------------------------------------
# Add a method for measuring the mixing scale.
#-------------------------------------------------------------------------------
def mixingScale(cycle, t, dt):
    si = []
    ci = []
    di = []
    ke = []
    for nodeL in nodeSet:
     xprof = mpi.reduce([x.x for x in nodeL.positions().internalValues()], mpi.SUM)
     yprof = mpi.reduce([x.y for x in nodeL.positions().internalValues()], mpi.SUM)
     vely = mpi.reduce([v.y for v in nodeL.velocity().internalValues()], mpi.SUM)
     hprof = mpi.reduce([1.0/sqrt(H.Determinant()) for H in nodeL.Hfield().internalValues()], mpi.SUM)
     rhoprof = mpi.reduce(nodeL.massDensity().internalValues(), mpi.SUM)
     if mpi.rank == 0:
      for j in range (len(xprof)):
        ke.append(0.5*rhoprof[j]*vely[j]*vely[j])
        if yprof[j] < 0.5:
          si.append(vely[j]*hprof[j]*hprof[j]*sin(4*pi*xprof[j])*exp(-4.0*pi*abs(yprof[j]-0.25)))
          ci.append(vely[j]*hprof[j]*hprof[j]*cos(4*pi*xprof[j])*exp(-4.0*pi*abs(yprof[j]-0.25)))
          di.append(hprof[j]*hprof[j]*exp(-4.0*pi*abs(yprof[j]-0.25)))
        else:
          si.append(vely[j]*hprof[j]*hprof[j]*sin(4*pi*xprof[j])*exp(-4.0*pi*abs((1.0-yprof[j])-0.25)))
          ci.append(vely[j]*hprof[j]*hprof[j]*cos(4*pi*xprof[j])*exp(-4.0*pi*abs((1.0-yprof[j])-0.25)))
          di.append(hprof[j]*hprof[j]*exp(-4.0*pi*abs((1.0-yprof[j])-0.25)))
    if mpi.rank == 0:
      S=sum(si)
      C=sum(ci)
      D=sum(di)
      M=sqrt((S/D)*(S/D)+(C/D)*(C/D))*2.0
      KE = max(ke)
      print("At time t = %s, Mixing Amp M = %s \n" % (t,M))
      with open(os.path.join(dataDir, mixFile), "a") as myfile:
        myfile.write("%s\t %s\t %s\n" % (t, M, KE))

if graphMixing:
    control.appendPeriodicTimeWork(mixingScale, mixInterval)
    myfile = open(os.path.join(dataDir, mixFile), "w")
    myfile.write("# time           mixamp                     KEMax\n")
    myfile.close()

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
if not steps is None:
    control.step(steps)

else:
    control.advance(goalTime, maxSteps)
    control.updateViz(control.totalSteps, integrator.currentTime, 0.0)
    control.dropRestartFile()

#    nsteps=int(goalTime/vizTime)
#    print "NSTEPS=",nsteps," dt=",vizTime
#    t=0.0
#    times=[]
#    amps=[]
#    for i in range(nsteps):
#      t=t+vizTime
#      print "TIME=",t
#      control.advance(t, t)
#      control.updateViz(control.totalSteps, integrator.currentTime, 0.0)
#      control.dropRestartFile()
# 
#      sums=0.0
#      sumc=0.0
#      sumd=0.0
#      for (nodes, vx) in ((nodes1, vx1),
#                            (nodes2, vx2)):
#         
#          pos = nodes.positions()
#          vel = nodes.velocity()
#          rho = nodes.massDensity()
#          si=numpy.zeros(nodes.numInternalNodes)
#          ci=numpy.zeros(nodes.numInternalNodes)
#          di=numpy.zeros(nodes.numInternalNodes)
#          for i in xrange(nodes.numInternalNodes):
#            hi2=pow(2.0/(nodes.Hfield()[i].Trace()),2.0)
#            posx=pos[i].x
#            posy=pos[i].y
#            if posy < 0.5:
#              si[i]=vel[i].y*hi2*sin(4*pi*posx)*exp(-4*pi*abs(posy-0.25))
#              ci[i]=vel[i].y*hi2*cos(4*pi*posx)*exp(-4*pi*abs(posy-0.25))
#              di[i]=vel[i].y*hi2*exp(-4*pi*abs(posy-0.25))
#            else:
#              si[i]=vel[i].y*hi2*sin(4*pi*posx)*exp(-4*pi*abs((1-posy)-0.25))
#              ci[i]=vel[i].y*hi2*cos(4*pi*posx)*exp(-4*pi*abs((1-posy)-0.25))
#              di[i]=vel[i].y*hi2*exp(-4*pi*abs((1-posy)-0.25))
#          sums+=numpy.sum(si)
#          sumc+=numpy.sum(ci)
#          sumd+=numpy.sum(di)
#          #sums = mpi.allreduce(sums,op=mpi.SUM)
#          #sumc = mpi.allreduce(sums,op=mpi.SUM)
#          #sumd = mpi.allreduce(sums,op=mpi.SUM)
#      amps.append(2.0*sqrt(pow(sums/sumd,2.0)+pow(sumc/sumd,2.0)))
#      times.append(t)
#      print "AMPS = ",2.0*sqrt(pow(sums/sumd,2.0)+pow(sumc/sumd,2.0)), " t=", t
#           
#rank = mpi.rank
#if rank == 0:
#  numpy.savetxt("mixing.txt",(times,amps))

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
    f = open(os.path.join(dataDir, "./serialDump.ascii"),'w')
    for i in range(len(serialData)):
      f.write("{0} {1} {2} {3} {4} {5} {6} {7}\n".format(i,serialData[i][0][0],serialData[i][0][1],0.0,serialData[i][1],serialData[i][2],serialData[i][3],serialData[i][4]))
    f.close()
