import shutil
from math import *
from Spheral3d import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from findLastRestart import *
from GenerateNodeDistribution3d import *
from CompositeNodeDistribution import *
from CentroidalVoronoiRelaxation import *

import mpi
import DistributeNodes

title("Kelvin-Helmholtz test problem in 3D")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(nx = 200,
            ny = 100,
            nz = 200,
            
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

            gamma = 5.0/3.0,
            mu = 1.0,

            nPerh = 1.51,

            SVPH = False,
            CRKSPH = False,
            PSPH = False,
            ASPH = False,  # This just chooses the H algorithm -- you can use this with CRKSPH for instance.
            filter = 0.0,   # CRKSPH filtering
            Qconstructor = MonaghanGingoldViscosity,
            KernelConstructor = BSplineKernel,
            order = 5,
            #Qconstructor = TensorMonaghanGingoldViscosity,
            linearConsistent = False,
            fcentroidal = 0.0,
            fcellPressure = 0.0,
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
            Cl = 1.0, 
            Cq = 1.0,
            linearInExpansion = False,
            Qlimiter = False,
            balsaraCorrection = False,
            epsilon2 = 1e-2,
            hmin = 0.0001, 
            hmax = 0.5,
            hminratio = 0.1,
            cfl = 0.5,
            useVelocityMagnitudeForDt = False,
            XSPH = False,
            epsilonTensile = 0.0,
            nTensile = 8,

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

            densityUpdate = RigorousSumDensity, # VolumeScaledDensity,
            compatibleEnergy = True,            # <--- Important!  rigorousBoundaries does not work with the compatibleEnergy algorithm currently.
            gradhCorrection = False,

            useVoronoiOutput = False,
            clearDirectories = False,
            restoreCycle = None,
            restartStep = 100,
            redistributeStep = 500,
            checkRestart = False,
            dataDir = "dumps-KelvinHelmholtz-3d",

            bArtificialConduction = False,
            arCondAlpha = 0.5,
            )
assert not(boolReduceViscosity and boolCullenViscosity)


# Decide on our hydro algorithm.
if SVPH:
    if ASPH:
        HydroConstructor = ASVPHFacetedHydro
    else:
        HydroConstructor = SVPHFacetedHydro
elif CRKSPH:
    Qconstructor = LimitedMonaghanGingoldViscosity
    if ASPH:
        HydroConstructor = ACRKSPHHydro
    else:
        HydroConstructor = CRKSPHHydro
elif PSPH:
    if ASPH:
        HydroConstructor = APSPHHydro
    else:
        HydroConstructor = PSPHHydro
else:
    if ASPH:
        HydroConstructor = ASPHHydro
    else:
        HydroConstructor = SPHHydro

dataDir = os.path.join(dataDir,
                       "rho1=%g-rho2=%g" % (rho1, rho2),
                       "vx1=%g-vx2=%g" % (abs(vx1), abs(vx2)),
                       "vxboost=%g-vyboost=%g" % (vxboost, vyboost),
                       str(HydroConstructor).split("'")[1].split(".")[-1],
                       "densityUpdate=%s" % (densityUpdate),
                       "compatibleEnergy=%s" % (compatibleEnergy),
                       "XSPH=%s" % XSPH,
                       "filter=%s" % filter,
                       "%s-Cl=%g-Cq=%g" % (str(Qconstructor).split("'")[1].split(".")[-1], Cl, Cq),
                       "%ix%ix%i" % (nx, ny*2,nz),
                       "nPerh=%g-Qhmult=%g" % (nPerh, Qhmult))
restartDir = os.path.join(dataDir, "restarts")
vizDir = os.path.join(dataDir, "visit")
restartBaseName = os.path.join(restartDir, "KelvinHelmholtz-3d")
vizBaseName = "KelvinHelmholtz-3d"

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
if KernelConstructor==NBSplineKernel:
  WT = TableKernel(NBSplineKernel(order), 1000)
  WTPi = TableKernel(NBSplineKernel(order), 1000, Qhmult)
else:
  WT = TableKernel(KernelConstructor(), 1000)
  WTPi = TableKernel(KernelConstructor(), 1000, Qhmult)
Wfbase = NBSplineKernel(9)
WTf = TableKernel(Wfbase, 1000, hmult=1.0/(nPerh*Wfbase.kernelExtent))
output("WT")
output("WTPi")
kernelExtent = WT.kernelExtent

#-------------------------------------------------------------------------------
# Make the NodeList.
#-------------------------------------------------------------------------------
nodesLow = makeFluidNodeList("Low density gas", eos,
                           hmin = hmin,
                           hmax = hmax,
                           hminratio = hminratio,
                           kernelExtent = kernelExtent,
                           nPerh = nPerh)
nodesHigh = makeFluidNodeList("High density gas", eos,
                           hmin = hmin,
                           hmax = hmax,
                           kernelExtent = kernelExtent,
                           hminratio = hminratio,
                           nPerh = nPerh)
nodeSet = [nodesLow, nodesHigh]
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
    generatorHigh = GenerateNodeDistribution3d(nx, ny, nz,
                                            rho = rho2,
                                            distributionType = "lattice",
                                            xmin = (0.0, 0.0, 0.0),
                                            xmax = (1.0, 0.5, 1.0),
                                            nNodePerh = nPerh,
                                            SPH = (not ASPH))
    generatorLow = GenerateNodeDistribution3d(nx, ny, nz,
                                             rho = rho1,
                                             distributionType = "lattice",
                                             xmin = (0.0, 0.5, 0.0),
                                             xmax = (1.0, 1.0, 1.0),
                                             nNodePerh = nPerh,
                                             SPH = (not ASPH))
    if mpi.procs > 1:
        from VoronoiDistributeNodes import distributeNodes3d
    else:
        from DistributeNodes import distributeNodes3d

    distributeNodes3d((nodesLow, generatorLow),
                      (nodesHigh, generatorHigh))


    # Finish initial conditions.
    rhom = 0.5*(rho1-rho2)
    vxm = 0.5*(vx1-vx2)
    for (nodes, vx) in ((nodesLow, vx1),
                        (nodesHigh, vx2)):
        pos = nodes.positions()
        vel = nodes.velocity()
        rho = nodes.massDensity()
        eps = nodes.specificThermalEnergy()
        mass = nodes.mass()
        for i in range(nodes.numInternalNodes):
            yval = pos[i].y
            xval = pos[i].x
            zval = pos[i].z
            velx = 0.0
            rho[i] = 0.0
            vely = delta*sin(4*pi*xval)*sin(4*pi*zval)
            if yval >= 0 and yval < 0.5:
               rho[i]=rho1 + rhom*exp((yval-0.5)/smooth)
               mass[i] *= rho[i]/rho1
               velx = vx1 + vxm*exp((yval-0.5)/smooth)
            elif yval >= 0.5 and yval < 1:
               rho[i]=rho2 - rhom*exp((0.5-yval)/smooth)
               mass[i] *= rho[i]/rho2
               velx = vx2 - vxm*exp((0.5-yval)/smooth)
            vel[i] = Vector(velx + vxboost, vely + vyboost,0)
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
    hydro = HydroConstructor(WT, q,
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
                             xmin = Vector(-2.0, -2.0, -2.0),
                             xmax = Vector(3.0, 3.0,3.0))
                             # xmin = Vector(x0 - 0.5*(x2 - x0), y0 - 0.5*(y2 - y0)),
                             # xmax = Vector(x2 + 0.5*(x2 - x0), y2 + 0.5*(y2 - y0)))
elif CRKSPH:
    Wf = NBSplineKernel(9)
    hydro = HydroConstructor(WT, WTPi, q,
                             Wfilter = WTf,
                             filter = filter,
                             cfl = cfl,
                             useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                             compatibleEnergyEvolution = compatibleEnergy,
                             XSPH = XSPH,
                             densityUpdate = densityUpdate,
                             HUpdate = HUpdate)
else:
    hydro = HydroConstructor(WT,
                             WTPi,
                             q,
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
elif boolCullenViscosity:
    evolveCullenViscosityMultiplier = CullenDehnenViscosity(WTPi,alphMax,alphMin,betaC,betaD,betaE,fKern,boolHopkinsCorrection)
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
xPlane0 = Plane(Vector(0.0, 0.0, 0.0), Vector( 1.0,  0.0, 0.0))
xPlane1 = Plane(Vector(1.0, 0.0, 0.0), Vector(-1.0,  0.0, 0.0))
yPlane0 = Plane(Vector(0.0, 0.0, 0.0), Vector( 0.0,  1.0, 0.0))
yPlane1 = Plane(Vector(0.0, 1.0, 0.0), Vector( 0.0, -1.0, 0.0))
zPlane0 = Plane(Vector(0.0, 0.0, 0.0), Vector( 0.0,  0.0, 1.0))
zPlane1 = Plane(Vector(0.0, 0.0, 1.0), Vector( 0.0,  0.0, -1.0))

xbc = PeriodicBoundary(xPlane0, xPlane1)
#ybc = PeriodicBoundary(yPlane0, yPlane1)
zbc = PeriodicBoundary(zPlane0, zPlane1)

ybc1 = ReflectingBoundary(yPlane0)
ybc2 = ReflectingBoundary(yPlane1)

bcSet = [xbc, zbc, ybc1, ybc2]

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
                            SPH = (not ASPH))
output("control")

#-------------------------------------------------------------------------------
# Add a method for measuring the mixing scale.
#-------------------------------------------------------------------------------
def mixingScale(cycle, t, dt):
    si = []
    ci = []
    di = []
    for nodeL in nodeSet:
     xprof = mpi.reduce([x.x for x in nodeL.positions().internalValues()], mpi.SUM)
     yprof = mpi.reduce([x.y for x in nodeL.positions().internalValues()], mpi.SUM)
     vely = mpi.reduce([v.y for v in nodeL.velocity().internalValues()], mpi.SUM)
     hprof = mpi.reduce([1.0/sqrt(H.Determinant()) for H in nodeL.Hfield().internalValues()], mpi.SUM)
     if mpi.rank == 0:
      for j in range (len(xprof)):
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
      print("At time t = %s, Mixing Amp M = %s \n" % (t,M))
      with open(mixFile, "a") as myfile:
        myfile.write("%s\t %s\n" % (t, M))

if graphMixing:
    control.appendPeriodicTimeWork(mixingScale, mixInterval)
    myfile = open(mixFile, "w")
    myfile.write("# time           mixamp\n")
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
