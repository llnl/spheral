#-------------------------------------------------------------------------------
# The Blob Test
# Based on the test presented in 
# 1.	Agertz O, Moore B, Stadel J, Potter D, Miniati F, Read J, et al. 
#       Fundamental differences between SPH and grid methods. Monthly Notices of
#       the Royal Astronomical Society. 2007; 380(3):963-978.
#       doi:10.1111/j.1365-2966.2007.12183.x.
#-------------------------------------------------------------------------------
import shutil
from math import *
from Spheral2d import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from findLastRestart import *
from GenerateNodeDistribution2d import *
from CloudMassFraction import *

import mpi
import DistributeNodes

title("2-D integrated hydro test --  Blob Test")

#-------------------------------------------------------------------------------
# Rejecter to help establishing initial conditions.
# This is used by GenerateNodeDistrubtion2d to carve out a spherical region.
#-------------------------------------------------------------------------------
class SphericalRejecterBlob:
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
                (y[i] - self.origin[1])**2 < R2):
                xnew.append(x[i])
                ynew.append(y[i])
                mnew.append(m[i])
                Hnew.append(H[i])
        return xnew, ynew, mnew, Hnew

#-------------------------------------------------------------------------------
# Rejecter to help establishing initial conditions.
# This is used by GenerateNodeDistrubtion2d to cut out a spherical region cavity
#-------------------------------------------------------------------------------
class SphericalRejecter:
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
                (y[i] - self.origin[1])**2 >= R2):
                xnew.append(x[i])
                ynew.append(y[i])
                mnew.append(m[i])
                Hnew.append(H[i])
        return xnew, ynew, mnew, Hnew

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(

    #Set external state (Internal state is dependent so dont set it)
    gamma = 5.0/3.0, # specific heat ratio
    rhoext = 1.0,    # density of external medium
    Pequi = 1.0,     # initial pressure
    mach = 2.7,      # mach number 
    chi = 10.0,      # Ratio of rhoblob/rhoext

    # Geometry Ambient medium box
    xb0 = 0.0,
    xb1 = 40.0,
    yb0 = 0.0,
    yb1 = 10.0,

    #Blob radius and central location
    br = 1.0,
    bx = 5.0,
    by = 5.0, 

    # Resolution and node seeding.
    nx1 = 256,            # num nodes spanning x
    ny1 = 64,             # num nodes spanning y
    massMatch = True,     # If False, match spatial resolution in blob  

    # kernel
    HUpdate = IdealH,
    nPerh = 1.35,
    KernelConstructor = NBSplineKernel,
    order = 5,
    hmin = 1e-5,
    hmax = 0.5,
    hminratio = 0.1,

    # hydro types
    svph = False,
    crksph = False,
    psph = False,
    fsisph = False,
    gsph = False, 

    # hydro options
    solid = False,                      # use solid node lists (fluid limit of solid hydro)
    asph = False,                       # This just chooses the H algorithm -- you can use this with CRKSPH for instance.
    xsph = False,                       # integrate position w/ smoothed velocity
    filter = 0.0,                       # For CRKSPH
    epsilonTensile = 0.0,               # tensile correction
    nTensile = 8,                       # tensile correction
    densityUpdate = RigorousSumDensity, # density method (RigorousSumDensity,IntegrateDensity)
    compatibleEnergy = True,            # activated 2nd loop to force energy conservation
    evolveTotalEnergy = False,          # formulation based on total energy instead of specific thermal energy
    gradhCorrection = True,             # correct for the temporal variation in h
    correctVelocityGradient = True,     # linear corrected velocity gradient

    # FSISPH parameters
    fsiSurfaceCoefficient = 0.00,           # adds additional repulsive force to material interfaces)
    fsiRhoStabilizeCoeff = 0.0,             # coefficient that smooths the density field
    fsiEpsDiffuseCoeff = 0.0,               # explicit diiffusion of the thermal energy
    fsiXSPHCoeff = 0.00,                    # fsi uses multiplier for XSPH instead of binary switch
    fsiInterfaceMethod = ModulusInterface,  # (HLLCInterface, ModulusInterface)
    fsiKernelMethod  = NeverAverageKernels, # (NeverAverageKernels, AlwaysAverageKernels, AverageInterfaceKernels)
    
    # GSPH parameters
    gsphEpsDiffuseCoeff = 0.0,
    gsphLinearCorrect = True,

    # artificial viscosity
    Cl = 1.0, 
    Cq = 1.0,
    Qconstructor = MonaghanGingoldViscosity,
    boolReduceViscosity = False,
    nhQ = 5.0,
    nhL = 10.0,
    boolCullenViscosity = False,
    alphMax = 2.0,
    alphMin = 0.02,
    betaC = 0.7,
    betaD = 0.05,
    betaE = 1.0,
    fKern = 1.0/3.0,
    boolHopkinsCorrection = True,
    aMin = 0.1,
    aMax = 2.0,
    linearConsistent = False,
    fcentroidal = 0.0,
    fcellPressure = 0.0,
    Qlimiter = False,
    balsaraCorrection = False,
    epsilon2 = 1e-2,

    # integrator
    cfl = 0.25,
    goalTKH = 2.5,  # Goal time in units of t_KH
    IntegratorConstructor = CheapSynchronousRK2Integrator,
    steps = None,
    vizCycle = 20,
    vizTime = 0.1,
    dt = 0.0001,
    dtMin = 1.0e-5, 
    dtMax = 0.1,
    dtGrowth = 2.0,
    maxSteps = None,
    statsStep = 10,
    domainIndependent = False,
    rigorousBoundaries = False,
    dtverbose = False,

    # outputs
    clearDirectories = False,
    restoreCycle = -1,
    restartStep = 200,
    dataDir = "dumps-blobtest-2d",
    histfilename = "cloud_mass_history.gnu",
    rhoThresholdFrac = 0.64,
    epsThresholdFrac = 0.9,
    massFracFreq = 10,
    )

# Check the input.
assert not (boolReduceViscosity and boolCullenViscosity)
assert not (compatibleEnergy and evolveTotalEnergy)
assert not svph 
assert sum([fsisph,psph,gsph,crksph,svph])<=1
assert not (fsisph and not solid)
# Decide on our hydro algorithm.
hydroname = 'SPH'
if svph:
    hydroname = "SVPH"
elif crksph:
    hydroname = "CRK"+hydroname
elif psph:
    hydroname = "P"+hydroname
elif fsisph:
    hydroname = "FSI"+hydroname
elif gsph:
    hydroname = "G"+hydroname
if asph: 
    hydorname = "A"+hydroname
if solid: 
    hydroname = "solid"+hydroname


# Build our directory paths.
densityUpdateLabel = {IntegrateDensity : "IntegrateDensity",
                      SumDensity : "SumDensity",
                      RigorousSumDensity : "RigorousSumDensity",
                      SumVoronoiCellDensity : "SumVoronoiCellDensity"}
baseDir = os.path.join(dataDir,
                       hydroname,
                       Qconstructor.__name__,
                       KernelConstructor.__name__,
                       densityUpdateLabel[densityUpdate],
                       "compatibleEnergy=%s" % compatibleEnergy,
                       "evolveTotalEnergy=%s" % evolveTotalEnergy,
                       "Cullen=%s" % boolCullenViscosity,
                       "XSPH=%s" % xsph,
                       "nPerh=%3.1f" % nPerh,
                       "fcentroidal=%1.3f" % fcentroidal,
                       "fcellPressure = %1.3f" % fcellPressure,
                       "massMatch=%s" % massMatch,
                       "%ix%i" % (nx1, ny1))
restartDir = os.path.join(baseDir, "restarts")
restartBaseName = os.path.join(restartDir, "blob-2d-%ix%i" % (nx1, ny1))

vizDir = os.path.join(baseDir, "visit")
if vizTime is None and vizCycle is None:
    vizBaseName = None
else:
    vizBaseName = "blobtest-2d-%ix%i" % (nx1, ny1)

# Figure out the goal time.
csext = sqrt(gamma*Pequi/rhoext)
vext = mach*csext
tCrush = 2.0*br*sqrt(chi)/vext
tKH = 1.6*tCrush
goalTime = goalTKH * tKH

print("Computed times (tCrush, tKH, goalTime) = (%g, %g, %g)" % (tCrush, tKH, goalTime))

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
import os, sys
if mpi.rank == 0:
    if clearDirectories and os.path.exists(baseDir):
        shutil.rmtree(baseDir)
    if not os.path.exists(restartDir):
        os.makedirs(restartDir)
    if not os.path.exists(vizDir):
        os.makedirs(vizDir)
mpi.barrier()

#-------------------------------------------------------------------------------
# Material properties.
#-------------------------------------------------------------------------------
mu = 1.0
eos1 = GammaLawGasMKS(gamma, mu)
eos2 = GammaLawGasMKS(gamma, mu)

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
# Make the NodeLists.
#-------------------------------------------------------------------------------
#Ambient Medium nodes
if solid:
    nodeListConstructor = makeSolidNodeList
else:
    nodeListConstructor = makeFluidNodeList
outerNodes = nodeListConstructor("outer", eos1,
                                 hmin = hmin,
                                 hmax = hmax,
                                 hminratio = hminratio,
                                 nPerh = nPerh)
#Blob nodes
innerNodes = nodeListConstructor("inner", eos2,
                                 hmin = hmin,
                                 hmax = hmax,
                                 hminratio = hminratio,
                                 nPerh = nPerh)
nodeSet = (outerNodes, innerNodes)
for nodes in nodeSet:
    output("nodes.name")
    output("    nodes.hmin")
    output("    nodes.hmax")
    output("    nodes.hminratio")
    output("    nodes.nodesPerSmoothingScale")

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
#Blob density is defined to be 10 times the external density
rhoblob=rhoext*chi
generatorOuter = GenerateNodeDistribution2d(nx1, ny1, rhoext,
                                            distributionType = "lattice",
                                            xmin = (xb0, yb0),
                                            xmax = (xb1, yb1),
                                            rejecter = SphericalRejecter(origin = (bx, by),
                                                                         radius = br),
                                            nNodePerh = nPerh,
                                            SPH = (not ASPH))
#generatorOuter = GenerateNodeDistribution2d(nx1, ny1, nz1, rhoext,
#                                            distributionType = "lattice",
#                                            xmin = (bx-br, by-br, bz-br),
#                                            xmax = (bx+br, by+br, bz+br),
#       				        rmin = br,
#                                            nNodePerh = nPerh,
#                                            SPH = (not ASPH))
#generatorInner = GenerateNodeDistribution2d(nx1, ny1, nz1, rhoblob,
#                                            distributionType = "lattice",
#                                            xmin = (xb0, yb0, zb0),
#                                            xmax = (xb1, yb1, zb1),
#                                            rejecter = SphericalRejecterBlob(origin = (bx, by, bz),
#                                                                         radius = br),
#                                            nNodePerh = nPerh,
#                                            SPH = (not ASPH))

if massMatch:
    # Figure out a mass matched resolution for the blob.
    mouter = (xb1 - xb0)*(yb1 - yb0)*rhoext/(nx1*ny1)
    nxinner = max(2, int(((2*br)**2*rhoblob/mouter)**(1.0/2.0) + 0.5))
    generatorInner = GenerateNodeDistribution2d(nxinner, nxinner, rhoblob,
                                                distributionType = "lattice",
                                                xmin = (bx-br, by-br),
                                                xmax = (bx+br, by+br),
                                                originreject = (bx, by),
                                                rreject = br,
                                                nNodePerh = nPerh,
                                                SPH = (not ASPH))
else:
    generatorInner = GenerateNodeDistribution2d(nx1, ny1, rhoblob,
                                                distributionType = "lattice",
                                                xmin = (xb0, yb0),
                                                xmax = (xb1, yb1),
                                                originreject = (bx, by),
                                                rreject = br,
                                                nNodePerh = nPerh,
                                                SPH = (not ASPH))

if mpi.procs > 1:
    from VoronoiDistributeNodes import distributeNodes2d
else:
    from DistributeNodes import distributeNodes2d

distributeNodes2d((outerNodes, generatorOuter),
                  (innerNodes, generatorInner))
for nodes in nodeSet:
    print(nodes.name, ":")
    output("    mpi.reduce(nodes.numInternalNodes, mpi.MIN)")
    output("    mpi.reduce(nodes.numInternalNodes, mpi.MAX)")
    output("    mpi.reduce(nodes.numInternalNodes, mpi.SUM)")
del nodes

# Set node specific thermal energies
for (nodes, rho) in ((outerNodes, rhoext),
                     (innerNodes, rhoblob)):
    eps0 = Pequi/((gamma - 1.0)*rho)
    nodes.specificThermalEnergy(ScalarField("tmp", nodes, eps0))
del nodes

vel = outerNodes.velocity() #wind velocity
for i in range(outerNodes.numInternalNodes):
    vel[i].x = vext

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node lists
#-------------------------------------------------------------------------------
db = DataBase()
output("db")
for nodes in nodeSet:
    db.appendNodeList(nodes)
del nodes
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Construct the artificial viscosity.
#-------------------------------------------------------------------------------
if not gsph:
    q = Qconstructor(Cl, Cq)
    q.epsilon2 = epsilon2
    q.limiter = Qlimiter
    q.balsaraShearCorrection = balsaraCorrection
    output("q")
    output("q.Cl")
    output("q.Cq")
    output("q.epsilon2")
    output("q.limiter")
    output("q.balsaraShearCorrection")

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
'''if svph:
    hydro = SVPH(W = WT, 
                 Q = q,
                 cfl = cfl,
                 compatibleEnergyEvolution = compatibleEnergy,
                 densityUpdate = densityUpdate,
                 XSVPH = XSPH,
                 linearConsistent = linearConsistent,
                 generateVoid = False,
                 HUpdate = HUpdate,
                 fcentroidal = fcentroidal,
                 fcellPressure = fcellPressure,
                 xmin = Vector(xb0 - (xb1 - xb0), yb0 - (yb1 - yb0)),
                xmax = Vector(xb1 + (xb1 - xb0), yb1 + (yb1 - yb0)))
'''
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
    sumDensityNodeListSwitch =[outerNodes,innerNodes]  
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
    limiter = VanLeerLimiter()
    waveSpeed = DavisWaveSpeed()
    solver = HLLC(limiter,waveSpeed,gsphLinearCorrect)
    hydro = GSPH(dataBase = db,
                riemannSolver = solver,
                W = WT,
                cfl=cfl,
                specificThermalEnergyDiffusionCoefficient = gsphEpsDiffuseCoeff,
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
output("hydro.cfl")
output("hydro.compatibleEnergyEvolution")
output("hydro.densityUpdate")
output("hydro.HEvolution")

packages = [hydro]

#-------------------------------------------------------------------------------
# Construct the MMRV physics object.
#-------------------------------------------------------------------------------
if boolReduceViscosity:
    evolveReducingViscosityMultiplier = MorrisMonaghanReducingViscosity(nhQ,nhL,aMin,aMax)
    packages.append(evolveReducingViscosityMultiplier)
elif boolCullenViscosity:
    evolveCullenViscosityMultiplier = CullenDehnenViscosity(WT,alphMax,alphMin,betaC,betaD,betaE,fKern,boolHopkinsCorrection)
    packages.append(evolveCullenViscosityMultiplier)

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
xPlane0 = Plane(Vector(xb0, yb0), Vector( 1.0,  0.0))
xPlane1 = Plane(Vector(xb1, yb0), Vector(-1.0,  0.0))
yPlane0 = Plane(Vector(xb0, yb0), Vector( 0.0,  1.0))
yPlane1 = Plane(Vector(xb0, yb1), Vector( 0.0, -1.0))

xbc = PeriodicBoundary(xPlane0, xPlane1)
ybc = PeriodicBoundary(yPlane0, yPlane1)

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
integrator.cullGhostNodes = False
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
# Build the history object to track the cloud mass fraction.
#-------------------------------------------------------------------------------
epsExt = Pequi/((gamma - 1.0)*rhoext)
massFracHistory = CloudMassFraction(r0 = br,
                                    rhoThreshold = rhoThresholdFrac * rhoblob,
                                    epsThreshold = epsThresholdFrac * epsExt,
                                    nodes = innerNodes,
                                    filename = os.path.join(baseDir, histfilename))

#-------------------------------------------------------------------------------
# Make the problem controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle,
                            vizBaseName = vizBaseName,
                            vizDir = vizDir,
                            vizStep = vizCycle,
                            vizTime = vizTime,
                            skipInitialPeriodicWork = svph,
                            SPH = (not ASPH))
output("control")

control.appendPeriodicWork(massFracHistory.sample, massFracFreq)
massFracHistory.flushHistory()

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
if not steps is None:
    control.step(steps)

else:
    control.advance(goalTime, maxSteps)
    control.updateViz(control.totalSteps, integrator.currentTime, 0.0)
    control.dropRestartFile()
