#ATS:test(SELF, "--CRKSPH True --cfl 0.25 --Cl 1.0 --Cq 1.0 --clearDirectories True --filter 0.0 --goalTime 3.0", label="Gresho CRK, filter=0.0", np=10)
#ATS:test(SELF, "--CRKSPH True --cfl 0.25 --Cl 1.0 --Cq 1.0 --clearDirectories True --filter 0.01 --goalTime 3.0", label="Gresho CRK, filter=0.01", np=10)
#ATS:test(SELF, "--CRKSPH True --cfl 0.25 --Cl 1.0 --Cq 1.0 --clearDirectories True --filter 0.1 --goalTime 3.0", label="Gresho CRK, filter=0.1", np=10)
#ATS:test(SELF, "--CRKSPH True --cfl 0.25 --Cl 1.0 --Cq 1.0 --clearDirectories True --filter 0.2 --goalTime 3.0", label="Gresho CRK, filter=0.2", np=10)
#-------------------------------------------------------------------------------
# The Gresho-Vortex Test
#-------------------------------------------------------------------------------
import shutil
from math import *
from Spheral2d import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from findLastRestart import *
from GenerateNodeDistribution2d import *
from CubicNodeGenerator import GenerateSquareNodeDistribution
from CentroidalVoronoiRelaxation import *
from siloPointmeshDump import siloPointmeshDump
from SpheralConservation import SpheralConservation

import mpi
import DistributeNodes

title("2-D integrated hydro test --  Gresho-Vortex Test")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(
    rho = 1.0,
    gamma = 5.0/3.0,

    # Translation
    velTx=0.0,
    velTy=0.0,

    # Geometry of Box
    x0 = -0.5,
    x1 =  0.5,
    y0 = -0.5,
    y1 =  0.5,
   
    #Center of Vortex
    xc=0.0,
    yc=0.0,

    # Resolution and node seeding.
    nx1 = 64,
    ny1 = 64,
    seed = "lattice",

    nPerh = 1.51,

    SVPH = False,
    CRKSPH = False,
    PSPH = False,
    ASPH = False,   # This just chooses the H algorithm -- you can use this with CRKSPH for instance.
    filter = 0.0,  # For CRKSPH
    KernelConstructor = NBSplineKernel,
    order = 5,
    Qconstructor = MonaghanGingoldViscosity,
    #Qconstructor = TensorMonaghanGingoldViscosity,
    boolReduceViscosity = False,
    nhQ = 5.0,
    nhL = 10.0,
    aMin = 0.1,
    aMax = 2.0,
    boolCullenViscosity = False,
    alphMax = 2.0,
    alphMin = 0.02,
    betaC = 0.7,
    betaD = 0.05,
    betaE = 1.0,
    fKern = 1.0/3.0,
    boolHopkinsCorrection = True,
    linearConsistent = False,
    fcentroidal = 0.0,
    fcellPressure = 0.0,
    Cl = 1.0, 
    Cq = 0.75,
    etaCritFrac = 1.0,
    etaFoldFrac = 0.2,
    linearInExpansion = False,
    Qlimiter = False,
    balsaraCorrection = False,
    epsilon2 = 1e-2,
    hmin = 1e-5,
    hmax = 0.5,
    hminratio = 0.1,
    cfl = 0.5,
    XSPH = False,
    epsilonTensile = 0.0,
    nTensile = 8,

    IntegratorConstructor = CheapSynchronousRK2Integrator,
    goalTime = 1.0,
    steps = None,
    vizCycle = 20,
    vizTime = 0.1,
    dt = 0.0001,
    dtMin = 1.0e-5, 
    dtMax = 0.1,
    dtGrowth = 2.0,
    maxSteps = None,
    statsStep = 10,
    HUpdate = IdealH,
    domainIndependent = False,
    rigorousBoundaries = False,
    dtverbose = False,

    densityUpdate = RigorousSumDensity, # VolumeScaledDensity,
    correctionOrder = LinearOrder,
    volumeType = RKSumVolume,
    compatibleEnergy = True,
    gradhCorrection = True,
    correctVelocityGradient = True,
    HopkinsConductivity = False,     # For PSPH
    evolveTotalEnergy = False,       # Only for SPH variants -- evolve total rather than specific energy

    useVoronoiOutput = False,
    clearDirectories = False,
    restoreCycle = None,
    restartStep = 200,
    dataDir = "dumps-greshovortex-xy",
    graphics = True,
    smooth = None,
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

#-------------------------------------------------------------------------------
# Build our directory paths.
#-------------------------------------------------------------------------------
densityUpdateLabel = {IntegrateDensity : "IntegrateDensity",
                      SumDensity : "SumDensity",
                      RigorousSumDensity : "RigorousSumDensity",
                      SumVoronoiCellDensity : "SumVoronoiCellDensity",
                      VoronoiCellDensity : "VoronoiCellDensity"}
baseDir = os.path.join(dataDir,
                       HydroConstructor.__name__,
                       Qconstructor.__name__,
                       KernelConstructor.__name__,
                       "Cl=%g_Cq=%g" % (Cl, Cq),
                       densityUpdateLabel[densityUpdate],
                       "compatibleEnergy=%s" % compatibleEnergy,
                       "Cullen=%s" % boolCullenViscosity,
                       "nPerh=%3.1f" % nPerh,
                       "fcentroidal=%f" % max(fcentroidal, filter),
                       "fcellPressure=%f" % fcellPressure,
                       "%ix%i" % (nx1, ny1))
restartDir = os.path.join(baseDir, "restarts")
restartBaseName = os.path.join(restartDir, "greshovortex-xy-%ix%i" % (nx1, ny1))

vizDir = os.path.join(baseDir, "visit")
if vizTime is None and vizCycle is None:
    vizBaseName = None
else:
    vizBaseName = "greshovortex-xy-%ix%i" % (nx1, ny1)

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
# Make the NodeLists.
#-------------------------------------------------------------------------------
nodes = makeFluidNodeList("fluid", eos,
                               hmin = hmin,
                               hmax = hmax,
                               hminratio = hminratio,
                               kernelExtent = kernelExtent,
                               nPerh = nPerh)
output("nodes.name")
output("    nodes.hmin")
output("    nodes.hmax")
output("    nodes.hminratio")
output("    nodes.nodesPerSmoothingScale")

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
rmin = 0.0
rmax = sqrt(2.0)*(x1-x0)

if seed=="latticeCylindrical":
    rmin = x1-8.0*nPerh/nx1
    rmax = x1-2.0*nPerh/nx1

generator = GenerateNodeDistribution2d(nx1, ny1, rho,
                                       distributionType = seed,
                                       xmin = (x0, y0),
                                       xmax = (x1, y1),
                                       #rmin = 0.0,
                                       theta = 2.0*pi,
                                       #rmax = sqrt(2.0)*(x1 - x0),
                                       rmax = rmax,
                                       rmin = rmin,
                                       nNodePerh = nPerh,
                                       SPH = (not ASPH))

if mpi.procs > 1:
    from VoronoiDistributeNodes import distributeNodes2d
else:
    from DistributeNodes import distributeNodes2d

distributeNodes2d((nodes, generator))
print(nodes.name, ":")
output("    mpi.reduce(nodes.numInternalNodes, mpi.MIN)")
output("    mpi.reduce(nodes.numInternalNodes, mpi.MAX)")
output("    mpi.reduce(nodes.numInternalNodes, mpi.SUM)")

#Set IC
vel = nodes.velocity()
eps = nodes.specificThermalEnergy()
pos = nodes.positions()
for i in range(nodes.numInternalNodes):
    xi, yi = pos[i]
    r2=(xi-xc)*(xi-xc)+(yi-yc)*(yi-yc)
    ri=sqrt(r2)
    vphi=0.0
    sinPhi=(yi-yc)/ri
    cosPhi=(xi-xc)/ri
    Pi=0.0
    if ri < 0.2:
       Pi=5.0+12.5*r2
       vphi=5*ri
    elif ri < 0.4 and ri >= 0.2:
       Pi=9.0+12.5*r2-20.0*ri+4.0*log(5.0*ri)
       vphi=2.0-5.0*ri
    else:
       Pi=3.0+4*log(2.0)
       #vphi is zero
    velx=velTx-vphi*sinPhi #translation velocity + azimuthal velocity 
    vely=velTy+vphi*cosPhi
    vel[i]=Vector(velx,vely)
    eps0 = Pi/((gamma - 1.0)*rho)
    eps[i]=eps0

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node lists
#-------------------------------------------------------------------------------
db = DataBase()
output("db")
db.appendNodeList(nodes)
output("db.numNodeLists")
output("db.numFluidNodeLists")

#------------------------------------------------------------------------------
# Build the conservation check.
#------------------------------------------------------------------------------
conserve = SpheralConservation(db)
print("Initial conservation check:")
for lab, q in (("M", conserve.massHistory), 
               ("pmom.x", [x.x for x in conserve.pmomHistory]), 
               ("pmom.y", [x.y for x in conserve.pmomHistory]), 
               ("KE", conserve.KEHistory),
               ("TE", conserve.TEHistory), 
               ("E", conserve.EHistory)):
    print("%10s: %g %s" % (lab, (q[-1] - q[0])/q[0], q))

#------------------------------------------------------------------------------
# Write the starting conditions.
#------------------------------------------------------------------------------
siloPointmeshDump(baseName = "GreshoOverlay_initial_nodes",
                  fields = [nodes.mass(), nodes.Hfield(),
                            nodes.velocity(), nodes.specificThermalEnergy()],
                  baseDirectory = "GreshoOverlay_initial_nodes")

#------------------------------------------------------------------------------
# Make a different generator to map to.
#------------------------------------------------------------------------------
generator2 = GenerateNodeDistribution2d(int(2*nx1), int(2*ny1), rho,
                                        distributionType = seed,
                                        xmin = (x0, y0),
                                        xmax = (x1, y1),
                                        #rmin = 0.0,
                                        theta = 2.0*pi,
                                        #rmax = sqrt(2.0)*(x1 - x0),
                                        rmax = rmax,
                                        rmin = rmin,
                                        nNodePerh = nPerh,
                                        SPH = (not ASPH))

# Do the mapping.
from overlayNodeList import overlayNodeList
overlayNodeList(nodes, generator2,
                removeUnusedNodes = False)

#------------------------------------------------------------------------------
# Conservation check.
#------------------------------------------------------------------------------
conserve.updateHistory()
print("Final conservation check:")
for lab, q in (("M", conserve.massHistory), 
               ("pmom.x", [x.x for x in conserve.pmomHistory]), 
               ("pmom.y", [x.y for x in conserve.pmomHistory]), 
               ("KE", conserve.KEHistory),
               ("TE", conserve.TEHistory), 
               ("E", conserve.EHistory)):
    print("%10s: %g %s" % (lab, (q[-1] - q[0])/q[0], q))

#------------------------------------------------------------------------------
# Write the overlayed results.
#------------------------------------------------------------------------------
siloPointmeshDump(baseName = "GreshoOverlay_overlay_nodes",
                  fields = [nodes.mass(), nodes.Hfield(),
                            nodes.velocity(), nodes.specificThermalEnergy()],
                  baseDirectory = "GreshoOverlay_overlay_nodes")
