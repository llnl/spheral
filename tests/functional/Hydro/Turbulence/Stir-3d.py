#-------------------------------------------------------------------------------
# 3-D stirred box test
#
#
#-------------------------------------------------------------------------------
from math import *
from Spheral3d import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from findLastRestart import *
from SpheralVisitDump import dumpPhysicsState
import mpi

from GenerateNodeDistribution3d import *

title("3-D stirred box turbulence test")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(seed = "lattice",

            nsize = 100,
            
            f_solenoidal = 0.66,
            kmin = 1,
            kmax = 50, #must be <= nsize/2
            perturbEach = 1,

            rho1 = 1.0,
            eps1 = 0.0,
            nPerh = 2.01,

            gamma = 5.0/3.0,
            mu = 1.0,
            
            SVPH = False,
            CRKSPH = False,
            ASPH = False,
            SPH = True,   # This just chooses the H algorithm -- you can use this with CRKSPH for instance.
            filter = 0.0,   # CRKSPH filtering
            Qconstructor = MonaghanGingoldViscosity,
            #Qconstructor = TensorMonaghanGingoldViscosity,
            linearConsistent = False,
            fcentroidal = 0.0,
            fcellPressure = 0.0,
            boolReduceViscosity = False,
            nh = 5.0,
            aMin = 0.1,
            aMax = 2.0,
            Qhmult = 1.0,
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
            dataDir = "stir-3d",
            outputFile = None,
            comparisonFile = None,
            
            serialDump = False, #whether to dump a serial ascii file at the end for viz
            )

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
else:
  if ASPH:
    HydroConstructor = ASPHHydro
  else:
    HydroConstructor = SPHHydro

xmin = (0.0, 0.0, 0.0)
xmax = (1.0, 1.0, 1.0)
nx = nsize
ny = nsize
nz = nsize

n = [nx,ny,nz]

dataDir = os.path.join(dataDir,
                       "rho1=%g" % rho1,
                       str(HydroConstructor).split("'")[1].split(".")[-1],
                       "densityUpdate=%s" % (densityUpdate),
                       "compatibleEnergy=%s" % (compatibleEnergy),
                       "XSPH=%s" % XSPH,
                       "filter=%s" % filter,
                       "%s-Cl=%g-Cq=%g" % (str(Qconstructor).split("'")[1].split(".")[-1], Cl, Cq),
                       "%ix%ix%i" % (nx, ny, nz),
                       "nPerh=%g-Qhmult=%g" % (nPerh, Qhmult))
restartDir = os.path.join(dataDir, "restarts")
vizDir = os.path.join(dataDir, "visit")
restartBaseName = os.path.join(restartDir, "stir-3d")
vizBaseName = "stir-3d"

#-------------------------------------------------------------------------------
# Helper functions for the generation of perturbations
#-------------------------------------------------------------------------------
def init_perturbations(dtype):
  kx = np.zeros(n, dtype=dtype)
  ky = np.zeros(n, dtype=dtype)
  kz = np.zeros(n, dtype=dtype)
  # perform fft k-ordering convention shifts
  for j in range(0,n[1]):
    for k in range(0,n[2]):
      kx[:,j,k] = n[0]*np.fft.fftfreq(n[0])
  for i in range(0,n[0]):
    for k in range(0,n[2]):
      ky[i,:,k] = n[1]*np.fft.fftfreq(n[1])
  for i in range(0,n[0]):
    for j in range(0,n[1]):
      kz[i,j,:] = n[2]*np.fft.fftfreq(n[2])
  
  kx = np.array(kx, dtype=dtype)
  ky = np.array(ky, dtype=dtype)
  kz = np.array(kz, dtype=dtype)
  k = np.sqrt(np.array(kx**2+ky**2+kz**2, dtype=dtype))
  
  # only use the positive frequencies
  inds = np.where(np.logical_and(k**2 >= kmin**2, k**2 < (kmax+1)**2))
  nr = len(inds[0])
  
  phasex = np.zeros(n, dtype=dtype)
  phasex[inds] = 2.*pi*np.random.uniform(size=nr)
  fx = np.zeros(n, dtype=dtype)
  fx[inds] = np.random.normal(size=nr)
  
  phasey = np.zeros(n, dtype=dtype)
  phasey[inds] = 2.*pi*np.random.uniform(size=nr)
  fy = np.zeros(n, dtype=dtype)
  fy[inds] = np.random.normal(size=nr)
  
  phasez = np.zeros(n, dtype=dtype)
  phasez[inds] = 2.*pi*np.random.uniform(size=nr)
  fz = np.zeros(n, dtype=dtype)
  fz[inds] = np.random.normal(size=nr)
  
  # rescale perturbation amplitude so that low number statistics
  # at low k do not throw off the desired power law scaling.
  for i in range(kmin, kmax+1):
    slice_inds = np.where(np.logical_and(k >= i, k < i+1))
    rescale = sqrt(np.sum(np.abs(fx[slice_inds])**2 + np.abs(fy[slice_inds])**2 + np.abs(fz[slice_inds])**2))
    fx[slice_inds] = fx[slice_inds]/rescale
    fy[slice_inds] = fy[slice_inds]/rescale
    fz[slice_inds] = fz[slice_inds]/rescale
  
  # set the power law behavior
  # wave number bins
  fx[inds] = fx[inds]*k[inds]**-(0.5*alpha)
  fy[inds] = fy[inds]*k[inds]**-(0.5*alpha)
  fz[inds] = fz[inds]*k[inds]**-(0.5*alpha)
  
  # add in phases
  fx = np.cos(phasex)*fx + 1j*np.sin(phasex)*fx
  fy = np.cos(phasey)*fy + 1j*np.sin(phasey)*fy
  fz = np.cos(phasez)*fz + 1j*np.sin(phasez)*fz
  
  return fx, fy, fz, kx, ky, kz

def normalize(fx, fy, fz):
  norm = np.sqrt(np.sum(fx**2 + fy**2 + fz**2)/np.product(n))
  fx = fx/norm
  fy = fy/norm
  fz = fz/norm
  return fx, fy, fz

def make_perturbations():
  fx, fy, fz, kx, ky, kz = init_perturbations(n, kmin, kmax, dtype)
  if f_solenoidal != None:
    k2 = kx**2+ky**2+kz**2
    # solenoidal part
    fxs = 0.; fys =0.; fzs = 0.
    if f_solenoidal != 0.0:
      fxs = np.real(fx - kx*(kx*fx+ky*fy+kz*fz)/np.maximum(k2,1e-16))
      fys = np.real(fy - ky*(kx*fx+ky*fy+kz*fz)/np.maximum(k2,1e-16))
      fzs = np.real(fz - kz*(kx*fx+ky*fy+kz*fz)/np.maximum(k2,1e-16))
      ind = np.where(k2 == 0)
      fxs[ind] = 0.; fys[ind] = 0.; fzs[ind] = 0.
      # need to normalize this before applying relative weighting of solenoidal / compressive components
      norm = np.sqrt(np.sum(fxs**2+fys**2+fzs**2))
      fxs = fxs/norm
      fys = fys/norm
      fzs = fzs/norm
    # compressive part
    # get a different random cube for the compressive part
    # so that we can target the RMS solenoidal fraction,
    # instead of setting a constant solenoidal fraction everywhere.
    fx, fy, fz, kx, ky, kz = init_perturbations(dtype)
    fxc = 0.; fyc =0.; fzc = 0.
    if f_solenoidal != 1.0:
      fxc = np.real(kx*(kx*fx+ky*fy+kz*fz)/np.maximum(k2,1e-16))
      fyc = np.real(ky*(kx*fx+ky*fy+kz*fz)/np.maximum(k2,1e-16))
      fzc = np.real(kz*(kx*fx+ky*fy+kz*fz)/np.maximum(k2,1e-16))
      ind = np.where(k2 == 0)
      fxc[ind] = 0.; fyc[ind] = 0.; fzc[ind] = 0.
      # need to normalize this before applying relative weighting of solenoidal / compressive components
      norm = np.sqrt(np.sum(fxc**2+fyc**2+fzc**2))
      fxc = fxc/norm
      fyc = fyc/norm
      fzc = fzc/norm
    # back to real space
    pertx = np.real(np.fft.ifftn(f_solenoidal*fxs + (1.-f_solenoidal)*fxc))
    perty = np.real(np.fft.ifftn(f_solenoidal*fys + (1.-f_solenoidal)*fyc))
    pertz = np.real(np.fft.ifftn(f_solenoidal*fzs + (1.-f_solenoidal)*fzc))
  else:
    # just convert to real space
    pertx = np.real(np.fft.ifftn(fx))
    perty = np.real(np.fft.ifftn(fy))
    pertz = np.real(np.fft.ifftn(fz))
  
  # subtract off COM (assuming uniform density)
  pertx = pertx-np.average(pertx)
  perty = perty-np.average(perty)
  pertz = pertz-np.average(pertz)
  # scale RMS of perturbation cube to unity
  pertx, perty, pertz = normalize(pertx, perty, pertz)
  return pertx, perty, pertz

#-------------------------------------------------------------------------------
# Periodic work function
#-------------------------------------------------------------------------------
class perturb(object):
  def __init__(self,nodeSet,directory):
    self.nodeSet = nodeSet
    self.directory = directory
  def __call__(self, cycle, time, dt):
    pertx, perty, pertz = make_perturbations()


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
output("WT")
kernelExtent = WT.kernelExtent

#-------------------------------------------------------------------------------
# Make the NodeList.
#-------------------------------------------------------------------------------
nodes1 = makeFluidNodeList("High density gas", eos,
                           hmin = hmin,
                           hmax = hmax,
                           hminratio = hminratio,
                           nPerh = nPerh,
                           kernelExtent = kernelExtent)
output("nodes1.nodesPerSmoothingScale")
output("nodes1.hmin")
output("nodes1.hmax")
output("nodes1.hminratio")

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
if restoreCycle is None:
    generator1 = GenerateNodeDistribution3d(nx, ny, nz, rho1, seed,
                                            xmin = xmin,
                                            xmax = xmax,
                                            nNodePerh = nPerh)

    if mpi.procs > 1:
      from VoronoiDistributeNodes import distributeNodes3d
    else:
      from DistributeNodes import distributeNodes3d
    distributeNodes3d((nodes1, generator1))
    output("mpi.reduce(nodes1.numInternalNodes, mpi.MIN)")
    output("mpi.reduce(nodes1.numInternalNodes, mpi.MAX)")
    output("mpi.reduce(nodes1.numInternalNodes, mpi.SUM)")

    # Set node specific thermal energies
    nodes1.specificThermalEnergy(ScalarField("tmp", nodes1, eps1))

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase3d()
output("db")
db.appendNodeList(nodes1)
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
                           xmin = Vector(-2.0, -2.0, -2.0),
                           xmax = Vector(3.0, 3.0, 3.0))
# xmin = Vector(x0 - 0.5*(x2 - x0), y0 - 0.5*(y2 - y0)),
# xmax = Vector(x2 + 0.5*(x2 - x0), y2 + 0.5*(y2 - y0)))
elif CRKSPH:
  hydro = HydroConstructor(W = WT,
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
# Create boundary conditions.
#-------------------------------------------------------------------------------
xPlane0 = Plane(Vector(*xmin), Vector(1.0, 0.0, 0.0))
xPlane1 = Plane(Vector(*xmax), Vector(-1.0,0.0, 0.0))

yPlane0 = Plane(Vector(*xmin), Vector(0.0, 1.0, 0.0))
yPlane1 = Plane(Vector(*xmax), Vector(0.0, -1.0, 0.0))

zPlane0 = Plane(Vector(*xmin), Vector(0.0, 0.0, 1.0))
zPlane1 = Plane(Vector(*xmax), Vector(0.0, 0.0, -1.0))

xbc = PeriodicBoundary(xPlane0, xPlane1)
ybc = PeriodicBoundary(yPlane0, yPlane1)
zbc = PeriodicBoundary(zPlane0, zPlane1)

for p in packages:
    p.appendBoundary(xbc)
    p.appendBoundary(ybc)
    p.appendBoundary(zbc)

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
pert = perturb([nodes1],dataDir)
control.appendPeriodicWork(pert,perturbEach)
output("control")

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
control.advance(goalTime)
