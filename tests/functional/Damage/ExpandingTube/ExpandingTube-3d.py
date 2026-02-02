#-------------------------------------------------------------------------------
# A tube of stainless steel undergoing expansion due to a projectile entering
# and smacking into a stop.
#
# See Vogler et al., 2003, International Journal of Impact Engineering, 29, 735
#-------------------------------------------------------------------------------
import sys
sys.path.append("../Utilities")
from Spheral3d import *
from SpheralTestUtilities import *
from SpheralController import *
from math import *

from bevelTubeEntrance import *
from NodeHistory import NodeHistory
from AverageStrain import AverageStrain
import mpi

#-------------------------------------------------------------------------------
# We'll work in CGuS (cm-gm-microsec) units
#-------------------------------------------------------------------------------
units = CGuS()

#-------------------------------------------------------------------------------
# Create a physics package to dump excessive velocities into thermal energy.
#-------------------------------------------------------------------------------
class VelocityDiffuser(Physics):

    def __init__(self, vmax):
        self.vmax = vmax
        self.vmax2 = vmax*vmax
        Physics.__init__(self)
        return

    def evaluateDerivatives(self, time, dt,dataBase, state, derivs):
        return

    def dt(self, dataBase, state, derivs, time):
        return pair_double_string(1e30, "No vote")

    def registerState(self, dataBase, state):
        return

    def registerDerivatives(self, dataBase, derivs):
        return

    def finalize(self, time, dt, dataBase, state, derivs):
        velocities = state.vectorFields(HydroFieldNames.velocity)
        specificEnergies = state.scalarFields(HydroFieldNames.specificThermalEnergy)
        assert len(velocities) == len(specificEnergies)
        for vel, eps in zip(velocities, specificEnergies):
            assert vel.numInternalElements == eps.numInternalElements
            for i in range(vel.numInternalElements):
                vmag2 = vel[i].magnitude2()
                if vmag2 > self.vmax2:
                    deps = 0.5*(vmag2 - self.vmax2)
                    assert deps > 0.0
                    vel[i] *= self.vmax/sqrt(vmag2)
                    eps[i] += deps
        return
        
#-------------------------------------------------------------------------------
# Identify ourselves!
#-------------------------------------------------------------------------------
title("3-D expanding tube impact strength/damage model test")

#-------------------------------------------------------------------------------
# Generic problem parameters
# All CGuS units.
#-------------------------------------------------------------------------------
commandLine(
    # How much of the 2 Pi geometry are we doing?
    phiFactor = 0.5,

    # Should we start the projectile at the end of the tube or in contact with
    # the plug?
    projectileOutside = True,

    # How much should we compress the projectile to allow it to slide into the tube?
    compressProjectile = 0.3,

    # Geometry
    tubeThickness = 0.3,       # cm
    rtubeInner = 0.5*1.27,     # cm
    ltube = 5.08,              # cm
    lCuAnvil = 1.0,            # cm
    lFoamAnvil = 1.0,          # cm
    lSteelAnvil = 1.0,         # cm

    # Numbers of nodes.
    nrtube = 15,
    nltube = 170, # 254,
    nrAnvil = 100,
    nrSteelAnvilCap = 10,

    # VISAR sampling parameters.
    dzVISAR = 0.04,            # cm
    drVISAR = 0.04,            # cm
    VISARsampleFrequency = 10,

    # Frequency for measuring the strain in the tube.
    strainFrequency = 5,

    # Inital velocity of the projectile.
    vzproj = -1.75e-1, # -1.92e-1,    # cm/usec

    # Maximum speed we're going to allow.
    vmax = None,   # 3e-1,            # cm/usec

    # Mass densities.
    rho0Steel = 7.94, # 7.85          # gm/cc
    rho0Lexan = 1.196,                # gm/cc
    rho0Air = 0.001205,               # gm/cc
    rho0Foam = 1.046,                 # gm/cc
    rho0Cu = 8.93,                    # gm/cc

    # Parameters for the damage model.
    DamageModelConstructor = ProbabilisticDamageModel,
    coupleType = PairMaxDamage,
    strainType = PseudoPlasticStrain,
    kWeibullSteelFactor = 1.0,
    mWeibullSteelFactor = 1.0,
    randomSeedSteel = 109482993,
    criticalDamageThreshold = 0.5,

    # Node seeding stuff.
    nPerh = 3.01,

    # Material specific bounds on the mass density.
    etaMinSteel = 0.9,
    etaMaxSteel = 1.1,
    etaMinLexan = 0.5,
    etaMaxLexan = 1.5,
    etaMinCu = 0.7,
    etaMaxCu = 1.3,
    etaMinFoam = 0.5,
    etaMaxFoam = 1.5,

    # Material specific limits on the smoothing scale.
    hminratio = 0.1,
    hmin = None,
    hmax = None,

    # Hydro parameters.
    hydroType = "FSISPH",
    cfl = 0.5,
    epsilonTensile = 0.3,
    nTensile = 4,
    HUpdate = IdealH,
    correctionOrder = LinearOrder,
    volumeType = RKSumVolume,
    asph = True,
    densityUpdate = IntegrateDensity,
    compatibleEnergy = True,
    evolveTotalEnergy = False,
    XSPH = True,
    gradhCorrection = False,
    useVelocityMagnitudeForDt = True,

    # Times, and simulation control.
    steps = None,
    goalTime = 100.0,       # usec
    dtSample = 0.5,         # usec
    dt = 5e-3,
    dtMin = None,
    dtMax = None,
    dtGrowth = 10.0,
    maxSteps = 200,
    statsStep = 10,
    redistributeStep = 10000,
    smoothIters = 0,
    dtverbose = False,

    # Restart and output files.
    restoreCycle = -1,
    restartStep = 200,
    baseDir = "dumps-expandingTube-3d",
    vizBaseName = "ExpandingTube-3d",
    vizStep = 1000,
    vizTime = 1.0,
    )

# Derived geometry.
assert phiFactor in (0.5, 1.0, 2.0)
phi = phiFactor * pi
rtubeOuter = rtubeInner + tubeThickness
rplug, lplug = rtubeInner, 0.5*ltube
rproj, lproj = rplug, lplug
rAnvil = 2.0*rtubeOuter
lAnvil = lCuAnvil + lFoamAnvil + lSteelAnvil

hydroType = hydroType.upper()

# Use the above geometry to define enclosing points of the materials for the
# node generators.
rminTube = rtubeInner
rmaxTube = rtubeOuter
zminTube = lAnvil
zmaxTube = lAnvil + ltube

rminPlug = 0.0
rmaxPlug = rplug
zminPlug = lAnvil
zmaxPlug = lAnvil + lplug

rminProj = 0.0
rmaxProj = rproj
if projectileOutside:
    zminProj = lAnvil + ltube
else:
    zminProj = zmaxPlug
zmaxProj = zminProj + lproj

rminInnerAir = 0.0
rmaxInnerAir = rtubeInner
zminInnerAir = lAnvil + lplug
zmaxInnerAir = lAnvil + ltube

rminSteelAnvil = 0.0
rmaxSteelAnvil = rAnvil
zminSteelAnvil = 0.0
zmaxSteelAnvil = lSteelAnvil

rminFoamAnvil = 0.0
rmaxFoamAnvil = rAnvil
zminFoamAnvil = lSteelAnvil
zmaxFoamAnvil = lSteelAnvil + lFoamAnvil

rminCuAnvil = 0.0
rmaxCuAnvil = rAnvil
zminCuAnvil = lSteelAnvil + lFoamAnvil
zmaxCuAnvil = lAnvil

rminSteelAnvilCap = rAnvil
rmaxSteelAnvilCap = rAnvil + float(nrSteelAnvilCap)/float(nrAnvil)*rAnvil
zminSteelAnvilCap = 0.0
zmaxSteelAnvilCap = lAnvil

# The geometry of the bevel on the inner tube opening surface.
tubeOpeningAngle = 1.8 * pi/180.0 # radians
zBevelBegin = lAnvil + ltube - 0.6

# Define the VISAR sampling points.
zVISARa = lAnvil + 2.5
zVISARb = lAnvil + 2.0
zVISARc = lAnvil + 1.5

rVISARa = rtubeOuter
rVISARb = rtubeOuter
rVISARc = rtubeOuter

# Derived numbers of nodes.
nrplug, nlplug = int(rplug/lplug * nltube/4), nltube//4
nrproj, nlproj = nrplug, nlplug
nlAnvil = int(nrAnvil * lAnvil/rAnvil)
nlSteelAnvil = int(nlAnvil * 0.5*lSteelAnvil/lAnvil + 0.5)
nlFoamAnvil = int(nlAnvil * 0.25*lFoamAnvil/lAnvil + 0.5)
nlCuAnvil = int(nlAnvil * 0.5*lCuAnvil/lAnvil + 0.5)
nrSteelAnvil = int(nlSteelAnvil * rAnvil/lSteelAnvil + 0.5)
nrFoamAnvil = int(nlFoamAnvil * rAnvil/lFoamAnvil + 0.5)
nrCuAnvil = int(nlCuAnvil * rAnvil/lCuAnvil + 0.5)
nlSteelAnvilCap = nlSteelAnvil + nlFoamAnvil + nlCuAnvil

# Inital velocity of the projectile.
v0proj = Vector(0.0, 0.0, vzproj)

# Damage model parameters.
kWeibullSteel = 8.8e5 * kWeibullSteelFactor
mWeibullSteel = 2.63  * mWeibullSteelFactor
volumeSteel = pi*(rtubeOuter**2 - rtubeInner**2)*ltube

# Restart and output files.
dataDir = os.path.join(baseDir,
                       hydroType,
                       "{}_{}_{}".format(DamageModelConstructor.__name__, str(strainType), str(coupleType)),
                       #str(DamageModelConstructor).split("'")[1],
                       "vzproj=%4.2e" % abs(vzproj),
                       "rho0Steel=%8.6f" % rho0Steel,
                       "k=%4.2f_m=%4.2f" % (kWeibullSteel, mWeibullSteel))
restartDir = os.path.join(dataDir, "restarts", "proc-%04i" % mpi.rank)
vizDir = os.path.join(dataDir, "visit")
restartBaseName = os.path.join(restartDir, "ExpandingTube")

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
import os, sys
if mpi.rank == 0:
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
# Stainless steel material properties.
# This attempts to mock the hardened steel properties.
#-------------------------------------------------------------------------------
## eosSteel = LinearPolynomialEquationOfStateCGS(rho0Steel,    # reference density  
##                                               etaMinSteel,  # etamin             
##                                               etaMaxSteel,  # etamax             
##                                               0.0,          # A0
##                                               1.649901e12,  # A1
##                                               1.674656e12,  # A2
##                                               0.832543e12,  # A3
##                                               1.93,         # B1
##                                               0.5,          # B2
##                                               0.0,          # B3
##                                               55.350)       # atomic weight

eosSteel = GruneisenEquationOfState(rho0Steel,    # reference density  
                                    etaMinSteel,  # etamin             
                                    etaMaxSteel,  # etamax             
                                    0.4529,       # C0                 
                                    1.50,         # S1                 
                                    0.0,          # S2                 
                                    0.0,          # S3                 
                                    1.84,         # gamma0             
                                    0.302,        # b                  
                                    55.350,       # atomic weight
                                    constants = units)

coldFitSteel = NinthOrderPolynomialFit(-1.06797724e-2,
                                       -2.06872020e-2,
                                        8.24893246e-1,
                                       -2.39505843e-2,
                                       -2.44522017e-2,
                                        5.38030101e-2,
                                        0.0,
                                        0.0,
                                        0.0,
                                        0.0)

meltFitSteel = NinthOrderPolynomialFit( 7.40464217e-2,
                                        2.49802214e-1,
                                        1.00445029,
                                       -1.36451475e-1,
                                        7.72897829e-3,
                                        5.06390305e-2,
                                        0.0,
                                        0.0,
                                        0.0,
                                        0.0)

strengthSteel = SteinbergGuinanStrength(eosSteel,
                                        7.700000e-1,        # G0
                                        2.2600,             # A
                                        4.5500e-04,         # B
                                        3.4000e-3,          # Y0
                                        2.5e-2,             # Ymax
                                        1.0e-3,             # Yp
                                        43.0000,            # beta
                                        0.0,                # gamma0
                                        0.35,               # nhard
                                        coldFitSteel,
                                        meltFitSteel)

#-------------------------------------------------------------------------------
# Lexan material properties.
# Note for lack of strength information about this material, I'm subsituting in
# the strength paramters for lucite here.  :)
#-------------------------------------------------------------------------------
eosLexan = GruneisenEquationOfState(rho0Lexan,    # reference density  
                                    etaMinLexan,  # etamin             
                                    etaMaxLexan,  # etamax             
                                    0.1933,       # C0                 
                                    3.49,         # S1                 
                                   -8.2,          # S2                 
                                    9.6,          # S3                 
                                    0.61,         # gamma0             
                                    0.0,          # b                  
                                    28423.0,      # atomic weight
                                    constants = units)

## strengthLexan = ConstantStrength(0.05e12,   # G0
##                                  0.0005e12) # Y0

coldFitLexan = NinthOrderPolynomialFit(-5.19191852e-3,
                                       -4.41500192e-3,
                                        2.84720528e-2,
                                        2.14093899e-2,
                                       -4.46412259e-3,
                                        1.24495222e-3,
                                        0.0,
                                        0.0,
                                        0.0,
                                        0.0)
meltFitLexan = NinthOrderPolynomialFit( 5.24383771e-4,
                                        1.49188457e-3,
                                        2.85704428e-2,
                                        2.13783662e-2,
                                       -4.45135120e-3,
                                        1.24138074e-3,
                                        0.0,
                                        0.0,
                                        0.0,
                                        0.0)
strengthLexan = SteinbergGuinanStrength(eosLexan,
                                        2.320000e-2,        # G0
                                        0.0,                # A
                                        0.0,                # B
                                        4.2000e-3,          # Y0
                                        1.0,                # Ymax
                                        1.0e-3,             # Yp
                                        0.0,                # beta
                                        0.0,                # gamma0
                                        0.0,                # nhard
                                        coldFitLexan,
                                        meltFitLexan)

#-------------------------------------------------------------------------------
# Copper material properties.
#-------------------------------------------------------------------------------
eosCu = GruneisenEquationOfState(rho0Cu,      # reference density  
                                 etaMinCu,    # etamin             
                                 etaMaxCu,    # etamax             
                                 0.394,       # C0                 
                                 1.489,       # S1                 
                                 0.0,         # S2                 
                                 0.0,         # S3                 
                                 2.02,        # gamma0             
                                 0.47,        # b                  
                                 63.57,       # atomic weight
                                 constants = units)
coldFitCu = NinthOrderPolynomialFit(-1.05111874e-2,
                                    -2.13429672e-2,
                                     6.92768584e-1,
                                    -2.45626513e-2,
                                    -2.48677403e-2,
                                     4.35373677e-2,
                                     0.0,
                                     0.0,
                                     0.0,
                                     0.0)
meltFitCu = NinthOrderPolynomialFit( 5.22055639e-2,
                                     1.90143176e-1,
                                     8.51351901e-1,
                                    -1.12049022e-1,
                                    -6.11436674e-3,
                                     4.36007831e-2,
                                     0.0,
                                     0.0,
                                     0.0,
                                     0.0)
strengthCu = SteinbergGuinanStrength(eosCu,
                                     4.770000e-1,        # G0
                                     2.8300,             # A
                                     3.7700e-04,         # B
                                     1.2000e-3,          # Y0
                                     6.4000e-3,          # Ymax
                                     1.0e-3,             # Yp
                                     36.0000,            # beta
                                     0.0,                # gamma0
                                     0.45,               # nhard
                                     coldFitCu,
                                     meltFitCu)

#-------------------------------------------------------------------------------
# Foam material properties. (Polystyrene CH)
#-------------------------------------------------------------------------------
eosFoam = GruneisenEquationOfState(rho0Foam,   # reference density  
                                   etaMinFoam, # etamin             
                                   etaMaxFoam, # etamax             
                                   0.189,      # C0                 
                                   2.965,      # S1                 
                                  -4.069,      # S2                 
                                   2.328,      # S3                 
                                   0.67,       # gamma0             
                                   0.0,        # b                  
                                   6.982,      # atomic weight
                                   constants = units)
strengthFoam = NullStrength()

#-------------------------------------------------------------------------------
# Air material properties.
#-------------------------------------------------------------------------------
gammaAir = 1.4
molecularWeightAir = 30.0
eosAir = GammaLawGas(gammaAir, molecularWeightAir,
                     constants = units)

#-------------------------------------------------------------------------------
# Create our interpolation kernels -- one for normal hydro interactions, and
# one for use with the artificial viscosity
#-------------------------------------------------------------------------------
WT = TableKernel(WendlandC4Kernel(), 100)
output("WT")

#-------------------------------------------------------------------------------
# Create the NodeLists.
#-------------------------------------------------------------------------------
nodesSteel = makeSolidNodeList("Stainless steel", eosSteel, strengthSteel,
                               kernelExtent = WT.kernelExtent)
nodesPlug = makeSolidNodeList("Lexan plug", eosLexan, strengthLexan, 
                              kernelExtent = WT.kernelExtent)
nodesProj = makeSolidNodeList("Lexan projectile", eosLexan, strengthLexan,
                              kernelExtent = WT.kernelExtent)
nodesSteelAnvil = makeSolidNodeList("Anvil (Steel)", eosSteel, strengthSteel, 
                                    kernelExtent = WT.kernelExtent)
nodesFoamAnvil = makeSolidNodeList("Anvil (Foam)", eosFoam, strengthFoam, 
                                   kernelExtent = WT.kernelExtent)
nodesCuAnvil = makeSolidNodeList("Anvil (Copper)", eosCu, strengthCu,
                                 kernelExtent = WT.kernelExtent)
nodeSet = [nodesSteel, nodesPlug, nodesProj,
           nodesSteelAnvil, nodesFoamAnvil, nodesCuAnvil]
for (n, etamin, etamax, rho0) in zip(nodeSet,
                                     [etaMinSteel, etaMinLexan, etaMinLexan, etaMinSteel, etaMinFoam, etaMinCu],
                                     [etaMaxSteel, etaMaxLexan, etaMaxLexan, etaMaxSteel, etaMaxFoam, etaMaxCu],
                                     [rho0Steel, rho0Lexan, rho0Lexan, rho0Steel, rho0Foam, rho0Cu]):
    n.nodesPerSmoothingScale = nPerh
    if hmin:
        n.hmin = hmin
    if hmax:
        n.hmax = hmax
    if hminratio:
        n.hminratio = hminratio
    n.rhoMin = etamin*rho0
    n.rhoMax = etamax*rho0
    output("n.name")
    output("  n.nodesPerSmoothingScale")
    output("  n.hmin")
    output("  n.hmax")
    output("  n.hminratio")
    output("  n.rhoMin")
    output("  n.rhoMax")

#-------------------------------------------------------------------------------
# Set node properties (positions, velocites, etc.)
#-------------------------------------------------------------------------------
print("Generating node distribution.")
from GenerateNodeDistribution3d import *
from CompositeNodeDistribution import *
from PeanoHilbertDistributeNodes import distributeNodes3d
generatorTube = GenerateNodeDistribution3d(nrtube,
                                           nltube,
                                           0,
                                           rho0Steel,
                                           "cylindrical",
                                           rmin = rminTube,
                                           rmax = rmaxTube,
                                           thetamin = 0,
                                           thetamax = phi,
                                           zmin = zminTube,
                                           zmax = zmaxTube,
                                           nNodePerh = nPerh,
                                           SPH = not asph)
generatorPlug = GenerateNodeDistribution3d(nrplug,
                                           nlplug,
                                           0,
                                           rho0Lexan,
                                           "cylindrical",
                                           rmin = rminPlug,
                                           rmax = rmaxPlug,
                                           thetamin = 0,
                                           thetamax = phi,
                                           zmin = zminPlug,
                                           zmax = zmaxPlug,
                                           nNodePerh = nPerh,
                                           SPH = not asph)
generatorProj = GenerateNodeDistribution3d(nrproj,
                                           nlproj,
                                           0,
                                           rho0Lexan,
                                           "cylindrical",
                                           rmin = rminProj,
                                           rmax = rmaxProj,
                                           thetamin = 0,
                                           thetamax = phi,
                                           zmin = zminProj,
                                           zmax = zmaxProj,
                                           nNodePerh = nPerh,
                                           SPH = not asph)
generatorSteelAnvil1 = GenerateNodeDistribution3d(nrSteelAnvil,
                                                  nlSteelAnvil,
                                                  0,
                                                  rho0Steel,
                                                  "cylindrical",
                                                  rmin = rminSteelAnvil,
                                                  rmax = rmaxSteelAnvil,
                                                  thetamin = 0,
                                                  thetamax = phi,
                                                  zmin = zminSteelAnvil,
                                                  zmax = zmaxSteelAnvil,
                                                  nNodePerh = nPerh,
                                                  SPH = not asph)
generatorSteelAnvil2 = GenerateNodeDistribution3d(nrSteelAnvilCap,
                                                  nlSteelAnvilCap,
                                                  0,
                                                  rho0Steel,
                                                  "cylindrical",
                                                  rmin = rminSteelAnvilCap,
                                                  rmax = rmaxSteelAnvilCap,
                                                  thetamin = 0,
                                                  thetamax = phi,
                                                  zmin = zminSteelAnvilCap,
                                                  zmax = zmaxSteelAnvilCap,
                                                  nNodePerh = nPerh,
                                                  SPH = not asph)
generatorSteelAnvil = CompositeNodeDistribution(generatorSteelAnvil1, generatorSteelAnvil2)
generatorFoamAnvil = GenerateNodeDistribution3d(nrFoamAnvil,
                                                nlFoamAnvil,
                                                0,
                                                rho0Foam,
                                                "cylindrical",
                                                rmin = rminFoamAnvil,
                                                rmax = rmaxFoamAnvil,
                                                thetamin = 0,
                                                thetamax = phi,
                                                zmin = zminFoamAnvil,
                                                zmax = zmaxFoamAnvil,
                                                nNodePerh = nPerh,
                                                SPH = not asph)
generatorCuAnvil = GenerateNodeDistribution3d(nrCuAnvil,
                                              nlCuAnvil,
                                              0,
                                              rho0Cu,
                                              "cylindrical",
                                              rmin = rminCuAnvil,
                                              rmax = rmaxCuAnvil,
                                              thetamin = 0,
                                              thetamax = phi,
                                              zmin = zminCuAnvil,
                                              zmax = zmaxCuAnvil,
                                              nNodePerh = nPerh,
                                              SPH = not asph)

distributeNodes3d((nodesSteel, generatorTube),
                  (nodesPlug, generatorPlug),
                  (nodesProj, generatorProj),
                  (nodesSteelAnvil, generatorSteelAnvil),
                  (nodesFoamAnvil, generatorFoamAnvil),
                  (nodesCuAnvil, generatorCuAnvil))
nGlobalNodes = 0
for n in nodeSet:
    print("Generator info for %s" % n.name)
    output("    mpi.allreduce(n.numInternalNodes, mpi.MIN)")
    output("    mpi.allreduce(n.numInternalNodes, mpi.MAX)")
    output("    mpi.allreduce(n.numInternalNodes, mpi.SUM)")
    nGlobalNodes += mpi.allreduce(n.numInternalNodes, mpi.SUM)
print("Total number of (internal) nodes in simulation: ", nGlobalNodes)

# Bevel the inner opening surface of the target tube.
numNodesBeveled = bevelTubeEntrance(nodesSteel,
                                    3,
                                    tubeOpeningAngle,
                                    rtubeInner,
                                    tubeThickness,
                                    zBevelBegin)
print("Beveled %i nodes in the tube opening." % mpi.allreduce(numNodesBeveled,
                                                              mpi.SUM))

# Adjust the diameter of the projectile inward a bit, so it will slide
# into the tube properly.
drProj = compressProjectile*nPerh*rproj/nrproj
projMultiplier = (rproj - drProj)/rproj
for i in range(nodesProj.numInternalNodes):
    nodesProj.positions()[i].x *= projMultiplier
    nodesProj.positions()[i].y *= projMultiplier

# Adjust the plug to match.
for i in range(nodesPlug.numInternalNodes):
    nodesPlug.positions()[i].x *= projMultiplier
    nodesPlug.positions()[i].y *= projMultiplier

# Iterate over the NodeLists and set some initial conditions.
for n, rho0 in [(nodesSteel, rho0Steel),
                (nodesPlug, rho0Lexan),
                (nodesProj, rho0Lexan),
                (nodesSteelAnvil, rho0Steel),
                (nodesFoamAnvil, rho0Foam),
                (nodesCuAnvil, rho0Cu)]:

    # Set node specific thermal energies
    u0 = 0.0 # n.equationOfState().specificThermalEnergy(rho0, 300.0)
    n.specificThermalEnergy(ScalarField("tmp", n, u0))
    print("Initial pressure for %s: %g" % (n.name,
                                           n.equationOfState().pressure(rho0, u0)))

# Set the projectile velocities.
nodesProj.velocity(VectorField("tmp", nodesProj, v0proj))

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node lists.
#-------------------------------------------------------------------------------
db = DataBase()
for n in nodeSet:
    db.appendNodeList(n)
output("db")
output("db.numNodeLists")
output("db.numFluidNodeLists")
output("db.numSolidNodeLists")

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
if hydroType == "CRKSPH":
    hydro = CRKSPH(dataBase = db,
                   W = WT,
                   order = correctionOrder,
                   cfl = cfl,
                   useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                   compatibleEnergyEvolution = compatibleEnergy,
                   evolveTotalEnergy = evolveTotalEnergy,
                   XSPH = XSPH,
                   densityUpdate = densityUpdate,
                   HUpdate = HUpdate,
                   ASPH = asph)

elif hydroType == "PSPH":
    hydro = PSPH(dataBase = db,
                 W = WT,
                 cfl = cfl,
                 useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                 compatibleEnergyEvolution = compatibleEnergy,
                 evolveTotalEnergy = evolveTotalEnergy,
                 correctVelocityGradient = True,
                 densityUpdate = densityUpdate,
                 HUpdate = HUpdate,
                 XSPH = XSPH,
                 ASPH = asph)

elif hydroType == "FSISPH":
    hydro = FSISPH(dataBase = db,
                   W = WT,
                   cfl = cfl,
                   useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                   compatibleEnergyEvolution = compatibleEnergy,
                   evolveTotalEnergy = evolveTotalEnergy,
                   linearCorrectGradients = True,
                   HUpdate = HUpdate,
                   ASPH = asph)

else:
    assert hydroType == "SPH"
    hydro = SPH(dataBase = db,
                W = WT,
                cfl = cfl,
                useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                compatibleEnergyEvolution = compatibleEnergy,
                evolveTotalEnergy = evolveTotalEnergy,
                gradhCorrection = gradhCorrection,
                correctVelocityGradient = True,
                densityUpdate = densityUpdate,
                HUpdate = HUpdate,
                XSPH = XSPH,
                epsTensile = epsilonTensile,
                nTensile = nTensile,
                ASPH = asph)

output("hydro")
try:
    output("hydro.kernel")
    output("hydro.PiKernel")
    output("hydro.XSPH")
except:
    pass
output("hydro.cfl")
output("hydro.compatibleEnergyEvolution")
output("hydro.evolveTotalEnergy")
output("hydro.densityUpdate")

packages = [hydro]

#-------------------------------------------------------------------------------
# Construct a damage model.
#-------------------------------------------------------------------------------
vizFields = []
if DamageModelConstructor is GradyKippTensorDamage:
    damageModel = DamageModelConstructor(nodeList = nodeSteel,
                                         kWeibull = kWeibullSteel,
                                         mWeibull = mWeibullSteel,
                                         volume = volumeSteel,
                                         kernel = WT,
                                         seed = randomSeedSteel,
                                         strainAlgorithm = strainType)

elif DamageModelConstructor is GradyKippTensorDamageOwen:
    damageModel = DamageModelConstructor(nodeList = nodesSteel,
                                         kWeibull = kWeibullSteel,
                                         mWeibull = mWeibullSteel,
                                         kernel = WT,
                                         seed = randomSeedSteel,
                                         strainAlgorithm = strainType)

elif DamageModelConstructor is JohnsonCookDamageWeibull:
    damageModel = DamageModelConstructor(nodes,
                                         D1 = D1,
                                         D2 = D2,
                                         D3 = D3,
                                         D4 = D4,
                                         D5 = D5,
                                         epsilondot0 = epsilondot0,
                                         Tcrit = Tcrit,
                                         sigmamax = sigmamax,
                                         efailmin = efailmin,
                                         aD1 = aD1,
                                         bD1 = bD1,
                                         eps0D1 = eps0D1,
                                         aD2 = aD2,
                                         bD2 = bD2,
                                         eps0D2 = eps0D2,
                                         seed = randomSeedSteel,
                                         domainIndependent = domainIndependent)
    vizFields += [damageModel.D1(), damageModel.D2()]

elif DamageModelConstructor is JohnsonCookDamageGaussian:
    damageModel = DamageModelConstructor(nodes,
                                         D1 = D1,
                                         D2 = D2,
                                         D3 = D3,
                                         D4 = D4,
                                         D5 = D5,
                                         epsilondot0 = epsilondot0,
                                         Tcrit = Tcrit,
                                         sigmamax = sigmamax,
                                         efailmin = efailmin,
                                         sigmaD1 = sigmaD1,
                                         sigmaD2 = sigmaD2,
                                         seed = randomSeed,
                                         domainIndependent = domainIndependent)

elif DamageModelConstructor is ProbabilisticDamageModel:
    damageModel = DamageModelConstructor(nodeList = nodesSteel,
                                         kernel = WT,
                                         kWeibull = kWeibullSteel,
                                         mWeibull = mWeibullSteel,
                                         seed = randomSeedSteel,
                                         strainAlgorithm = strainType)

output("damageModel")
packages.append(damageModel)

#-------------------------------------------------------------------------------
# Construct the velocity diffuser to keep the velocities under control.
#-------------------------------------------------------------------------------
if vmax:
    CHP = VelocityDiffuser(vmax)
    output("CHP")
    packages.append(CHP)

#-------------------------------------------------------------------------------
# Construct an integrator.
#-------------------------------------------------------------------------------
integrator = CheapSynchronousRK2Integrator(db, packages)
integrator.lastDt = dt
if dtMin:
    integrator.dtMin = dtMin
if dtMax:
    integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
integrator.verbose = dtverbose
output("integrator")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.dtGrowth")
output("integrator.verbose")

#-------------------------------------------------------------------------------
# Construct boundary conditions, and add them to our physics packages.
#-------------------------------------------------------------------------------
xbcPlane = Plane(Vector(0.0, 0.0, 0.0), Vector(1.0, 0.0, 0.0))
ybcPlane = Plane(Vector(0.0, 0.0, 0.0), Vector(0.0, 1.0, 0.0))
zbcPlane = Plane(Vector(0.0, 0.0, 0.0), Vector(0.0, 0.0, 1.0))
xbc = ReflectingBoundary(xbcPlane)
ybc = ReflectingBoundary(ybcPlane)
zbc = ReflectingBoundary(zbcPlane)

if phiFactor == 0.5:
    bcs = [xbc, ybc, zbc]
elif phiFactor == 1.0:
    bcs = [ybc, zbc]
else:
    bcs = [zbc]
for package in integrator.physicsPackages():
    for bc in bcs:
        package.appendBoundary(bc)

#-------------------------------------------------------------------------------
# Monitor the evolution of the mass averaged strain.
#-------------------------------------------------------------------------------
strainHistory = AverageStrain(damageModel,
                              os.path.join(dataDir, "strainhistory.txt"))

#-------------------------------------------------------------------------------
# Select the nodes for the VISAR sampling.
#-------------------------------------------------------------------------------
def xymagnitude(vec):
    return sqrt((vec.x)**2 + (vec.y)**2)

class SelectVISARNodes:
    def __init__(self, z0, dz, dr, label):
        self.z0 = z0
        self.dz = dz
        self.dr = dr
        self.label = label
        return
    def __call__(self, nodes):
        r = nodes.positions()
        potentials = [i for i in range(nodes.numInternalNodes)
                      if abs(r[i].z - self.z0) < self.dz]
        rxymax = mpi.allreduce(max([xymagnitude(r[i]) for i in potentials] + [-1e50]), mpi.MAX)
        result = [i for i in potentials if xymagnitude(r[i]) > rxymax - self.dr]
        print("Selected %i %s velocimetry test points." % (mpi.allreduce(len(result), mpi.SUM),
                                                           self.label))
        return result

#-------------------------------------------------------------------------------
# Sampling function to measure the average velocity at the VISAR probe sites.
#-------------------------------------------------------------------------------
def averageCylindricalRadialVelocity(nodes, indicies):
    r = nodes.positions()
    m = nodes.mass()
    v = nodes.velocity()

    massSum = 1e-30
    result = 0.0
    for i in indicies:
        assert i >= 0 and i < nodes.numInternalNodes
        massSum += m[i]
        runit = Vector(r[i].x, r[i].y, 0.0).unitVector()
        result += m[i] * (v[i].dot(runit))

    globalMassSum = mpi.allreduce(massSum, mpi.SUM)
    globalResult = mpi.allreduce(result, mpi.SUM)
    assert globalMassSum > 0.0
    globalResult /= globalMassSum
    return globalResult

#-------------------------------------------------------------------------------
# Build the history objects to simulate the VISAR velocity probes.
#-------------------------------------------------------------------------------
nodesA = SelectVISARNodes(zVISARa, dzVISAR, drVISAR, "A")
nodesB = SelectVISARNodes(zVISARb, dzVISAR, drVISAR, "B")
nodesC = SelectVISARNodes(zVISARc, dzVISAR, drVISAR, "C")
VISARa = NodeHistory(nodesSteel, nodesA, averageCylindricalRadialVelocity,
                     dataDir + "/VISAR-a")
VISARb = NodeHistory(nodesSteel, nodesB, averageCylindricalRadialVelocity,
                     dataDir + "/VISAR-b")
VISARc = NodeHistory(nodesSteel, nodesC, averageCylindricalRadialVelocity,
                     dataDir + "/VISAR-c")
VISARa.nodeFlags.name = "VISAR a points"
VISARb.nodeFlags.name = "VISAR b points"
VISARc.nodeFlags.name = "VISAR c points"

vizFields += [VISARa.nodeFlags,
              VISARb.nodeFlags,
              VISARc.nodeFlags]

#-------------------------------------------------------------------------------
# Build the controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            volumeType = volumeType,
                            restoreCycle = restoreCycle,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            redistributeStep = redistributeStep,
                            restartBaseName = restartBaseName,
                            vizFields = vizFields,
                            vizBaseName = vizBaseName,
                            vizDir = vizDir,
                            vizStep = vizStep,
                            vizTime = vizTime,
                            periodicWork = [(strainHistory.sample, strainFrequency),
                                            (VISARa.sample, VISARsampleFrequency),
                                            (VISARb.sample, VISARsampleFrequency),
                                            (VISARc.sample, VISARsampleFrequency)],
                            SPH = not ASPH)

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
if not steps is None:
    control.step(steps)
    raise ValueError("Completed %i steps." % steps)
else:
    control.advance(goalTime)

#-------------------------------------------------------------------------------
# Now we can collect some info on the fragment population.
#-------------------------------------------------------------------------------
from identifyFragments import *
linkRadius = 1.0
damageThreshold = 0.9
fragIDs = identifyFragments(nodesSteel,
                            linkRadius,
                            damageModel.damage(),
                            damageThreshold,
                            True)
fragProps = fragmentProperties(nodesSteel, fragIDs, damageModel.strain())
viz(fields = [fragIDs],
    filename = "ExpandingTube-3d-frags-nodust-linkRadius=%f-damageThreshold=%f" %
    (linkRadius, damageThreshold))
if mpi.rank == 0:
    f = open("frag-props-t=%8.3fusec-lr=%f-dt=%f" % (control.time()/1e-6, linkRadius, damageThreshold), "w")
    import pickle
    pickle.dump(fragProps, f)
    f.close()
