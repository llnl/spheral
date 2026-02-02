#-------------------------------------------------------------------------------
# The Planar Noh test case run in 1-D.
#
# W.F. Noh 1987, JCP, 72, 78-120.
#-------------------------------------------------------------------------------
#
# Ordinary SPH
#
#ATS:t0 = test(      SELF, "--graphics None --clearDirectories True  --checkError True   --restartStep 20", label="Planar Noh problem -- 1-D (serial)")
#ATS:t1 = testif(t0, SELF, "--graphics None --clearDirectories False --checkError False  --restartStep 20 --restoreCycle 20 --steps 20 --checkRestart True", label="Planar Noh problem -- 1-D (serial) RESTART CHECK")
#ATS:t2 = test(      SELF, "--graphics None --clearDirectories True  --checkError True  --dataDirBase 'dumps-planar-restartcheck' --restartStep 20", np=2, label="Planar Noh problem -- 1-D (parallel)")
#ATS:t3 = testif(t2, SELF, "--graphics None --clearDirectories False --checkError False --dataDirBase 'dumps-planar-restartcheck' --restartStep 20 --restoreCycle 20 --steps 20 --checkRestart True", np=2, label="Planar Noh problem -- 1-D (parallel) RESTART CHECK")
#ATS:t4 = test(      SELF, "--graphics None --clearDirectories True  --checkError True  --dataDirBase 'dumps-planar-reproducing' --domainIndependent True --outputFile 'Noh-planar-1proc-reproducing.txt'", label="Planar Noh problem -- 1-D (serial reproducing test setup)")
#ATS:t5 = testif(t4, SELF, "--graphics None --clearDirectories False  --checkError True  --dataDirBase 'dumps-planar-reproducing' --domainIndependent True --outputFile 'Noh-planar-4proc-reproducing.txt' --comparisonFile 'Noh-planar-1proc-reproducing.txt'", np=4, label="Planar Noh problem -- 1-D (4 proc reproducing test)")
#
# Ordinary SPH restart check for SidreFileIO
#
#ATS:t10 = test(       SELF, "--graphics None --clearDirectories True  --checkError True   --dataDir 'dumps-planar-sidre' --restartStep 20 --restartFileConstructor SidreFileIO", label="Planar Noh problem -- 1-D (serial) with Sidre")
#ATS:t11 = testif(t10, SELF, "--graphics None --clearDirectories False --checkError False  --dataDir 'dumps-planar-sidre' --restartStep 20 --restartFileConstructor SidreFileIO --restoreCycle 20 --steps 20 --checkRestart True", label="Planar Noh problem -- 1-D (serial) RESTART CHECK with Sidre")
#ATS:t12 = test(       SELF, "--graphics None --clearDirectories True  --checkError True  --dataDir 'dumps-planar-sidre-parallel' --restartStep 20 --restartFileConstructor SidreFileIO", np=2, label="Planar Noh problem -- 1-D (parallel) with Sidre")
#ATS:t13 = testif(t12, SELF, "--graphics None --clearDirectories False --checkError False --dataDir 'dumps-planar-sidre-parallel' --restartStep 20 --restartFileConstructor SidreFileIO --restoreCycle 20 --steps 20 --checkRestart True", np=2, label="Planar Noh problem -- 1-D (parallel) RESTART CHECK with Sidre")
#ATS:t14 = test(       SELF, "--graphics None --clearDirectories True  --checkError True  --dataDir 'dumps-planar-spio' --restartStep 20 --restartFileConstructor SidreFileIO --SPIOFileCountPerTimeslice 1", np=6, label="Planar Noh problem -- 1-D (parallel) with Sidre (SPIO check)")
#ATS:t15 = testif(t14, SELF, "--graphics None --clearDirectories False --checkError False --dataDir 'dumps-planar-spio' --restartStep 20 --restartFileConstructor SidreFileIO --SPIOFileCountPerTimeslice 1 --restoreCycle 20 --steps 20 --checkRestart True", np=6, label="Planar Noh problem -- 1-D (parallel) RESTART CHECK with Sidre (SPIO check)")
#
# Ordinary solid SPH
#
#ATS:t100 = test(        SELF, "--solid True --graphics None --clearDirectories True  --checkError True   --restartStep 20", label="Planar Noh problem with solid SPH -- 1-D (serial)")
#ATS:t101 = testif(t100, SELF, "--solid True --graphics None --clearDirectories False --checkError False  --restartStep 20 --restoreCycle 20 --steps 20 --checkRestart True", label="Planar Noh problem with solid SPH -- 1-D (serial) RESTART CHECK")
#ATS:t102 = test(        SELF, "--solid True --graphics None --clearDirectories True  --checkError True  --dataDirBase 'dumps-planar-restartcheck' --restartStep 20", np=2, label="Planar Noh problem with solid SPH -- 1-D (parallel)")
#ATS:t103 = testif(t102, SELF, "--solid True --graphics None --clearDirectories False --checkError False --dataDirBase 'dumps-planar-restartcheck' --restartStep 20 --restoreCycle 20 --steps 20 --checkRestart True", np=2, label="Planar Noh problem with solid SPH -- 1-D (parallel) RESTART CHECK")
#ATS:t104 = test(        SELF, "--solid True --graphics None --clearDirectories True  --checkError True  --dataDirBase 'dumps-planar-reproducing' --domainIndependent True --outputFile 'Noh-planar-1proc-reproducing.txt'", label="Planar Noh problem with solid SPH -- 1-D (serial reproducing test setup)")
#ATS:t105 = testif(t104, SELF, "--solid True --graphics None --clearDirectories False  --checkError True  --dataDirBase 'dumps-planar-reproducing' --domainIndependent True --outputFile 'Noh-planar-4proc-reproducing.txt' --comparisonFile 'Noh-planar-1proc-reproducing.txt'", np=4, label="Planar Noh  problem with solid SPH -- 1-D (4 proc reproducing test)")
#
# CRK
#
#ATS:t200 = test(        SELF, "--hydroType CRKSPH --cfl 0.25 --KernelConstructor NBSplineKernel --order 7 --nPerh 1.01 --Cl 2.0 --Cq 1.0 --graphics None --clearDirectories True --checkError True --restartStep 20", label="Planar Noh problem with CRK -- 1-D (serial)")
#ATS:t201 = testif(t200, SELF, "--hydroType CRKSPH --cfl 0.25 --KernelConstructor NBSplineKernel --order 7 --nPerh 1.01 --Cl 2.0 --Cq 1.0 --graphics None --clearDirectories False --checkError False --restartStep 20 --restoreCycle 20 --steps 20 --checkRestart True", label="Planar Noh problem with CRK -- 1-D (serial) RESTART CHECK")
#ATS:t202 = test(        SELF, "--hydroType CRKSPH --cfl 0.25 --KernelConstructor NBSplineKernel --order 7 --nPerh 1.01 --Cl 2.0 --Cq 1.0 --graphics None --clearDirectories True  --checkError False  --dataDirBase 'dumps-planar-CRK-reproducing' --domainIndependent True --outputFile 'Noh-planar-1proc-reproducing.txt' --steps 100", label="Planar Noh problem with CRK -- 1-D (serial reproducing test setup)")
#ATS:t203 = testif(t202, SELF, "--hydroType CRKSPH --cfl 0.25 --KernelConstructor NBSplineKernel --order 7 --nPerh 1.01 --Cl 2.0 --Cq 1.0 --graphics None --clearDirectories False  --checkError False  --dataDirBase 'dumps-planar-CRK-reproducing' --domainIndependent True --outputFile 'Noh-planar-4proc-reproducing.txt' --steps 100 --comparisonFile 'Noh-planar-1proc-reproducing.txt'", np=4, label="Planar Noh problem with CRK -- 1-D (4 proc reproducing test)")
#
# PSPH
#
#ATS:t300 = test(        SELF, "--hydroType PSPH --graphics None --clearDirectories True --checkError True --restartStep 20", label="Planar Noh problem with PSPH -- 1-D (serial)")
#ATS:t301 = testif(t300, SELF, "--hydroType PSPH --graphics None --clearDirectories False --checkError False --restartStep 20 --restoreCycle 20 --steps 20 --checkRestart True", label="Planar Noh problem with PSPH -- 1-D (serial) RESTART CHECK")
#
# Solid FSISPH
#
#ATS:t400 = test(        SELF, "--hydroType FSISPH --solid True --graphics None --clearDirectories True --checkError True --restartStep 20", label="Planar Noh problem with FSISPH -- 1-D (serial)", fsisph=True)
#ATS:t401 = testif(t400, SELF, "--hydroType FSISPH --solid True --graphics None --clearDirectories False --checkError False --restartStep 20 --restoreCycle 20 --steps 20 --checkRestart True", label="Planar Noh problem with FSISPH -- 1-D (serial) RESTART CHECK", fsisph=True)
#
# GSPH
#
#ATS:t500 = test(        SELF, "--hydroType GSPH --gsphReconstructionGradient RiemannGradient --graphics None --clearDirectories True --checkError True --restartStep 20", label="Planar Noh problem with GSPH and RiemannGradient -- 1-D (serial)", gsph=True)
#ATS:t501 = testif(t500, SELF, "--hydroType GSPH --gsphReconstructionGradient RiemannGradient --graphics None --clearDirectories False --checkError False --restartStep 20 --restoreCycle 20 --steps 20 --checkRestart True", label="Planar Noh problem with GSPH and RiemannGradient -- 1-D (serial) RESTART CHECK", gsph=True)
#ATS:t502 = test(        SELF, "--hydroType GSPH --gsphReconstructionGradient HydroAccelerationGradient --graphics None --clearDirectories True --checkError True --tol 5e-2 --restartStep 20", label="Planar Noh problem with GSPH and and HydroAccelerationGradient -- 1-D (serial)", gsph=True)
#ATS:t503 = testif(t502, SELF, "--hydroType GSPH --gsphReconstructionGradient HydroAccelerationGradient --graphics None --clearDirectories False --checkError False --restartStep 20 --restoreCycle 20 --steps 20 --checkRestart True", label="Planar Noh problem with GSPH and HydroAccelerationGradient -- 1-D (serial) RESTART CHECK", gsph=True)
#ATS:t504 = test(        SELF, "--hydroType GSPH --gsphReconstructionGradient SPHGradient --graphics None --clearDirectories True --checkError True --tol 0.1 --restartStep 20", label="Planar Noh problem with GSPH and SPHGradient -- 1-D (serial)", gsph=True)
#ATS:t505 = testif(t504, SELF, "--hydroType GSPH --gsphReconstructionGradient SPHGradient --graphics None --clearDirectories False --checkError False --restartStep 20 --restoreCycle 20 --steps 20 --checkRestart True", label="Planar Noh problem with GSPH and SPHGradient -- 1-D (serial) RESTART CHECK", gsph=True)
#
# MFM
#
#ATS:t600 = test(        SELF, "--hydroType MFM --gsphReconstructionGradient RiemannGradient --graphics None --clearDirectories True --checkError False --restartStep 20", label="Planar Noh problem with MFM  -- 1-D (serial)")
#ATS:t601 = testif(t600, SELF, "--hydroType MFM --gsphReconstructionGradient RiemannGradient --graphics None --clearDirectories False --checkError False --restartStep 20 --restoreCycle 20 --steps 20 --checkRestart True", label="Planar Noh problem with MFM -- 1-D (serial) RESTART CHECK")

import os, shutil, sys
import numpy as np
from SolidSpheral1d import *
from SpheralTestUtilities import *

from GenerateNodeDistribution1d import GenerateNodeDistribution1d
from SortAndDivideRedistributeNodes import distributeNodes1d

title("1-D integrated hydro test -- planar Noh problem")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(KernelConstructor = NBSplineKernel,
            order = 5,

            nx1 = 100,
            rho1 = 1.0,
            eps1 = 0.0,
            smallPressure = False, #If set to True eps is not zero but small. 
            x0 = 0.0,
            x1 = 1.0,
            xwall = 0.0,
            nPerh = 1.35,

            vr0 = -1.0, 
            vrSlope = 0.0,

            gamma = 5.0/3.0,
            mu = 1.0,

            solid = False,                     # If true, use the fluid limit of the solid hydro option
            inflow = False,                    # Should we impose inflow boundaries?

            hydroType = "SPH",                 # one of (SPH, SVPH, CRKSPH, PSPH, FSISPH, GSPH, MFM, MFV)
            crktype = "default",               # one of ("default", "variant")
            asph = False,                      # For H update algorithm, applies to all hydros
            gsphReconstructionGradient = RiemannGradient, #one of (RiemannGradient, HydroAccelerationGradient, SPHGradient, MixedGradient, OnlyDvDxGradient)
            evolveTotalEnergy = False,         # Only for SPH variants -- evolve total rather than specific energy
            boolReduceViscosity = False,
            HopkinsConductivity = False,       # For PSPH
            nhQ = 5.0,
            nhL = 10.0,
            aMin = 0.1,
            aMax = 2.0,
            boolCullenViscosity = False,
            cullenUseHydroDerivatives = True,  # Reuse the hydro calculation of DvDx.
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
            Qhmult = 1.0,
            Cl = None, 
            Cq = None,
            etaCritFrac = None,
            linearInExpansion = None,
            quadraticInExpansion = None,
            balsaraCorrection = None,
            epsilon2 = None,
            QcorrectionOrder = None,
            hmin = 0.0001, 
            hmax = 0.1,
            cfl = 0.5,
            useVelocityMagnitudeForDt = False,
            XSPH = False,
            epsilonTensile = 0.0,
            nTensile = 4.0,
            fhourglass = 0.0,

            IntegratorConstructor = CheapSynchronousRK2Integrator,
            goalTime = 0.6,
            steps = None,
            dt = 0.0001,
            dtMin = 1.0e-5, 
            dtMax = 0.1,
            dtGrowth = 2.0,
            dtverbose = False,
            updateBoundaryFrequency = 1,
            maxSteps = None,
            statsStep = 1,
            smoothIters = 0,
            HUpdate = IdealH,
            correctionOrder = LinearOrder,
            volumeType = RKSumVolume,
            densityUpdate = RigorousSumDensity, # VolumeScaledDensity,
            compatibleEnergy = True,
            gradhCorrection = True,
            correctVelocityGradient = True,
            domainIndependent = True,
            cullGhostNodes = True,
            
            bArtificialConduction = False,
            arCondAlpha = 0.5,

            clearDirectories = True,
            checkError = False,
            checkRestart = False,
            restoreCycle = None,
            restartStep = 10000,
            dataDirBase = "dumps-planar-Noh",
            restartBaseName = "Noh-planar-1d",
            restartFileConstructor = SiloFileIO,
            SPIOFileCountPerTimeslice = None,
            outputFile = None,
            comparisonFile = None,
            normOutputFile = None,
            writeOutputLabel = True,
            doCompare = True,

            # Parameters for the test acceptance.,
            tol = 1.0e-5,

            graphics = True,
            )

hydroType = hydroType.upper()

assert not(boolReduceViscosity and boolCullenViscosity)
assert not(hydroType == "GSPH" and (boolReduceViscosity or boolCullenViscosity))
assert not(hydroType == "FSISPH" and not solid)
if smallPressure:
    P0 = 1.0e-6
    eps1 = P0/((gamma - 1.0)*rho1)
   
# Build a path spec that varies a bit based on the hydro choice
hydroPath = hydroType
if hydroType == "CRKSPH":
    hydroPath = os.path.join(hydroPath,
                             str(volumeType),
                             str(correctionOrder))
elif hydroType in ("GSPH", "MFM", "MFV"):
    hydroPath = os.path.join(hydroPath, str(gsphReconstructionGradient))

if solid:
    hydroPath = "Solid" + hydroPath

if dataDirBase:
    dataDir = os.path.join(dataDirBase,
                           hydroPath,
                           "nPerh=%f" % nPerh,
                           "compatibleEnergy=%s" % compatibleEnergy,
                           "fhourglass=%s" % fhourglass,
                           "Cullen=%s" % boolCullenViscosity)
    restartDir = os.path.join(dataDir, "restarts")
    restartBaseName = os.path.join(restartDir, "Noh-planar-1d-%i" % nx1)
else:
    dataDir = None
    restartDir = None
    restartBaseName = None

dx = (x1 - x0)/nx1

#-------------------------------------------------------------------------------
# The reference values for error norms checking for pass/fail
#-------------------------------------------------------------------------------
LnormRef = {"SPH": {"Mass density" : {"L1"   : 0.05376370586899846,   
                                      "L2"   : 0.01472935554844709,   
                                      "Linf" : 1.6558627223391338},    
                    "Pressure    " : {"L1"   : 0.018082144742236444,  
                                      "L2"   : 0.005431965024057943,  
                                      "Linf" : 0.6289446614721329},   
                    "Velocity    " : {"L1"   : 0.024463871274705278,  
                                      "L2"   : 0.008419536302504558,  
                                      "Linf" : 0.8561295316236415},   
                    "Spec Therm E" : {"L1"   : 0.010557215425476638,  
                                      "L2"   : 0.0033659949510588386, 
                                      "Linf" : 0.355220540682649},    
                    "h           " : {"L1"   : 0.00043625606815746957,
                                      "L2"   : 0.00012010712699702793,
                                      "Linf" : 0.008480811209733824}},
            "CRKSPH": {"Mass density" : {"L1"   : 0.05064113393844768,   
                                         "L2"   : 0.015297215762507312,  
                                         "Linf" : 1.6768360873659973},    
                       "Pressure    " : {"L1"   : 0.01687133903296828,   
                                         "L2"   : 0.006328924998534429,  
                                         "Linf" : 0.7822574604543725},   
                       "Velocity    " : {"L1"   : 0.007746512026971996,  
                                         "L2"   : 0.0029862099280521903, 
                                         "Linf" : 0.20321897008372736},  
                       "Spec Therm E" : {"L1"   : 0.005051748111924938,  
                                         "L2"   : 0.0015094950940911932, 
                                         "Linf" : 0.14418618728403587},  
                       "h           " : {"L1"   : 0.00019175169182527455,
                                         "L2"   : 6.850786936014129e-05, 
                                         "Linf" : 0.004376337346566557}},
            "FSISPH": {"Mass density" : {"L1"   : 0.080317978098225,     
                                         "L2"   : 0.01731304525219977,   
                                         "Linf" : 1.803838026495527},     
                       "Pressure    " : {"L1"   : 0.020655347086429156,  
                                         "L2"   : 0.005363359619262443,  
                                         "Linf" : 0.6139872258598699},   
                       "Velocity    " : {"L1"   : 0.0260229353392815,    
                                         "L2"   : 0.008514799936559351,  
                                         "Linf" : 0.8539687338327078},   
                       "Spec Therm E" : {"L1"   : 0.012595057574721798,  
                                         "L2"   : 0.003152955099116865,  
                                         "Linf" : 0.32847148447477337},  
                       "h           " : {"L1"   : 0.0004623051617065188, 
                                         "L2"   : 0.00012256543877716504,
                                         "Linf" : 0.008641461274901919}},
            "PSPH": {"Mass density" : {"L1"   : 0.06067866550282133,   
                                       "L2"   : 0.015430245737443435,  
                                       "Linf" : 1.707010689252927},     
                     "Pressure    " : {"L1"   : 0.022914718927884124,  
                                       "L2"   : 0.005975425367223863,  
                                       "Linf" : 0.6676159215788076},   
                     "Velocity    " : {"L1"   : 0.02585232906638296,   
                                       "L2"   : 0.008736748109104572,  
                                       "Linf" : 0.8666199732590047},   
                     "Spec Therm E" : {"L1"   : 0.011826327206410471,  
                                       "L2"   : 0.0036116821761209244, 
                                       "Linf" : 0.36993717326050946},  
                     "h           " : {"L1"   : 0.00044462158158787294,
                                       "L2"   : 0.00011990796335122118,
                                       "Linf" : 0.00843114543207368}}, 
            "GSPH": {"Mass density" : {"L1"   : 0.04835000721781902,   
                                       "L2"   : 0.014738218758202467,  
                                       "Linf" : 1.680911284343678},     
                     "Pressure    " : {"L1"   : 0.020488554417202706,  
                                       "L2"   : 0.006258439801775443,  
                                       "Linf" : 0.7290690753591823},   
                     "Velocity    " : {"L1"   : 0.022636006592029485,  
                                       "L2"   : 0.00780524881260205,   
                                       "Linf" : 0.8751592871887031},   
                     "Spec Therm E" : {"L1"   : 0.012494783387052708,  
                                       "L2"   : 0.0040445381100901465, 
                                       "Linf" : 0.4079223432067062},   
                     "h           " : {"L1"   : 0.0004272365926225631, 
                                       "L2"   : 0.00012051038774533798,
                                       "Linf" : 0.00840917801325954}}, 
            "MFM": {"Mass density" : {"L1"   : 0.0873630138456682,    
                                      "L2"   : 0.02097262837445441,   
                                      "Linf" : 2.259098555266673},     
                    "Pressure    " : {"L1"   : 0.029892808080209176,  
                                      "L2"   : 0.007283817451547816,  
                                      "Linf" : 0.8859532047529353},   
                    "Velocity    " : {"L1"   : 0.03650211122120961,   
                                      "L2"   : 0.009996952619228788,  
                                      "Linf" : 0.9499807356516067},   
                    "Spec Therm E" : {"L1"   : 0.015673619503887928,  
                                      "L2"   : 0.004023662041147362,  
                                      "Linf" : 0.40775993256842896},  
                    "h           " : {"L1"   : 0.0005398542934066652, 
                                      "L2"   : 0.0001315019397588108, 
                                      "Linf" : 0.009141770037360247}},
}

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
if mpi.rank == 0:
    if clearDirectories and os.path.exists(dataDir):
        shutil.rmtree(dataDir)
    if not os.path.exists(restartDir):
        os.makedirs(restartDir)
mpi.barrier()

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
WT = TableKernel(Wbase, 1000)
kernelExtent = WT.kernelExtent
output("WT")

#-------------------------------------------------------------------------------
# Make the NodeList.
#-------------------------------------------------------------------------------
if solid:
    nodes1 = makeSolidNodeList("nodes1", eos, 
                               hmin = hmin,
                               hmax = hmax,
                               nPerh = nPerh,
                               kernelExtent = kernelExtent)
else:
    nodes1 = makeFluidNodeList("nodes1", eos, 
                               hmin = hmin,
                               hmax = hmax,
                               nPerh = nPerh,
                               kernelExtent = kernelExtent)
    
output("nodes1")
output("nodes1.hmin")
output("nodes1.hmax")
output("nodes1.nodesPerSmoothingScale")

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
gen = GenerateNodeDistribution1d(n = nx1,
                                 rho = rho1,
                                 xmin = x0,
                                 xmax = x1,
                                 nNodePerh = nPerh)
distributeNodes1d((nodes1, gen))
output("nodes1.numNodes")

# Set node specific thermal energies
nodes1.specificThermalEnergy(ScalarField("tmp", nodes1, eps1))
nodes1.massDensity(ScalarField("tmp", nodes1, rho1))

# Set node velocities
pos = nodes1.positions()
vel = nodes1.velocity()
for ix in range(nodes1.numNodes):
    if pos[ix].x > xwall:
        vel[ix].x = vr0 + vrSlope*pos[ix].x
    else:
        vel[ix].x = -vr0 + vrSlope*pos[ix].x

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase()
output("db")
output("db.appendNodeList(nodes1)")
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
if hydroType == "SVPH":
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
                 xmin = Vector(-100.0),
                 xmax = Vector( 100.0))
elif hydroType == "CRKSPH":
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
                   crktype = crktype)
elif hydroType == "PSPH":
    hydro = PSPH(dataBase = db,
                 W = WT,
                 cfl = cfl,
                 useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                 compatibleEnergyEvolution = compatibleEnergy,
                 evolveTotalEnergy = evolveTotalEnergy,
                 correctVelocityGradient = correctVelocityGradient,
                 densityUpdate = densityUpdate,
                 HUpdate = HUpdate,
                 XSPH = XSPH)

elif hydroType == "FSISPH":
    hydro = FSISPH(dataBase = db,
                   W = WT,
                   cfl = cfl,
                   interfaceMethod = HLLCInterface,
                   sumDensityNodeLists=[nodes1],                       
                   densityStabilizationCoefficient = 0.00,
                   useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                   compatibleEnergyEvolution = compatibleEnergy,
                   evolveTotalEnergy = evolveTotalEnergy,
                   linearCorrectGradients = correctVelocityGradient,
                   HUpdate = HUpdate)
elif hydroType == "GSPH":
    limiter = VanLeerLimiter()
    waveSpeed = DavisWaveSpeed()
    solver = HLLC(limiter,
                  waveSpeed,
                  True)
    hydro = GSPH(dataBase = db,
                riemannSolver = solver,
                W = WT,
                cfl=cfl,
                useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                compatibleEnergyEvolution = compatibleEnergy,
                correctVelocityGradient=correctVelocityGradient,
                evolveTotalEnergy = evolveTotalEnergy,
                XSPH = XSPH,
                gradientType = gsphReconstructionGradient,
                densityUpdate=densityUpdate,
                HUpdate = IdealH,
                epsTensile = epsilonTensile,
                nTensile = nTensile)
elif hydroType == "MFM":
    limiter = VanLeerLimiter()
    waveSpeed = DavisWaveSpeed()
    solver = HLLC(limiter,
                  waveSpeed,
                  True)
    hydro = MFM(dataBase = db,
                riemannSolver = solver,
                W = WT,
                cfl=cfl,
                useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                compatibleEnergyEvolution = compatibleEnergy,
                correctVelocityGradient=correctVelocityGradient,
                evolveTotalEnergy = evolveTotalEnergy,
                XSPH = XSPH,
                gradientType = gsphReconstructionGradient,
                densityUpdate=densityUpdate,
                HUpdate = IdealH,
                epsTensile = epsilonTensile,
                nTensile = nTensile)
elif hydroType == "MFV":
    limiter = VanLeerLimiter()
    waveSpeed = DavisWaveSpeed()
    solver = HLLC(limiter,
                  waveSpeed,
                  True)
    hydro = MFV(dataBase = db,
                riemannSolver = solver,
                W = WT,
                cfl=cfl,
                nodeMotionType=NodeMotionType.Lagrangian,
                useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                compatibleEnergyEvolution = compatibleEnergy,
                correctVelocityGradient=correctVelocityGradient,
                evolveTotalEnergy = evolveTotalEnergy,
                XSPH = XSPH,
                gradientType = gsphReconstructionGradient,
                densityUpdate=densityUpdate,
                HUpdate = IdealH,
                epsTensile = epsilonTensile,
                nTensile = nTensile)
else:
    assert hydroType == "SPH"
    hydro = SPH(dataBase = db,
                W = WT,
                cfl = cfl,
                useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                compatibleEnergyEvolution = compatibleEnergy,
                evolveTotalEnergy = evolveTotalEnergy,
                gradhCorrection = gradhCorrection,
                XSPH = XSPH,
                correctVelocityGradient = correctVelocityGradient,
                densityUpdate = densityUpdate,
                epsTensile = epsilonTensile,
                nTensile = nTensile,
                HUpdate = HUpdate,
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
output("hydro.densityUpdate")

packages = [hydro]

#-------------------------------------------------------------------------------
# Set the artificial viscosity parameters.
#-------------------------------------------------------------------------------
if not hydroType in ("GSPH", "MFM", "MFV"):
    q = hydro.Q
    if not Cl is None:
        q.Cl = Cl
    if not Cq is None:
        q.Cq = Cq
    if not epsilon2 is None:
        q.epsilon2 = epsilon2
    if not balsaraCorrection is None:
        q.balsaraShearCorrection = balsaraCorrection
    if not QcorrectionOrder is None:
        q.QcorrectionOrder = QcorrectionOrder
    output("q")
    output("q.Cl")
    output("q.Cq")
    output("q.epsilon2")
    output("q.balsaraShearCorrection")
    output("q.requireVelocityGradient()")
    output("q.rigorousVelocityGradient")
    if hasattr(q, "linearInExpansion") and not linearInExpansion is None:
        q.linearInExpansion = linearInExpansion
        output("q.linearInExpansion")
    if hasattr(q, "quadraticInExpansion") and not quadraticInExpansion is None:
        q.quadraticInExpansion = quadraticInExpansion
        output("q.quadraticInExpansion")
    if hasattr(q, "etaCritFrac") and not etaCritFrac is None:
        q.etaCritFrac = etaCritFrac
        output("q.etaCritFrac")

#-------------------------------------------------------------------------------
# Construct the MMRV physics object.
#-------------------------------------------------------------------------------
if boolReduceViscosity:
    evolveReducingViscosityMultiplier = MorrisMonaghanReducingViscosity(nhQ,nhL,aMin,aMax)
    packages.append(evolveReducingViscosityMultiplier)
elif boolCullenViscosity:
    evolveCullenViscosityMultiplier = CullenDehnenViscosity(WT,alphMax,alphMin,betaC,betaD,betaE,fKern,boolHopkinsCorrection,cullenUseHydroDerivatives)
    packages.append(evolveCullenViscosityMultiplier)

#-------------------------------------------------------------------------------
# Construct the Artificial Conduction physics object.
#-------------------------------------------------------------------------------
if bArtificialConduction:
    #q.reducingViscosityCorrection = True
    ArtyCond = ArtificialConduction(WT,arCondAlpha)
    
    packages.append(ArtyCond)

#-------------------------------------------------------------------------------
# Optionally construct an hourglass control object.
#-------------------------------------------------------------------------------
if fhourglass > 0.0:
    hg = SubPointPressureHourglassControl(fhourglass)
    output("hg")
    output("hg.fHG")
    packages.append(hg)

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
bcs = []
if x0 == xwall:
    xPlane0 = Plane(Vector(0.0), Vector(1.0))
    xbc0 = ReflectingBoundary(xPlane0)
    bcs = [xbc0]

if inflow:
    right_inflow = InflowOutflowBoundary(db, Plane(Vector(x1), Vector(-1)))
    bcs.append(right_inflow)
    packages.append(right_inflow)
    if x0 != xwall:
        left_inflow = InflowOutflowBoundary(db, Plane(Vector(x0), Vector(1)))
        bcs.append(left_inflow)
        packages.append(left_inflow)

for p in packages:
    for bc in bcs:
        p.appendBoundary(bc)

#-------------------------------------------------------------------------------
# Construct an integrator.
#-------------------------------------------------------------------------------
integrator = IntegratorConstructor(db)
for p in packages:
    integrator.appendPhysicsPackage(p)
del p
integrator.lastDt = dt
integrator.dtMin = dtMin
integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
integrator.updateBoundaryFrequency = updateBoundaryFrequency
integrator.domainDecompositionIndependent = domainIndependent
integrator.cullGhostNodes = cullGhostNodes
integrator.verbose = dtverbose
output("integrator")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.dtGrowth")
output("integrator.updateBoundaryFrequency")
output("integrator.domainDecompositionIndependent")
output("integrator.cullGhostNodes")
output("integrator.verbose")

#-------------------------------------------------------------------------------
# Make the problem controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator,
                            kernel = WT,
                            volumeType = volumeType,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            restartFileConstructor = restartFileConstructor,
                            SPIOFileCountPerTimeslice = SPIOFileCountPerTimeslice,
                            restoreCycle = restoreCycle,
                            SPH = not asph,
                            )
output("control")

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
if not steps is None:
    if checkRestart:
        control.setRestartBaseName(restartBaseName + "_CHECK")
    control.step(steps)
    if checkRestart:
        control.setRestartBaseName(restartBaseName)

    # Are we doing the restart test?
    if checkRestart:
        state0 = State(db, integrator.physicsPackages())
        state0.copyState()
        print(control.totalSteps)
        control.loadRestartFile(control.totalSteps)
        state1 = State(db, integrator.physicsPackages())
        if not state1 == state0:
            raise ValueError("The restarted state does not match!")
        else:
            print("Restart check PASSED.")

else:
    if control.time() < goalTime:
        control.step(5)
        control.advance(goalTime, maxSteps)

# If running the performance test, stop here
if not doCompare:
    sys.exit(0)

#-------------------------------------------------------------------------------
# Compute the analytic answer.
#-------------------------------------------------------------------------------
import mpi
import NohAnalyticSolution
rlocal = [pos.x for pos in nodes1.positions().internalValues()]
r = mpi.reduce(rlocal, mpi.SUM)
h1 = 1.0/(nPerh*dx)
answer = NohAnalyticSolution.NohSolution(1,
                                         r = r,
                                         v0 = -1.0,
                                         h0 = 1.0/h1)

# Compute the simulated specific entropy.
rho = mpi.allreduce(nodes1.massDensity().internalValues(), mpi.SUM)
Pf = ScalarField("pressure", nodes1)
nodes1.pressure(Pf)
P = mpi.allreduce(Pf.internalValues(), mpi.SUM)
A = [Pi/rhoi**gamma for (Pi, rhoi) in zip(P, rho)]

# The analytic solution for the simulated entropy.
xprof = mpi.allreduce([x.x for x in nodes1.positions().internalValues()], mpi.SUM)
xans, vans, uans, rhoans, Pans, hans = answer.solution(control.time(), xprof)
Aans = [Pi/rhoi**gamma for (Pi, rhoi) in zip(Pans,  rhoans)]
L1 = 0.0
for i in range(len(rho)):
    L1 = L1 + abs(rho[i]-rhoans[i])
L1_tot = L1 / len(rho)
if mpi.rank == 0 and outputFile:
    print("L1=",L1_tot,"\n")
    with open("Converge.txt", "a") as myfile:
        myfile.write("%s %s\n" % (nx1, L1_tot))

#-------------------------------------------------------------------------------
# Plot the final state.
#-------------------------------------------------------------------------------
if graphics:
    from SpheralMatplotlib import *
    rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotState(db)
    plotAnswer(answer, control.time(), rhoPlot, velPlot, epsPlot, PPlot, HPlot = HPlot)
    EPlot = plotEHistory(control.conserve)
    plots = [(rhoPlot, "Noh-planar-rho.png"),
             (velPlot, "Noh-planar-vel.png"),
             (epsPlot, "Noh-planar-eps.png"),
             (PPlot, "Noh-planar-P.png"),
             (HPlot, "Noh-planar-h.png")]

    # Plot the specific entropy.
    Aplot = newFigure()
    Aplot.plot(xprof, A, "ro", label="Simulation")
    Aplot.plot(xprof, Aans, "b-", label="Solution")
    Aplot.set_title("Specific entropy")
    plots.append((Aplot, "Noh-planar-A.png"))
    
    if hydroType == "CRKSPH":
        volPlot = plotFieldList(control.RKCorrections.volume, 
                                winTitle = "volume",
                                colorNodeLists = False, plotGhosts = False)
        splot = plotFieldList(control.RKCorrections.surfacePoint,
                              winTitle = "surface point",
                              colorNodeLists = False)
        plots += [(volPlot, "Noh-planar-vol.png"),
                  (splot, "Noh-planar-surfacePoint.png")]

    if boolCullenViscosity:
        cullAlphaPlot = plotFieldList(q.ClMultiplier(),
                                      winTitle = "Cullen alpha")
        cullDalphaPlot = plotFieldList(evolveCullenViscosityMultiplier.DalphaDt(),
                                       winTitle = "Cullen DalphaDt")
        plots += [(cullAlphaPlot, "Noh-planar-Cullen-alpha.png"),
                  (cullDalphaPlot, "Noh-planar-Cullen-DalphaDt.png")]

    if boolReduceViscosity:
        alphaPlotQ = plotFieldList(q.reducingViscosityMultiplierQ(),
                                  winTitle = "rvAlphaQ",
                                  colorNodeLists = False, plotGhosts = False)
        alphaPlotL = plotFieldList(q.reducingViscosityMultiplierL(),
                                   winTitle = "rvAlphaL",
                                   colorNodeLists = False, plotGhosts = False)

    # # Plot the grad h correction term (omega)
    # omegaPlot = plotFieldList(hydro.omegaGradh(),
    #                           winTitle = "grad h correction",
    #                           colorNodeLists = False)

    # Make hardcopies of the plots.
    for p, filename in plots:
        p.figure.savefig(os.path.join(dataDir, filename))

#-------------------------------------------------------------------------------
# Measure the difference between the simulation and analytic answer.
#-------------------------------------------------------------------------------
rmin, rmax = 0.05, 0.35   # Throw away anything with r < rwall to avoid wall heating.
rhoprof = mpi.reduce(nodes1.massDensity().internalValues(), mpi.SUM)
P = ScalarField("pressure", nodes1)
nodes1.pressure(P)
Pprof = mpi.reduce(P.internalValues(), mpi.SUM)
vprof = mpi.reduce([v.x for v in nodes1.velocity().internalValues()], mpi.SUM)
epsprof = mpi.reduce(nodes1.specificThermalEnergy().internalValues(), mpi.SUM)
hprof = mpi.reduce([1.0/H.xx for H in nodes1.Hfield().internalValues()], mpi.SUM)
xprof = mpi.reduce([x.x for x in nodes1.positions().internalValues()], mpi.SUM)

#-------------------------------------------------------------------------------
# If requested, write out the state in a global ordering to a file.
#-------------------------------------------------------------------------------
if outputFile:
    outputFile = os.path.join(dataDir, outputFile)
    from SpheralTestUtilities import multiSort
    mof = mortonOrderIndices(db)
    mo = mpi.reduce(mof[0].internalValues(), mpi.SUM)
    mprof = mpi.reduce(nodes1.mass().internalValues(), mpi.SUM)
    rhoprof = mpi.reduce(nodes1.massDensity().internalValues(), mpi.SUM)
    P = ScalarField("pressure", nodes1)
    nodes1.pressure(P)
    Pprof = mpi.reduce(P.internalValues(), mpi.SUM)
    vprof = mpi.reduce([v.x for v in nodes1.velocity().internalValues()], mpi.SUM)
    epsprof = mpi.reduce(nodes1.specificThermalEnergy().internalValues(), mpi.SUM)
    hprof = mpi.reduce([1.0/H.xx for H in nodes1.Hfield().internalValues()], mpi.SUM)
    if mpi.rank == 0:
        multiSort(mo, xprof, rhoprof, Pprof, vprof, epsprof, hprof)
        f = open(outputFile, "w")
        f.write(("#  " + 20*"'%s' " + "\n") % ("x", "m", "rho", "P", "v", "eps", "h", "mo",
                                               "rhoans", "Pans", "vans", "epsans", "hans",
                                               "x_UU", "m_UU", "rho_UU", "P_UU", "v_UU", "eps_UU", "h_UU"))
        for (xi, mi, rhoi, Pi, vi, epsi, hi, moi,
             rhoansi, Pansi, vansi, uansi, hansi) in zip(xprof, mprof, rhoprof, Pprof, vprof, epsprof, hprof, mo,
                                                         rhoans, Pans, vans, uans, hans):
            f.write((7*"%16.12e " + "%i " + 5*"%16.12e " + 7*"%i " + '\n') % 
                    (xi, mi, rhoi, Pi, vi, epsi, hi, moi,
                     rhoansi, Pansi, vansi, uansi, hansi,
                     unpackElementUL(packElementDouble(xi)),
                     unpackElementUL(packElementDouble(mi)),
                     unpackElementUL(packElementDouble(rhoi)),
                     unpackElementUL(packElementDouble(Pi)),
                     unpackElementUL(packElementDouble(vi)),
                     unpackElementUL(packElementDouble(epsi)),
                     unpackElementUL(packElementDouble(hi))))
        f.close()

        #---------------------------------------------------------------------------
        # Also we can optionally compare the current results with another file.
        #---------------------------------------------------------------------------
        if comparisonFile:
            comparisonFile = os.path.join(dataDir, comparisonFile)
            import filecmp
            assert filecmp.cmp(outputFile, comparisonFile)

#------------------------------------------------------------------------------
# Compute the error.
#------------------------------------------------------------------------------
failure = False
if mpi.rank == 0 :
    xans, vans, epsans, rhoans, Pans, hans = answer.solution(control.time(), xprof)
    import Pnorm
    print("Quantity \t\tL1 \t\t\t\tL2 \t\t\t\tLinf")
    failure = False

    if normOutputFile:
       f = open(normOutputFile, "a")
       if writeOutputLabel:
          f.write(("#" + 13*"%17s " + "\n") % ('"nx"',
                                               '"rho L1"', '"rho L2"', '"rho Linf"',
                                               '"P L1"',   '"P L2"',   '"P Linf"',
                                               '"vel L1"', '"vel L2"', '"vel Linf"',
                                               '"E L1"', '"E L2"', '"E Linf"',
                                               '"h L1"',   '"h L2"',   '"h Linf"'))
       f.write("%5i " % nx1)
    for (name, data, ans) in [("Mass density", rhoprof, rhoans),
                              ("Pressure    ", Pprof, Pans),
                              ("Velocity    ", vprof, vans),
                              ("Spec Therm E", epsprof, epsans),
                              ("h           ", hprof, hans)]:
        assert len(data) == len(ans)
        error = [data[i] - ans[i] for i in range(len(data))]
        Pn = Pnorm.Pnorm(error, xprof)
        L1 = Pn.gridpnorm(1, rmin, rmax)
        L2 = Pn.gridpnorm(2, rmin, rmax)
        Linf = Pn.gridpnorm("inf", rmin, rmax)
        print(f"{name}\t\t{L1} \t\t{L2} \t\t{Linf}")
        if normOutputFile:
           f.write((3*"%16.12e ") % (L1, L2, Linf))

        if checkError and not (np.allclose(L1, LnormRef[hydroType][name]["L1"], tol, tol) and
                               np.allclose(L2, LnormRef[hydroType][name]["L2"], tol, tol) and
                               np.allclose(Linf, LnormRef[hydroType][name]["Linf"], tol, tol)):
            print("Failing Lnorm tolerance for ", name, (L1, L2, Linf), LnormRef[hydroType][name])
            failure = True
  
    if normOutputFile:
       f.write("\n")

    if checkError and failure:
        raise ValueError("Error bounds violated.")

Eerror = (control.conserve.EHistory[-1] - control.conserve.EHistory[0])/control.conserve.EHistory[0]
print("Total energy error: %g" % Eerror)
if compatibleEnergy and abs(Eerror) > 1e-13:
    raise ValueError("Energy error outside allowed bounds.")

# Check that SPIO is writing the expected amount of files also need to check if mpi is enabled to see if we are using Spio
if (control.restartFileConstructor is SidreFileIO) and (mpi.rank == 0) and (not mpi.is_fake_mpi()) and (control.SPIOFileCountPerTimeslice is not None):
    if not control.SPIOFileCountPerTimeslice is len(os.listdir(os.path.join(os.getcwd(), control.restartBaseName + "_cycle%i" % control.totalSteps))):
        raise ValueError("The amount of restart files written does not match the amount expected based on input!")
