import os, shutil
from Spheral2d import *
from SpheralTestUtilities import *
from GenerateNodeDistribution2d import *

title("Simple 2 node interpolation test")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(KernelConstructor = BSplineKernel,
            rho0 = 1.0,
            eps = 0.0,
            x0 = 0.0,
            x1 = 1.0,
            nPerh = 1.25,

            gamma = 5.0/3.0,
            mu = 1.0,

            SVPH = False,
            CRKSPH = False,
            #Qconstructor = MonaghanGingoldViscosity,
            Qconstructor = TensorMonaghanGingoldViscosity,
            boolReduceViscosity = False,
            nhQ = 5.0,
            nhL = 10.0,
            aMin = 0.1,
            aMax = 2.0,
            linearConsistent = False,
            fcentroidal = 0.0,
            fcellPressure = 0.0,
            Qhmult = 1.0,
            Cl = 1.0, 
            Cq = 1.0,
            Qlimiter = False,
            epsilon2 = 1e-2,
            hmin = 0.0001, 
            hmax = 2.0,
            cfl = 0.5,
            XSPH = False,
            epsilonTensile = 0.3,
            nTensile = 4.0,
            hourglass = None,
            hourglassOrder = 0,
            hourglassLimiter = 0,
            hourglassFraction = 0.5,
            filter = 0.0,
            momentumConserving = True, # For CRKSPH
            
            vizCycle = None,
            vizTime = 0.1,

            IntegratorConstructor = CheapSynchronousRK2Integrator,
            goalTime = 0.6,
            steps = None,
            dt = 0.0001,
            dtMin = 1.0e-5, 
            dtMax = 0.1,
            dtGrowth = 2.0,
            dtverbose = False,
            rigorousBoundaries = False,
            updateBoundaryFrequency = 1,
            maxSteps = None,
            statsStep = 1,
            smoothIters = 0,
            HUpdate = IdealH,
            densityUpdate = RigorousSumDensity, # VolumeScaledDensity,
            compatibleEnergy = True,
            gradhCorrection = True,
            domainIndependent = True,
            cullGhostNodes = True,
            
            bArtificialConduction = False,
            arCondAlpha = 0.5,

            clearDirectories = True,
            checkError = False,
            checkRestart = False,
            checkEnergy = True,
            restoreCycle = None,
            restartStep = 10000,
            dataDir = "dumps-2p",
            restartBaseName = "2p",
            outputFile = None,
            comparisonFile = None,

            graphics = True,
            serialDump = True #whether to dump a serial ascii file at the end for viz
            )


vizDir = os.path.join(dataDir, "visit")
vizBaseName = "two-particles"
restartDir = os.path.join(dataDir, "restarts")
restartBaseName = os.path.join(restartDir, "2p-%i" % 1)

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
mpi.barrier()

#-------------------------------------------------------------------------------
# Material properties.
#-------------------------------------------------------------------------------
eos = GammaLawGasMKS(gamma, mu)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
WT = TableKernel(KernelConstructor(), 1000)
WTPi = TableKernel(KernelConstructor(), 1000, Qhmult)
output("WT")
output("WTPi")

#-------------------------------------------------------------------------------
# Make the NodeList.
#-------------------------------------------------------------------------------
nodes = makeFluidNodeList("nodes", eos,
                          numInternal = 2,
                          hmin = hmin,
                          hmax = hmax,
                          nPerh = nPerh)
output("nodes")
output("nodes.hmin")
output("nodes.hmax")
output("nodes.nodesPerSmoothingScale")

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------

mass    = nodes.mass()
rho     = nodes.massDensity()
H       = nodes.Hfield()
pos     = nodes.positions()

H[0].xx = 1
H[0].yy = 1
H[1].xx = 1
H[1].yy = 1
mass[0] = 1
mass[1] = 1
rho[0]  = rho0
rho[1]  = rho0
pos[0]  = Vector(x0,0.5)
pos[1]  = Vector(x1,0.5)

#-------------------------------------------------------------------------------
# DataBase
#-------------------------------------------------------------------------------
db = DataBase()
db.appendNodeList(nodes)

from siloPointmeshDump import siloPointmeshDump
siloPointmeshDump("one-node-test",
                  fields = [mass, rho, H],
                  fieldLists = [])


if serialDump:
    serialData = []
    i,j = 0,0
    
    f = open(dataDir + "/one-node.ascii",'w')
    #&ident,&x,&y,&z,&h,&mass,&rho,&temp
    for j in range(nodes.numInternalNodes):
        f.write("{0} {1} {2} {3} {4} {5} {6} {7}\n".format(j,nodes.positions()[j][0],
                                                                   nodes.positions()[j][1],
                                                                   0.0,
                                                                   1.0/(nodes.Hfield()[j].xx),
                                                                   nodes.mass()[j],
                                                                   nodes.massDensity()[j],
                                                                   0.0))
    f.close()





