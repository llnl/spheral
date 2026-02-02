#ATS:test(SELF, "--CRKSPH=False --n=50 --cfl=0.25 --Cl=1.0 --Cq=2.0 --filter=0 --nPerh=1.01 --balsaraCorrection=True --fractionPressureSupport=0.25 --serialDump=True --compatibleEnergy=False --goalTime=50", label="Kepler SPH balsara no-compat, nPerh=1.5 fp=0.05", np=20)
#ATS:test(SELF, "--CRKSPH=False --n=50 --cfl=0.25 --Cl=1.0 --Cq=2.0 --filter=0 --nPerh=1.01 --balsaraCorrection=True --fractionPressureSupport=0.25 --serialDump=True --compatibleEnergy=True --goalTime=50", label="Kepler SPH balsara w-compat, nPerh=1.5 fp=0.05", np=20)
#ATS:test(SELF, "--CRKSPH=False --n=50 --cfl=0.25 --Cl=1.0 --Cq=2.0 --filter=0 --nPerh=1.01 --balsaraCorrection=False --fractionPressureSupport=0.25 --serialDump=True --compatibleEnergy=False --goalTime=50", label="Kepler SPH no-compat, nPerh=1.5 fp=0.05", np=20)
#ATS:test(SELF, "--CRKSPH=False --n=50 --cfl=0.25 --Cl=1.0 --Cq=2.0 --filter=0 --nPerh=1.01 --balsaraCorrection=False --fractionPressureSupport=0.25 --serialDump=True --compatibleEnergy=True --goalTime=50", label="Kepler SPH w-compat, nPerh=1.5 fp=0.05", np=20)


#-------------------------------------------------------------------------------
# This test problem sets up a gas disk in a fixed potential from a softened
# point mass.  The fractionPressureSupport parameter selects the ratio of
# pressure and rotational support in the disk.  If all is working properly,
# the disk should be stable as initialized and should just rotate without any
# radial evolution.
#-------------------------------------------------------------------------------
from Spheral2d import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from findLastRestart import *
from math import *
import SpheralPointmeshSiloDump

# Load the mpi module if we're parallel.
import mpi
#mpi, rank, procs = mpi.loadmpi()

from GenerateNodeDistribution2d import *

title("2-D Keplerian disk with arbitrary pressure support.")

# serialDump thing for external viz
class sDump(object):
    def __init__(self,nodeSet,directory):
        self.nodeSet = nodeSet
        self.directory = directory
    def __call__(self, cycle, time, dt):
        procs = mpi.procs
        rank = mpi.rank
        serialData = []
        i,j = 0,0
        for i in range(procs):
            for nodeL in self.nodeSet:
                if rank == i:
                    for j in range(nodeL.numInternalNodes):
                        serialData.append([nodeL.positions()[j],
                                           3.0/(nodeL.Hfield()[j].Trace()),
                                           nodeL.mass()[j],nodeL.massDensity()[j],
                                           nodeL.specificThermalEnergy()[j]])

        serialData = mpi.reduce(serialData,mpi.SUM)
        if rank == 0:
            f = open(self.directory + "/serialDump" + str(cycle) + ".ascii",'w')
            for i in range(len(serialData)):
                f.write("{0} {1} {2} {3} {4} {5} {6} {7}\n".format(i,serialData[i][0][0],serialData[i][0][1],0.0,serialData[i][1],serialData[i][2],serialData[i][3],serialData[i][4]))
            f.close()

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(asph = False,

            n = 100,
            thetaMin = 0.0,
            thetaMax = 2.0*pi,
            rmin = 0.0,
            rmax = 3.0,
            nPerh = 1.51,

            # Properties of the central gravitating particle.
            G0 = 1.0,
            M0 = 1.0,
            Rc = 0.5,
            R0 = Vector(0.0, 0.0),

            # Properties of the gas disk.
            rho0  = 1.0,
            rd0   = 10.0,
            sig   = 2.5,
            Rcutoff = 0.5,

            # Material properties of the gas.
            polytropicIndex = 2.0,
            mu = 1.0,

            SVPH = False,
            CRKSPH = False,
            ASPH = False,
            SPH = True,   # This just chooses the H algorithm -- you can use this with CRKSPH for instance.
            
            XSPH = False,
            
            epsilonTensile = 0.0,
            nTensile = 8,
            
            # Hydro
            Qconstructor = MonaghanGingoldViscosity2d,
            #Qconstructor = TensorMonaghanGingoldViscosity2d,
            KernelConstructor = NBSplineKernel,
            order = 5,
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
            Cq = 0.75,
            Qlimiter = False,
            balsaraCorrection = False,
            epsilon2 = 1e-4,
            negligibleSoundSpeed = 1e-5,
            csMultiplier = 0.1,
            correctionOrder = LinearOrder,
            hmin = 0.004,
            hmax = 0.5,
            hminratio = 0.1,
            compatibleEnergy = True,
            gradhCorrection = False,
            HEvolution = IdealH,
            sumForMassDensity = RigorousSumDensity,
            densityUpdate = RigorousSumDensity, # VolumeScaledDensity,
            HUpdate = IdealH,
            filter = 0.0,
            volumeType = RKSumVolume,

            # Timestep constraints
            cfl = 0.5,
            deltaPhi = 0.01,
            domainIndependent = False,

            # Integrator and run time.
            IntegratorConstructor = CheapSynchronousRK2Integrator,
            steps = None,
            goalTime = 10.0,
            dt = 0.0001,
            dtMin = 1.0e-5,
            dtMax = 0.1,
            dtGrowth = 2.0,
            maxSteps = None,
            statsStep = 10,
            redistributeStep = None,
            restartStep = 500,
            restoreCycle = -1,
            smoothIters = 0,
            rigorousBoundaries = True,
            dtverbose = False,
            
            serialDump = False,
            serialDumpEach = 100,

            histFile = "history.ascii",
            writeHistory = False,
            historyInterval = 2.0,
            clearDirectories = False,

            dataDir = "twomat-%i",

            outputFile = None,
            comparisonFile = None,
            
            vizCycle = None,
            vizTime = 1.0,
            vizMethod = SpheralPointmeshSiloDump.dumpPhysicsState
            )

polytropicConstant1 = G0*M0/(3.0*Rc*sqrt(rho0))
polytropicConstant2 = G0*M0/(3.0*Rc*sqrt(rho0*0.5))

# Decide on our hydro algorithm.
if SVPH:
    if ASPH:
        HydroConstructor = ASVPHFacetedHydro
    else:
        HydroConstructor = SVPHFacetedHydro
elif CRKSPH:
    if ASPH:
        HydroConstructor = ACRKSPHHydro
    else:
        HydroConstructor = CRKSPHHydro
    Qconstructor = LimitedMonaghanGingoldViscosity
else:
    if ASPH:
        HydroConstructor = ASPHHydro
    else:
        HydroConstructor = SPHHydro

# Data output info.
dataDir = dataDir % n
viscString = "MG"
if balsaraCorrection:
    viscString = "Balsara"
elif boolCullenViscosity:
    viscString = "Cullen"
dataDir = os.path.join(dataDir, "CRK=%s-Visc=%s-nPerh=%f-compatible=%s-volume=%s" % (CRKSPH,viscString,nPerh,compatibleEnergy,volumeType))
dataDir = os.path.join(dataDir, "Cl=%f-Cq=%f" % (Cl,Cq))
restartBaseName = "%s/KeplerianDisk-n=%i" % (dataDir,
                                                  n)

vizDir = os.path.join(dataDir, "visit")
vizBaseName = "Kepler-disk-2d"

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
import os, sys
if mpi.rank == 0:
    if not os.path.exists(dataDir):
        os.makedirs(dataDir)
    if not os.path.exists(vizDir):
        os.makedirs(vizDir)
mpi.barrier()

#-------------------------------------------------------------------------------
# Define a helper class that knows how to specify our requested radial profiles
# for rho, v, and eps.
#-------------------------------------------------------------------------------
class KeplerianPressureDiskProfile:
    def __init__(self,G,M,n,rc,rho0):
        self.G = G
        self.M = M
        self.GM = G*M
        self.gamma = (n+1.0)/n
        self.rc = rc
        self.rho0 = rho0
        
        self.K = G*M/(3.0*rc*sqrt(rho0))
        return

    def rho(self,r):
        a = self.GM*(self.gamma-1.0)/(self.K*self.gamma*sqrt(r**2+self.rc**2))
        return pow(a,1.0/(self.gamma-1.0))

    def pressure(self,r):
        return self.K*self.rho(r)**self.gamma
    
    def __call__(self,r):
        return self.rho(r)


#-------------------------------------------------------------------------------
# Create a polytrope for the equation of state.
#-------------------------------------------------------------------------------
eos1 = PolytropicEquationOfStateMKS(polytropicConstant1,
                                     polytropicIndex, mu)

eos2 = PolytropicEquationOfStateMKS(polytropicConstant2,
                                    polytropicIndex, mu)


#-------------------------------------------------------------------------------
# Create our interpolation kernels -- one for normal hydro interactions, and
# one for use with the artificial viscosity
#-------------------------------------------------------------------------------
WT = TableKernel(KernelConstructor(order), 1000)
WTPi = TableKernel(KernelConstructor(order), 1000)
output('WT')
output('WTPi')

#-------------------------------------------------------------------------------
# Create the NodeList and distribute it's nodes.
#-------------------------------------------------------------------------------

diskNodes1 = makeFluidNodeList("diskNodes1", eos1,
                               hmin = hmin,
                               hmax = hmax,
                               hminratio = hminratio,
                               nPerh = nPerh)

diskNodes2 = makeFluidNodeList("diskNodes2", eos2,
                               hmin = hmin,
                               hmax = hmax,
                               hminratio = hminratio,
                               nPerh = nPerh)

output("diskNodes1")
output("diskNodes1.hmin")
output("diskNodes1.hmax")
output("diskNodes1.hminratio")
output("diskNodes1.nodesPerSmoothingScale")
#output("diskNodes.epsilonTensile")
#output("diskNodes.nTensile")
#output("diskNodes.XSPH")

# Construct the neighbor object and associate it with the node list.
#neighbor1 = TreeNeighbor(diskNodes1,
#                               kernelExtent = WT.kernelExtent)
#diskNodes1.registerNeighbor(neighbor1)
#diskNodes2.registerNeighbor(neighbor2)

# Build the radial profile object that knows how to create the keplerian disk
# profile.
diskProfile1 = KeplerianPressureDiskProfile(G0, M0, polytropicIndex, Rc, rho0)
diskProfile2 = KeplerianPressureDiskProfile(G0, M0, polytropicIndex, Rc, rho0*0.5)

# Set node positions, masses, and H's for this domain.
from VoronoiDistributeNodes import distributeNodes2d as distributeNodes
print("Generating node distribution.")
generator1 = GenerateNodesMatchingProfile2d(n*0.25, diskProfile1,
                                            rmin = rmin,
                                            rmax = rmax*0.25,
                                            thetaMin = thetaMin,
                                            thetaMax = thetaMax,
                                            nNodePerh = nPerh)
n1 = generator1.globalNumNodes()
generator2 = GenerateNodesMatchingProfile2d(n*0.75, diskProfile2,
                                            rmin = rmax*0.27,
                                            rmax = rmax,
                                            thetaMin = thetaMin,
                                            thetaMax = thetaMax,
                                            nNodePerh = nPerh,
                                            m0 = generator1.m0)
n1 = generator1.globalNumNodes()
n2 = generator2.globalNumNodes()

print("Distributing nodes amongst processors.")
distributeNodes((diskNodes1, generator1),(diskNodes2,generator2))
output('mpi.reduce(diskNodes1.numInternalNodes, mpi.MIN)')
output('mpi.reduce(diskNodes1.numInternalNodes, mpi.MAX)')
output('mpi.reduce(diskNodes1.numInternalNodes, mpi.SUM)')

# Loop over the nodes, and set the specific energies and velocities.
for nodes in [diskNodes1,diskNodes2]:
    for i in range(nodes.numInternalNodes):
        r = nodes.positions()[i].magnitude()
        #nodes.specificThermalEnergy()[i] = diskProfile.eps(r)

#-------------------------------------------------------------------------------
# Set an external pressure on the disk equivalent to the pressure at the
# cutoff radius.
#-------------------------------------------------------------------------------
externalPressure = eos2.polytropicConstant*diskProfile2.rho(1.01*rmax)**eos2.gamma_
eos2.externalPressure = externalPressure

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase()
output('db')
output('db.appendNodeList(diskNodes1)')
output('db.appendNodeList(diskNodes2)')
output('db.numNodeLists')
output('db.numFluidNodeLists')

#-------------------------------------------------------------------------------
# Construct the artificial viscosity.
#-------------------------------------------------------------------------------
q = Qconstructor(Cl, Cq)
q.limiter = Qlimiter
q.balsaraShearCorrection = balsaraCorrection
q.epsilon2 = epsilon2
q.negligibleSoundSpeed = negligibleSoundSpeed
q.csMultiplier = csMultiplier
output('q')
output('q.Cl')
output('q.Cq')
output('q.limiter')
output('q.epsilon2')
output('q.negligibleSoundSpeed')
output('q.csMultiplier')
output('q.balsaraShearCorrection')

#-------------------------------------------------------------------------------
# Create the gravity physics object.
#-------------------------------------------------------------------------------
gravity = PointPotential(G0, M0, Rc, R0)
gravity.deltaPotentialFraction = deltaPhi
output("gravity.G")
output("gravity.mass")
output("gravity.coreRadius")
output("gravity.origin")
output("gravity.deltaPotentialFraction")

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
if SVPH:
    hydro = HydroConstructor(W = WT, 
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
                             xmin = Vector(-2.0, -2.0),
                             xmax = Vector(3.0, 3.0))
# xmin = Vector(x0 - 0.5*(x2 - x0), y0 - 0.5*(y2 - y0)),
# xmax = Vector(x2 + 0.5*(x2 - x0), y2 + 0.5*(y2 - y0)))
elif CRKSPH:
    hydro = HydroConstructor(W = WT,
                             WPi = WTPi, 
                             Q = q,
                             filter = filter,
                             cfl = cfl,
                             compatibleEnergyEvolution = compatibleEnergy,
                             XSPH = XSPH,
                             densityUpdate = densityUpdate,
                             correctionOrder = correctionOrder,
                             volumeType = volumeType,
                             HUpdate = HUpdate)
else:
    hydro = HydroConstructor(W = WT,
                             WPi = WTPi,
                             Q = q,
                             cfl = cfl,
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
# Construct a time integrator, and add the physics packages.
#-------------------------------------------------------------------------------
integrator = IntegratorConstructor(db)
for p in packages:
    integrator.appendPhysicsPackage(gravity)
    integrator.appendPhysicsPackage(p)
integrator.lastDt = dt
integrator.dtMin = dtMin
integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
integrator.domainDecompositionIndependent = domainIndependent
integrator.verbose = dtverbose
integrator.rigorousBoundaries = rigorousBoundaries

# Blago!  Currently a problem with periodic boundaries.
integrator.cullGhostNodes = False

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
# Build the controller to run the simulation.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            redistributeStep = redistributeStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            vizMethod = vizMethod,
                            vizBaseName = vizBaseName,
                            vizDir = vizDir,
                            vizStep = vizCycle,
                            vizTime = vizTime,
                            vizDerivs = True,
                            restoreCycle = restoreCycle)


if serialDump:
    dump = sDump([diskNodes1,diskNodes2],dataDir)
    control.appendPeriodicWork(dump,serialDumpEach)
output('control')

#-------------------------------------------------------------------------------
# Function to measure the angular momentum and radial coordinate.
#-------------------------------------------------------------------------------

        

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
if steps is None:
    control.advance(goalTime)
else:
    control.step(steps)

if outputFile:
    outputFile = os.path.join(dataDir, outputFile)
    from SpheralTestUtilities import multiSort
    P1 = ScalarField("pressure",diskNodes1)
    P2 = ScalarField("pressure",diskNodes2)
    diskNodes1.pressure(P1)
    diskNodes2.pressure(P2)

    xprof1 = mpi.reduce([x.x for x in diskNodes1.positions().internalValues()], mpi.SUM)
    yprof1 = mpi.reduce([y.y for y in diskNodes1.positions().internalValues()], mpi.SUM)
    rhoprof1 = mpi.reduce(diskNodes1.massDensity().internalValues(), mpi.SUM)
    Pprof1 = mpi.reduce(P1.internalValues(), mpi.SUM)
    rprof1 = mpi.reduce([ri.magnitude() for ri in diskNodes1.positions().internalValues()], mpi.SUM)
    vx1 = mpi.reduce([v.x for v in diskNodes1.velocity().internalValues()], mpi.SUM)
    vy1 = mpi.reduce([v.y for v in diskNodes1.velocity().internalValues()], mpi.SUM)

    xprof2 = mpi.reduce([x.x for x in diskNodes2.positions().internalValues()], mpi.SUM)
    yprof2 = mpi.reduce([y.y for y in diskNodes2.positions().internalValues()], mpi.SUM)
    rhoprof2 = mpi.reduce(diskNodes2.massDensity().internalValues(), mpi.SUM)
    Pprof2 = mpi.reduce(P2.internalValues(), mpi.SUM)
    rprof2 = mpi.reduce([ri.magnitude() for ri in diskNodes2.positions().internalValues()], mpi.SUM)
    vx2 = mpi.reduce([v.x for v in diskNodes2.velocity().internalValues()], mpi.SUM)
    vy2 = mpi.reduce([v.y for v in diskNodes2.velocity().internalValues()], mpi.SUM)

    np1 = int(diskNodes1.numInternalNodes)
    np2 = int(diskNodes2.numInternalNodes)
    if np1 is None:
        np1 = 0
    np1 = mpi.reduce(np1,mpi.SUM)
    if np2 is None:
        np2 = 0
    np2 = mpi.reduce(np2,mpi.SUM)

    vprof1 = []
    vprof2 = []
    if mpi.rank == 0:
        for i in range(np1):
            vprof1.append(xprof1[i]*vx1[i]/rprof1[i]+yprof1[i]*vy1[i]/rprof1[i])
        for i in range(np2):
            vprof2.append(xprof2[i]*vx2[i]/rprof2[i]+yprof2[i]*vy2[i]/rprof2[i])

    mof = mortonOrderIndices(db)
    mo1 = mpi.reduce(mof[0].internalValues(),mpi.SUM)
    mo2 = mpi.reduce(mof[1].internalValues(),mpi.SUM)
    if mpi.rank == 0:
        multiSort(rprof1,mo1,xprof1,yprof1,rhoprof1,Pprof1,vprof1)
        multiSort(rprof2,mo2,xprof2,yprof2,rhoprof2,Pprof2,vprof2)
        f = open(outputFile, "w")
        f.write("r x y rho P v mortonOrder\n")
        for (ri, xi, yi, rhoi, Pi, vi, mi) in zip(rprof1,xprof1,yprof1,rhoprof1,Pprof1,vprof1,mo1):
            f.write((7*"%16.12e "+"\n") % (ri,xi,yi,rhoi,Pi,vi,mi))
        for (ri, xi, yi, rhoi, Pi, vi, mi) in zip(rprof2,xprof2,yprof2,rhoprof2,Pprof2,vprof2,mo2):
            f.write((7*"%16.12e "+"\n") % (ri,xi,yi,rhoi,Pi,vi,mi))
        
        f.close()
        if comparisonFile:
            comparisonFile = os.path.join(dataDir, comparisonFile)
            import filecmp
            assert filecmp.cmp(outputFile,comparisonFile)
