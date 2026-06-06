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
            fractionPressureSupport = 0.5,
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
            restoreCycle = None,
            smoothIters = 0,
            rigorousBoundaries = True,
            dtverbose = False,
            
            serialDump = False,
            serialDumpEach = 100,

            histFile = "history.ascii",
            writeHistory = True,
            historyInterval = 2.0,
            clearDirectories = False,
            
            vizCycle = None,
            vizTime = 1.0,
            vizMethod = SpheralPointmeshSiloDump.dumpPhysicsState
            )

polytropicConstant = G0*M0/(3.0*Rc*sqrt(rho0))

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
dataDir = "owen-%i" % n
dataDir = os.path.join(dataDir, "fp=%f" % (fractionPressureSupport))
viscString = "MG"
if balsaraCorrection:
    viscString = "Balsara"
elif boolCullenViscosity:
    viscString = "Cullen"
dataDir = os.path.join(dataDir, "CRK=%s-Visc=%s-nPerh=%f-compatible=%s-volume=%s" % (CRKSPH,viscString,nPerh,compatibleEnergy,volumeType))
dataDir = os.path.join(dataDir, "Cl=%f-Cq=%f" % (Cl,Cq))
restartBaseName = "%s/KeplerianDisk-f=%f-n=%i" % (dataDir,
                                                  fractionPressureSupport,
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
# If we're restarting, find the set of most recent restart files.
#-------------------------------------------------------------------------------
if restoreCycle is None:
    restoreCycle = findLastRestart(restartBaseName)

#-------------------------------------------------------------------------------
# Define a helper class that knows how to specify our requested radial profiles
# for rho, v, and eps.
#-------------------------------------------------------------------------------
class KeplerianPressureDiskProfile:
    def __init__(self, G, M, f, n, rc,rho0):
        self.G      = G
        self.M      = M
        self.GM     = G*M
        self.f      = f
        self.gamma  = (n+1.0)/n
        self.rc     = rc
        self.rho0   = rho0
        
        self.K      = G*M/(3.0*rc*sqrt(rho0))
        
        return

    def rho(self, r):
        a = self.GM*(self.gamma-1.0)/(self.K*self.gamma*sqrt(r**2+self.rc**2))
        return pow(a,1.0/(self.gamma-1.0))

    def velocity(self, r):
        return sqrt((1.0-self.f)*self.GM*r**2/pow(r**2+self.rc**2,3.0/2.0))

    def pressure(self, r):
        return self.f*self.K*self.rho(r)**self.gamma

    def __call__(self, r):
        return self.rho(r)

#-------------------------------------------------------------------------------
# Create a polytrope for the equation of state.
#-------------------------------------------------------------------------------
eos = PolytropicEquationOfStateMKS(fractionPressureSupport*polytropicConstant,
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
diskNodes = None
if asph:
    diskNodes = AsphNodeList("diskNodes", eos, WT, WTPi)
    diskNodes.hmin = hmin
    diskNodes.hmax = hmax
    diskNodes.hminratio = hminratio
    diskNodes.nodesPerSmoothingScale = nPerh
else:
    diskNodes = makeFluidNodeList("diskNodes", eos,
                                  hmin = hmin,
                                  hmax = hmax,
                                  hminratio = hminratio,
                                  nPerh = nPerh)

output("diskNodes")
output("diskNodes.hmin")
output("diskNodes.hmax")
output("diskNodes.hminratio")
output("diskNodes.nodesPerSmoothingScale")
#output("diskNodes.epsilonTensile")
#output("diskNodes.nTensile")
#output("diskNodes.XSPH")

# Construct the neighbor object and associate it with the node list.
neighbor1 = TreeNeighbor(diskNodes,
                               kernelExtent = WT.kernelExtent)
diskNodes.registerNeighbor(neighbor1)

# Build the radial profile object that knows how to create the keplerian disk
# profile.
diskProfile = KeplerianPressureDiskProfile(G0, M0, fractionPressureSupport, polytropicIndex, Rc, rho0)

# Set node positions, masses, and H's for this domain.
if restoreCycle is None:
    from VoronoiDistributeNodes import distributeNodes2d as distributeNodes
    print("Generating node distribution.")
    generator = GenerateNodesMatchingProfile2d(n, diskProfile,
                                               rmin = rmin,
                                               rmax = rmax,
                                               thetaMin = thetaMin,
                                               thetaMax = thetaMax,
                                               nNodePerh = nPerh)
    n1 = generator.globalNumNodes()

    print("Distributing nodes amongst processors.")
    distributeNodes((diskNodes, generator))
    output('mpi.reduce(diskNodes.numInternalNodes, mpi.MIN)')
    output('mpi.reduce(diskNodes.numInternalNodes, mpi.MAX)')
    output('mpi.reduce(diskNodes.numInternalNodes, mpi.SUM)')

    # Loop over the nodes, and set the specific energies and velocities.
    for nodes in [diskNodes]:
        for i in range(nodes.numInternalNodes):
            r = nodes.positions()[i].magnitude()
            runit = nodes.positions()[i].unitVector()
            vunit = Vector(-runit.y, runit.x)
            vt = diskProfile.velocity(r)
            #nodes.specificThermalEnergy()[i] = diskProfile.eps(r)
            nodes.velocity()[i] = Vector(vt*vunit.x, vt*vunit.y)

#-------------------------------------------------------------------------------
# Set an external pressure on the disk equivalent to the pressure at the
# cutoff radius.
#-------------------------------------------------------------------------------
externalPressure = eos.polytropicConstant*diskProfile.rho(1.01*rmax)**eos.gamma_
eos.externalPressure = externalPressure

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase()
output('db')
output('db.appendNodeList(diskNodes)')
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
    dump = sDump([diskNodes],dataDir)
    control.appendPeriodicWork(dump,serialDumpEach)
output('control')

#-------------------------------------------------------------------------------
# Function to measure the angular momentum and radial coordinate.
#-------------------------------------------------------------------------------
def ParticleHistory(cycle, t, dt):
    proc = mpi.rank
    m = diskNodes.mass().internalValues()
    r = diskNodes.positions().internalValues()
    v = diskNodes.velocity().internalValues()
    for i in range(len(m)):
        if i == 100 and proc == 5:
            rm = r[i].magnitude()
            vp = (-1.0*r[i].y*v[i].x + r[i].x*v[i].y)/rm
            am = m[i]*vp
            with open(os.path.join(dataDir,histFile), "a") as myFile:
                myFile.write("%g\t%g\t%g\t%g\n" % (t,vp,am,rm))

if writeHistory:
    control.appendPeriodicTimeWork(ParticleHistory, historyInterval)
    if os.path.isfile(os.path.join(dataDir,histFile))== False:
        myFile = open(os.path.join(dataDir,histFile), "w")
        myFile.write("t\tvp\tam\trm\n")
        myFile.close()
        

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
if steps is None:
    control.advance(goalTime)
else:
    control.step(steps)

#rPlot = plotNodePositions2d(db, colorNodeLists=0, colorDomains=1)
#
#for localGoalTime in xrange(1.0, goalTime + 1.0, 1.0):
#    if control.time() < localGoalTime:
#        print "Advancing to time %f" % localGoalTime
#        control.step(5)
#        control.advance(localGoalTime)
#        control.dropRestartFile()
#
#        # Report the cumulative angular momentum error.
#        lz = [x.z for x in control.conserve.amomHistory]
#        print "Lz error: ", (max(lz) - min(lz))/lz[0]
#
#        # Plot the current state.
#        rPlot = plotNodePositions2d(db, colorNodeLists=0, colorDomains=1)
#
#        # Plot the final state.
#        rhoPlot, vrPlot, epsPlot, PPlot, HPlot = plotRadialState(db)
#        va = azimuthalVelocityFieldList(db.fluidPosition, db.fluidVelocity)
#        velPlot = plotFieldList(va,
#                                xFunction = "%s.magnitude()",
#                                plotStyle = "points",
#                                winTitle = "Velocity",
#                                lineTitle = "Simulation")
#
#        # Overplot the analytic solution.
#        import KeplerianPressureDiskSolution
#        answer = KeplerianPressureDiskSolution.KeplerianPressureDiskSolution(fractionPressureSupport,
#                                                                             polytropicConstant,
#                                                                             polytropicIndex,
#                                                                             G0,
#                                                                             M0,
#                                                                             Rc,
#                                                                             rmin = rmin,
#                                                                             rmax = rmax)
#        plotAnswer(answer, control.time(), rhoPlot, velPlot, epsPlot, PPlot)
