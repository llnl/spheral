#ATS:t0 = test(SELF, "--clearDirectories True --domainIndependent True --outputFile 'PlateImpact-1proc-reproducing.txt'", np=1, label="5 material flyer plate impact test -- 1-D (1 processor reproducing)")
#ATS:t1 = testif(t0, SELF, "--clearDirectories False --domainIndependent True --comparisonFile 'PlateImpact-1proc-reproducing.txt'--outputFile 'PlateImpact-2proc-reproducing.txt'", np=2, label="5 material flyer plate impact test -- 1-D (2 processor reproducing)")
#-------------------------------------------------------------------------------
# A five material flyer plate test, based on an experiment.
#-------------------------------------------------------------------------------
from SolidSpheral1d import *
from SpheralTestUtilities import *
from findLastRestart import *
from SpheralVisitDump import SpheralVisitDump
from SpheralGnuPlotUtilities import *
from math import *
import shutil
import mpi

#-------------------------------------------------------------------------------
# Generic problem parameters
# All CGS units.
#-------------------------------------------------------------------------------
commandLine(
    seed = 'lattice',

    Sapphire1Thickness = 1.15,
    TantalumThickness = 0.50225,
    Sapphire2Thickness = 0.3171,
    TungstenCarbideThickness = 0.16,
    PMMAThickness = 0.616,

    nxSapphire1 = 181,
    nxTantalum = 251,
    nxSapphire2 = 50,
    nxTungstenCarbide = 60,
    nxPMMA = 30,

    nPerh = 1.25,

    rhoSapphire = 3.985,
    rhoTantalum = 16.69,
    rhoTungstenCarbide = 14.9,
    rhoPMMA = 1.182,

    etaminSapphire = 0.1,
    etamaxSapphire = 10.0,
    etaminTantalum = 0.5,
    etamaxTantalum = 1.5,
    etaminTungstenCarbide = 0.5,
    etamaxTungstenCarbide = 1.5,
    etaminPMMA = 0.5,
    etamaxPMMA = 1.5,

    v1 = Vector(-1.33e4),

    HydroConstructor = SolidSPHHydro,
    Qconstructor = MonaghanGingoldViscosity,
    Cl = 0.5,
    Cq = 1.0,
    Qlimiter = False,
    balsaraCorrection = False,
    epsilon2 = 1e-2,
    negligibleSoundSpeed = 1e-5,
    csMultiplier = 1e-4,
    hmin = 1e-10,
    hmax = 0.5,
    cfl = 0.5,
    useVelocityMagnitudeForDt = True,
    XSPH = False,
    epsilonTensile = 0.3,
    nTensile = 4,

    goalTime = 5.0e-6,
    steps = None,
    dtSample = 1e-7,
    dt = 1e-10,
    dtMin = 1e-12,
    dtMax = 1e-5,
    dtGrowth = 2.0,
    maxSteps = None,
    statsStep = 10,
    smoothIters = 0,
    HEvolution = IdealH,
    sumForMassDensity = IntegrateDensity,
    compatibleEnergy = True,
    gradhCorrection = True,

    restartStep = 1000,
    restartDir = "dumps",
    readRestartFile = False,
    restoreCycle = None,
    clearDirectories = False,

    # Parameters for the test acceptance.
    vtol = 1.0e-3,

    # Should we run in domain independent mode, and if so should we check
    # for domain independence?
    domainIndependent = False,
    outputFile = None,
    comparisonFile = None,
    )

Sapphire1Range = (0.0, Sapphire1Thickness)
TantalumRange = (Sapphire1Range[1], Sapphire1Range[1] + TantalumThickness)
Sapphire2Range = (TantalumRange[1], TantalumRange[1] + Sapphire2Thickness)
TungstenCarbideRange = (Sapphire2Range[1], Sapphire2Range[1] + TungstenCarbideThickness)
PMMARange = (TungstenCarbideRange[1], TungstenCarbideRange[1] + PMMAThickness)

xmin = 0.0
xmax = PMMARange[1]

restartBaseName = restartDir + "/PlateImpact-%i-%i-%i-%i-%i" % (nxSapphire1,
                                                                nxTantalum,
                                                                nxSapphire2,
                                                                nxTungstenCarbide,
                                                                nxPMMA)

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
import os, sys
if mpi.rank == 0:
    if clearDirectories and os.path.exists(restartDir):
        shutil.rmtree(restartDir)
    if not os.path.exists(restartDir):
        os.makedirs(restartDir)
mpi.barrier()

#-------------------------------------------------------------------------------
# If we're restarting, find the set of most recent restart files.
#-------------------------------------------------------------------------------
if readRestartFile and (restoreCycle is None):
    restoreCycle = findLastRestart(restartBaseName)

#-------------------------------------------------------------------------------
# Define a class to track the history of the interface.
#-------------------------------------------------------------------------------
class InterfaceHistory:

    def __init__(self,
                 sapphireIndex,
                 tantalumIndex,
                 sapphireProc,
                 tantalumProc,
                 sapphireNodes,
                 tantalumNodes,
                 filename):
        self.restart = RestartableObject(self)
        self.sapphireIndex = sapphireIndex
        self.tantalumIndex = tantalumIndex
        self.sapphireProc = sapphireProc
        self.tantalumProc = tantalumProc
        self.sapphireNodes = sapphireNodes
        self.tantalumNodes = tantalumNodes
        self.filename = filename
        self.file = None
        self.timeHistory = []
        self.sapphireVelocity = []
        self.tantalumVelocity = []
        return

    def sample(self, t):
        sapphireVelocity0 = 1e40
        tantalumVelocity0 = 1e40
        if mpi.rank == self.sapphireProc:
            sapphireVelocity0 = self.sapphireNodes.velocity()[self.sapphireIndex].x
        if mpi.rank == self.tantalumProc:
            tantalumVelocity0 = self.tantalumNodes.velocity()[self.tantalumIndex].x
        sapphireVelocity = mpi.allreduce(sapphireVelocity0, mpi.MIN)
        tantalumVelocity = mpi.allreduce(tantalumVelocity0, mpi.MIN)
        self.timeHistory.append(t)
        self.sapphireVelocity.append(sapphireVelocity)
        self.tantalumVelocity.append(tantalumVelocity)
        if mpi.rank == 0:
            if self.file is None:
                self.file = open(self.filename, "w")
            self.file.write("%g \t%g \t%g\n" % (t, sapphireVelocity, tantalumVelocity))
            self.file.flush()
        return

    def label(self):
        return "InterfaceHistory"

    def dumpState(self, file, path):
        file.writeObject(self.sapphireIndex, path + "/sapphireIndex")
        file.writeObject(self.tantalumIndex, path + "/tantalumIndex")
        file.writeObject(self.sapphireProc, path + "/sapphireProc")
        file.writeObject(self.tantalumProc, path + "/tantalumProc")
        file.writeObject(self.filename, path + "/filename")
        file.writeObject(self.timeHistory, path + "/timeHistory")
        file.writeObject(self.sapphireVelocity, path + "/sapphireVelocity")
        file.writeObject(self.tantalumVelocity, path + "/tantalumVelocity")
        return

    def restoreState(self, file, path):
        self.sapphireIndex = file.readObject(path + "/sapphireIndex")
        self.tantalumIndex = file.readObject(path + "/tantalumIndex")
        self.sapphireProc = file.readObject(path + "/sapphireProc")
        self.tantalumProc = file.readObject(path + "/tantalumProc")
        self.filename = file.readObject(path + "/filename")
        self.timeHistory = file.readObject(path + "/timeHistory")
        self.sapphireVelocity = file.readObject(path + "/sapphireVelocity")
        self.tantalumVelocity = file.readObject(path + "/tantalumVelocity")
        return

#-------------------------------------------------------------------------------
# Identify ourselves!
#-------------------------------------------------------------------------------
title("1-D Plate impact strength test")
print("Sapphire 1 in range [%f, %f]" % Sapphire1Range)
print("Tantalum   in range [%f, %f]" % TantalumRange)
print("Sapphire 2 in range [%f, %f]" % Sapphire2Range)
print("TungCarbid in range [%f, %f]" % TungstenCarbideRange)
print("PMMA       in range [%f, %f]" % PMMARange)

#-------------------------------------------------------------------------------
# Sapphire material properties.
#-------------------------------------------------------------------------------
eosSapphire  = GruneisenEquationOfStateCGS(rhoSapphire,    # reference density  
                                           etaminSapphire, # etamin             
                                           etamaxSapphire, # etamax             
                                           1.119e6,        # C0
                                           1.0,            # S1                 
                                           0.0,            # S2                 
                                           0.0,            # S3                 
                                           0.0,            # gamma0             
                                           0.0,            # b                  
                                           20.39)          # atomic weight      

#-------------------------------------------------------------------------------
# Tantalum material properties.
#-------------------------------------------------------------------------------
eosTantalum = GruneisenEquationOfStateCGS(rhoTantalum,    # reference density  
                                          etaminTantalum, # etamin             
                                          etamaxTantalum, # etamax             
                                          0.341e6, # C0                 
                                          1.2,     # S1                 
                                          0.0,     # S2                 
                                          0.0,     # S3                 
                                          1.67,    # gamma0             
                                          0.42,    # b                  
                                          180.948) # atomic weight
coldFitTantalum = NinthOrderPolynomialFit(-6.86446714e9,
                                          -1.17070812e10,
                                           9.70252276e11,
                                          -4.12557402e11,
                                           1.13401823e11,
                                          -1.86584799e10,
                                           0.0,
                                           0.0,
                                           0.0,
                                           0.0)
meltFitTantalum = NinthOrderPolynomialFit(9.24414908e10,
                                          2.53949977e11,
                                          1.06113848e12,
                                         -5.28947636e11,
                                          1.67906438e11,
                                         -2.92459765e10,
                                          0.0,
                                          0.0,
                                          0.0,
                                          0.0)
strengthTantalum = SteinbergGuinanStrengthCGS(eosTantalum,
                                              6.900000e11,        # G0
                                              1.4500e-12,         # A
                                              1.3000e-04,         # B
                                              7.7000e9,           # Y0
                                              1.1e10,             # Ymax
                                              1.0e-3,             # Yp
                                              10.0000,            # beta
                                              0.0,                # gamma0
                                              0.1,                # nhard
                                              coldFitTantalum,
                                              meltFitTantalum)

# Suspending this until I believe the Steinberg-Guinan-Lund strength model
# is working correctly.
# strengthTantalum = SteinbergGuinanLundStrengthCGS(eosTantalum,
#                                                     6.900000e11,        # G0
#                                                     1.4500e-12,         # A
#                                                     1.3000e-04,         # B
#                                                     7.7000e9,           # Y0
#                                                     1.1e10,             # Ymax
#                                                     1.0e-3,             # Yp
#                                                     10.0000,            # beta
#                                                     0.0,                # gamma0
#                                                     0.1,                # nhard
#                                                     7.10e-7,            # C1
#                                                     1.20e17,            # C2
#                                                     7.195e3,            # UK
#                                                     8.2e9,              # YP
#                                                     4.5e9,              # YTmax
#                                                     coldFitTantalum,
#                                                     meltFitTantalum)

#-------------------------------------------------------------------------------
# Tungsten carbide material properties.
#-------------------------------------------------------------------------------
eosTungstenCarbide = GruneisenEquationOfStateCGS(rhoTungstenCarbide, # reference density  
                                                 etaminTungstenCarbide, 
                                                 etamaxTungstenCarbide, 
                                                 0.519e6,# C0                 
                                                 1.16,   # S1                 
                                                 0.0,    # S2                 
                                                 0.0,    # S3                 
                                                 1.5,    # gamma0             
                                                 0.0,    # b                  
                                                 155.89) # atomic weight      
coldFitTungstenCarbide = NinthOrderPolynomialFit(0.0,
                                                 0.0,
                                                 0.0,
                                                 0.0,
                                                 0.0,
                                                 0.0,
                                                 0.0,
                                                 0.0,
                                                 0.0,
                                                 0.0)
meltFitTungstenCarbide = NinthOrderPolynomialFit(0.0,
                                                 0.0,
                                                 0.0,
                                                 0.0,
                                                 0.0,
                                                 0.0,
                                                 0.0,
                                                 0.0,
                                                 0.0,
                                                 0.0)
strengthTungstenCarbide = SteinbergGuinanStrengthCGS(eosTungstenCarbide,
                                                     2.540000e12,        # G0
                                                     0.0,                # A
                                                     0.0,                # B
                                                     5.1000e10,          # Y0
                                                     1.0e12,             # Ymax
                                                     1.0e-3,             # Yp
                                                     0.0,                # beta
                                                     0.0,                # gamma0
                                                     0.0,                # nhard
                                                     coldFitTungstenCarbide,
                                                     meltFitTungstenCarbide)

#-------------------------------------------------------------------------------
# PMMA material properties.
#-------------------------------------------------------------------------------
eosPMMA = GruneisenEquationOfStateCGS(rhoPMMA,    # reference density  
                                      etaminPMMA, # etamin             
                                      etamaxPMMA, # etamax             
                                      0.218e6,# C0                 
                                      2.088,  # S1                 
                                     -1.124, # S2                 
                                      0.0,    # S3                 
                                      0.85,   # gamma0             
                                      0.0,    # b                  
                                      17.035) # atomic weight      
coldFitPMMA = NinthOrderPolynomialFit(-5.19191852e9,
                                      -4.41500192e9,
                                       2.84720528e10,
                                       2.14093899e10,
                                      -4.46412259e9,
                                       1.24495222e9,
                                       0.0,
                                       0.0,
                                       0.0,
                                       0.0)
meltFitPMMA = NinthOrderPolynomialFit(5.24383771e8,
                                      1.49188457e9,
                                      2.85704428e10,
                                      2.13783662e10,
                                     -4.45135120e9,
                                      1.24138074e9,
                                      0.0,
                                      0.0,
                                      0.0,
                                      0.0)
strengthPMMA = SteinbergGuinanStrengthCGS(eosPMMA,
                                          2.320000e10,        # G0
                                          0.0,                # A
                                          0.0,                # B
                                          4.2000e9,           # Y0
                                          1.0e12,             # Ymax
                                          1.0e-3,             # Yp
                                          0.0,                # beta
                                          0.0,                # gamma0
                                          0.0,                # nhard
                                          coldFitPMMA,
                                          meltFitPMMA)

#-------------------------------------------------------------------------------
# Create our interpolation kernels -- one for normal hydro interactions, and
# one for use with the artificial viscosity
#-------------------------------------------------------------------------------
WT = TableKernel(BSplineKernel(), 1000)
output('WT')

#-------------------------------------------------------------------------------
# Create the NodeLists.
#-------------------------------------------------------------------------------
nodesSapphire1 = makeSolidNodeList("Sapphire 1", eosSapphire,
                                   nPerh = nPerh,
                                   hmin = hmin,
                                   hmax = hmax,
                                   rhoMin = etaminSapphire*rhoSapphire,
                                   rhoMax = etamaxSapphire*rhoSapphire)
nodesTantalum = makeSolidNodeList("Tantalum", eosTantalum, strengthTantalum,
                                   nPerh = nPerh,
                                   hmin = hmin,
                                   hmax = hmax,
                                   rhoMin = etaminTantalum*rhoTantalum,
                                   rhoMax = etamaxTantalum*rhoTantalum)
nodesSapphire2 = makeSolidNodeList("Sapphire 2", eosSapphire,
                                   nPerh = nPerh,
                                   hmin = hmin,
                                   hmax = hmax,
                                   rhoMin = etaminSapphire*rhoSapphire,
                                   rhoMax = etamaxSapphire*rhoSapphire)
nodesTungstenCarbide = makeSolidNodeList("Tungsten Carbide", eosTungstenCarbide, strengthTungstenCarbide,
                                   nPerh = nPerh,
                                   hmin = hmin,
                                   hmax = hmax,
                                   rhoMin = etaminTungstenCarbide*rhoTungstenCarbide,
                                   rhoMax = etamaxTungstenCarbide*rhoTungstenCarbide)
nodesPMMA = makeSolidNodeList("PMMA", eosPMMA, strengthPMMA,
                              nPerh = nPerh,
                              hmin = hmin,
                              hmax = hmax,
                              rhoMin = etaminPMMA*rhoPMMA,
                              rhoMax = etamaxPMMA*rhoPMMA)

nodeLists = [nodesSapphire1, nodesTantalum, nodesSapphire2, nodesTungstenCarbide, nodesPMMA]

#-------------------------------------------------------------------------------
# Create an instance of our history object to find the interface history.
#-------------------------------------------------------------------------------
interfaceHistory = InterfaceHistory(None, None, None, None,
                                    nodesSapphire1, nodesTantalum,
                                    "PlateImpact-interface-history.txt")

#-------------------------------------------------------------------------------
# Set node properties (positions, masses, H's, etc.)
#-------------------------------------------------------------------------------
if restoreCycle is None:
    from DistributeNodes import distributeNodesInRange1d
    print("Generating node distribution.")
    distributeNodesInRange1d([(nodesSapphire1, nxSapphire1, rhoSapphire, Sapphire1Range),
                              (nodesTantalum, nxTantalum, rhoTantalum, TantalumRange),
                              (nodesSapphire2, nxSapphire2, rhoSapphire, Sapphire2Range),
                              (nodesTungstenCarbide, nxTungstenCarbide, rhoTungstenCarbide, TungstenCarbideRange),
                              (nodesPMMA, nxPMMA, rhoPMMA, PMMARange)])

    # Explicitly set the mass densities again, 'cause there's some roundoff
    # issue python->C++ with doing it from the NodeGenerators.  This is
    # necessary to get the initial pressures to be zero to roundoff.
    for nodes, rho0 in [(nodesSapphire1, rhoSapphire),
                        (nodesTantalum, rhoTantalum),
                        (nodesSapphire2,  rhoSapphire),
                        (nodesTungstenCarbide, rhoTungstenCarbide),
                        (nodesPMMA, rhoPMMA)]:
        nodes.massDensity(ScalarField("tmp", nodes, rho0))

    # Set the node velocities.
    nodesSapphire2.velocity(VectorField("tmp", nodesSapphire2, v1))
    nodesTungstenCarbide.velocity(VectorField("tmp", nodesTungstenCarbide, v1))
    nodesPMMA.velocity(VectorField("tmp", nodesPMMA, v1))

    # Find the initial thermal energy to give us a zero pressure.
    for nodes in nodeLists:
        Pf = ScalarField("pressure", nodes)
        nodes.pressure(Pf)
        P = list(Pf.internalValues())
        Pmin = mpi.allreduce(min(P + [1e100]), mpi.MIN)
        Pmax = mpi.allreduce(max(P + [-1e100]), mpi.MAX)
        print("Initial pressures for %s : [%g : %g]" % (nodes.name, Pmin, Pmax))
    del nodes

    # Find the interface nodes between the Sapphire 2 and Tungsten.
    sapphirelist = list(zip([pos.x for pos in nodesSapphire1.positions().internalValues()],
                       list(range(nodesSapphire1.numInternalNodes)),
                       [mpi.rank]*nodesSapphire1.numInternalNodes)) + [(-1e20, None, None)]
    tantalumlist = list(zip([pos.x for pos in nodesTantalum.positions().internalValues()],
                        list(range(nodesTantalum.numInternalNodes)),
                        [mpi.rank]*nodesTantalum.numInternalNodes)) + [(-1e20, None, None)]
    maxsapphire = mpi.allreduce(max(sapphirelist), mpi.MAX)
    maxtantalum = mpi.allreduce(max(tantalumlist), mpi.MAX)
    assert maxsapphire[2] is not None
    assert maxtantalum[2] is not None
    interfaceHistory.sapphireIndex = maxsapphire[1]
    interfaceHistory.sapphireProc = maxsapphire[2]
    interfaceHistory.tantalumIndex = maxtantalum[1]
    interfaceHistory.tantalumProc = maxtantalum[2]

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase()
for nodes in nodeLists:
    db.appendNodeList(nodes)
del nodes
output("db")
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Construct the artificial viscosities for the problem.
#-------------------------------------------------------------------------------
q = Qconstructor(Cl, Cq)
q.limiter = Qlimiter
q.balsaraShearCorrection = balsaraCorrection
q.epsilon2 = epsilon2
q.negligibleSoundSpeed = negligibleSoundSpeed
q.csMultiplier = csMultiplier
output("q")
output("q.Cl")
output("q.Cq")
output("q.limiter")
output("q.epsilon2")
output("q.negligibleSoundSpeed")
output("q.csMultiplier")
output("q.balsaraShearCorrection")

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
hydro = HydroConstructor(W = WT,
                         Q = q,
                         cfl = cfl,
                         useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                         compatibleEnergyEvolution = compatibleEnergy,
                         gradhCorrection = gradhCorrection,
                         densityUpdate = sumForMassDensity,
                         HUpdate = HEvolution,
                         epsTensile = epsilonTensile,
                         nTensile = nTensile)
output("hydro")
output("hydro.kernel()")
output("hydro.PiKernel()")
output("hydro.cfl")
output("hydro.useVelocityMagnitudeForDt")
output("hydro.compatibleEnergyEvolution")
output("hydro.gradhCorrection")
output("hydro.XSPH")
output("hydro.sumForMassDensity")
output("hydro.HEvolution")
output("hydro.epsilonTensile")
output("hydro.nTensile")

#-------------------------------------------------------------------------------
# Create boundary conditions, and give them to the physis packages.
#-------------------------------------------------------------------------------
xPlane0 = Plane(Vector(xmin), Vector(1.0))
xbc0 = ReflectingBoundary(xPlane0)
##xPlane1 = Plane(Vector(xmax), Vector(-1.0))
##xbc1 = ReflectingBoundary(xPlane1)

hydro.appendBoundary(xbc0)

#-------------------------------------------------------------------------------
# Construct a predictor corrector integrator, and add the physics packages.
#-------------------------------------------------------------------------------
integrator = CheapSynchronousRK2Integrator(db)
integrator.appendPhysicsPackage(hydro)
integrator.lastDt = dt
if dtMin:
    integrator.dtMin = dtMin
if dtMax:
    integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
integrator.domainDecompositionIndependent = domainIndependent
output("integrator")
output("integrator.havePhysicsPackage(hydro)")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.dtGrowth")
output("integrator.domainDecompositionIndependent")

#-------------------------------------------------------------------------------
# Build the controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName)
output("control")

#-------------------------------------------------------------------------------
# Smooth the initial conditions/restore state.
#-------------------------------------------------------------------------------
if restoreCycle is not None:
    control.loadRestartFile(restoreCycle)
else:
    control.iterateIdealH(hydro)
    control.smoothState(smoothIters)

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
print("Tantalum diagnositic node is: ", interfaceHistory.tantalumIndex, interfaceHistory.tantalumProc)
print("Sapphire diagnositic node is: ", interfaceHistory.sapphireIndex, interfaceHistory.sapphireProc)
if not steps is None:
    control.step(steps)
    raise RuntimeError("Completed %i steps." % steps)

while control.time() < goalTime:
    nextGoalTime = min(control.time() + dtSample, goalTime)
    control.advance(nextGoalTime, maxSteps)
    interfaceHistory.sample(control.time())
    control.dropRestartFile()

    tantalumS = 1e40
    tantalumDSDt = 1e40
    tantalumPS = 1e40
    if mpi.rank == interfaceHistory.tantalumProc:
        tantalumS = nodesTantalum.deviatoricStress()[interfaceHistory.tantalumIndex].xx
        tantalumPS = nodesTantalum.plasticStrain()[interfaceHistory.tantalumIndex]
    print("Tantalum deviatoric stress and plastic strain:  %g %g" % (mpi.allreduce(tantalumS, mpi.MIN),
                                                                     mpi.allreduce(tantalumPS, mpi.MIN)))

#-------------------------------------------------------------------------------
# Read in the reference sapphire velocity history.
#-------------------------------------------------------------------------------
f = open("IntegrateDensity/PlateImpact-interface-history.txt", "r")
refthist = []
refvhist = []
for line in f:
    if line[0] != "#":
        t, vs, vt = (eval(x) for x in line.split())
        refthist.append(t)
        refvhist.append(vs)

#-------------------------------------------------------------------------------
# Did the test pass?
#-------------------------------------------------------------------------------
EError = (control.conserve.EHistory[-1] - control.conserve.EHistory[0])/control.conserve.EHistory[0]
print("Energy discrepency: %g" % EError)

assert len(interfaceHistory.sapphireVelocity) == len(refvhist)
for vtest, vref in zip(interfaceHistory.sapphireVelocity, refvhist):
    assert fuzzyEqual(vtest, vref, vtol)

print("Test passed.")

#-------------------------------------------------------------------------------
# If requested, write out the state in a global ordering to a file.
#-------------------------------------------------------------------------------
if outputFile:
    outputFile = os.path.join(restartDir, outputFile)
    pos = db.fluidPosition
    rho = db.fluidMassDensity
    P = db.newFluidScalarFieldList(0.0, "pressure")
    db.fluidPressure(P)
    vel = db.fluidVelocity
    eps = db.fluidSpecificThermalEnergy
    Hfield = db.fluidHfield
    xprof = mpi.reduce([x.x for x in internalValues(pos)], mpi.SUM)
    rhoprof = mpi.reduce(internalValues(rho), mpi.SUM)
    Pprof = mpi.reduce(internalValues(P), mpi.SUM)
    vprof = mpi.reduce([v.x for v in internalValues(vel)], mpi.SUM)
    epsprof = mpi.reduce(internalValues(eps), mpi.SUM)
    hprof = mpi.reduce([1.0/sqrt(H.Determinant()) for H in internalValues(Hfield)], mpi.SUM)
    mof = mortonOrderIndices(db)
    mo = mpi.reduce(internalValues(mof), mpi.SUM)
    if mpi.rank == 0:
        multiSort(mo, xprof, rhoprof, Pprof, vprof, epsprof, hprof)
        f = open(outputFile, "w")
        for xi, rhoi, Pi, vi, epsi, hi, mi in zip(xprof, rhoprof, Pprof, vprof, epsprof, hprof, mo):
            f.write((6*"%16.12e " + "%i\n") % (xi, rhoi, Pi, vi, epsi, hi, mi))
            f.write((6*"%i " + "\n") % (unpackElementUL(packElementDouble(xi)),
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
            comparisonFile = os.path.join(restartDir, comparisonFile)
            import filecmp
            assert filecmp.cmp(outputFile, comparisonFile)
