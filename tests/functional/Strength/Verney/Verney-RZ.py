#-------------------------------------------------------------------------------
# Spherical (2D RZ) version.
# An idealized strength test test where an imploding shell is ultimately stopped
# at a known radius due to plastic work dissipation.
#
# See LA-14379, Howell & Ball 2002 JCP
#-------------------------------------------------------------------------------
from math import *
import shutil
import mpi

from SpheralRZ import *
from SpheralTestUtilities import *
from findLastRestart import *
from NodeHistory import NodeHistory
from GenerateNodeDistribution2d import *
from VoronoiDistributeNodes import distributeNodes2d

#-------------------------------------------------------------------------------
# Identify ourselves!
#-------------------------------------------------------------------------------
title("RZ Verney imploding shell test.")

#-------------------------------------------------------------------------------
# Use (cm, gm, usec) units as in paper.
#-------------------------------------------------------------------------------
units = PhysicalConstants(0.01,  # Unit length in m
                          0.001, # Unit mass in kg
                          1e-6)  # Unit length in sec

#-------------------------------------------------------------------------------
# The analytic function paramterizing the shell evolution.
#-------------------------------------------------------------------------------
def F(alpha, lamb, R0, R1, n):

    class integfunc(ScalarScalarFunctor):
        def __init__(self, R0, R1):
            ScalarScalarFunctor.__init__(self)
            self.alpha = R1/R0 - 1.0
            self.thpt = 3.0*(self.alpha +
                             self.alpha*self.alpha +
                             self.alpha*self.alpha*self.alpha)
            return
        def __call__(self, x):
            return x*x*log(1.0 + self.thpt/(x*x*x))

    func = integfunc(R0, R1)
    return simpsonsIntegrationDouble(func, lamb, 1.0, n)

#-------------------------------------------------------------------------------
# Generic problem parameters
# All (cm, gm, usec) units.
#-------------------------------------------------------------------------------
commandLine(nr = 10,                     # Radial resolution of the shell in points
            seed = "constantDTheta",     # "lattice" or "constantDTheta"
            geometry = "quadrant",       # choose ("quadrant", "full").

            kernelOrder = 5,
            nPerh = 1.35,

            # Material specific bounds on the mass density.
            etamin = 1e-3,
            etamax = 1e3,

            # How many shells should we create tracer histories for?
            nshells = 10,

            # Hydro parameters.
            hydroType = "SPH",  # (SPH, CRKSPH, FSISPH)
            sph = False,   # This just chooses the H algorithm -- you can use this with CRKSPH for instance.
            Cl = None,
            Cq = None,
            linearInExpansion = None,
            Qlimiter = None,
            balsaraCorrection = None,
            epsilon2 = None,
            negligibleSoundSpeed = None,
            csMultiplier = None,
            hmin = 1e-5,
            hmax = 10.0,
            cfl = 0.25,
            useVelocityMagnitudeForDt = False,
            XSPH = True,
            epsilonTensile = 0.0,
            nTensile = 4,
            filter = 0.0,
            HUpdate = IdealH,
            densityUpdate = IntegrateDensity,
            compatibleEnergy = False,
            evolveTotalEnergy = True,
            gradhCorrection = False,

            # Time integration
            IntegratorConstructor = CheapSynchronousRK2Integrator,
            goalTime = 150.0,
            steps = None,
            dt = 1e-6,
            dtMin = 1e-6,
            dtMax = 10.0,
            dtGrowth = 2.0,
            maxSteps = None,
            statsStep = 10,
            domainIndependent = False,
            dtverbose = False,
            restoreCycle = -1,
            restartStep = 500,
            sampleFreq = 10,
            vizTime = 1.0,
            vizStep = 100,

            graphics = True,

            clearDirectories = False,
            dataDirBase = "dumps-Verney-Be-RZ",
            outputFile = "Verney-Be-RZ.gnu",
        )

assert seed in ("lattice", "constantDTheta")
assert geometry in ("quadrant", "full")

hydroType = hydroType.upper()
assert hydroType in ("SPH", "CRKSPH", "FSISPH")

# Material parameters for this test problem.
rho0Be = 1.845
R0, R1 = 8.0, 10.0        # Inner, outer initial radius
G0, Y0 = 1.510, 0.0033    # Shear modulus and yield strength
r0 = 4.0                  # Expected final inner radius
delta = R1 - R0
lamb = r0/R0
alpha = delta/R0
Fval = F(alpha, lamb, R0, R1, 1000)
u0 = sqrt(4.0*Y0*R1*Fval/(3.0*rho0Be*delta))
print("  lambda = %s\n  alpha = %s\n  F = %s\n  u0 = %s\n" % (lamb, alpha, Fval, u0))

hydroPath = "Solid" + hydroType

# Directories.
dataDir = os.path.join(dataDirBase,
                       hydroPath,
                       "densityUpdate=%s" % densityUpdate,
                       "compatibleEnergy=%s_totalEnergy=%s" % (compatibleEnergy, evolveTotalEnergy),
                       seed,
                       geometry,
                       "nr=%i" % nr)
restartDir = os.path.join(dataDir, "restarts")
vizDir = os.path.join(dataDir, "visit")
vizBaseName = "Verney-Be-%i" % nr
restartBaseName = os.path.join(restartDir, "Verney-%i" % nr)
historyOutputName = os.path.join(dataDir, "Verney-shell%i.dat")

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
# Material properties.
#-------------------------------------------------------------------------------
eosBe = OsborneEquationOfState( 1.85, etamin, etamax,     # Parameters from Howell & Ball 2002
                                0.951168,   # a1
                                0.345301,   # a2pos
                               -0.345301,   # asneg
                                0.926914,   # b0
                                2.948420,   # b1
                                0.507979,   # b2pos
                                0.507979,   # b2neg
                                0.564362,   # c0
                                0.620422,   # c1
                                0.0,        # c2pos
                                0.0,        # c2neg
                                0.8,
                                9.015,
                                units,
                                minimumPressure = -0.1)
strengthModelBe = ConstantStrength(G0, Y0)

#-------------------------------------------------------------------------------
# Create our interpolation kernels -- one for normal hydro interactions, and
# one for use with the artificial viscosity
#-------------------------------------------------------------------------------
WT = TableKernel(NBSplineKernel(kernelOrder), 1000)
output("WT")

#-------------------------------------------------------------------------------
# Create the NodeLists.
#-------------------------------------------------------------------------------
nodesBe = makeSolidNodeList("Beryllium", eosBe, strengthModelBe,
                            nPerh = nPerh,
                            hmin = hmin,
                            hmax = hmax,
                            rhoMin = etamin*rho0Be,
                            rhoMax = etamax*rho0Be,
                            xmin = -100.0*Vector.one,
                            xmax =  100.0*Vector.one)
nodeSet = [nodesBe]

#-------------------------------------------------------------------------------
# Set node properties (positions, masses, H's, etc.)
#-------------------------------------------------------------------------------
print("Generating node distribution.")
nx = int(R1/(R1 - R0) * nr + 0.5)
if geometry == "quadrant":
    xmin = (0.0, 0.0)
    xmax = (R1, R1)
    thetamax = 0.5*pi
else:
    xmin = (-R1, 0.0)
    xmax = ( R1, R1)
    nx = 2*nx
    thetamax = pi
if seed == "constantDTheta":
    nx = nr

gen = RZGenerator(GenerateNodeDistribution2d(nx, nx, rho0Be,
                                             distributionType = seed,
                                             xmin = xmin,
                                             xmax = xmax,
                                             theta = thetamax,
                                             rmin = R0, 
                                             rmax = R1))

distributeNodes2d((nodesBe, gen))
output("mpi.reduce(nodesBe.numInternalNodes, mpi.MIN)")
output("mpi.reduce(nodesBe.numInternalNodes, mpi.MAX)")
output("mpi.reduce(nodesBe.numInternalNodes, mpi.SUM)")

# Set node velocites.
pos = nodesBe.positions()
vel = nodesBe.velocity()
for i in range(nodesBe.numInternalNodes):
    ri = pos[i].magnitude()
    rhat = pos[i].unitVector()
    vel[i] = -u0 * (R0/ri)**2 * rhat

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase()
for n in nodeSet:
    db.appendNodeList(n)
del n
output("db")
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
if hydroType == "CRKSPH":
    hydro = CRKSPH(dataBase = db,
                   W = WT,
                   filter = filter,
                   cfl = cfl,
                   compatibleEnergyEvolution = compatibleEnergy,
                   evolveTotalEnergy = evolveTotalEnergy,
                   XSPH = XSPH,
                   densityUpdate = densityUpdate,
                   HUpdate = HUpdate)
elif hydroType == "FSISPH":
    hydro = FSISPH(dataBase = db,
                   W = WT,
                   # filter = filter, # TODO: add filter to FSISPH too
                   cfl = cfl,
                   interfaceMethod = HLLCInterface,
                   densityStabilizationCoefficient = 0.00,
                   compatibleEnergyEvolution = compatibleEnergy,
                   HUpdate = HUpdate,
                   epsTensile = epsilonTensile,
                   nTensile = nTensile)
else:
    assert hydroType == "SPH"
    hydro = SPH(dataBase = db,
                W = WT,
                filter = filter,
                cfl = cfl,
                compatibleEnergyEvolution = compatibleEnergy,
                evolveTotalEnergy = evolveTotalEnergy,
                gradhCorrection = gradhCorrection,
                densityUpdate = densityUpdate,
                HUpdate = HUpdate,
                XSPH = XSPH,
                epsTensile = epsilonTensile,
                nTensile = nTensile)
output("hydro")
output("hydro.cfl")
output("hydro.useVelocityMagnitudeForDt")
output("hydro._smoothingScaleMethod.HEvolution")
output("hydro.densityUpdate")
output("hydro.compatibleEnergyEvolution")
output("hydro.kernel")
if hydroType != "FSISPH":
    output("hydro.PiKernel")
    # TODO: add this output to FSISPH too

#-------------------------------------------------------------------------------
# Tweak the artificial viscosity.
#-------------------------------------------------------------------------------
q = hydro.Q
if not Cl is None:
    q.Cl = Cl
if not Cq is None:
    q.Cq = Cq
if not Qlimiter is None:
    q.limiter = Qlimiter
if not epsilon2 is None:
    q.epsilon2 = epsilon2
if not linearInExpansion is None:
    q.linearInExpansion = linearInExpansion
output("q")
output("q.Cl")
output("q.Cq")
output("q.limiter")
output("q.epsilon2")
output("q.linearInExpansion")
output("q.quadraticInExpansion")

#-------------------------------------------------------------------------------
# Boundary conditions.
#-------------------------------------------------------------------------------
if geometry == "quadrant":
    xPlane = Plane(Vector(0.0, 0.0), Vector(1.0, 0.0))
    xbc = ReflectingBoundary(xPlane)
    hydro.appendBoundary(xbc)

#-------------------------------------------------------------------------------
# Construct a time integrator.
#-------------------------------------------------------------------------------
integrator = IntegratorConstructor(db)
integrator.appendPhysicsPackage(hydro)
integrator.lastDt = dt
if dtMin:
    integrator.dtMin = dtMin
if dtMax:
    integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
integrator.domainDecompositionIndependent = domainIndependent
integrator.verbose = dtverbose
output("integrator")
output("integrator.havePhysicsPackage(hydro)")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.dtGrowth")
output("integrator.domainDecompositionIndependent")

#-------------------------------------------------------------------------------
# Select points to follow histories.
#-------------------------------------------------------------------------------
def verneySample(nodes, indices):
    mass = nodes.mass()
    pos = nodes.positions()
    vel = nodes.velocity()
    rho = nodes.massDensity()
    eps = nodes.specificThermalEnergy()
    P = ScalarField("pressure", nodes)
    nodes.pressure(P)
    ps = nodes.plasticStrain()
    H = nodes.Hfield()
    S = nodes.deviatoricStress()
    mshell = mpi.allreduce(sum([mass[i] for i in indices] + [0.0]), mpi.SUM)
    assert mshell > 0.0
    Srr = []
    for i in indices:
        Si = SymTensor(S[i])
        Si.rotationalTransform(rotationMatrix(pos[i].unitVector()))
        Srr.append(Si.xx)
    rshell = mpi.allreduce(sum([mass[i]*pos[i].magnitude() for i in indices] + [0.0]), mpi.SUM)
    vshell = mpi.allreduce(sum([mass[i]*vel[i].dot(pos[i].unitVector()) for i in indices] + [0.0]), mpi.SUM)
    rhoshell = mpi.allreduce(sum([mass[i]*rho[i] for i in indices] + [0.0]), mpi.SUM)
    epsshell = mpi.allreduce(sum([mass[i]*eps[i] for i in indices] + [0.0]), mpi.SUM)
    Pshell = mpi.allreduce(sum([mass[i]*P[i] for i in indices] + [0.0]), mpi.SUM)
    Sshell = mpi.allreduce(sum([mass[i]*Srri for i, Srri in zip(indices, Srr)] + [0.0]), mpi.SUM)
    strainShell = mpi.allreduce(sum([mass[i]*ps[i] for i in indices] + [0.0]), mpi.SUM)
    hshell = mpi.allreduce(sum([mass[i]*2.0/(H[i].Trace()) for i in indices] + [0.0]), mpi.SUM)
    rshell /= mshell
    vshell /= mshell
    rhoshell /= mshell
    epsshell /= mshell
    Pshell /= mshell
    Sshell /= mshell
    strainShell /= mshell
    hshell /= mshell
    return rshell, vshell, rhoshell, epsshell, Pshell, Sshell, strainShell, hshell

# Find shells of points in binned radii
histories = []
dr = (R1 - R0)/nshells
pos = nodesBe.positions()
shellIndices = [[] for i in range(nr)]
for i in range(nodesBe.numInternalNodes):
    ishell = min(nr - 1, int((pos[i].magnitude() - R0)/dr + 0.5))
    shellIndices[ishell].append(i)
for ishell in range(nshells):
    n = mpi.allreduce(len(shellIndices[ishell]), mpi.SUM)
    print("Selected %i nodes for shell %i." % (n, ishell))
    if n > 0:
        histories.append(NodeHistory(nodesBe, shellIndices[ishell], verneySample, historyOutputName % ishell, 
                                     labels = ("r", "vel", "rho", "eps", "pressure", "Srr", "plastic strain", "h")))

#-------------------------------------------------------------------------------
# Build the controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle,
                            vizBaseName = vizBaseName,
                            vizDir = vizDir,
                            vizTime = vizTime,
                            vizStep = vizStep,
                            periodicWork = [(hist.sample, sampleFreq) for hist in histories])
output("control")

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
if not steps is None:
    control.step(steps)
else:
    control.advance(goalTime)
    control.dropRestartFile()
    control.conserve.writeHistory(os.path.join(dataDir, "conservation_history.gnu"))

    # Check the tracer info for the innermost-shell, reporting the discrepancy
    # from the expected final stopping radius.  Not quite correct since we're sampling
    # the interior of the shell, but it gives us something quantitative.
    histories[0].sample(control.totalSteps, control.time(), control.lastDt())
    rsim0 = histories[0].sampleHistory[-1][0]
    print("Simulation measured final inner shell radius of %g cm: errror %g cm." % (rsim0,
                                                                                    rsim0 - r0))
Eerror = (control.conserve.EHistory[-1] - control.conserve.EHistory[0])/control.conserve.EHistory[0]
print("Total energy error: %g" % Eerror)

#-------------------------------------------------------------------------------
# If requested, write out the state in a global ordering to a file.
#-------------------------------------------------------------------------------
if outputFile:
    from SpheralTestUtilities import multiSort
    state = State(db, integrator.physicsPackages())
    outputFile = os.path.join(dataDir, outputFile)
    pos = state.vectorFields(HydroFieldNames.position)
    rho = state.scalarFields(HydroFieldNames.massDensity)
    P = state.scalarFields(HydroFieldNames.pressure)
    vel = state.vectorFields(HydroFieldNames.velocity)
    eps = state.scalarFields(HydroFieldNames.specificThermalEnergy)
    Hfield = state.symTensorFields(HydroFieldNames.H)
    S = state.symTensorFields(SolidFieldNames.deviatoricStress)
    ps = state.scalarFields(SolidFieldNames.plasticStrain)
    rprof = mpi.reduce([x.magnitude() for x in internalValues(pos)], mpi.SUM)
    rhoprof = mpi.reduce(internalValues(rho), mpi.SUM)
    Pprof = mpi.reduce(internalValues(P), mpi.SUM)
    vprof = mpi.reduce([v.x for v in internalValues(vel)], mpi.SUM)
    epsprof = mpi.reduce(internalValues(eps), mpi.SUM)
    hprof = mpi.reduce([2.0/(H.Trace()) for H in internalValues(Hfield)], mpi.SUM)
    sprof = mpi.reduce([x.xx for x in internalValues(S)], mpi.SUM)
    psprof = mpi.reduce(internalValues(ps), mpi.SUM)
    mof = mortonOrderIndices(db)
    mo = mpi.reduce(internalValues(mof), mpi.SUM)
    if mpi.rank == 0:
        multiSort(mo, rprof, rhoprof, Pprof, vprof, epsprof, hprof, sprof, psprof)
        f = open(outputFile, "w")
        f.write(("#" + 17*' "%16s"' + "\n") % ("r", "rho", "P", "v", "eps", "h", "S", "plastic strain", "m", 
                                               "int(r)", "int(rho)", "int(P)", "int(v)", "int(eps)", "int(h)", "int(S)", "int(ps)"))
        for (ri, rhoi, Pi, vi, epsi, hi, si, psi, mi) in zip(rprof, rhoprof, Pprof, vprof, epsprof, hprof, sprof, psprof, mo):
            f.write((8*"%16.12e " + 9*"%i " + "\n") %
                    (ri, rhoi, Pi, vi, epsi, hi, si, psi, mi,
                     unpackElementUL(packElementDouble(ri)),
                     unpackElementUL(packElementDouble(rhoi)),
                     unpackElementUL(packElementDouble(Pi)),
                     unpackElementUL(packElementDouble(vi)),
                     unpackElementUL(packElementDouble(epsi)),
                     unpackElementUL(packElementDouble(hi)),
                     unpackElementUL(packElementDouble(si)),
                     unpackElementUL(packElementDouble(psi))))
        f.close()

#-------------------------------------------------------------------------------
# Plot the state.
#-------------------------------------------------------------------------------
if graphics:
    from SpheralGnuPlotUtilities import *
    state = State(db, integrator.physicsPackages())
    rhoPlot = plotFieldList(state.scalarFields("mass density"),
                            xFunction = "%s.magnitude()",
                            plotStyle="points",
                            winTitle="rho @ %g" % (control.time()))
    radialVelocity = radialVelocityFieldList(db.fluidPosition,
                                             db.fluidVelocity)
    velPlot = plotFieldList(radialVelocity,
                            xFunction = "%s.magnitude()",
                            plotStyle="points",
                            winTitle="vel @ %g" % (control.time()))
    mPlot = plotFieldList(state.scalarFields("mass"),
                          xFunction = "%s.magnitude()",
                          plotStyle="points",
                          winTitle="mass @ %g" % (control.time()))
    PPlot = plotFieldList(state.scalarFields("pressure"),
                          xFunction = "%s.magnitude()",
                          plotStyle="points",
                          winTitle="pressure @ %g" % (control.time()))
    hPlot = plotFieldList(state.symTensorFields("H"),
                          xFunction = "%s.magnitude()",
                          yFunction = "2.0/%s.Trace()",
                          plotStyle="points",
                          winTitle="h @ %g" % (control.time()))
    psPlot = plotFieldList(state.scalarFields(SolidFieldNames.plasticStrain),
                           xFunction = "%s.magnitude()",
                           plotStyle="points",
                           winTitle="plastic strain @ %g" % (control.time()))
