#-------------------------------------------------------------------------------
# A rod of aluminum undergoing compaction by a piston.  This test is based on the
# first example test from
#
# Jutzi, M., Benz, W., & Michel, P. (2008). Numerical simulations of impacts
# involving porous bodies.I. Implementing sub-resolution porosity in a 3D SPH
# hydrocode. Icarus, 198(1), 242–255.
#-------------------------------------------------------------------------------
#
# Ordinary SPH
#
#ATS:t0 = test(      SELF, "--graphics False --clearDirectories True  --checkError True  --dataDirBase dumps-PlanarCompaction-1d-sph --restartStep 100000 --postCleanup True", np=4, label="Planar porous aluminum compaction problem -- 1-D (4 proc)")
#ATS:t1 = test(      SELF, "--graphics False --clearDirectories True  --checkError False --dataDirBase dumps-PlanarCompaction-1d-sph-restart --restartStep 100 --steps 200", label="Planar porous aluminum compaction problem -- 1-D (serial, restart test step 1)")
#ATS:t2 = testif(t1, SELF, "--graphics False --clearDirectories False --checkError False --dataDirBase dumps-PlanarCompaction-1d-sph-restart --restartStep 100 --steps 100 --checkRestart True --restoreCycle 100 --postCleanup True", label="Planar porous aluminum compaction problem -- 1-D (serial, restart test step 2)")
#
# Ordinary SPH with no initial porosity
#
#ATS:t3 = testif(t0, SELF, "--alpha0 1.0 --useDamage False --goalTime 1.0 --graphics False --clearDirectories True  --checkError True  --dataDirBase dumps-PlanarCompaction-1d-sph --restartStep 100000 --postCleanup True", np=4, label="Planar porous aluminum compaction problem with no initial porosity -- 1-D (4 proc, no damage model)")
#ATS:t4 = testif(t3, SELF, "--alpha0 1.0 --useDamage True  --goalTime 1.0 --graphics False --clearDirectories True  --checkError True  --dataDirBase dumps-PlanarCompaction-1d-sph --restartStep 100000 --postCleanup True", np=4, label="Planar porous aluminum compaction problem with no initial porosity -- 1-D (4 proc, with damage model)")
#
# FSISPH
#
#ATS:t10 = test(       SELF, "--graphics False --clearDirectories True  --checkError True  --hydroType FSISPH --dataDirBase dumps-PlanarCompaction-1d-fsisph --restartStep 100000 --postCleanup True", np=4, label="Planar porous aluminum compaction problem -- 1-D (FSISPH, 4 proc)", fsisph=True)
#ATS:t11 = test(       SELF, "--graphics False --clearDirectories True  --checkError False --hydroType FSISPH --dataDirBase dumps-PlanarCompaction-1d-fsisph-restart --restartStep 100 --steps 200", label="Planar porous aluminum compaction problem -- 1-D (FSISPH, serial, restart test step 1)", fsisph=True)
#ATS:t12 = testif(t11, SELF, "--graphics False --clearDirectories False --checkError False --hydroType FSISPH --dataDirBase dumps-PlanarCompaction-1d-fsisph-restart --restartStep 100 --steps 100 --checkRestart True --restoreCycle 100 --postCleanup True", label="Planar porous aluminum compaction problem -- 1-D (FSISPH, serial, restart test step 2)", fsisph=True)
#
# CRKSPH
#
#ATS:t20 = test(       SELF, "--graphics False --clearDirectories True  --checkError True  --hydroType CRKSPH --dataDirBase dumps-PlanarCompaction-1d-crksph --restartStep 100000 --postCleanup True", np=4, label="Planar porous aluminum compaction problem -- 1-D (CRKSPH, 4 proc)")
#ATS:t21 = test(       SELF, "--graphics False --clearDirectories True  --checkError False --hydroType CRKSPH --dataDirBase dumps-PlanarCompaction-1d-crksph-restart --restartStep 100 --steps 200", label="Planar porous aluminum compaction problem -- 1-D (CRKSPH, serial, restart test step 1)")
#ATS:t22 = testif(t21, SELF, "--graphics False --clearDirectories False --checkError False --hydroType CRKSPH --dataDirBase dumps-PlanarCompaction-1d-crksph-restart --restartStep 100 --steps 100 --checkRestart True --restoreCycle 100 --postCleanup True", label="Planar porous aluminum compaction problem -- 1-D (CRKSPH, serial, restart test step 2)")

from SolidSpheral1d import *
from SpheralTestUtilities import *
from SpheralMatplotlib import *
from PlanarCompactionSolution import *
from math import *
import os, sys, shutil
import mpi
import Pnorm

#-------------------------------------------------------------------------------
# Identify ourselves!
#-------------------------------------------------------------------------------
title("1-D planar compaction of a poroous Al rod")

#-------------------------------------------------------------------------------
# Build our units (cm, gm, microsec)
#-------------------------------------------------------------------------------
units = CGuS()

#-------------------------------------------------------------------------------
# Generic problem parameters
# All (cm, gm, usec) units.
#-------------------------------------------------------------------------------
commandLine(nx = 500,                          # Number of internal free points
            vpiston = -45.8e-3,                # Fixed velocity on right-end piston
            nxpiston = 10,                     # How many points should we put in the piston boundary 
            nPerh = 4.01,

            # Select the frame of the fluid so we're either pushing a piston or running into a wall
            boundary = "wall",                 # (wall, piston) 

            # Material models
            material = "aluminum melosh89",    # One of valid materials in our MaterialLibrary
            EOSConstructor = TillotsonEquationOfState,

            # Porosity
            PorousModel = PalphaPorosity,
            alpha0 = 1.275,
            jutziStateUpdate = True,
            fdt = 0.5,                         # Timestep control fractional change in alpha

            # Hydro
            hydroType = "SPH",                 # SPH, CRKSPH, FSISPH
            cfl = 0.25,
            useVelocityMagnitudeForDt = True,
            XSPH = False,
            epsilonTensile = 0.3,
            nTensile = 4,
            volumeType = RKSumVolume,
            Cl = None,                         # Artificial viscosity linear coefficient
            Cq = None,                         # Artificial viscosity quadratic coefficient

            # Optionally allow damage
            useDamage = False,

            # Time integration
            IntegratorConstructor = VerletIntegrator,
            goalTime = 3.5,
            steps = None,
            dt = 1e-10,
            dtMin = 1e-10,
            dtMax = 10.0,
            dtGrowth = 2.0,
            maxSteps = None,
            statsStep = 10,
            densityUpdate = IntegrateDensity,
            compatibleEnergy = True,
            domainIndependent = True,
            dtverbose = False,

            # Options for implicit integrators
            ftol = None,
            convergenceTolerance = None,
            maxIterations = None,
            maxAllowedDtMultiplier = None,
            beta = None,

            # Problem control
            restoreCycle = -1,
            restartStep = 500,

            # Output
            graphics = True,
            clearDirectories = False,
            postCleanup = False,
            dataDirBase = "dumps-PlanarCompaction-1d",
            checkError = False,
            checkRestart = False,
            outputFile = None,
            comparisonFile = None,

            # Parameters for the test acceptance.,
            tol = 1.0e-5,
            )

hydroType = hydroType.upper()

boundary = boundary.upper()
assert boundary in ("PISTON", "WALL")

dataDir = os.path.join(dataDirBase,
                       PorousModel.__name__,
                       material.replace(" ", "_") + "_" + EOSConstructor.__name__,
                       hydroType,
                       boundary,
                       "nx=%i" % nx)

restartDir = os.path.join(dataDir, "restarts")
restartBaseName = os.path.join(restartDir, "PlanarCompaction-%i" % nx)

#-------------------------------------------------------------------------------
# The reference values for error norms checking for pass/fail
#-------------------------------------------------------------------------------
LnormRef = {("SPH", 1.275) : {"Mass density" : {"L1"   : 0.06784186300927694,
                                                "L2"   : 0.012774373437299643,
                                                "Linf" : 0.6245679444354701},
                              "Spec Therm E" : {"L1"   : 0.0001200742460791407,
                                                "L2"   : 2.2616105613742583e-05,
                                                "Linf" : 0.0010923440797387786},
                              "velocity    " : {"L1"   : 0.004921931042558655,
                                                "L2"   : 0.0009173594117158436,
                                                "Linf" : 0.0448725433453345},
                              "pressure    " : {"L1"   : 0.0022217375280911347,
                                                "L2"   : 0.00039479153550769805,
                                                "Linf" : 0.018793913196205617},
                              "alpha       " : {"L1"   : 0.0590391542763204,
                                                "L2"   : 0.007963583413760916,
                                                "Linf" : 0.2738180402369801},
                              "h           " : {"L1"   : 0.00043261838627472803,
                                                "L2"   : 8.062946952637553e-05,
                                                "Linf" : 0.014201309070925212}},

            ("SPH", 1.0) : {"Mass density" : {"L1"   : 0.005416588903113661,  
                                              "L2"   : 0.001994539770945366,  
                                              "Linf" : 0.22293554904196977},  
                            "Spec Therm E" : {"L1"   : 2.9911432798797925e-05,
                                              "L2"   : 1.0394781088376176e-05,
                                              "Linf" : 0.001048750342669405}, 
                            "velocity    " : {"L1"   : 0.0010892623667423551, 
                                              "L2"   : 0.0004049477641523086, 
                                              "Linf" : 0.04538028008374117},  
                            "pressure    " : {"L1"   : 0.0018120800940308072, 
                                              "L2"   : 0.000663492813279845,  
                                              "Linf" : 0.07305139125993931},  
                            "alpha       " : {"L1"   : 0.0,                   
                                              "L2"   : 0.0,                   
                                              "Linf" : 0.0},                  
                            "h           " : {"L1"   : 6.11935370574284e-05,  
                                              "L2"   : 3.08688837519722e-05,  
                                              "Linf" : 0.014204020975568329}},

            ("FSISPH", 1.275) : {"Mass density" : {"L1"   : 0.06781314493410028,
                                                   "L2"   : 0.012767602580471844,
                                                   "Linf" : 0.6245698198195724},
                                 "Spec Therm E" : {"L1"   : 0.00011965457154529887,
                                                   "L2"   : 2.254867764655585e-05,
                                                   "Linf" : 0.0010923441579020216},
                                 "velocity    " : {"L1"   : 0.004913235403786584,
                                                   "L2"   : 0.0009163181868064007,
                                                   "Linf" : 0.044872645505280966},
                                 "pressure    " : {"L1"   : 0.002210222289968727,
                                                   "L2"   : 0.00039387994606202237,
                                                   "Linf" : 0.018793847854203908},
                                 "alpha       " : {"L1"   : 0.05903530352469833,
                                                   "L2"   : 0.007960855343380124,
                                                   "Linf" : 0.2738175996911776},
                                 "h           " : {"L1"   : 0.0004319674641541397,
                                                   "L2"   : 8.05536967933465e-05,
                                                   "Linf" : 0.014201309071287523}},

            ("CRKSPH", 1.275) : {"Mass density" : {"L1"   : 0.0679707220773017,
                                                   "L2"   : 0.012783621322270034,
                                                   "Linf" : 0.6245687018640083},
                                 "Spec Therm E" : {"L1"   : 0.0001201303520771047,
                                                   "L2"   : 2.2622550331357265e-05,
                                                   "Linf" : 0.00109234444748564},
                                 "velocity    " : {"L1"   : 0.004926121368235753,
                                                   "L2"   : 0.0009173629106343205,
                                                   "Linf" : 0.044872623201390446},
                                 "pressure    " : {"L1"   : 0.0022258225868044238,
                                                   "L2"   : 0.00039496818856079414,
                                                   "Linf" : 0.018793890216913627},
                                 "alpha       " : {"L1"   : 0.05909030710035451,
                                                   "L2"   : 0.007965976598237856,
                                                   "Linf" : 0.2738173870953471},
                                 "h           " : {"L1"   : 0.00043274827065262975,
                                                   "L2"   : 8.062817025125132e-05,
                                                   "Linf" : 0.014201309070451522}},
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
# Al material parameters
# Note, in the test as presented in the reference there is no strength model,
# just pressure response.
#-------------------------------------------------------------------------------
eosS = EOSConstructor(material,
                      etamin = 0.1,
                      etamax = 10.0,
                      units = units)
strengthModelS = NullStrength()
rhoS0 = eosS.referenceDensity
eps0 = 0.0
Ps0 = eosS.pressure(rhoS0, eps0)
cS0 = eosS.soundSpeed(rhoS0, eps0)
print("Aluminum equation of state initial properties:")
output("  eosS")
output(" rhoS0")
output("  eps0")
output("   Ps0")
output("   cS0")

#-------------------------------------------------------------------------------
# Create our interpolation kernels
#-------------------------------------------------------------------------------
WT = TableKernel(WendlandC4Kernel(), 200)
output("WT")

#-------------------------------------------------------------------------------
# Create the NodeLists.
#-------------------------------------------------------------------------------
nodes = makeSolidNodeList("aluminium", eosS, strengthModelS,
                          nPerh = nPerh,
                          xmin = -10.0*Vector.one,
                          xmax =  10.0*Vector.one)
nodeSet = [nodes]

#-------------------------------------------------------------------------------
# Set node properties (positions, masses, H's, etc.)
#-------------------------------------------------------------------------------
print("Generating node distribution.")
from DistributeNodes import distributeNodesInRange1d
xmin = -1.0
xmax = 1.0
nxtot = nx
dx = (xmax - xmin)/nx
if boundary == "PISTON":
    xmax += nxpiston*dx
    nxtot += nxpiston
rho0 = rhoS0/alpha0
distributeNodesInRange1d([(nodes, nxtot, rho0, (xmin, xmax))])
output("mpi.reduce(nodes.numInternalNodes, mpi.MIN)")
output("mpi.reduce(nodes.numInternalNodes, mpi.MAX)")
output("mpi.reduce(nodes.numInternalNodes, mpi.SUM)")

# Set node velocites and find the boundary nodes
pos = nodes.positions()
vel = nodes.velocity()
for i in range(nodes.numInternalNodes):
    if boundary == "WALL":
        vel[i].x = -vpiston
    elif pos[i].x > 1.0:
        vel[i].x = vpiston

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
output("db.numSolidNodeLists")

#-------------------------------------------------------------------------------
# Boundary conditions
#-------------------------------------------------------------------------------
if boundary == "WALL":
    xbc0 = InflowOutflowBoundary(db, Plane(Vector(xmin), Vector(1.0)))
    xbc1 = ReflectingBoundary(Plane(Vector(xmax), Vector(-1.0)))
    bcs = [xbc1]
else:
    indices = [i for i in range(nodes.numInternalNodes) if pos[i].x > 1.0]
    xbc0 = ReflectingBoundary(Plane(Vector(xmin), Vector(1.0)))
    xbc1 = ConstantVelocityBoundary(nodes, indices)
    bcs = [xbc0, xbc1]

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
if hydroType == "CRKSPH":
    hydro = CRKSPH(dataBase = db,
                   W = WT,
                   cfl = cfl,
                   compatibleEnergyEvolution = compatibleEnergy,
                   XSPH = XSPH,
                   densityUpdate = densityUpdate,
                   useVelocityMagnitudeForDt = useVelocityMagnitudeForDt)
elif hydroType == "SPH":
    hydro = SPH(dataBase = db,
                W = WT,
                cfl = cfl,
                compatibleEnergyEvolution = compatibleEnergy,
                densityUpdate = densityUpdate,
                XSPH = XSPH,
                epsTensile = epsilonTensile,
                nTensile = nTensile,
                useVelocityMagnitudeForDt = useVelocityMagnitudeForDt)
elif hydroType == "FSISPH":
    hydro = FSISPH(dataBase = db,
                   W = WT,
                   cfl = cfl,
                   compatibleEnergyEvolution = compatibleEnergy,
                   useVelocityMagnitudeForDt = useVelocityMagnitudeForDt)
if not Cl is None:
    hydro.Q.Cl = Cl
if not Cq is None:
    hydro.Q.Cq = Cq
output("hydro")
output("  hydro.cfl")
output("  hydro.useVelocityMagnitudeForDt")
output("  hydro._smoothingScaleMethod.HEvolution")
output("  hydro.Q")
output("  hydro.Q.Cl")
output("  hydro.Q.Cq")
if hydro != "FSISPH":
    output("  hydro.densityUpdate")
    output("  hydro.compatibleEnergyEvolution")

packages = [hydro]

#-------------------------------------------------------------------------------
# Construct a porosity model
#-------------------------------------------------------------------------------
# First we need some conversion factors from CGS to our units
cgs = CGS()
mCGSconv = cgs.unitMassKg/units.unitMassKg
lCGSconv = cgs.unitLengthMeters/units.unitLengthMeters
tCGSconv = cgs.unitTimeSec/units.unitTimeSec
PCGSconv = mCGSconv/(tCGSconv*tCGSconv*lCGSconv)
vCGSconv = lCGSconv/tCGSconv

# Parameters for P-alpha
Pe = 8e8      # dynes/cm^2
Ps = 7e9      # dynes/cm^2
cS0 = 5.35e5  # cm/sec
ce = 4.11e5   # cm/sec
phi0 = 1.0 - 1.0/alpha0

# Parameters for Strain-alpha
def epsX_from_alphaX(alphaX,
                     phi0,
                     epsE,
                     kappa):
    if alphaX < 1.001:
        epsX = -1e10
    else:
        alpha0 = 1.0/(1.0 - phi0)
        epsX = log(alphaX/alpha0)/kappa + epsE
    return epsX
epsE = -1.88e-4                 # Elastic compaction limit (porosity)
kappa = 0.999                   # Exponential factor for distention (porosity)
alphaX = 1.0                    # transition distension exponential compation -> powerlaw
cS0 = 5.35e5                    # cm/sec
ce = 4.11e5   # cm/sec
epsX = epsX_from_alphaX(1.0, phi0, epsE, kappa)   # Transition from exponential to power-law distention (porosity)
gammaS0 = eosS.gamma(rhoS0, eps0)

if PorousModel is PalphaPorosity:
    print("P-alpha porosity model parameters:")
    output("  alpha0")
    output("      Pe")
    output("      Ps")
    output("     cS0")
    output("      ce")
    print("Computing cS0 from solid EOS yields ", eosS.soundSpeed(rhoS0, 0.0))
    porosityAl = PalphaPorosity(nodes,
                                phi0 = phi0,
                                Pe = Pe * PCGSconv,
                                Pt = Pe * PCGSconv,
                                Ps = Ps * PCGSconv,
                                alphat = None,
                                n1 = 0.0,
                                n2 = 2.0,
                                cS0 = cS0 * vCGSconv,
                                c0 = ce * vCGSconv,
                                jutziStateUpdate = jutziStateUpdate)
    porosityAl.fdt = fdt
    output("porosityAl")
    output("  porosityAl.Pe")
    output("  porosityAl.Pt")
    output("  porosityAl.Ps")
    output("  porosityAl.alphae")
    output("  porosityAl.alphat")
    output("  porosityAl.n1")
    output("  porosityAl.n2")
    output("  porosityAl.cS0")
    output("  porosityAl.fdt")
    output("  porosityAl.jutziStateUpdate")

elif PorousModel is StrainPorosity:
    print("Strain-alpha porosity model parameters:")
    output("  alpha0")
    output("    epsE")
    output("    epsX")
    output("   kappa")
    output("  alphaX")
    output("     cS0")
    output("      ce")
    porosityAl = StrainPorosity(nodes,           # Nodes we will make porous
                                phi0,            # Initial porosity
                                epsE,            # Elastic compaction limit
                                epsX,            # Exponetial->power-law transition
                                kappa,           # Compaction rate
                                gammaS0,         # Reference gamma at full density
                                cS0 * vCGSconv,  # Reference sound speed at full density
                                ce * vCGSconv,   # Reference sound speed at initial porosity
                                jutziStateUpdate = jutziStateUpdate)
    porosityAl.fdt = fdt
    output("porosityAl")
    output("  porosityAl.epsE")
    output("  porosityAl.epsX")
    output("  porosityAl.kappa")
    output("  porosityAl.gammaS0")
    output("  porosityAl.fdt")
    output("  porosityAl.jutziStateUpdate")

else:
    raise(RuntimeError, "Unknown porosity model")

packages.append(porosityAl)

#-------------------------------------------------------------------------------
# Optionally test damage as well
#-------------------------------------------------------------------------------
if useDamage:
    damage = ProbabilisticDamageModel(materialName = "aluminum",
                                      nodeList = nodes,
                                      units = units,
                                      kernel = WT)
    packages.append(damage)

#-------------------------------------------------------------------------------
# Construct a time integrator.
#-------------------------------------------------------------------------------
integrator = IntegratorConstructor(db, packages)
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
output("integrator.havePhysicsPackage(porosityAl)")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.dtGrowth")
output("integrator.domainDecompositionIndependent")

# Special stuff for implicit integrators
if isinstance(integrator, ImplicitIntegrator):
    if beta:
        integrator.beta = beta
    if convergenceTolerance:
        integrator.convergenceTolerance = convergenceTolerance
    if maxIterations:
        integrator.maxIterations = maxIterations
    if maxAllowedDtMultiplier:
        integrator.maxAllowedDtMultiplier = maxAllowedDtMultiplier
    output("integrator.beta")
    output("integrator.convergenceTolerance")
    output("integrator.maxIterations")
    output("integrator.maxAllowedDtMultiplier")
try:   # This will only work for BackwardEuler currently
    if ftol:
        integrator.ftol = ftol
    output("integrator.ftol")
except:
    pass

#-------------------------------------------------------------------------------
# Add the boundary conditions.
#-------------------------------------------------------------------------------
for package in integrator.physicsPackages():
    for bc in bcs:
        package.appendBoundary(bc)

#-------------------------------------------------------------------------------
# Build the controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            volumeType = volumeType,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle)
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

    # Optionally check restart state
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
    control.advance(goalTime)
    control.dropRestartFile()

Eerror = (control.conserve.EHistory[-1] - control.conserve.EHistory[0])/control.conserve.EHistory[0]
print("Total energy change: %g" % Eerror)

#-------------------------------------------------------------------------------
# Compute the analytic solution
#-------------------------------------------------------------------------------
# Need the smoothing scale for the solution if we want to plot and compare it
state = State(db, integrator.physicsPackages())
H = state.symTensorFields(HydroFieldNames.H)
h = db.newFluidScalarFieldList(0.0, "h")
for i in range(nodes.numInternalNodes):
    h[0][i] = 1.0/H[0][i].xx

solution = PlanarCompactionSolution(eos = eosS,
                                    vpiston = vpiston,
                                    eps0 = eps0,
                                    alpha0 = alpha0,
                                    alphat = None,
                                    Pe = Pe * PCGSconv,
                                    Pt = Pe * PCGSconv,
                                    Ps = Ps * PCGSconv,
                                    n1 = 0.0,
                                    n2 = 2.0,
                                    cS0 = cS0 * vCGSconv,
                                    c0 = ce * vCGSconv,
                                    h0 = nPerh * dx,
                                    nPoints = 1000,
                                    pistonFrame = (boundary == "WALL"))

#-------------------------------------------------------------------------------
# Compute the error norms
# We clip off the ends of the simulation to avoid wall heating effects in the
# model.
#-------------------------------------------------------------------------------
up = abs(vpiston)
t = control.time()
if boundary == "WALL":
    xmin, xmax = -1.0 + up*t, 1.0
else:
    xmin, xmax = -1.0, 1.0 - up*t
dxbound = 10*dx
xmin += dxbound
xmax -= dxbound
xprof = np.array(mpi.allreduce([x.x for x in pos.internalValues()], mpi.SUM))
rhoprof = np.array(mpi.allreduce(state.scalarFields(HydroFieldNames.massDensity)[0].internalValues(), mpi.SUM))
epsprof = np.array(mpi.allreduce(state.scalarFields(HydroFieldNames.specificThermalEnergy)[0].internalValues(), mpi.SUM))
vprof = np.array(mpi.allreduce([x.x for x in state.vectorFields(HydroFieldNames.velocity)[0].internalValues()], mpi.SUM))
Pprof = np.array(mpi.allreduce(state.scalarFields(HydroFieldNames.pressure)[0].internalValues(), mpi.SUM))
hprof = np.array(mpi.allreduce(h[0].internalValues(), mpi.SUM))
alphaprof = np.array(mpi.allreduce(state.scalarFields(SolidFieldNames.porosityAlpha)[0].internalValues(), mpi.SUM))
multiSort(xprof, rhoprof, epsprof, vprof, Pprof, hprof, alphaprof)
xans, vans, epsans, rhoans, Pans, hans = solution.solution(t, xprof)
xans, alphaans = solution.alpha_solution(t, xprof)
failure = False
if mpi.rank == 0:
    print("Quantity \t\tL1 \t\t\t\tL2 \t\t\t\tLinf")
    for (name, data, ans) in [("Mass density", rhoprof, rhoans),
                              ("Spec Therm E", epsprof, epsans),
                              ("velocity    ", vprof, vans),
                              ("pressure    ", Pprof, Pans),
                              ("alpha       ", alphaprof, alphaans),
                              ("h           ", hprof, hans)]:
        assert len(data) == len(ans)
        error = data - ans
        Pn = Pnorm.Pnorm(error, xprof)
        L1 = Pn.gridpnorm(1, xmin, xmax)
        L2 = Pn.gridpnorm(2, xmin, xmax)
        Linf = Pn.gridpnorm("inf", xmin, xmax)
        print(f"{name}\t\t{L1} \t\t{L2} \t\t{Linf}")

        if checkError and not (np.allclose(L1, LnormRef[(hydroType, alpha0)][name]["L1"], tol, tol) and
                               np.allclose(L2, LnormRef[(hydroType, alpha0)][name]["L2"], tol, tol) and
                               np.allclose(Linf, LnormRef[(hydroType, alpha0)][name]["Linf"], tol, tol)):
            print("Failing Lnorm tolerance for ", name, (L1, L2, Linf), LnormRef[(hydroType, alpha0)][name])
            failure = True
sys.stdout.flush()

failure = mpi.allreduce(failure, mpi.MAX)
if checkError and failure:
    raise ValueError("Error bounds violated")

#-------------------------------------------------------------------------------
# Plot the state.
#-------------------------------------------------------------------------------
if graphics:
    rhoPlot = plotFieldList(state.scalarFields("mass density"),
                            plotStyle = "o-",
                            lineTitle = "Simulation",
                            xlabel = r"$x$",
                            ylabel = r"$\rho$ (gm/cm$^3$)",
                            winTitle = r"$\rho$ @ %g %i" % (control.time(), mpi.procs))
    velPlot = plotFieldList(state.vectorFields("velocity"),
                            yFunction = "%s.x",
                            plotStyle = "o-",
                            lineTitle = "Simulation",
                            xlabel = r"$x$",
                            ylabel = r"$vel$ cm/$\mu$sec",
                            winTitle = "vel @ %g %i" % (control.time(), mpi.procs))
    PPlot = plotFieldList(state.scalarFields("pressure"),
                          plotStyle = "o-",
                          lineTitle = "Simulation",
                          xlabel = r"$x$",
                          ylabel = r"$P$ (Mbar)",
                          winTitle = "pressure @ %g %i" % (control.time(), mpi.procs),
                          semilogy = True)
    uPlot = plotFieldList(state.scalarFields("specific thermal energy"),
                          plotStyle = "o-",
                          lineTitle = "Simulation",
                          xlabel = r"$x$",
                          ylabel = r"$\varepsilon$ (Mbar cm$^3$/gm)",
                          winTitle = "specific thermal energy @ %g %i" % (control.time(), mpi.procs))
    SPlot = plotFieldList(state.symTensorFields("deviatoric stress"),
                          yFunction = "%s.xx",
                          plotStyle = "o-",
                          lineTitle = "Simulation",
                          xlabel = r"$x$",
                          ylabel = r"$S_{xx}$ (Mbar)",
                          winTitle = "$S_{xx}$ @ %g %i" % (control.time(), mpi.procs))
    hPlot = plotFieldList(h,
                          plotStyle = "o-",
                          lineTitle = "Simulation",
                          xlabel = r"$x$",
                          ylabel = r"$h$ (cm)",
                          winTitle = "h @ %g %i" % (control.time(), mpi.procs))
    csPlot = plotFieldList(state.scalarFields(HydroFieldNames.soundSpeed),
                           plotStyle = "o-",
                           lineTitle = "Simulation",
                           xlabel = r"$x$",
                           ylabel = r"$c_S$",
                           winTitle = "Sound speed @ %g %i" %  (control.time(), mpi.procs))
    alphaPlot = plotField(porosityAl.alpha,
                          plotStyle = "o-",
                          lineTitle = "Simulation",
                          xlabel = r"$x$",
                          ylabel = r"$\alpha$",
                          winTitle = r"$\alpha$ @ %g %i" %  (control.time(), mpi.procs))
    DalphaDtPlot = plotField(porosityAl.DalphaDt,
                             plotStyle = "o-",
                             lineTitle = "Simulation",
                             xlabel = r"$x$",
                             ylabel = r"$D\alpha/Dt$",
                             winTitle = r"$D\alpha/Dt$ @ %g %i" %  (control.time(), mpi.procs))
    phiPlot = plotField(porosityAl.phi(),
                        plotStyle = "o-",
                        lineTitle = "Simulation",
                        xlabel = r"$x$",
                        ylabel = r"$\phi$",
                        winTitle = r"$\phi$ @ %g %i" %  (control.time(), mpi.procs))
    fDSplot = plotField(porosityAl.fDS,
                         plotStyle = "o-",
                         lineTitle = "Simulation",
                         xlabel = r"$x$",
                         ylabel = r"$f_{DS}$",
                         winTitle = r"$f_{DS}$ @ %g %i" %  (control.time(), mpi.procs))
    Dplot = plotField(nodes.damage(),
                      yFunction = "%s.xx",
                      plotStyle = "o-",
                      lineTitle = "Simulation",
                      xlabel = r"$x$",
                      ylabel = r"$D_{xx}$",
                      winTitle = "Damage @ %g %i" %  (control.time(), mpi.procs))

    if PorousModel is PalphaPorosity:
        dPdRplot = plotField(porosityAl.partialPpartialRho,
                             plotStyle = "o-",
                             lineTitle = "Simulation",
                             xlabel = r"$x$",
                             ylabel = r"$\partial P/\partial \rho$",
                             winTitle = r"$\partial P/\partial \rho$ @ %g %i" %  (control.time(), mpi.procs))
        dPdUplot = plotField(porosityAl.partialPpartialEps,
                             plotStyle = "o-",
                             lineTitle = "Simulation",
                             xlabel = r"$x$",
                             ylabel = r"$\partial P/\partial \varepsilon$",
                             winTitle = r"$\partial P/\partial \varepsilon$ @ %g %i" %  (control.time(), mpi.procs))

    # Add the solution
    plotAnswer(solution, control.time(),
               rhoPlot = rhoPlot,
               velPlot = velPlot,
               epsPlot = uPlot,
               PPlot = PPlot,
               HPlot = hPlot)
    x, cs_solution = solution.soundSpeed_solution(control.time())
    csPlot.plot(x, cs_solution,
                "k-",
                label = "Solution")
    csPlot.axes.legend()
    x, alpha_solution = solution.alpha_solution(control.time())
    alphaPlot.plot(x, alpha_solution,
                   "k-",
                   label = "Solution")
    alphaPlot.axes.legend()
    phiPlot.plot(x, 1.0 - 1.0/alpha_solution,
                   "k-",
                   label = "Solution")
    phiPlot.axes.legend()

    plots = [(rhoPlot, "rho.png"),
             (velPlot, "vel.png"),
             (PPlot, "pressure.png"),
             (SPlot, "devstress.png"),
             (uPlot, "u.png"),
             (hPlot, "h.png"),
             (alphaPlot, "alpha.png"),
             (phiPlot, "phi.png"),
             #(dPdRplot, "dPdR.png"),
             #(dPdUplot, "dPdU.png"),
             (fDSplot, "fDS.png"),
             (Dplot, "damage.png")]

    # Save the figures.
    for p, fname in plots:
        savefig(p, os.path.join(dataDir, fname))

#-------------------------------------------------------------------------------
# Final cleanup
#-------------------------------------------------------------------------------
if postCleanup and mpi.rank == 0 and os.path.exists(dataDirBase):
    print("Removing output path ", dataDirBase)
    shutil.rmtree(dataDirBase)
