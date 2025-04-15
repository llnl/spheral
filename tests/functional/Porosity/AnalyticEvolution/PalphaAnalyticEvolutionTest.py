#-------------------------------------------------------------------------------
# Use forced terms for the inputs to the P-alpha model (pressure, derivatives,
# hydro time derivatives, etc.) to force the P-alpha distension evolution to
# follow a known analytic solution.
#-------------------------------------------------------------------------------
from Spheral1d import *
from SpheralTestUtilities import *
from NodeHistory import *
from math import *
import os, shutil
import mpi

#-------------------------------------------------------------------------------
# Identify ourselves!
#-------------------------------------------------------------------------------
title("Analytic evolution test of P-alpha porosity model")

#-------------------------------------------------------------------------------
# Build our units (cm, gm, microsec)
#-------------------------------------------------------------------------------
units = CGuS()

#-------------------------------------------------------------------------------
# Generic problem parameters
# All (cm, gm, usec) units.
#-------------------------------------------------------------------------------
commandLine(

    # Porosity
    alpha0 = 1.5,                      # Initial distention
    alphae = 1.4,                      # Elastic limiting distention
    alphat = 1.2,                      # Plastic limiting distention
    c0frac = 0.5,                      # Fraction of initial solid sound speed for initial porous sound speed
    Pe = 0.1,                          # Elastic limiting pressure
    Pt = 0.5,                          # Transition pressure in plastic regime
    Ps = 2.0,                          # A solid pressure fitting parameter, where alpha should go to 1
    n1 = 1.0,
    n2 = 1.0,

    # Fake hydro
    evolution = "sinusoidal",          # linear, sinusoidal
    DrhoDt0 = 1.0,
    DuDt0 = 2.0,

    # Initial conditions (initial density set by alpha0 and rhoS0)
    eps0 = 0.0,

    # Linear-polynomial EOS
    rhoS0 = 2.0,
    a0 = 0.0,
    a1 = 1.5,
    a2 = 0.0,
    a3 = 0.0,
    b0 = 0.75,
    b1 = 0.0,
    b2 = 0.0,
    mu = 2.0,

    IntegratorConstructor = ForwardEulerIntegrator,
    goalTime = 5.0,
    steps = None,
    sampleCycle = 1,
    dt = 1e-2,
    dtMin = 1e-2,
    dtMax = 1e-2,
    dtverbose = False,
    graphics = True,
    clearDirectories = False,
    dataDir = "dumps-PalphaAnalyticEvolution",
    outputFileName = "PorosityTimeHistory.gnu",
)

evolution = evolution.lower()
assert evolution in ("linear", "sinusoidal")

assert 1.0 <= alpha0
assert 1.0 <= alphat <= alpha0
assert c0frac <= 1.0
assert Pe <= Pt

phi0 = 1.0 - 1.0/alpha0

outputFileName = os.path.join(dataDir, outputFileName)
restartDir = os.path.join(dataDir, "restarts")
restartBaseName = os.path.join(restartDir, "PalphaAnalyticEvolution")

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
import os, sys
if mpi.rank == 0:
    if clearDirectories and os.path.exists(dataDir):
        shutil.rmtree(dataDir)
    if not os.path.exists(dataDir):
        os.makedirs(dataDir)
    if not os.path.exists(restartDir):
        os.makedirs(restartDir)
mpi.barrier()

#-------------------------------------------------------------------------------
# Material models
#-------------------------------------------------------------------------------
eosS = LinearPolynomialEquationOfState(referenceDensity = rhoS0,
                                       etamin = 0.0,
                                       etamax = 1e5,
                                       a0 = a0,
                                       a1 = a1,
                                       a2 = a2,
                                       a3 = a3,
                                       b0 = b0,
                                       b1 = b1,
                                       b2 = b2,
                                       atomicWeight = mu,
                                       constants = units)

strengthModelS = NullStrength()
eps0 = 0.0
Ps0 = eosS.pressure(rhoS0, eps0)
cS0 = eosS.soundSpeed(rhoS0, eps0)
print("Material initial properties:")
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
                          numInternal = 1,
                          nPerh = 4.01,
                          xmin = -10.0*Vector.one,
                          xmax =  10.0*Vector.one)
nodeSet = [nodes]

#-------------------------------------------------------------------------------
# Set node properties (positions, masses, H's, etc.)
#-------------------------------------------------------------------------------
pos = nodes.positions()
vel = nodes.velocity()
rho = nodes.massDensity()
eps = nodes.specificThermalEnergy()
H = nodes.Hfield()
rho[0] = rhoS0 * (1.0 - phi0)
eps[0] = eps0
H[0].xx = 1.0

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
# We create a fake hydro class to force known time derviative values
#-------------------------------------------------------------------------------
class FakeHydro(Physics):
    def __init__(self,
                 db,
                 DrhoDt0,
                 DuDt0,
                 evolution):
        Physics.__init__(self)
        assert evolution in ("linear", "sinusoidal")
        self.DrhoDt0 = DrhoDt0
        self.DuDt0 = DuDt0
        self.evolution = evolution
        self.pressure = db.newFluidScalarFieldList(0.0, HydroFieldNames.pressure)
        self.DrhoDt = db.newFluidScalarFieldList(0.0, "delta " + HydroFieldNames.massDensity)
        self.DuDt = db.newFluidScalarFieldList(0.0, "delta " + HydroFieldNames.specificThermalEnergy)
        return

    def evaluateDerivatives(self, t, dt, db, state, derivs):
        DrhoDt = derivs.scalarFields("delta " + HydroFieldNames.massDensity)
        DuDt = derivs.scalarFields("delta " + HydroFieldNames.specificThermalEnergy)
        assert (DrhoDt.numFields == 1 and DuDt.numFields == 1 and
                DrhoDt[0].numElements == 1 and DuDt[0].numElements == 1)
        if self.evolution == "linear":
            ft = 1.0
        elif self.evolution == "sinusoidal":
            ft = sin(2.0*pi*t)
        DrhoDt[0][0] = ft * self.DrhoDt0
        DuDt[0][0] = ft * self.DuDt0
        return

    def dt(self, db, state, derivs, t):
        rho = state.scalarFields(HydroFieldNames.massDensity)
        return pair_double_string(1.0, "Fake timestep")

    def registerState(self, db, state):
        rho = db.fluidMassDensity
        state.enroll(rho, ScalarIncrementBoundedState(0.0))
        state.enroll(db.fluidSpecificThermalEnergy, ScalarIncrementState())
        state.enroll(self.pressure, PressurePolicy())
        return

    def registerDerivatives(self, db, derivs):
        derivs.enroll(self.DrhoDt)
        derivs.enroll(self.DuDt)
        return

    def label(self):
        return "FakeHydro"

#-------------------------------------------------------------------------------
# Construct the fake hydro physics object.  We force known derivatives for
# the hydro state in order to drive and test the porosity evolution.
#-------------------------------------------------------------------------------
hydro = FakeHydro(db, DrhoDt0, DuDt0, evolution)
output("hydro")
output("hydro.DrhoDt0")
output("hydro.DuDt0")

packages = [hydro]

#-------------------------------------------------------------------------------
# Construct a porosity model
#-------------------------------------------------------------------------------
phi0 = 1.0 - 1.0/alpha0
c0 = c0frac * cS0

print("P-alpha porosity model parameters:")
output("  alpha0")
output("    phi0")
porosityAl = PalphaPorosity(nodes,
                            phi0 = phi0,
                            Pe = Pe,
                            Pt = Pt,
                            Ps = Ps,
                            alphae = alphae,
                            alphat = alphat,
                            n1 = n1,
                            n2 = n2,
                            cS0 = cS0,
                            c0 = c0)

output("porosityAl")
output("  porosityAl.Pe")
output("  porosityAl.Pt")
output("  porosityAl.alphae")
output("  porosityAl.alphat")
output("  porosityAl.n1")
output("  porosityAl.n2")
output("  porosityAl.cS0")
output("  porosityAl.phi().internalValues()")
packages.append(porosityAl)

#-------------------------------------------------------------------------------
# Construct a time integrator.
#-------------------------------------------------------------------------------
integrator = IntegratorConstructor(db)
for package in packages:
    integrator.appendPhysicsPackage(package)
integrator.lastDt = dt
if dtMin:
    integrator.dtMin = dtMin
if dtMax:
    integrator.dtMax = dtMax
integrator.verbose = dtverbose
output("integrator")
output("integrator.havePhysicsPackage(hydro)")
output("integrator.havePhysicsPackage(porosityAl)")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.dtGrowth")
output("integrator.domainDecompositionIndependent")

#-------------------------------------------------------------------------------
# Use a NodeHistory object to track the time evolution of the porosity and other
# state.
#-------------------------------------------------------------------------------
def sampleState(nodes, nodeIndices):
    rho = nodes.massDensity()
    eps = nodes.specificThermalEnergy()
    alpha = porosityAl.alpha
    dPdR = porosityAl.partialPpartialRho
    dPdu = porosityAl.partialPpartialEps
    return alpha[0], 1.0 - 1.0/alpha[0], rho[0], eps[0], hydro.pressure(0,0), dPdR[0], dPdu[0]
    
history = NodeHistory(nodes,
                      nodeIndices = [0],
                      sampleMethod = sampleState,
                      filename = outputFileName,
                      labels = ("alpha", "phi", "rho", "eps", "P", "dPdR", "dPdu"))

#-------------------------------------------------------------------------------
# Build the controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            restartStep = None,
                            iterateInitialH = False,
                            periodicWork = [(history, sampleCycle)],
                            restartBaseName = restartBaseName)
output("control")

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
if not steps is None:
    control.step(steps)
else:
    control.advance(goalTime)
    control.dropRestartFile()

#-------------------------------------------------------------------------------
# Plot the state.
#-------------------------------------------------------------------------------
if graphics:
    from matplotlib import pyplot as plt
    def makePlot(x, y, marker, xlabel=None, ylabel=None, title=None):
        result = plt.figure().add_subplot(111)
        plt.plot(x, y, marker)
        if title:    plt.title(title)
        if xlabel:   plt.xlabel(xlabel)
        if ylabel:   plt.ylabel(ylabel)
        return result

    # def alpha_of_P(P):
    #     "The analytic expectation for how alpha(P) should vary"
    #     if P <= Pe:
    #         return 1.0 + (alphae - 1.0)*((Ps - P)/(

    t = history.timeHistory
    alphaPlot = makePlot(history.timeHistory, [x[0] for x in history.sampleHistory], "r-", r"$t$", r"$\alpha$", r"$\alpha$ history")
    phiPlot = makePlot(history.timeHistory, [x[1] for x in history.sampleHistory], "r-", r"$t$", r"$\phi$", r"$\phi$ history")
    rhoPlot = makePlot(history.timeHistory, [x[2] for x in history.sampleHistory], "r-", r"$t$", r"$\rho$", r"Bulk $\rho$ history")
    epsPlot = makePlot(history.timeHistory, [x[3] for x in history.sampleHistory], "r-", r"$t$", r"$\varepsilon$", r"$\varepsilon$ history")
    PPlot = makePlot(history.timeHistory, [x[4] for x in history.sampleHistory], "r-", r"$t$", r"$P$", r"Bulk $P$ history")
    dPdRPlot = makePlot(history.timeHistory, [x[5] for x in history.sampleHistory], "r-", r"$t$", r"$\partial P/\partial \rho$", r"Solid $\partial P/\partial \rho$ history")
    dPduPlot = makePlot(history.timeHistory, [x[6] for x in history.sampleHistory], "r-", r"$t$", r"$\partial P/\partial \varepsilon$", r"Solid $\partial P/\partial \varepsilon$ history")
    alphaVSPplot = makePlot([x[4] for x in history.sampleHistory],
                            [x[0] for x in history.sampleHistory], "ro", r"$P$", r"$\alpha$", r"$\alpha(P)$")

    plots = [(alphaPlot, "alpha.png"),
             (phiPlot, "phi.png"),
             (rhoPlot, "rho.png"),
             (epsPlot, "eps.png"),
             (PPlot, "pressure.png"),
             (dPdRPlot, "partialPpartialRho.png"),
             (dPduPlot, "partialPpartialEps.png"),
             (alphaVSPplot, "alpha_vs_pressure.png")]

    # Save the figures.
    for p, fname in plots:
        p.figure.savefig(os.path.join(dataDir, fname))

    plt.show()
