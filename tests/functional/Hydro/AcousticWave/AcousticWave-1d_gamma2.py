#-------------------------------------------------------------------------------
# A 1-D acoustic wave test.  Here we propogate a simple sound wave round and
# round in a periodic box.  This specific example is based on the test case
# described in D.J. Price's dissertation as an example of the effect of the
# grad h terms.
#-------------------------------------------------------------------------------
import os, shutil
from math import *
from Spheral1d import *
from SpheralTestUtilities import *
import mpi
import numpy as np
#import matplotlib.pyplot as plt

def smooth(x,window_len=11,window='hanning'):
    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")
    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")
    if window_len<3:
        return x
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")
    s=np.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')
    y=np.convolve(w/w.sum(),s,mode='same')
    return y[window_len:-window_len+1]

title("Acoustic wave propagation test.")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(nx1 = 100,
            x0 = 0.0,
            x1 = 1.0,

            rho1 = 1.0,
            eps1 = 1.0,
            A = 1.0e-6,
            kfreq = 1.0,

            cs2 = 1.0,
            mu = 1.0,
            gamma = 5.0/3.0,
            P1 = 3.0/5.0,

            nPerh = 1.51,

            Qconstructor = MonaghanGingoldViscosity,
            #Qconstructor = TensorMonaghanGingoldViscosity,
            Cl = 0.0,
            Cq = 0.0,
            linearInExpansion = False,
            Qlimiter = False,
            epsilon2 = 1e-30,
            hmin = 1.0e-10,
            hmax = 0.1,
            cfl = 0.5,
            XSPH = False,
            epsilonTensile = 0.0,
            nTensile = 4,
            filter = 0.0,
            KernelConstructor = BSplineKernel,
            order = 5,

            SVPH = False,
            CRKSPH = False,
            TSPH = False,
            PSPH = False,
            IntegratorConstructor = CheapSynchronousRK2Integrator,
            steps = None,
            goalTime = 5.0,
            dt = 1.0e-10,
            dtMin = 1.0e-10, 
            dtMax = 0.1,
            dtGrowth = 2.0,
            dtverbose = False,
            rigorousBoundaries = False,
            maxSteps = None,
            statsStep = 1,
            smoothIters = 0,
            HUpdate = IntegrateH,
            correctionOrder = LinearOrder,
            densityUpdate = RigorousSumDensity,
            compatibleEnergy = True,
            gradhCorrection = True,
            correctVelocityGradient = True,
            linearConsistent = False,

            restoreCycle = None,
            restartStep = 10000,
            clearDirectories = True,
            dataDirBase = "dumps-planar-AcousticWave-1d",
            outputFile = "AcousticWave-planar-1d.gnu",
            normOutputFile = None,
            writeOutputLabel = True,

            graphics = "gnu",

            checkReversibility = False,
            )

if SVPH:
    HydroConstructor = SVPHFacetedHydro
elif CRKSPH:
    HydroConstructor = CRKSPHHydro
    Qconstructor = LimitedMonaghanGingoldViscosity
elif TSPH:
    HydroConstructor = TaylorSPHHydro
elif PSPH:
    HydroConstructor = PSPHHydro
else:
    HydroConstructor = SPHHydro

dataDir = os.path.join(dataDirBase,
                       str(HydroConstructor).split("'")[1].split(".")[-1],
                       str(Qconstructor).split("'")[1].split(".")[-1],
                       "nx=%i" % nx1)
restartDir = os.path.join(dataDir, "restarts")
restartBaseName = os.path.join(restartDir, "AcousticWave-planar-1d-%i" % nx1)

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
#eos = IsothermalEquationOfStateMKS(cs2, mu)
eos = GammaLawGasMKS(gamma, mu)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
if KernelConstructor==NBSplineKernel:
  WT = TableKernel(NBSplineKernel(order), 1000)
else:
  WT = TableKernel(KernelConstructor(), 1000)
output("WT")
kernelExtent = WT.kernelExtent
output("WT")

#-------------------------------------------------------------------------------
# Make the NodeList.
#-------------------------------------------------------------------------------
nodes1 = makeFluidNodeList("nodes1", eos,
                           hmin = hmin,
                           hmax = hmax,
                           kernelExtent = kernelExtent,
                           nPerh = nPerh)
output("nodes1")
output("nodes1.hmin")
output("nodes1.hmax")
output("nodes1.nodesPerSmoothingScale")

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
from DistributeNodes import distributeNodesInRange1d
distributeNodesInRange1d([(nodes1, nx1, rho1, (x0, x1))],
                         nPerh = nPerh)
nNodesThisDomain1 = nodes1.numInternalNodes
output("nodes1.numNodes")

# Find the cumulative mass at each point.
Mi = ScalarField("Cumulative mass", nodes1)
positions = mpi.allreduce([(nodes1.positions()[i].x, i, mpi.rank)
                           for i in range(nodes1.numInternalNodes)], mpi.SUM)
assert len(positions) == nx1
positions.sort()
Msum = 0.0
mi = rho1/nx1
for (x, i, proc) in positions:
    Msum += mi
    if proc == mpi.rank:
        assert i < nodes1.numInternalNodes
        Mi[i] = Msum
assert fuzzyEqual(Msum, rho1)

# Define the function which we are going to solve for the node positions.
twopi = 2.0*pi
class MassFunctor(PairScalarFunctor):
    def __init__(self, Mcumulative):
        PairScalarFunctor.__init__(self)
        self.Mcumulative = Mcumulative
        return
    def __call__(self, x):
        return pair_double_double(self.Mcumulative - rho1*(x + A/(twopi*kfreq)*(1.0 - cos(twopi*kfreq*x))),
                                  -rho1*(1.0 + A*sin(twopi*kfreq*x)))

# Set the node positions, velocities, and densities.
from newtonRaphson import *
cs = sqrt(cs2)
pos = nodes1.positions()
vel = nodes1.velocity()
rho = nodes1.massDensity()
eps = nodes1.specificThermalEnergy()
mass = nodes1.mass()
H = nodes1.Hfield()
dx = (x1 - x0)/nx1
xi = x0
for i in range(nodes1.numInternalNodes):
    func0 = MassFunctor(max(0.0, Mi[i] - mi))
    func1 = MassFunctor(Mi[i])
    xi0 = newtonRaphsonFindRoot(func0, xi, xi + 2.0*dx, 1.0e-18, 1.0e-18)
    xi1 = newtonRaphsonFindRoot(func1, xi, xi + 2.0*dx, 1.0e-18, 1.0e-18)
    rhoi0 = rho1*(1.0 + A*sin(twopi*kfreq*(xi0 - x0)/(x1 - x0)))
    rhoi1 = rho1*(1.0 + A*sin(twopi*kfreq*(xi1 - x0)/(x1 - x0)))
    xi = x0 + (x1 - x0)*(rhoi0*xi0 + rhoi1*xi1)/(rhoi0 + rhoi1)
    pos[i].x = xi
    vel[i].x = A*cs*sin(twopi*kfreq*(xi - x0)/(x1 - x0))
    rho[i] = rho1*(1.0 + A*sin(twopi*kfreq*(xi - x0)/(x1 - x0)))
    mass[i] = rho1*((xi1 - xi0) - A/(twopi*kfreq)*(cos(twopi*kfreq*xi1) - cos(twopi*kfreq*xi0)))
    Pi = P1+A*sin(twopi*kfreq*(xi - x0)/(x1 - x0))
    eps[i] = Pi/(((5.0/3.0) - 1.0)*rho[i])
    H[i] *= rho[i]/rho1
# xi0 = 0.0
# dx0 = (x1 - x0)/nx1
# for i in xrange(nodes1.numInternalNodes):
#     dxi0 = dx0*(1.0 - A*sin(twopi*kfreq*(xi0 - x0)/(x1 - x0)))
#     xi = xi0 + 0.5*dxi0
#     pos[i].x = xi
#     vel[i].x = A*cs*sin(twopi*kfreq*(xi - x0)/(x1 - x0))
#     rho[i] = rho1*(1.0 + A*sin(twopi*kfreq*(xi - x0)/(x1 - x0)))
#     xi0 += dxi0

# Compute the summation correction for the density, and apply it to the mass per point.
m0 = rho1*dx
Hdet0 = 1.0/(nPerh*dx)
rhoscale = m0*WT.kernelValue(0.0, Hdet0)
deta = 1.0/nPerh
for i in range(1, int(WT.kernelExtent * (nPerh + 1))):
    rhoscale += 2.0*m0*WT.kernelValue(i*deta, Hdet0)
rhoscale = rho1/rhoscale
print("Compute analytic rho scaling of %16.12e." % rhoscale)
for i in range(nodes1.numInternalNodes):
    mass[i] *= rhoscale

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase()
output("db")
output("db.appendNodeList(nodes1)")
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Construct the artificial viscosity.
#-------------------------------------------------------------------------------
q = Qconstructor(Cl, Cq, linearInExpansion = linearInExpansion)
q.epsilon2 = epsilon2
q.limiter = Qlimiter
output("q")
output("q.Cl")
output("q.Cq")
output("q.epsilon2")
output("q.limiter")

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
if SVPH:
    hydro = HydroConstructor(W = WT, 
                             Q = q,
                             cfl = cfl,
                             compatibleEnergyEvolution = compatibleEnergy,
                             XSVPH = XSPH,
                             linearConsistent = linearConsistent,
                             densityUpdate = densityUpdate,
                             HUpdate = HUpdate,
                             xmin = Vector(-100.0),
                             xmax = Vector( 100.0))
elif CRKSPH:
    hydro = HydroConstructor(W = WT, 
                             Q = q,
                             filter = filter,
                             cfl = cfl,
                             compatibleEnergyEvolution = compatibleEnergy,
                             correctionOrder = correctionOrder,
                             XSPH = XSPH,
                             densityUpdate = densityUpdate,
                             HUpdate = HUpdate)

elif TSPH:
    hydro = HydroConstructor(W = WT, 
                             Q = q,
                             cfl = cfl,
                             compatibleEnergyEvolution = compatibleEnergy,
                             XSPH = XSPH,
                             HUpdate = HUpdate)
else:
    hydro = HydroConstructor(W = WT, 
                             Q = q,
                             cfl = cfl,
                             compatibleEnergyEvolution = compatibleEnergy,
                             gradhCorrection = gradhCorrection,
                             correctVelocityGradient = correctVelocityGradient,
                             XSPH = XSPH,
                             densityUpdate = densityUpdate,
                             HUpdate = HUpdate,
                             epsTensile = epsilonTensile,
                             nTensile = nTensile)
output("hydro")

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
xPlane0 = Plane(Vector(x0), Vector( 1.0))
xPlane1 = Plane(Vector(x1), Vector(-1.0))
xbc = PeriodicBoundary(xPlane0, xPlane1)
hydro.appendBoundary(xbc)

#-------------------------------------------------------------------------------
# Construct a time integrator.
#-------------------------------------------------------------------------------
integrator = IntegratorConstructor(db)
integrator.appendPhysicsPackage(hydro)
integrator.lastDt = dt
integrator.dtMin = dtMin
integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
integrator.rigorousBoundaries = rigorousBoundaries
integrator.verbose = dtverbose
output("integrator")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.dtGrowth")
output("integrator.rigorousBoundaries")

#-------------------------------------------------------------------------------
# Make the problem controller.
#-------------------------------------------------------------------------------
print("Making controller.")
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle)
output("control")

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
if steps is None:
    if control.time() < goalTime:
        control.advance(goalTime, maxSteps)
    if checkReversibility:
        for i in range(nodes1.numNodes):
            vel[i] = -vel[i]
        control.advance(2*goalTime, maxSteps)
else:
    control.step(steps)

#-------------------------------------------------------------------------------
# Compute the analytic answer.
#-------------------------------------------------------------------------------
import AcousticWaveSolution
xlocal = [pos.x for pos in nodes1.positions().internalValues()]
xprof = mpi.reduce(xlocal, mpi.SUM)
dx = (x1 - x0)/nx1
h1 = nPerh*dx
answer = AcousticWaveSolution.AcousticWaveSolution(eos, cs, rho1, x0, x1, A, twopi*kfreq, h1)
#print "\n\nPERIOD=",1.0/(kfreq*cs)

### Compute the simulated specific entropy.
##rho = mpi.allreduce(nodes1.massDensity().internalValues(), mpi.SUM)
##P = mpi.allreduce(nodes1.pressure().internalValues(), mpi.SUM)
##A = [Pi/rhoi**gamma for (Pi, rhoi) in zip(P, rho)]

### The analytic solution for the simulated entropy.
##xans, vans, uans, rhoans, Pans, hans = answer.solution(control.time(), xprof)
##Aans = [Pi/rhoi**gamma for (Pi, rhoi) in zip(Pans,  rhoans)]

#-------------------------------------------------------------------------------
# Plot the final state.
#-------------------------------------------------------------------------------
if graphics == "gnu":
    from SpheralGnuPlotUtilities import *
    state = State(db, integrator.physicsPackages())
    rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotState(state)
    if mpi.rank == 0:
        plotAnswer(answer, control.time(), rhoPlot, velPlot, epsPlot, PPlot, HPlot, xprof)
    cs = state.scalarFields(HydroFieldNames.soundSpeed)
    csPlot = plotFieldList(cs, winTitle="Sound speed", colorNodeLists=False)
    EPlot = plotEHistory(control.conserve)

    # Plot the correction terms.

    # Plot the grad h correction term (omega)

    if SVPH:
        volPlot = plotFieldList(hydro.volume(),
                                winTitle = "volume",
                                colorNodeLists = False)
    elif CRKSPH:
        APlot = plotFieldList(hydro.A(),
                              winTitle = "A",
                              colorNodeLists = False)
        BPlot = plotFieldList(hydro.B(),
                              yFunction = "%s.x",
                              winTitle = "B",
                              colorNodeLists = False)

    else:
        omegaPlot = plotFieldList(hydro.omegaGradh(),
                                  winTitle = "grad h correction",
                                  colorNodeLists = False)
Eerror = (control.conserve.EHistory[-1] - control.conserve.EHistory[0])/control.conserve.EHistory[0]
print("Total energy error: %g" % Eerror)

#-------------------------------------------------------------------------------
# If requested, write out the state in a global ordering to a file.
#-------------------------------------------------------------------------------
if outputFile:
    outputFile = os.path.join(dataDir, outputFile)
    from SpheralTestUtilities import multiSort
    mprof = mpi.reduce(nodes1.mass().internalValues(), mpi.SUM)
    rhoprof = mpi.reduce(nodes1.massDensity().internalValues(), mpi.SUM)
    P = ScalarField("pressure", nodes1)
    nodes1.pressure(P)
    Pprof = mpi.reduce(P.internalValues(), mpi.SUM)
    vprof = mpi.reduce([v.x for v in nodes1.velocity().internalValues()], mpi.SUM)
    epsprof = mpi.reduce(nodes1.specificThermalEnergy().internalValues(), mpi.SUM)
    hprof = mpi.reduce([1.0/H.xx for H in nodes1.Hfield().internalValues()], mpi.SUM)
    xans, vans, uans, rhoans, Pans, hans = answer.solution(control.time(), xprof)

    labels = ["x", "m", "rho", "P", "v", "eps", "h", 
              "rhoans", "Pans", "vans", "epsans", "hans"]
    stuff = [xprof, mprof, rhoprof, Pprof, vprof, epsprof, hprof, 
             rhoans, Pans, vans, uans, hans]
    if CRKSPH:
        Aprof = mpi.reduce(hydro.A()[0].internalValues(), mpi.SUM)
        Bprof = mpi.reduce([x.x for x in hydro.B()[0].internalValues()], mpi.SUM)
        labels += ["A", "B"]
        stuff += [Aprof, Bprof]

    if mpi.rank == 0:
        multiSort(*tuple(stuff))
        f = open(outputFile, "w")
        f.write(("#  " + len(labels)*"'%s' " + "\n") % tuple(labels))
        for tup in zip(*tuple(stuff)):
            assert len(tup) == len(labels)
            f.write((len(tup)*"%16.12e " + "\n") % tup)
        f.close()

        # While we're at it compute and report the error norms.
        import Pnorm
        print("\tQuantity \t\tL1 \t\t\tL2 \t\t\tLinf")
        if normOutputFile:
            f = open(normOutputFile, "a")
            if writeOutputLabel:
                f.write(("#" + 13*"%17s " + "\n") % ('"nx"',
                                                     '"rho L1"', '"rho L2"', '"rho Linf"',
                                                     '"P L1"',   '"P L2"',   '"P Linf"',
                                                     '"vel L1"', '"vel L2"', '"vel Linf"',
                                                     '"h L1"',   '"h L2"',   '"h Linf"'))
            f.write("%16i " % nx1)
        xmin, xmax = x0, x1
        for (name, data, ans) in [("Mass Density", rhoprof, rhoans),
                                  ("Pressure", Pprof, Pans),
                                  ("Velocity", vprof, vans),
                                  ("h       ", hprof, hans)]:
            assert len(data) == len(ans)
            error = [data[i] - ans[i] for i in range(len(data))]
            Pn = Pnorm.Pnorm(error, xprof)
            L1 = Pn.pnormAverage(1, xmin, xmax)
            L2 = Pn.pnormAverage(2, xmin, xmax)
            Linf = Pn.pnormAverage("inf", xmin, xmax)
            print("\t%s \t\t%g \t\t%g \t\t%g" % (name, L1, L2, Linf))
            if normOutputFile:
                f.write((3*"%16.12e ") % (L1, L2, Linf))
        if normOutputFile:
            f.write("\n")
            f.close()

if compatibleEnergy and abs(Eerror) > 1e-5:
    raise ValueError("Energy error outside allowed bounds.")
