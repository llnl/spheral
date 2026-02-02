#-------------------------------------------------------------------------------
# A 1-D acoustic wave test.  Here we propogate a simple sound wave round and
# round in a periodic box.  
# This version has been adapted as a 1D test bed for some hourglass control
# ideas for use with CRKSPH.
#-------------------------------------------------------------------------------
import sys, mpi
from math import *
from Spheral1d import *
from SpheralTestUtilities import *
import numpy as np

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
            A = 0.005,
            kfreq = 1.0,

            cs2 = 1.0,
            mu = 1.0,

            nPerh = 2.01,

            ranfrac = 0.01,  # How much to randomly perturb intial condition

            Qconstructor = MonaghanGingoldViscosity,
            #Qconstructor = TensorMonaghanGingoldViscosity,
            Cl = 0.0,
            Cq = 0.0,
            Qlimiter = False,
            epsilon2 = 1e-2,
            hmin = 0.0001, 
            hmax = 0.1,
            cfl = 0.5,
            XSPH = False,
            epsilonTensile = 0.0,
            nTensile = 4,
            filter = 0.01,

            SVPH = False,
            CRKSPH = True,
            TSPH = False,
            IntegratorConstructor = CheapSynchronousRK2Integrator,
            steps = None,
            goalTime = 5.0,
            dt = 0.0001,
            dtMin = 1.0e-5, 
            dtMax = 0.1,
            dtGrowth = 2.0,
            rigorousBoundaries = False,
            maxSteps = None,
            statsStep = 1,
            smoothIters = 0,
            HEvolution = IdealH,
            densityUpdate = IntegrateDensity,
            compatibleEnergy = False,
            gradhCorrection = True,
            linearConsistent = False,

            hourglassIters = 100,

            restoreCycle = None,
            restartStep = 10000,
            restartBaseName = "dumps-AcousticWave-1d",

            graphics = "gnu",
            )

#-------------------------------------------------------------------------------
# Material properties.
#-------------------------------------------------------------------------------
eos = IsothermalEquationOfStateMKS(cs2, mu)

##gamma = 5.0/3.0
##eps1 = cs2/(gamma*(gamma - 1.0))
##eos = GammaLawGasMKS(gamma, mu)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
WT = TableKernel(BSplineKernel(), 1000)
WTPi = TableKernel(BSplineKernel(), 1000)
output("WT")
output("WTPi")

#-------------------------------------------------------------------------------
# Make the NodeList.
#-------------------------------------------------------------------------------
nodes1 = makeFluidNodeList("nodes1", eos,
                           hmin = hmin,
                           hmax = hmax,
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
dx = 1.0/nx1
import random
from newtonRaphson import *
cs = sqrt(cs2)
pos = nodes1.positions()
vel = nodes1.velocity()
rho = nodes1.massDensity()
for i in range(nodes1.numInternalNodes):
    func0 = MassFunctor(max(0.0, Mi[i] - mi))
    func1 = MassFunctor(Mi[i])
    xi0 = newtonRaphsonFindRoot(func0, 0.0, 1.0, 1.0e-15, 1.0e-15)
    xi1 = newtonRaphsonFindRoot(func1, 0.0, 1.0, 1.0e-15, 1.0e-15)
    xi = x0 + (x1 - x0)*0.5*(xi0 + xi1)
    pos[i].x = xi + ranfrac*dx*random.uniform(-1.0, 1.0)
    vel[i].x = 0.5*(A*cs*sin(twopi*kfreq*(xi0 - x0)/(x1 - x0)) +
                    A*cs*sin(twopi*kfreq*(xi1 - x0)/(x1 - x0)))
    rho[i] = rho1*0.5*((1.0 + A*sin(twopi*kfreq*(xi0 - x0)/(x1 - x0))) +
                       (1.0 + A*sin(twopi*kfreq*(xi1 - x0)/(x1 - x0))))

##    # Set the specific thermal energy.
##    P1 = cs2*rhoi
##    nodes1.specificThermalEnergy()[i] = P1/((gamma - 1.0)*rhoi)

#nodes1.specificThermalEnergy(ScalarField("tmp", nodes1, eps1))

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
q = Qconstructor(Cl, Cq)
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
    hydro = SVPHFacetedHydro(WT, q,
                             cfl = cfl,
                             compatibleEnergyEvolution = compatibleEnergy,
                             XSVPH = XSPH,
                             linearConsistent = linearConsistent,
                             densityUpdate = densityUpdate,
                             HUpdate = HEvolution,
                             xmin = Vector(-100.0),
                             xmax = Vector( 100.0))
elif CRKSPH:
    hydro = CRKSPHHydro(WT, WTPi, q,
                      filter = filter,
                      cfl = cfl,
                      compatibleEnergyEvolution = compatibleEnergy,
                      XSPH = XSPH,
                      densityUpdate = densityUpdate,
                      HUpdate = HEvolution)

elif TSPH:
    hydro = TaylorSPHHydro(WT, q,
                           cfl = cfl,
                           compatibleEnergyEvolution = compatibleEnergy,
                           XSPH = XSPH,
                           HUpdate = HEvolution)
else:
    hydro = SPHHydro(WT, WTPi, q,
                     cfl = cfl,
                     compatibleEnergyEvolution = compatibleEnergy,
                     gradhCorrection = gradhCorrection,
                     XSPH = XSPH,
                     densityUpdate = densityUpdate,
                     HUpdate = HEvolution,
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
else:
    control.step(steps)

#-------------------------------------------------------------------------------
# Compute the analytic answer.
#-------------------------------------------------------------------------------
sys.path.append("../../../tests/Hydro/AcousticWave")
import AcousticWaveSolution
xlocal = [pos.x for pos in nodes1.positions().internalValues()]
xglobal = mpi.reduce(xlocal, mpi.SUM)
dx = (x1 - x0)/nx1
h1 = 1.0/(nPerh*dx)
answer = AcousticWaveSolution.AcousticWaveSolution(eos, cs, rho1, x0, x1, A, twopi*kfreq, h1)

### Compute the simulated specific entropy.
##rho = mpi.allreduce(nodes1.massDensity().internalValues(), mpi.SUM)
##P = mpi.allreduce(nodes1.pressure().internalValues(), mpi.SUM)
##A = [Pi/rhoi**gamma for (Pi, rhoi) in zip(P, rho)]

### The analytic solution for the simulated entropy.
##xans, vans, uans, rhoans, Pans, hans = answer.solution(control.time(), xglobal)
##Aans = [Pi/rhoi**gamma for (Pi, rhoi) in zip(Pans,  rhoans)]

#-------------------------------------------------------------------------------
# Plot the final state.
#-------------------------------------------------------------------------------
if graphics == "gnu":

    from momentHourglass import *
    from positionHourglass import *
    for i in range(hourglassIters):
        print("Iteration %i: %i" % (i, positionHourglass(db, WT, integrator.uniqueBoundaryConditions(), minfrac=0.0, maxfrac=1.0)))
        control.iterateIdealH()
    control.step()

    from SpheralGnuPlotUtilities import *
    state = State(db, integrator.physicsPackages())
    rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotState(state)
    if mpi.rank == 0:
        plotAnswer(answer, control.time(), rhoPlot, velPlot, epsPlot, PPlot, HPlot, xglobal)
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
        A0=hydro.A0()
	print("ARRAY LENGTH:")
        print((A0[0].__len__()))
        tmp=[]
        for i in range(A0[0].__len__()):
		tmp.append(A0[0][i])
        A=np.array(tmp)
        #ret=smooth(A,11,'hamming') 
        CoeffBx = Gnuplot.Data(A,
                               with_ = "points",
                               #with_ = "lines",
                               title = "Bx",
                               inline = True)
        p0 = generateNewGnuPlot()
        p0.plot(CoeffBx)
        p0.title("COEFF")
        p0.refresh()
        #print(A0.size())
        #A=np.array(A0)
        #ret=smooth(A,11,'hamming')
        volPlot = plotFieldList(hydro.volume(),
                                winTitle = "volume",
                                colorNodeLists = False)
        A0Plot = plotFieldList(hydro.A0(),
                               winTitle = "A0",
                               colorNodeLists = False)
        APlot = plotFieldList(hydro.A(),
                              winTitle = "A",
                              colorNodeLists = False)
        BPlot = plotFieldList(hydro.B(),
                              yFunction = "%s.x",
                              winTitle = "B",
                              colorNodeLists = False)

        db.updateConnectivityMap()
        cm = db.connectivityMap()
        m0 = hydro.m0()
        for bc in (xbc,):
            bc.applyFieldListGhostBoundary(m0)
            bc.finalizeGhostBoundary()
        m0smooth = interpolateCRKSPH(m0, db.fluidPosition, db.fluidMass, db.fluidHfield,
                                   hydro.A(), hydro.B(), cm, WT)
        m0var = db.newFluidScalarFieldList(0.0, "m0 variation")
        for i in range(nodes1.numInternalNodes):
            m0var[0][i] = abs(m0smooth[0][i] - m0[0][i])/m0[0][i]
        m0varPlot = plotFieldList(m0var,
                                  winTitle = "abs(<m0> - m0)/m0",
                                  plotStyle = "lines")

    else:
        omegaPlot = plotFieldList(hydro.omegaGradh(),
                                  winTitle = "grad h correction",
                                  colorNodeLists = False)

Eerror = (control.conserve.EHistory[-1] - control.conserve.EHistory[0])/control.conserve.EHistory[0]
print("Total energy error: %g" % Eerror)
if abs(Eerror) > 1e-13:
    raise ValueError("Energy error outside allowed bounds.")
