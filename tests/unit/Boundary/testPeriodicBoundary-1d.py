#ATS:t0 = test(SELF, "--graphics None",       label="Periodic boundary unit test -- 1-D (serial)")
#ATS:t1 = test(SELF, "--graphics None", np=2, label="Periodic boundary unit test -- 1-D (parallel)")
#-------------------------------------------------------------------------------
# 1D test of periodic boundaries -- we simply allow a pressureless fluid to
# cycle around a box and check the sum density 
#-------------------------------------------------------------------------------
from math import *
from Spheral1d import *
from SpheralTestUtilities import *
import mpi

title("1D periodic boundary test.")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(nx1 = 100,
            x0 = 0.0,
            x1 = 1.0,

            rho1 = 1.0,
            cs2 = 1.0,
            mu = 1.0,
            vx1 = 1.0,

            nPerh = 2.01,

            hmin = 0.0001, 
            hmax = 0.1,
            cfl = 0.5,

            tol = 1.0e-10,
            steps = 500,
            dt = 0.0001,
            dtMin = 1.0e-5, 
            dtMax = 0.1,
            dtGrowth = 2.0,
            dtverbose = False,
            rigorousBoundaries = False,
            maxSteps = None,
            statsStep = 1,
            smoothIters = 0,
            HEvolution = IdealH,
            densityUpdate = RigorousSumDensity,
            compatibleEnergy = True,
            gradhCorrection = True,
            linearConsistent = False,

            restoreCycle = None,
            restartStep = 10000,
            restartBaseName = "dumps-AcousticWave-1d",

            graphics = True,
            )

#-------------------------------------------------------------------------------
# Material properties.
#-------------------------------------------------------------------------------
eos = IsothermalEquationOfStateMKS(cs2, mu)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
WT = TableKernel(BSplineKernel(), 100)

#-------------------------------------------------------------------------------
# Make the NodeList.
#-------------------------------------------------------------------------------
nodes1 = makeFluidNodeList("nodes1", eos,
                           hmin = hmin,
                           hmax = hmax,
                           nPerh = nPerh)

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
from DistributeNodes import distributeNodesInRange1d
distributeNodesInRange1d([(nodes1, nx1, rho1, (x0, x1))],
                         nPerh = nPerh)
nNodesThisDomain1 = nodes1.numInternalNodes

# Set the node positions, velocities, and densities.
nodes1.velocity(VectorField("tmp velocity", nodes1, Vector(vx1)))

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase()
db.appendNodeList(nodes1)

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
hydro = SPH(dataBase = db,
            W = WT, 
            cfl = cfl,
            densityUpdate = densityUpdate,
            HUpdate = HEvolution)

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
integrator = CheapSynchronousRK2Integrator(db)
integrator.appendPhysicsPackage(hydro)
integrator.lastDt = dt
integrator.dtMin = dtMin
integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
integrator.rigorousBoundaries = rigorousBoundaries
integrator.verbose = dtverbose

#-------------------------------------------------------------------------------
# Function to track the density and raise the alarm if it doesn't stay within
# bounds.
#-------------------------------------------------------------------------------
def checkRho(steps, t, dt):
    rho = nodes1.massDensity()
    rhoMin = rho.min()
    rhoMax = rho.max()
    #print "Rho range : [%16.12e, %16.12e]" % (rhoMin, rhoMax)
    if abs(rhoMin/rhoMax - 1.0) > tol:
        if graphics:
            from SpheralMatplotlib import plotState
            state = State(db, integrator.physicsPackages())
            rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotState(state, plotGhosts=True)
        pos = nodes1.positions()
        for i in range(nodes1.numInternalNodes):
            if rho[i] == rhoMin:
                sys.stderr.write("rho min @ %i %s\n" % (i, pos[i]))
            if rho[i] == rhoMax:
                sys.stderr.write("rho max @ %i %s\n" % (i, pos[i]))
        raise ValueError("rho outside bounds : [%16.12e, %16.12e]" % (rhoMin, rhoMax))
    return

#-------------------------------------------------------------------------------
# Make the problem controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle)
control.appendPeriodicWork(checkRho, 1)

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
control.step(steps)
print("** PASS **")

#-------------------------------------------------------------------------------
# Plot the final state.
#-------------------------------------------------------------------------------
if graphics:
    from SpheralMatplotlib import *
    state = State(db, integrator.physicsPackages())
    rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotState(state, plotGhosts=True)
