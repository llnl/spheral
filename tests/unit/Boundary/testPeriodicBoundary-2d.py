#ATS:t0 = test(SELF, "",       label="Periodic boundary unit test -- 2-D (serial)")
#ATS:t1 = test(SELF, "", np=2, label="Periodic boundary unit test -- 2-D (parallel)")
#-------------------------------------------------------------------------------
# 2D test of periodic boundaries -- we simply allow a pressureless fluid to
# cycle around a box and check the sum density 
#-------------------------------------------------------------------------------
from math import *
from Spheral2d import *
from SpheralTestUtilities import *
from SpheralPointmeshSiloDump import dumpPhysicsState
import mpi

title("2D periodic boundary test.")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(nx1 = 20,
            ny1 = 20,
            x0 = 0.0,
            x1 = 1.0,
            y0 = 0.0,
            y1 = 1.0,

            rho1 = 1.0,
            cs2 = 1.0,
            mu = 1.0,
            vx1 = 1.0,
            vy1 = 1.0,

            nPerh = 2.01,

            hmin = 0.0001, 
            hmax = 0.5,
            cfl = 0.5,

            tol = 1.0e-3,
            steps = 300,
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
            domainIndependent = False,

            restoreCycle = None,
            restartStep = 10000,
            restartBaseName = "dumps-AcousticWave-2d",
            )

#-------------------------------------------------------------------------------
# Material properties.
#-------------------------------------------------------------------------------
eos = IsothermalEquationOfStateMKS(cs2, mu)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
WT = TableKernel(BSplineKernel(), 1000)
WTPi = TableKernel(BSplineKernel(), 1000)

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
from GenerateNodeDistribution2d import GenerateNodeDistribution2d
gen1 = GenerateNodeDistribution2d(nx1, ny1,
                                  rho = rho1,
                                  distributionType = "lattice",
                                  xmin = (x0, y0),
                                  xmax = (x1, y1),
                                  nNodePerh = nPerh,
                                  SPH = True)
if mpi.procs > 1:
    from PeanoHilbertDistributeNodes import distributeNodes2d
else:
    from DistributeNodes import distributeNodes2d
distributeNodes2d((nodes1, gen1))

# Set the node positions, velocities, and densities.
nodes1.velocity(VectorField("tmp velocity", nodes1, Vector(vx1, vy1)))

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
            densityUpdate = RigorousSumDensity,
            HUpdate = HEvolution)

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
xPlane0 = Plane(Vector(x0,y0), Vector( 1.0,  0.0))
xPlane1 = Plane(Vector(x1,y1), Vector(-1.0,  0.0))
yPlane0 = Plane(Vector(x0,y0), Vector( 0.0,  1.0))
yPlane1 = Plane(Vector(x1,y1), Vector( 0.0, -1.0))
xbc = PeriodicBoundary(xPlane0, xPlane1)
ybc = PeriodicBoundary(yPlane0, yPlane1)
hydro.appendBoundary(xbc)
hydro.appendBoundary(ybc)

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
integrator.domainDecompositionIndependent = domainIndependent
integrator.verbose = dtverbose

#-------------------------------------------------------------------------------
# Function to track the density and raise the alarm if it doesn't stay within
# bounds.
#-------------------------------------------------------------------------------
def checkRho(steps, t, dt):
    rho = nodes1.massDensity()
    rhoMin = rho.min()
    rhoMax = rho.max()
    fluc = abs(rhoMin/rhoMax - 1.0)
    print("rho bounds : [%g, %g], variation %g" % (rhoMin, rhoMax, fluc))
    if fluc > tol:
        dumpPhysicsState(integrator,
                         baseFileName = "periodic_boundary_2d",
                         fieldLists = [db.fluidMass, db.fluidMassDensity, db.fluidVelocity, db.fluidHfield],
                         currentTime = t,
                         currentCycle = steps,
                         dumpGhosts = False)
        raise ValueError("rho fluctuation outside bounds")
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
