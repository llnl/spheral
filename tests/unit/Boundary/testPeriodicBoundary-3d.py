#ATS:t0 = test(SELF, "", np=10, label="Periodic boundary unit test -- 3-D (parallel)")
#-------------------------------------------------------------------------------
# 3D test of periodic boundaries -- we simply allow a pressureless fluid to
# cycle around a box and check the sum density 
#-------------------------------------------------------------------------------
from math import *
from Spheral3d import *
from SpheralTestUtilities import *
from SpheralPointmeshSiloDump import dumpPhysicsState
import mpi

title("3D periodic boundary test.")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(nx = 20,
            ny = 20,
            nz = 20,
            x0 = 0.0,
            x1 = 1.0,
            y0 = 0.0,
            y1 = 1.0,
            z0 = 0.0,
            z1 = 1.0,

            rho1 = 1.0,
            cs2 = 1.0,
            mu = 1.0,
            vx1 = 1.0,
            vy1 = 1.0,
            vz1 = 1.0,

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
            restartBaseName = None
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
from GenerateNodeDistribution3d import GenerateNodeDistribution3d
gen1 = GenerateNodeDistribution3d(nx, ny, nz,
                                  rho = rho1,
                                  distributionType = "lattice",
                                  xmin = (x0, y0, z0),
                                  xmax = (x1, y1, z1),
                                  nNodePerh = nPerh,
                                  SPH = True)
if mpi.procs > 1:
    from PeanoHilbertDistributeNodes import distributeNodes3d
else:
    from DistributeNodes import distributeNodes3d
distributeNodes3d((nodes1, gen1))

# Set the node positions, velocities, and densities.
nodes1.velocity(VectorField("tmp velocity", nodes1, Vector(vx1, vy1, vz1)))

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
loVect = Vector(x0, y0, z0)
hiVect = Vector(x1, y1, z1)
bcs = []
for i in range(0,3):
    nVect = Vector(0., 0., 0.)
    nVect[i] = 1.
    plane0 = Plane(loVect, nVect)
    plane1 = Plane(hiVect, -nVect)
    bcs.append(PeriodicBoundary(plane0, plane1))
# Segfault occurs if hydro is append directly in previous loop
for i in bcs:
    hydro.appendBoundary(i)

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
# Make the problem controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle)

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
control.step(steps)
