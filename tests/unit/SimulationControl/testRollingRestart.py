#ATS:t0 = test(SELF, "--dataDir dumps-rolling-restart-serial", label="SpheralController rolling restarts (serial)")
#ATS:t1 = testif(t0, SELF, "--dataDir dumps-rolling-restart-serial --restoreCycle 12", label="SpheralController rolling restarts (serial) RESTART")
#ATS:t2 = test(SELF, "--dataDir dumps-rolling-restart-parallel", np=2, label="SpheralController rolling restarts (parallel)")
#ATS:t3 = testif(t2, SELF, "--dataDir dumps-rolling-restart-parallel --restoreCycle 12", np=2, label="SpheralController rolling restarts (parallel) RESTART")
#-------------------------------------------------------------------------------
# Check the SpheralController rolling restart bookkeeping: rolling restart
# files are written on their own cadence, only the most recent
# rollingRestartMax of them are kept, permanent (restartStep) restart files are
# never removed, and a restarted run picks up the existing rolling restarts.
#-------------------------------------------------------------------------------
import os, shutil, glob, mpi
from Spheral1d import *
from SpheralTestUtilities import *
from GenerateNodeDistribution1d import GenerateNodeDistribution1d
from SortAndDivideRedistributeNodes import distributeNodes1d

commandLine(nx = 50,
            restartStep = 6,
            rollingRestartStep = 2,
            rollingRestartMax = 2,
            steps = 13,
            restoreCycle = None,
            dataDir = "dumps-rolling-restart")

restartDir = os.path.join(dataDir, "restarts")
restartBaseName = os.path.join(restartDir, "rolling")
if mpi.rank == 0 and restoreCycle is None:
    if os.path.exists(dataDir):
        shutil.rmtree(dataDir)
    os.makedirs(restartDir)
mpi.barrier()

#-------------------------------------------------------------------------------
# A trivial 1D SPH problem.
#-------------------------------------------------------------------------------
eos = GammaLawGasMKS(5.0/3.0, 1.0)
WT = TableKernel(BSplineKernel(), 1000)
nodes = makeFluidNodeList("nodes", eos)
gen = GenerateNodeDistribution1d(n = nx,
                                 rho = 1.0,
                                 xmin = 0.0,
                                 xmax = 1.0)
distributeNodes1d((nodes, gen))
nodes.specificThermalEnergy(ScalarField("tmp", nodes, 1.0))

db = DataBase()
db.appendNodeList(nodes)

hydro = SPH(dataBase = db,
            W = WT)
# Keep references to the boundaries: the physics package only holds pointers.
bcs = [ReflectingBoundary(Plane(Vector(0.0), Vector( 1.0))),
       ReflectingBoundary(Plane(Vector(1.0), Vector(-1.0)))]
for bc in bcs:
    hydro.appendBoundary(bc)

integrator = CheapSynchronousRK2Integrator(db)
integrator.appendPhysicsPackage(hydro)

control = SpheralController(integrator,
                            restartStep = restartStep,
                            rollingRestartStep = rollingRestartStep,
                            rollingRestartMax = rollingRestartMax,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle)
firstCycle = control.totalSteps
assert firstCycle == (restoreCycle or 0)

#-------------------------------------------------------------------------------
# Which cycles currently have a restart file for this rank?
#-------------------------------------------------------------------------------
def restartCycles(control):
    result = set()
    prefix = os.path.basename(control.restartBaseName) + "_cycle"
    for path in glob.glob(control.restartBaseName + "_cycle*"):
        cycle = os.path.basename(path)[len(prefix):].split(".")[0]
        if cycle.isdigit():
            result.add(int(cycle))
    return result

def waitForDeletes(control):
    if control._rollingExecutor is not None:
        control._rollingExecutor.shutdown(wait = True)
        control._rollingExecutor = None
    mpi.barrier()

#-------------------------------------------------------------------------------
# Run, and check the restart files left behind.  On a restarted run the rolling
# restarts from the first run should be picked up and retired in turn.
#-------------------------------------------------------------------------------
control.step(steps)
waitForDeletes(control)
lastCycle = firstCycle + steps

# Permanent restarts: every restartStep, plus the one the controller forces at
# the end of each advance (so the end of the first run too, if this is a restart).
permanent = set(range(0, lastCycle + 1, restartStep)) | {lastCycle}
if restoreCycle is not None:
    permanent.add(steps)
rolling = sorted(set(range(0, lastCycle + 1, rollingRestartStep)) - permanent)
expectedRolling = set(rolling[-rollingRestartMax:])
found = restartCycles(control)
print("Restart cycles found: ", sorted(found))
assert permanent - {0} <= found, "Missing permanent restart files: %s" % sorted(permanent - found)
assert found - permanent == expectedRolling, \
    "Expected rolling restart cycles %s, found %s" % (sorted(expectedRolling), sorted(found - permanent))
assert control._rollingCycles == sorted(expectedRolling)
print("PASS")
