#ATS:t0 = test(SELF,       "--nr 10 --numViz 0 --graphics None --timeStepChoice AccelerationRatio --steps=40 --restartStep 20  --dataDir 'Collisionless_Sphere_Collapse_AccelerationRatio' --clearDirectories True --outputFile 'Collisionless_sphere_collapse_AccelerationRatio_data.gnu' --checkRef True", np=1, label="Collisionless sphere gravitational collapse restart test (serial, acceleration ratio) INITIAL RUN")
#
#ATS:t1 = testif(t0, SELF, "--nr 10 --numViz 0 --graphics None --timeStepChoice AccelerationRatio --steps 20 --restartStep 100 --dataDir 'Collisionless_Sphere_Collapse_AccelerationRatio' --clearDirectories False --outputFile 'Collisionless_sphere_collapse_AccelerationRatio_data.gnu' --checkRef True --restoreCycle 20 --checkRestart True", np=1, label="Collisionless sphere gravitational collapse restart test (serial, acceleration ratio) RESTARTED CHECK")
#
#ATS:t2 = test(SELF,       "--nr 10 --numViz 0 --graphics None --timeStepChoice DynamicalTime --steps=40 --restartStep 20   --dataDir 'Collisionless_Sphere_Collapse_DynamicalTime' --clearDirectories True --outputFile 'Collisionless_sphere_collapse_DynamicalTime_data.gnu' --checkRef True", np=1, label="Collisionless sphere gravitational collapse restart test (serial, dynamical time) INITIAL RUN")
#
#ATS:t3 = testif(t2, SELF, "--nr 10 --numViz 0 --graphics None --timeStepChoice DynamicalTime --steps 20 --restartStep 100  --dataDir 'Collisionless_Sphere_Collapse_DynamicalTime' --clearDirectories False --outputFile 'Collisionless_sphere_collapse_DynamicalTime_data.gnu' --checkRef True --restoreCycle 20 --checkRestart True", np=1, label="Collisionless sphere gravitational collapse restart test (serial, dynamical time) RESTARTED CHECK")

#-------------------------------------------------------------------------------
# Create a spherical distribution of collisionless points, which will of course 
# promptly collapse under their own self-gravity.
#-------------------------------------------------------------------------------
import os, shutil
from Spheral3d import *
from SpheralTestUtilities import *
from NodeHistory import *
from GenerateNodeDistribution3d import *
from PeanoHilbertDistributeNodes import distributeNodes3d
from math import *
import numpy as np

print("3-D N-Body Gravity test -- collisionless sphere.")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(

    # Initial particle stuff
    nr = 50,                       # radial number of particles to seed
    r0 = 1.0,                      # (AU) initial radius of the sphere
    M0 = 1.0,                      # (earth masses) total mass of the sphere
    plummerLength = 1.0e-2,        # (AU) Plummer softening scale
    opening = 1.0,                 # (dimensionless, OctTreeGravity) opening parameter for tree walk
    fdt = 0.1,                     # (dimensionless, OctTreeGravity) timestep multiplier

    # Problem control
    steps = None,                  # Optionally advance a fixed number of steps
    numCollapseTimes = 1.0,        # What time to advance to in units of the collapse time for the sphere

    # Which N-body method should we use?
    nbody = OctTreeGravity,
    timeStepChoice = AccelerationRatio,

    # Output
    clearDirectories = False,
    dataDir = "dumps_Collisionless_Sphere_Collapse",
    baseNameRoot = "sphere_collapse_%i",
    restoreCycle = -1,
    restartStep = 100,
    numViz = 200,
    verbosedt = False,

    graphics = True,

    # Parameters purely for test checking
    outputFile = "Collisionless_sphere_collapse.gnu",
    checkRestart = False,
    checkRef = False,
    tol = 5.0e-5,
    )

# Reference values for tests
if timeStepChoice == AccelerationRatio:
    coefsRef = np.array([ 8.33175998e+00,  1.24358171e-12, -2.83895427e-23])
    sigmaPhiRef = 9.25853628363647
elif timeStepChoice == DynamicalTime:
    coefsRef = np.array([ 8.29009280e+00,  1.20192051e-12, -2.60731358e-23])
    sigmaPhiRef = 8.968145544554178

# Convert to MKS units.
AU = 149597870700.0  # m
Mearth = 5.9722e24   # kg
r0 *= AU
M0 *= Mearth
plummerLength *= AU

# Compute the analytically expected collapse time.
G = MKS().G
rho0 = M0/(4.0/3.0*pi*r0*r0*r0)
tdyn = sqrt(3.0*pi/(16*G*rho0))
collapseTime = tdyn/sqrt(2.0)

# Miscellaneous problem control parameters.
dt = collapseTime / 100
goalTime = collapseTime * numCollapseTimes
dtMin, dtMax = 0.1*dt, 100.0*dt
dtGrowth = 2.0
maxSteps = None
statsStep = 10
smoothIters = 0
if numViz > 0:
    vizTime = goalTime / numViz
else:
    vizTime = goalTime

baseName = baseNameRoot % nr
restartDir = os.path.join(dataDir, "restarts")
visitDir = os.path.join(dataDir, "visit")
restartBaseName = os.path.join(restartDir, baseName + "_restart")

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
import os, sys
if mpi.rank == 0:
    if clearDirectories and os.path.exists(dataDir):
        shutil.rmtree(dataDir)
    if not os.path.exists(restartDir):
        os.makedirs(restartDir)
    if not os.path.exists(visitDir):
        os.makedirs(visitDir)
mpi.barrier()

#-------------------------------------------------------------------------------
# For now we have set up a fluid node list, even though this is collisionless
# problem.  Fix at some point!
# In the meantime, set up the hydro objects this script isn't really going to
# need.
#-------------------------------------------------------------------------------
WT = TableKernel(BSplineKernel(), 1000)
eos = GammaLawGasMKS3d(gamma = 5.0/3.0, mu = 1.0)

#-------------------------------------------------------------------------------
# Make the NodeList, and set our initial properties.
#-------------------------------------------------------------------------------
nodes = makeFluidNodeList("nodes", eos,
                          xmin = -100.0*r0 * Vector.one,
                          xmax =  100.0*r0 * Vector.one)

generator = GenerateNodeDistribution3d(2*nr, 2*nr, 2*nr, rho0,
                                       distributionType = "lattice",
                                       xmin = (-r0, -r0, -r0),
                                       xmax = ( r0,  r0,  r0),
                                       rmin = 0.0,
                                       rmax = r0)
distributeNodes3d((nodes, generator))
output("mpi.reduce(nodes.numInternalNodes, mpi.MIN)")
output("mpi.reduce(nodes.numInternalNodes, mpi.MAX)")
output("mpi.reduce(nodes.numInternalNodes, mpi.SUM)")

# Renormalize the node masses to get our desired total mass.
mass = nodes.mass()
msum = mpi.allreduce(sum(nodes.mass().internalValues()), mpi.SUM)
fmass = M0/msum
print("Renormalizing masses by ", fmass)
for i in range(nodes.numInternalNodes):
    mass[i] *= fmass

#-------------------------------------------------------------------------------
# DataBase
#-------------------------------------------------------------------------------
db = DataBase()
db.appendNodeList(nodes)

#-------------------------------------------------------------------------------
# Gimme gravity.
#-------------------------------------------------------------------------------
if nbody is NBodyGravity:
    gravity = NBodyGravity(plummerSofteningLength = plummerLength,
                           maxDeltaVelocity = 1e-2*r0/tdyn,
                           G = G)
# elif nbody is FractalGravity:
#     gravity = FractalGravity(G = G,
#                              xmin = Vector(-1.5*r0, -1.5*r0, -1.5*r0),
#                              xmax = Vector( 1.5*r0,  1.5*r0,  1.5*r0),
#                              periodic = False,
#                              ngrid = 64,
#                              nlevelmax = 1,
#                              minHighParticles = 10,
#                              padding = 0,
#                              maxDeltaVelocity = 1e-2*v0)
elif nbody is OctTreeGravity:
    gravity = OctTreeGravity(G = G,
                             softeningLength = plummerLength,
                             opening = opening,
                             ftimestep = fdt,
                             timeStepChoice = timeStepChoice)

#-------------------------------------------------------------------------------
# Construct a time integrator.
#-------------------------------------------------------------------------------
integrator = SynchronousRK2Integrator(db)
integrator.appendPhysicsPackage(gravity)
integrator.lastDt = dtMin
integrator.dtMin = dtMin
integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
integrator.verbose = verbosedt

#-------------------------------------------------------------------------------
# Build the problem controller to follow the problem evolution.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle,
                            vizBaseName = baseName,
                            vizDir = visitDir,
                            vizTime = vizTime)

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
if not steps is None:
    if checkRestart:
        control.setRestartBaseName(restartBaseName + "_CHECK")
    control.step(steps)
    if checkRestart:
        control.setRestartBaseName(restartBaseName)

    # Are we doing the restart test?
    if checkRestart:
        state0 = State(db, integrator.physicsPackages())
        state0.copyState()
        control.loadRestartFile(control.totalSteps)
        state1 = State(db, integrator.physicsPackages())
        if not state1 == state0:
            raise ValueError("The restarted state does not match!")
        else:
            print("Restart check PASSED.")

else:
    print("Advancing to %g sec = %g years" % (goalTime, goalTime/(365.24*24*3600)))
    control.advance(goalTime)

#-------------------------------------------------------------------------------
# Compute the radial profiles
#-------------------------------------------------------------------------------
from SpheralTestUtilities import multiSort
import numpy.polynomial.polynomial as poly
xprof = mpi.reduce(nodes.positions().internalValues(), mpi.SUM)
rprof = [x.magnitude() for x in xprof]
vprof = mpi.reduce(nodes.velocity().internalValues(), mpi.SUM)
Hprof = mpi.reduce(nodes.Hfield().internalValues(), mpi.SUM)
phi = gravity.potential
phiprof = mpi.reduce(phi[0].internalValues(), mpi.SUM)
mof = mortonOrderIndices(db)
mo = mpi.reduce(mof[0].internalValues(), mpi.SUM)
coefs, phifit = None, None
if mpi.rank == 0:
    from SpheralTestUtilities import multiSort
    multiSort(rprof, mo, xprof, vprof, Hprof, phiprof)

    # Fit phi(r)
    coefs = poly.polyfit(rprof, np.log(-np.array(phiprof)), 2)
    phifit = -np.exp(poly.polyval(rprof, coefs))
    print("Fit coefficients: ", coefs)
    sigphi = np.std(np.array(phiprof) - phifit)
    print("Standard deviation: ", sigphi)

coefs = mpi.bcast(coefs, root=0)
phifit = mpi.bcast(phifit, root=0)

if graphics:
    from SpheralMatplotlib import *
    EPlot = plotEHistory(control.conserve)
    phiPlot = plotFieldList(gravity.potential, xFunction="%s.magnitude()", plotStyle="bo", winTitle="Gravitational potential $\phi$")
    phiPlot.plot(rprof, phifit, "k-", label="Fit")
    phiPlot.axes.legend()

#-------------------------------------------------------------------------------
# If requested, write out the profiles
#-------------------------------------------------------------------------------
if outputFile and mpi.rank == 0:
    outputFile = os.path.join(dataDir, outputFile)
    f = open(outputFile, "w")
    f.write(("# " + 14*"%15s " + "\n") % ("r", "x", "y", "z", "vx", "vy", "vz", "Hxx", "Hxy", "Hxz", "Hyy", "Hyz", "Hzz", "phi"))
    for (ri, xi, vi, Hi, phii, moi) in zip(rprof, xprof, vprof, Hprof, phiprof, mo):
        f.write((14*" %16.12e" + "\n") % (ri, xi.x, xi.y, xi.z, vi.x, vi.y, vi.z, 
                                          Hi.xx, Hi.xy, Hi.xz, Hi.yy, Hi.yz, Hi.zz, phii))
    f.close()

#-------------------------------------------------------------------------------
# Check the answer againt references
#-------------------------------------------------------------------------------
if checkRef:
    print((coefs, coefsRef))
    print((sigphi, sigmaPhiRef))
    assert np.absolute(coefs - coefsRef).max() < tol
    assert abs(sigphi - sigmaPhiRef)/sigmaPhiRef < tol
    print("Pass")
