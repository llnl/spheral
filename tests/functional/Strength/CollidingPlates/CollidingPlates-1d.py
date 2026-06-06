# SPH non-compatible
#ATS:t10 = test(SELF,        "--graphics False --hydroType SPH --restartStep 50 --clearDirectories True  --domainIndependent True --compatibleEnergy False --outputFile 'CollidingPlates-1d-1proc-reproducing.txt' --checkReference True", np=1, label="Colliding plates domain independence test SERIAL SPH Non-compatible RUN")
#ATS:t11 = testif(t10, SELF, "--graphics False --hydroType SPH                  --clearDirectories False --domainIndependent True --compatibleEnergy False --outputFile 'CollidingPlates-1d-1proc-reproducing.txt' --checkReference True --restoreCycle 50", np=1, label="Colliding plates domain independence test SERIAL SPH RESTART Non-compatible RUN")
#ATS:t12 = testif(t10, SELF, "--graphics False --hydroType SPH                  --clearDirectories False --domainIndependent True --compatibleEnergy False --outputFile 'CollidingPlates-1d-4proc-reproducing.txt' --checkReference True --comparisonFile 'dumps-CollidingPlates-1d/SPH/compatibleEnergy=False/100/CollidingPlates-1d-1proc-reproducing.txt'", np=4, label="Colliding plates domain independence test 4 DOMAIN SPH Non-compatible RUN")
#
# SPH compatible
#
#ATS:t20 = test(SELF,        "--graphics False --hydroType SPH --restartStep 50 --clearDirectories True  --domainIndependent True --compatibleEnergy True --outputFile 'CollidingPlates-1d-1proc-reproducing.txt' --checkReference True", np=1, label="Colliding plates domain independence test SERIAL SPH Compatible RUN")
#ATS:t21 = testif(t20, SELF, "--graphics False --hydroType SPH                  --clearDirectories False --domainIndependent True --compatibleEnergy True --outputFile 'CollidingPlates-1d-1proc-reproducing.txt' --checkReference True --restoreCycle 50", np=1, label="Colliding plates domain independence test SERIAL SPH RESTART Compatible RUN")
#ATS:t22 = testif(t20, SELF, "--graphics False --hydroType SPH                  --clearDirectories False --domainIndependent True --compatibleEnergy True --outputFile 'CollidingPlates-1d-4proc-reproducing.txt' --checkReference True --comparisonFile 'dumps-CollidingPlates-1d/SPH/compatibleEnergy=True/100/CollidingPlates-1d-1proc-reproducing.txt'", np=4, label="Colliding plates domain independence test 4 DOMAIN SPH Compatible RUN")
#
# CRKSPH compatible
#
#ATS:t30 = test(SELF,        "--graphics False --hydroType CRKSPH --restartStep 50 --clearDirectories True  --domainIndependent True --compatibleEnergy True --outputFile 'CollidingPlates-1d-1proc-reproducing.txt' --checkReference True", np=1, label="Colliding plates domain independence test SERIAL CRKSPH Compatible RUN")
#ATS:t31 = testif(t30, SELF, "--graphics False --hydroType CRKSPH                  --clearDirectories False --domainIndependent True --compatibleEnergy True --outputFile 'CollidingPlates-1d-1proc-reproducing.txt' --checkReference True --restoreCycle 50", np=1, label="Colliding plates domain independence test SERIAL CRKSPH RESTART Compatible RUN")
#ATS:t32 = testif(t30, SELF, "--graphics False --hydroType CRKSPH                  --clearDirectories False --domainIndependent True --compatibleEnergy True --outputFile 'CollidingPlates-1d-4proc-reproducing.txt' --checkReference True --comparisonFile 'dumps-CollidingPlates-1d/CRKSPH/compatibleEnergy=True/100/CollidingPlates-1d-1proc-reproducing.txt'", np=4, label="Colliding plates domain independence test 4 DOMAIN CRKSPH Compatible RUN")
#
# FSISPH compatible    -- note: domain independence not working with FSI currently
#
#ATS:t40 = test(SELF,        "--graphics False --hydroType FSISPH --restartStep 50 --clearDirectories True  --domainIndependent True --compatibleEnergy True --outputFile 'CollidingPlates-1d-1proc-reproducing.txt' --checkReference True", np=1, label="Colliding plates domain independence test SERIAL FSISPH Compatible RUN")
#ATS:t41 = testif(t40, SELF, "--graphics False --hydroType FSISPH                  --clearDirectories False --domainIndependent True --compatibleEnergy True --outputFile 'CollidingPlates-1d-1proc-reproducing.txt' --checkReference True --restoreCycle 50", np=1, label="Colliding plates domain independence test SERIAL FSISPH RESTART Compatible RUN")
###ATS:t42 = testif(t40, SELF, "--graphics False --hydroType FSISPH                  --clearDirectories False --domainIndependent True --compatibleEnergy True --outputFile 'CollidingPlates-1d-4proc-reproducing.txt' --checkReference True --comparisonFile 'dumps-CollidingPlates-1d/CRKSPH/compatibleEnergy=True/100/CollidingPlates-1d-1proc-reproducing.txt'", np=4, label="Colliding plates domain independence test 4 DOMAIN FSISPH Compatible RUN")

#-------------------------------------------------------------------------------
# A pair of steel plates colliding at the origin.  This is a useful test of
# how our reflecting boundary condition handles problems with strength.
#-------------------------------------------------------------------------------
from SolidSpheral1d import *
from SpheralTestUtilities import *
from findLastRestart import *
from SpheralVisitDump import dumpPhysicsState
from identifyFragments import identifyFragments, fragmentProperties
from math import *
import shutil
import mpi

#-------------------------------------------------------------------------------
# Identify ourselves!
#-------------------------------------------------------------------------------
title("1-D colliding plates strength test")

#-------------------------------------------------------------------------------
# Pick our units
#-------------------------------------------------------------------------------
units = CGuS()

#-------------------------------------------------------------------------------
# Generic problem parameters
# All CGuS units.
#-------------------------------------------------------------------------------
commandLine(# Geometry
            length = 3.0,
            radius = 0.5,
            nx = 100,
            reflect = False,
            v0 = 1.0e-2,        # Initial collision velocity

            # Material specific bounds on the mass density.
            rho0 = 7.9,
            etamin = 0.5,
            etamax = 1.5,

            # Hydro
            hydroType = "SPH",  # (SPH, CRKSPH, FSISPH)
            nPerh = 1.25,
            hmin = 1e-5,
            hmax = 1.0,
            HUpdate = IdealH,
            densityUpdate = IntegrateDensity,
            compatibleEnergy = True,
            evolveTotalEnergy = False,
            gradhCorrection = False,
            correctVelocityGradient = True,

            # Time advancement.
            IntegratorConstructor = CheapSynchronousRK2Integrator,
            goalTime = 2.0,
            steps = None,
            dt = 1e-6,
            dtMin = 1e-6,
            dtMax = 0.1,
            dtGrowth = 2.0,
            dumpFrac = 1.0,
            maxSteps = None,
            statsStep = 10,
            smoothIters = 0,
            domainIndependent = False,

            restoreCycle = None,
            restartStep = 500,

            # Diagnostics
            graphics = True,
            testtol = 1.0e-3,
            clearDirectories = False,
            checkReference = True,
            dataDirBase = "dumps-CollidingPlates-1d",
            outputFile = None,
            comparisonFile = None,
            )

hydroType = hydroType.upper()
assert hydroType in ("SPH", "CRKSPH", "FSISPH")

dataDir = os.path.join(dataDirBase,
                       hydroType,
                       f"compatibleEnergy={compatibleEnergy}",
                       str(nx))
restartDir = os.path.join(dataDir, "restarts")
restartBaseName = os.path.join(restartDir, "CollidingPlates-%i" % nx)

if reflect:
    xmin = 0.0
else:
    xmin = -0.5*length
xmax =  0.5*length

dtSample = dumpFrac*goalTime

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
# Stainless steel material properties.
#-------------------------------------------------------------------------------
eos = GruneisenEquationOfState(rho0,    # reference density  
                               etamin,  # etamin             
                               etamax,  # etamax             
                               0.457,   # C0                 
                               1.49,    # S1                 
                               0.0,     # S2                 
                               0.0,     # S3                 
                               1.93,    # gamma0             
                               0.5,     # b                  
                               55.350,  # atomic weight
                               constants = units)

coldFit = NinthOrderPolynomialFit(-1.06797724e-2,
                                  -2.06872020e-2,
                                   8.24893246e-1,
                                  -2.39505843e-2,
                                  -2.44522017e-2,
                                   5.38030101e-2,
                                   0.0,
                                   0.0,
                                   0.0,
                                   0.0)
meltFit = NinthOrderPolynomialFit(7.40464217e-2,
                                  2.49802214e-1,
                                  1.00445029e00,
                                 -1.36451475e-1,
                                  7.72897829e-3,
                                  5.06390305e-2,
                                  0.0,
                                  0.0,
                                  0.0,
                                  0.0)
strengthModel = SteinbergGuinanStrength(eos,
                                        7.700000e-1,        # G0
                                        2.2600,             # A
                                        4.5500e-04,         # B
                                        3.4000e-3,          # Y0
                                        2.5e-2,             # Ymax
                                        1.0e-3,             # Yp
                                        43.0000,            # beta
                                        0.0,                # gamma0
                                        0.35,               # nhard
                                        coldFit,
                                        meltFit)

#-------------------------------------------------------------------------------
# Create our interpolation kernels -- one for normal hydro interactions, and
# one for use with the artificial viscosity
#-------------------------------------------------------------------------------
WT = TableKernel(BSplineKernel(), 1000)
output("WT")

#-------------------------------------------------------------------------------
# Create the NodeLists.
#-------------------------------------------------------------------------------
nodes = makeSolidNodeList("Stainless steel", eos, strengthModel,
                          nPerh = nPerh,
                          hmin = hmin,
                          hmax = hmax,
                          rhoMin = etamin*rho0,
                          rhoMax = etamax*rho0)
nodeSet = [nodes]

#-------------------------------------------------------------------------------
# Set node properties (positions, masses, H's, etc.)
#-------------------------------------------------------------------------------
eps0 = 0.0
print("Generating node distribution.")
from DistributeNodes import distributeNodesInRange1d
distributeNodesInRange1d([(nodes, nx, rho0, (xmin, xmax))])
output("mpi.reduce(nodes.numInternalNodes, mpi.MIN)")
output("mpi.reduce(nodes.numInternalNodes, mpi.MAX)")
output("mpi.reduce(nodes.numInternalNodes, mpi.SUM)")

# Set node specific thermal energies
T = ScalarField("temperature", nodes, 300.0)
rho = nodes.massDensity()
eps = nodes.specificThermalEnergy()
eos.setSpecificThermalEnergy(eps, rho, T)

# Set node velocites.
pos = nodes.positions()
vel = nodes.velocity()
for i in range(nodes.numInternalNodes):
    if pos[i].x < 0.0:
        vel[i].x = v0
    else:
        vel[i].x = -v0

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
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
if hydroType == "CRKSPH":
    hydro = CRKSPH(dataBase = db,
                   W = WT,
                   compatibleEnergyEvolution = compatibleEnergy,
                   densityUpdate = densityUpdate,
                   HUpdate = HUpdate)

elif hydroType == "FSISPH":
    hydro = FSISPH(dataBase = db,
                   W = WT,
                   interfaceMethod = HLLCInterface,
                   densityStabilizationCoefficient = 0.00,
                   compatibleEnergyEvolution = compatibleEnergy,
                   linearCorrectGradients = correctVelocityGradient,
                   HUpdate = HUpdate)

else:
    assert hydroType == "SPH"
    hydro = SPH(dataBase = db,
                W = WT,
                compatibleEnergyEvolution = compatibleEnergy,
                evolveTotalEnergy = evolveTotalEnergy,
                gradhCorrection = gradhCorrection,
                correctVelocityGradient = correctVelocityGradient,
                densityUpdate = densityUpdate,
                HUpdate = HUpdate)

output("hydro")
output("hydro.cfl")
output("hydro.useVelocityMagnitudeForDt")
output("hydro._smoothingScaleMethod.HEvolution")
output("hydro.densityUpdate")
output("hydro.compatibleEnergyEvolution")
output("hydro.Q")
output("hydro.Q.Cl")
output("hydro.Q.Cq")
output("hydro.Q.epsilon2")
output("hydro.Q.limiter")
output("hydro.Q.balsaraShearCorrection")

#-------------------------------------------------------------------------------
# Boundary conditions.
#-------------------------------------------------------------------------------
if reflect:
    plane0 = Plane(Vector(0.0), Vector(1.0))
    xbc0 = ReflectingBoundary(plane0)
    hydro.appendBoundary(xbc0)

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
output("integrator")
output("integrator.havePhysicsPackage(hydro)")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.dtGrowth")
output("integrator.domainDecompositionIndependent")

#-------------------------------------------------------------------------------
# Build the controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle)
output("control")

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
if steps is not None:
    control.step(steps)

    pos = nodes.positions()
    rho = nodes.massDensity()
    vel = nodes.velocity()
    eps = nodes.specificThermalEnergy()
    P = ScalarField("pressure", nodes)
    nodes.pressure(P)
    S = nodes.deviatoricStress()
    if reflect:
        points = [50, 0]
    else:
        points = [49, 50]
    for p in points:
        print(pos[p], rho[p], vel[p], eps[p], P[p], S[p])

else:
    while control.time() < goalTime:
        nextGoalTime = min(control.time() + dtSample, goalTime)
        control.advance(nextGoalTime, maxSteps)
        control.dropRestartFile()

Eerror = (control.conserve.EHistory[-1] - control.conserve.EHistory[0])/control.conserve.EHistory[0]
print("Total energy error: %g" % Eerror)
if compatibleEnergy and abs(Eerror) > 1e-13:
    raise ValueError("Energy error outside allowed bounds.")

#-------------------------------------------------------------------------------
# Plot the state.
#-------------------------------------------------------------------------------
if graphics:
    from SpheralMatplotlib import *
    state = State(db, integrator.physicsPackages())
    rhoPlot = plotFieldList(state.scalarFields("mass density"),
                            plotStyle = "r-",
                            winTitle = "rho @ %g %i" % (control.time(), mpi.procs))
    velPlot = plotFieldList(state.vectorFields("velocity"),
                            yFunction = "%s.x",
                            plotStyle = "r-",
                            winTitle = "vel @ %g %i" % (control.time(), mpi.procs))
    mPlot = plotFieldList(state.scalarFields("mass"),
                          plotStyle = "r-",
                          winTitle = "mass @ %g %i" % (control.time(), mpi.procs))
    PPlot = plotFieldList(state.scalarFields("pressure"),
                          plotStyle = "r-",
                          winTitle = "pressure @ %g %i" % (control.time(), mpi.procs))
    SPlot = plotFieldList(state.symTensorFields(SolidFieldNames.deviatoricStress),
                          yFunction = "%s.xx",
                          plotStyle = "r-",
                          winTitle = "Deviatoric stress @ %g %i" % (control.time(), mpi.procs))

#-------------------------------------------------------------------------------
# If requested, write out the state in a global ordering to a file.
#-------------------------------------------------------------------------------
if outputFile:
    state = State(db, integrator.physicsPackages())
    outputFile = os.path.join(dataDir, outputFile)
    pos = state.vectorFields(HydroFieldNames.position)
    rho = state.scalarFields(HydroFieldNames.massDensity)
    P = state.scalarFields(HydroFieldNames.pressure)
    vel = state.vectorFields(HydroFieldNames.velocity)
    eps = state.scalarFields(HydroFieldNames.specificThermalEnergy)
    Hfield = state.symTensorFields(HydroFieldNames.H)
    S = state.symTensorFields(SolidFieldNames.deviatoricStress)
    xprof = mpi.reduce([x.x for x in internalValues(pos)], mpi.SUM)
    rhoprof = mpi.reduce(internalValues(rho), mpi.SUM)
    Pprof = mpi.reduce(internalValues(P), mpi.SUM)
    vprof = mpi.reduce([v.x for v in internalValues(vel)], mpi.SUM)
    epsprof = mpi.reduce(internalValues(eps), mpi.SUM)
    hprof = mpi.reduce([1.0/sqrt(H.Determinant()) for H in internalValues(Hfield)], mpi.SUM)
    sprof = mpi.reduce([x.xx for x in internalValues(S)], mpi.SUM)
    mof = mortonOrderIndices(db)
    mo = mpi.reduce(internalValues(mof), mpi.SUM)
    if mpi.rank == 0:
        thpt = list(zip(mo, xprof, rhoprof, Pprof, vprof, epsprof, hprof, sprof))
        thpt = [thpt[i] for i in range(len(thpt)) if thpt[i][1] > 0.0]
        thpt.sort()
        mo, xprof, rhoprof, Pprof, vprof, epsprof, hprof, sprof = list(zip(*thpt))
        f = open(outputFile, "w")
        f.write(("#" + 7*" %16s" + "\n") % ("x", "rho", "P", "v", "eps", "h", "S"))
        for (xi, rhoi, Pi, vi, epsi, hi, si, mi) in zip(xprof, rhoprof, Pprof, vprof, epsprof, hprof, sprof, mo):
            f.write((7*"%16.12e " + "\n") %
                    (xi, rhoi, Pi, vi, epsi, hi, si))
        f.close()

        #---------------------------------------------------------------------------
        # Check the floating values for the state against reference data.
        #---------------------------------------------------------------------------
        if checkReference:
            referenceFile = os.path.join("Reference",
                                         hydroType,
                                         f"compatibleEnergy={compatibleEnergy}",
                                         str(nx),
                                         "CollidingPlates-1d-1proc-reproducing.txt")
            import filearraycmp as fcomp
            assert fcomp.filearraycmp(outputFile, referenceFile, testtol, testtol)
            print("Floating point comparison test passed.")

        #---------------------------------------------------------------------------
        # Also we can optionally compare the current results with another file for
        # bit level consistency.
        #---------------------------------------------------------------------------
        if comparisonFile and BuildData.cxx_compiler_id != "IntelLLVM":
            import filecmp
            print("Compare files : %s     <--->     %s" % (outputFile, comparisonFile))
            assert filecmp.cmp(outputFile, comparisonFile)
