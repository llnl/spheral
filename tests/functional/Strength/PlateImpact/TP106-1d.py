#-------------------------------------------------------------------------------
# An idealized strength test test mocking a flyer plate of Aluminum impacting 
# a fixed plate of Tantalum.
#
# See Weseloh & Najjar, LLNL-TR-484921
#-------------------------------------------------------------------------------
from Spheral1d import *
from SpheralTestUtilities import *
from findLastRestart import *
from NodeHistory import NodeHistory
from math import *
import shutil
import mpi

#-------------------------------------------------------------------------------
# Identify ourselves!
#-------------------------------------------------------------------------------
title("1-D TP106 flyer plate test.")

#-------------------------------------------------------------------------------
# Use (cm, gm, usec) units as in paper.
#-------------------------------------------------------------------------------
units = PhysicalConstants(0.01,  # Unit length in m
                          0.001, # Unit mass in kg
                          1e-6)  # Unit length in sec

#-------------------------------------------------------------------------------
# Generic problem parameters
# All (cm, gm, usec) units.
#-------------------------------------------------------------------------------
commandLine(nxAl = 100,
            nxTa = 0,       # 0 implies automatically mass match against Al

            # Material specific bounds on the mass density.
            etamin = 1e-3,
            etamax = 1e3,

            # Should we run with strength?
            useStrength = True,

            # Hydro parameters.
            crksph = False,
            hmin = 1e-5,
            hmax = 0.1,
            cfl = 0.25,
            useVelocityMagnitudeForDt = False,
            XSPH = False,
            epsilonTensile = 0.0,
            nTensile = 4,
            filter = 0.0,
            HUpdate = IdealH,
            densityUpdate = IntegrateDensity,
            compatibleEnergy = True,
            evolveTotalEnergy = False,
            correctionOrder = LinearOrder,                # CRKSPH
            volumeType = RKVoronoiVolume,                 # CRKSPH

            # Time integration
            IntegratorConstructor = CheapSynchronousRK2Integrator,
            goalTime = 3.0,
            steps = None,
            dt = 1e-10,
            dtMin = 1e-6,
            dtMax = 10.0,
            dtGrowth = 2.0,
            dumpFrac = 0.005,
            maxSteps = None,
            statsStep = 10,
            domainIndependent = False,
            dtverbose = False,
            restoreCycle = -1,
            restartStep = 500,
            sampleFreq = 10,

            graphics = True,

            clearDirectories = False,
            dataDirBase = "dumps-TP106-1d",
            outputFile = "TP106-profiles.gnu",
            )

# Figure out the mass matched resolution.
rho0Al = 2.785
rho0Ta = 16.654
if nxTa == 0:
    nxTa = int(rho0Ta/rho0Al * nxAl)
    print("Selected %i Ta points to mass match with %i Al points." % (nxTa, nxAl))

# Hydro constructor.
if crksph:
    hydroname = os.path.join("CRKSPH",
                             str(correctionOrder),
                             str(volumeType))
    kernelOrder = 7
    nPerh = 1.01
else:
    hydroname = "SPH"
    kernelOrder = 5
    nPerh = 1.51

# Directories.
dataDir = os.path.join(dataDirBase,
                       "strength=%s" % useStrength,
                       hydroname,
                       "densityUpdate=%s" % densityUpdate,
                       "nxAl=%i_nxTa=%i" % (nxAl, nxTa))
restartDir = os.path.join(dataDir, "restarts")
visitDir = os.path.join(dataDir, "visit")
restartBaseName = os.path.join(restartDir, "TP106-%i" % (nxAl + nxTa))
tracerOutputName = os.path.join(dataDir, "TP106-tracer%i.dat")

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
# Material properties.
#-------------------------------------------------------------------------------
# Al
eosAl = GruneisenEquationOfState(rho0Al,  # reference density  
                                 etamin,  # etamin             
                                 etamax,  # etamax             
                                 0.5328,  # C0                 
                                 1.338,   # S1                 
                                 0.0,     # S2                 
                                 0.0,     # S3                 
                                 2.0,     # gamma0             
                                 0.0,     # b                  
                                 26.98,  # atomic weight
                                 units)
eosAl.energyMultiplier = 0.0  # <-- Hack for this particular problem

coldFit = NinthOrderPolynomialFit(0.0,
                                  0.0,
                                  0.0,
                                  0.0,
                                  0.0,
                                  0.0,
                                  0.0,
                                  0.0,
                                  0.0,
                                  0.0)
meltFitAl = NinthOrderPolynomialFit(0.0796,
                                    0.0,
                                    0.0,
                                    0.0,
                                    0.0,
                                    0.0,
                                    0.0,
                                    0.0,
                                    0.0,
                                    0.0)
strengthModelAl = SteinbergGuinanStrength(eosAl,
                                          0.2860,        # G0
                                          0.6470,        # Gmax
                                          6.52,          # A
                                          0.0,           # B
                                          0.0026,        # Y0
                                          0.0076,        # Ymax
                                          1.0e-3,        # Yp (delta)
                                          310.0,         # beta
                                          0.0,           # gamma0
                                          0.185,         # nhard
                                          coldFit,
                                          meltFitAl)

# Ta
eosTa = GruneisenEquationOfState(rho0Ta,  # reference density  
                                 etamin,  # etamin             
                                 etamax,  # etamax             
                                 0.3414,  # C0                 
                                 1.201,   # S1                 
                                 0.0,     # S2                 
                                 0.0,     # S3                 
                                 2.0,     # gamma0             
                                 0.0,     # b                  
                                 180.95,  # atomic weight
                                 units)
eosTa.energyMultiplier = 0.0  # <-- Hack for this particular problem

meltFitTa = NinthOrderPolynomialFit(0.0546,
                                    0.0,
                                    0.0,
                                    0.0,
                                    0.0,
                                    0.0,
                                    0.0,
                                    0.0,
                                    0.0,
                                    0.0)
strengthModelTa = SteinbergGuinanStrength(eosTa,
                                          0.6900,        # G0
                                          0.9860,        # Gmax
                                          1.45,          # A
                                          0.0,           # B
                                          0.0077,        # Y0
                                          0.0110,        # Ymax
                                          1.0e-3,        # Yp (delta)
                                          10.0,          # beta
                                          0.0,           # gamma0
                                          0.1,           # nhard
                                          coldFit,
                                          meltFitTa)

#-------------------------------------------------------------------------------
# If we're not using strength, override the strength models.
#-------------------------------------------------------------------------------
if not useStrength:
    strengthModelAl = NullStrength()
    strengthModelTa = NullStrength()

#-------------------------------------------------------------------------------
# Create our interpolation kernels -- one for normal hydro interactions, and
# one for use with the artificial viscosity
#-------------------------------------------------------------------------------
WT = TableKernel(NBSplineKernel(kernelOrder), 1000)
output("WT")

#-------------------------------------------------------------------------------
# Create the NodeLists.
#-------------------------------------------------------------------------------
nodesAl = makeSolidNodeList("Aluminum", eosAl, strengthModelAl,
                            nPerh = nPerh,
                            hmin = hmin,
                            hmax = hmax,
                            rhoMin = etamin*rho0Al,
                            rhoMax = etamax*rho0Al,
                            xmin = -100.0*Vector.one,
                            xmax =  100.0*Vector.one)
nodesTa = makeSolidNodeList("Tantalum", eosTa, strengthModelTa,
                            nPerh = nPerh,
                            hmin = hmin,
                            hmax = hmax,
                            rhoMin = etamin*rho0Ta,
                            rhoMax = etamax*rho0Ta,
                            xmin = -100.0*Vector.one,
                            xmax =  100.0*Vector.one)
nodeSet = [nodesAl, nodesTa]

#-------------------------------------------------------------------------------
# Set node properties (positions, masses, H's, etc.)
#-------------------------------------------------------------------------------
print("Generating node distribution.")
from DistributeNodes import distributeNodesInRange1d
distributeNodesInRange1d([(nodesAl, nxAl, rho0Al, (-1.025, -0.025)),
                          (nodesTa, nxTa, rho0Ta, ( 0.025,  1.025))])

# Set node velocites.
vel = nodesAl.velocity()
for i in range(nodesAl.numInternalNodes):
    vel[i].x = 0.18

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

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
if crksph:
    hydro = CRKSPH(dataBase = db,
                   order = correctionOrder,
                   filter = filter,
                   cfl = cfl,
                   compatibleEnergyEvolution = compatibleEnergy,
                   evolveTotalEnergy = evolveTotalEnergy,
                   XSPH = XSPH,
                   densityUpdate = densityUpdate,
                   HUpdate = HUpdate)
else:
    hydro = SPH(dataBase = db,
                W = WT,
                filter = filter,
                cfl = cfl,
                compatibleEnergyEvolution = compatibleEnergy,
                evolveTotalEnergy = evolveTotalEnergy,
                densityUpdate = densityUpdate,
                HUpdate = HUpdate,
                XSPH = XSPH,
                epsTensile = epsilonTensile,
                nTensile = nTensile)
output("hydro")
output("hydro.cfl")
output("hydro.useVelocityMagnitudeForDt")
output("hydro.HEvolution")
output("hydro.densityUpdate")
output("hydro.compatibleEnergyEvolution")
output("hydro.evolveTotalEnergy")

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
integrator.verbose = dtverbose
output("integrator")
output("integrator.havePhysicsPackage(hydro)")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.dtGrowth")
output("integrator.domainDecompositionIndependent")

#-------------------------------------------------------------------------------
# Select points to follow histories.
#-------------------------------------------------------------------------------
def tp106tracersample(nodes, indices):
    assert mpi.allreduce(len(indices), mpi.SUM) == 1
    if len(indices) == 1:
        i = indices[0]
        assert i < nodes.numInternalNodes
        pos = nodes.positions()
        vel = nodes.velocity()
        rho = nodes.massDensity()
        eps = nodes.specificThermalEnergy()
        P = ScalarField("pressure", nodes)
        nodes.pressure(P)
        result = [pos[i].x, vel[i].x, rho[i], eps[i], P[i]]
    else:
        result = []
    result = mpi.allreduce(result, mpi.SUM)
    return tuple(result)

# Al sample points.
histories = []
tracerNumber = 1
for nodes, samplePositions in ((nodesAl, (-0.0375, -0.0625, -0.9875)),
                               (nodesTa, ( 0.0375,  0.0625,  0.9875))):
    pos = nodes.positions()
    for x0 in samplePositions:
        thpt = [(abs(pos[i].x - x0), i) for i in range(nodes.numInternalNodes)] + [(1e100, -1)]
        thpt.sort()
        dxmin = mpi.allreduce(thpt[0][0], mpi.MIN)
        if thpt[0][0] == dxmin:
            i = thpt[0][1]
            indices = [i]
            sys.stderr.write("Tracer %i is node %i @ %s.\n" % (tracerNumber, i, pos[i]))
        else:
            indices = []
        histories.append(NodeHistory(nodes, indices, tp106tracersample, tracerOutputName % tracerNumber, 
                                     labels = ("pos", "vel", "rho", "eps", "P")))
        tracerNumber += 1

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
# Add the diagnostics to the controller.
#-------------------------------------------------------------------------------
for hist in histories:
    control.appendPeriodicWork(hist.sample, sampleFreq)
    hist.flushHistory()

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
if not steps is None:
    control.step(steps)
else:
    control.advance(goalTime)
    control.dropRestartFile()

#-------------------------------------------------------------------------------
# If requested, write out the state in a global ordering to a file.
#-------------------------------------------------------------------------------
if outputFile:
    from SpheralTestUtilities import multiSort
    state = State(db, integrator.physicsPackages())
    outputFile = os.path.join(dataDir, outputFile)
    pos = state.vectorFields(HydroFieldNames.position)
    rho = state.scalarFields(HydroFieldNames.massDensity)
    P = state.scalarFields(HydroFieldNames.pressure)
    vel = state.vectorFields(HydroFieldNames.velocity)
    eps = state.scalarFields(HydroFieldNames.specificThermalEnergy)
    Hfield = state.symTensorFields(HydroFieldNames.H)
    S = state.symTensorFields(SolidFieldNames.deviatoricStress)
    ps = state.scalarFields("plastic strain")
    xprof = mpi.reduce([x.x for x in internalValues(pos)], mpi.SUM)
    rhoprof = mpi.reduce(internalValues(rho), mpi.SUM)
    Pprof = mpi.reduce(internalValues(P), mpi.SUM)
    vprof = mpi.reduce([v.x for v in internalValues(vel)], mpi.SUM)
    epsprof = mpi.reduce(internalValues(eps), mpi.SUM)
    hprof = mpi.reduce([1.0/sqrt(H.Determinant()) for H in internalValues(Hfield)], mpi.SUM)
    sprof = mpi.reduce([x.xx for x in internalValues(S)], mpi.SUM)
    psprof = mpi.reduce(internalValues(ps), mpi.SUM)
    mof = mortonOrderIndices(db)
    mo = mpi.reduce(internalValues(mof), mpi.SUM)
    if mpi.rank == 0:
        multiSort(mo, xprof, rhoprof, Pprof, vprof, epsprof, hprof, sprof, psprof)
        f = open(outputFile, "w")
        f.write(("#" + 17*" %16s" + "\n") % ("x", "rho", "P", "v", "eps", "h", "S", psprof, "m", 
                                             "int(x)", "int(rho)", "int(P)", "int(v)", "int(eps)", "int(h)", "int(S)", "int(ps)"))
        for (xi, rhoi, Pi, vi, epsi, hi, si, psi, mi) in zip(xprof, rhoprof, Pprof, vprof, epsprof, hprof, sprof, psprof, mo):
            f.write((8*"%16.12e " + 9*"%i " + "\n") %
                    (xi, rhoi, Pi, vi, epsi, hi, si, psi, mi,
                     unpackElementUL(packElementDouble(xi)),
                     unpackElementUL(packElementDouble(rhoi)),
                     unpackElementUL(packElementDouble(Pi)),
                     unpackElementUL(packElementDouble(vi)),
                     unpackElementUL(packElementDouble(epsi)),
                     unpackElementUL(packElementDouble(hi)),
                     unpackElementUL(packElementDouble(si)),
                     unpackElementUL(packElementDouble(psi))))
        f.close()

#-------------------------------------------------------------------------------
# Plot the state.
#-------------------------------------------------------------------------------
if graphics:
    from SpheralGnuPlotUtilities import *
    state = State(db, integrator.physicsPackages())
    rhoPlot = plotFieldList(state.scalarFields("mass density"),
                            plotStyle="linespoints",
                            winTitle="rho @ %g %i" % (control.time(), mpi.procs))
    velPlot = plotFieldList(state.vectorFields("velocity"),
                            yFunction = "%s.x",
                            plotStyle="linespoints",
                            winTitle="vel @ %g %i" % (control.time(), mpi.procs))
    mPlot = plotFieldList(state.scalarFields("mass"),
                          plotStyle="linespoints",
                          winTitle="mass @ %g %i" % (control.time(), mpi.procs))
    PPlot = plotFieldList(state.scalarFields("pressure"),
                          plotStyle="linespoints",
                          winTitle="pressure @ %g %i" % (control.time(), mpi.procs))
    hPlot = plotFieldList(state.symTensorFields("H"),
                          yFunction = "1.0/%s.xx",
                          plotStyle="linespoints",
                          winTitle="h @ %g %i" % (control.time(), mpi.procs))
    psPlot = plotFieldList(state.scalarFields("plastic strain"),
                           plotStyle="linespoints",
                           winTitle="plastic strain @ %g %i" % (control.time(), mpi.procs))
