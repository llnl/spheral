#-------------------------------------------------------------------------------
# A 1-D acoustic wave test.  In this version we establish a standing wave
# in a reflecting box.
#-------------------------------------------------------------------------------
import os, shutil
from math import *
from Spheral1d import *
from SpheralTestUtilities import *
import mpi
import numpy as np
#import matplotlib.pyplot as plt

title("Standing wave propagation test.")

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

            nPerh = 1.51,

            Qconstructor = MonaghanGingoldViscosity,
            #Qconstructor = TensorMonaghanGingoldViscosity,
            Cl = 0.0,
            Cq = 0.0,
            linearInExpansion = False,
            Qlimiter = False,
            epsilon2 = 1e-2,
            hmin = 0.0001, 
            hmax = 0.1,
            cfl = 0.5,
            XSPH = False,
            epsilonTensile = 0.0,
            nTensile = 4,
            filter = 0.0,

            SVPH = False,
            CRKSPH = False,
            TSPH = False,
            IntegratorConstructor = CheapSynchronousRK2Integrator,
            steps = None,
            goalTime = 5.0,
            dt = 0.0001,
            dtMin = 1.0e-5, 
            dtMax = 0.1,
            dtGrowth = 2.0,
            dtverbose = False,
            rigorousBoundaries = False,
            maxSteps = None,
            statsStep = 1,
            smoothIters = 0,
            HUpdate = IdealH,
            densityUpdate = RigorousSumDensity,
            compatibleEnergy = True,
            gradhCorrection = True,
            linearConsistent = False,

            restoreCycle = None,
            restartStep = 10000,
            clearDirectories = True,
            dataDirBase = "dumps-planar-StandingWave-1d",
            outputFile = "StandingWave-planar-1d.gnu",
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
else:
    HydroConstructor = SPHHydro

dataDir = os.path.join(dataDirBase,
                       str(HydroConstructor).split("'")[1].split(".")[-1],
                       str(Qconstructor).split("'")[1].split(".")[-1],
                       "nx=%i" % nx1)
restartDir = os.path.join(dataDir, "restarts")
restartBaseName = os.path.join(restartDir, "StandingWave-planar-1d-%i" % nx1)

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
eos = IsothermalEquationOfStateMKS(cs2, mu)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
WT = TableKernel(BSplineKernel(), 10000)
output("WT")

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

L = x1 - x0
def Minterval(xi0, xi1):
    return rho1*((xi1 - xi0) +
                 A*L/(pi*kfreq)*(sin(pi*kfreq/L*(xi1 - x0)) -
                                 sin(pi*kfreq/L*(xi0 - x0))))

# Grab the analytic answer.
import StandingWaveSolution
dx = (x1 - x0)/nx1
h1 = nPerh*dx
cs = sqrt(cs2)
answer = StandingWaveSolution.StandingWaveSolution(eos, cs, rho1, x0, x1, A, kfreq, h1)

# Set the node positions, velocities, and densities.
pos = nodes1.positions()
vel = nodes1.velocity()
rho = nodes1.massDensity()
mass = nodes1.mass()
eps = nodes1.specificThermalEnergy()
xvals = [p.x for p in pos.internalValues()]
xans, vans, uans, rhoans, Pans, hans = answer.solution(0.0, xvals)
for i in range(nodes1.numInternalNodes):
    vel[i].x = vans[i]
    eps[i] = uans[i]

# imin = int((pos[0].x - x0)/dx + 0.5)
# imax = int((pos[-1].x - x0)/dx + 0.5) + 1
# for i in xrange(nodes1.numInternalNodes):
#     xi0 = pos[i].x - 0.5*dx
#     xi1 = pos[i].x + 0.5*dx
#     mass[i] = Minterval(xi0, xi1)
#     rho[i] = mass[i]/dx

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
xbc0 = ReflectingBoundary(xPlane0)
xbc1 = ReflectingBoundary(xPlane1)
for bc in (xbc0, xbc1):
    hydro.appendBoundary(bc)

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
# Plot the final state.
#-------------------------------------------------------------------------------
xlocal = [pos.x for pos in nodes1.positions().internalValues()]
xprof = mpi.reduce(xlocal, mpi.SUM)
if graphics == "gnu":
    from SpheralGnuPlotUtilities import *
    state = State(db, integrator.physicsPackages())
    rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotState(state)
    if mpi.rank == 0:
        plotAnswer(answer, control.time(), rhoPlot, velPlot, epsPlot, PPlot, HPlot, xprof)
        #plotAnswer(answer.plus, control.time(), rhoPlot, velPlot, epsPlot, PPlot, HPlot, xprof)
        #plotAnswer(answer.minus, control.time(), rhoPlot, velPlot, epsPlot, PPlot, HPlot, xprof)
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
        A0prof = mpi.reduce(hydro.A0()[0].internalValues(), mpi.SUM)
        Aprof = mpi.reduce(hydro.A()[0].internalValues(), mpi.SUM)
        Bprof = mpi.reduce([x.x for x in hydro.B()[0].internalValues()], mpi.SUM)
        labels += ["A0", "A", "B"]
        stuff += [A0prof, Aprof, Bprof]

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
        xmin, xmax = 0.1, 0.9
        for (name, data, ans) in [("Mass Density", rhoprof, rhoans),
                                  ("Pressure", Pprof, Pans),
                                  ("Velocity", vprof, vans),
                                  ("h       ", hprof, hans)]:
            assert len(data) == len(ans)
            error = [data[i] - ans[i] for i in range(len(data))]
            Pn = Pnorm.Pnorm(error, xprof)
            L1 = Pn.gridpnorm(1, xmin, xmax)
            L2 = Pn.gridpnorm(2, xmin, xmax)
            Linf = Pn.gridpnorm("inf", xmin, xmax)
            print("\t%s \t\t%g \t\t%g \t\t%g" % (name, L1, L2, Linf))
            if normOutputFile:
                f.write((3*"%16.12e ") % (L1, L2, Linf))
        if normOutputFile:
            f.write("\n")
            f.close()

if compatibleEnergy and abs(Eerror) > 1e-8:
    raise ValueError("Energy error outside allowed bounds.")
