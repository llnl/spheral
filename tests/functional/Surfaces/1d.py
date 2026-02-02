#-------------------------------------------------------------------------------
# Free expansion of gas into a vacuum.
#-------------------------------------------------------------------------------
from math import *
from Spheral1d import *
from SpheralTestUtilities import *
import mpi
import os, shutil
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

title("1d boundary test")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(nx1 = 100,
            x0 = 0.0,
            x1 = 1.0,

            rho1 = 1.0,
            eps1 = 1.0,

            gamma = 5.0/3.0,
            mu = 1.0,

            nPerh = 1.51,

            Qconstructor = MonaghanGingoldViscosity,
            Cl = 1.0,
            Cq = 1.0,
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
            goalTime = 0.5,
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
            densityUpdate = IntegrateDensity, # RigorousSumDensity,
            compatibleEnergy = True,
            gradhCorrection = True,
            linearConsistent = False,

            restoreCycle = None,
            restartStep = 10000,
            clearDirectories = True,
            dataDirBase = "dumps-1d",
            outputFile = "planar-1d.gnu",

            graphics = "gnu",
            vizCycle = None,
            vizTime = None,
            )

if SVPH:
    HydroConstructor = SVPHFacetedHydro
elif CRKSPH:
    HydroConstructor = CRKSPHHydro
    #Qconstructor = LimitedMonaghanGingoldViscosity
    #linearInExpansion = True
elif TSPH:
    HydroConstructor = TaylorSPHHydro
else:
    HydroConstructor = SPHHydro

dataDir = os.path.join(dataDirBase,
                       str(HydroConstructor).split("'")[1].split(".")[-1],
                       str(Qconstructor).split("'")[1].split(".")[-1],
                       "nx=%i" % nx1)
restartDir = os.path.join(dataDir, "restarts")
restartBaseName = os.path.join(restartDir, "AcousticWave-planar-1d-%i" % nx1)
vizDir = os.path.join(dataDir, "visit")
if vizTime is None and vizCycle is None:
  vizBaseName = None
else:
  vizBaseName = "1d-test-%i" % (nx1)

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
nNodesThisDomain1 = nodes1.numInternalNodes
output("nodes1.numNodes")

eps = nodes1.specificThermalEnergy(ScalarField("tmp", nodes1, eps1))

# Compute the summation correction for the density, and apply it to the mass per point.
mass = nodes1.mass()
dx = (x1 - x0)/nx1
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
xbc = ReflectingBoundary(xPlane0)
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
import SpheralVisitDump
vizMethod = SpheralVisitDump.dumpPhysicsState

print("Making controller.")
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle,
                            vizMethod = vizMethod,
                            vizBaseName = vizBaseName,
                            vizDir = vizDir,
                            vizStep = vizCycle,
                            vizTime = vizTime)
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
# Plot the final state.
#-------------------------------------------------------------------------------
xlocal = [pos.x for pos in nodes1.positions().internalValues()]
xprof = mpi.reduce(xlocal, mpi.SUM)
if graphics == "gnu":
    from SpheralGnuPlotUtilities import *
    state = State(db, integrator.physicsPackages())
    rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotState(state)
    cs = state.scalarFields(HydroFieldNames.soundSpeed)
    csPlot = plotFieldList(cs, winTitle="Sound speed", colorNodeLists=False)
    sn = state.vectorFields("Surface Normal")
    snPlot = plotFieldList(sn, yFunction = "%s.x", winTitle="Surface Normal")
    EPlot = plotEHistory(control.conserve)

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

    labels = ["x", "m", "rho", "P", "v", "eps", "h"]
    stuff = [xprof, mprof, rhoprof, Pprof, vprof, epsprof, hprof]
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

#------------------------------------------------------------------------------
# Check the energy error.
#------------------------------------------------------------------------------
Eerror = (control.conserve.EHistory[-1] - control.conserve.EHistory[0])/control.conserve.EHistory[0]
print("Total energy error: %g" % Eerror)
if compatibleEnergy and abs(Eerror) > 1e-10:
    raise ValueError("Energy error outside allowed bounds.")

