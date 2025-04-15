#-------------------------------------------------------------------------------
# A 1-D acoustic wave test.  Here we propogate a simple sound wave round and
# round in a periodic box.  This specific example is based on the test case
# described in D.J. Price's dissertation as an example of the effect of the
# grad h terms.
#-------------------------------------------------------------------------------
import os, shutil, sys
from math import *
from Spheral1d import *
from Spheral import ScalarPairScalarFunctor as PairScalarFunctor
from SpheralTestUtilities import *
import mpi
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

title("Acoustic wave propagation test.")

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

            nPerh = 3.01,

            Cl = 1.0,
            Cq = 2.0,
            linearInExpansion = False,
            Qlimiter = False,
            epsilon2 = 1e-30,
            hmin = 1.0e-10,
            hmax = 0.1,
            cfl = 0.25,
            XSPH = False,
            epsilonTensile = 0.0,
            nTensile = 4,
            filter = 0.0,
            KernelConstructor = WendlandC2Kernel,
            order = 5,

            svph = False,
            crksph = False,
            psph = False,
            fsisph = False,
            solid = False,
            IntegratorConstructor = CheapSynchronousRK2Integrator,
            correctionOrder = LinearOrder,
            steps = None,
            goalTime = 5.0,
            dt = 1.0e-2,
            dtMin = 1.0e-10, 
            dtMax = 0.1,
            dtGrowth = 2.0,
            dtverbose = False,
            rigorousBoundaries = False,
            maxSteps = None,
            statsStep = 1,
            smoothIters = 0,
            HUpdate = IntegrateH,
            densityUpdate = RigorousSumDensity,
            compatibleEnergy = True,
            evolveTotalEnergy = False,
            gradhCorrection = True,
            linearConsistent = False,

            ftol = None,
            convergenceTolerance = None,
            maxIterations = None,
            maxAllowedDtMultiplier = None,
            beta = None,

            restoreCycle = None,
            restartStep = 10000,
            clearDirectories = True,
            dataDirBase = "dumps-planar-AcousticWave-1d",
            outputFile = "AcousticWave-planar-1d.gnu",
            normOutputFile = "Limited_asciiDump.dat",
            writeOutputLabel = True,

            graphics = True,

            checkReversibility = False,
            )

if svph:
    hydroname = "SVPH"
elif crksph:
    hydroname = "CRKSPH"
elif psph:
    hydroname = "PSPH"
elif fsisph:
    hydroname = "FSISPH"
else:
    hydroname = "SPH"
normOutputFile = hydroname+"_"+normOutputFile
dataDir = os.path.join(dataDirBase,
                       hydroname,
                       "nx=%i" % nx1)
restartDir = os.path.join(dataDir, "restarts")
restartBaseName = os.path.join(restartDir, "AcousticWave-planar-1d-%i" % nx1)

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
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
if KernelConstructor==NBSplineKernel:
  WT = TableKernel(NBSplineKernel(order), 10000)
else:
  WT = TableKernel(KernelConstructor(), 10000)
output("WT")
kernelExtent = WT.kernelExtent
output("WT")

#-------------------------------------------------------------------------------
# Make the NodeList.
#-------------------------------------------------------------------------------
if solid:
    makeNL = makeSolidNodeList
else:
    makeNL = makeFluidNodeList
nodes1 = makeNL("nodes1", eos,
                hmin = hmin,
                hmax = hmax,
                kernelExtent = kernelExtent,
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
from newtonRaphson import *
cs = sqrt(cs2)
pos = nodes1.positions()
vel = nodes1.velocity()
rho = nodes1.massDensity()
mass = nodes1.mass()
H = nodes1.Hfield()
dx = (x1 - x0)/nx1
xi = x0
for i in range(nodes1.numInternalNodes):
    func0 = MassFunctor(max(0.0, Mi[i] - mi))
    func1 = MassFunctor(Mi[i])
    xi0 = newtonRaphsonFindRoot(func0, xi, xi + 2.0*dx, 1.0e-18, 1.0e-18)
    xi1 = newtonRaphsonFindRoot(func1, xi, xi + 2.0*dx, 1.0e-18, 1.0e-18)
    rhoi0 = rho1*(1.0 + A*sin(twopi*kfreq*(xi0 - x0)/(x1 - x0)))
    rhoi1 = rho1*(1.0 + A*sin(twopi*kfreq*(xi1 - x0)/(x1 - x0)))
    xi = x0 + (x1 - x0)*(rhoi0*xi0 + rhoi1*xi1)/(rhoi0 + rhoi1)
    pos[i].x = xi
    vel[i].x = A*cs*sin(twopi*kfreq*(xi - x0)/(x1 - x0))
    rho[i] = rho1*(1.0 + A*sin(twopi*kfreq*(xi - x0)/(x1 - x0)))
    mass[i] = rho1*((xi1 - xi0) - A/(twopi*kfreq)*(cos(twopi*kfreq*xi1) - cos(twopi*kfreq*xi0)))
    H[i] *= rho[i]/rho1
# xi0 = 0.0
# dx0 = (x1 - x0)/nx1
# for i in xrange(nodes1.numInternalNodes):
#     dxi0 = dx0*(1.0 - A*sin(twopi*kfreq*(xi0 - x0)/(x1 - x0)))
#     xi = xi0 + 0.5*dxi0
#     pos[i].x = xi
#     vel[i].x = A*cs*sin(twopi*kfreq*(xi - x0)/(x1 - x0))
#     rho[i] = rho1*(1.0 + A*sin(twopi*kfreq*(xi - x0)/(x1 - x0)))
#     xi0 += dxi0

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
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
if svph:
    hydro = SVPH(dataBase = db,
                 W = WT,
                 cfl = cfl,
                 compatibleEnergyEvolution = compatibleEnergy,
                 evolveTotalEnergy = evolveTotalEnergy,
                 XSVPH = XSPH,
                 linearConsistent = linearConsistent,
                 densityUpdate = densityUpdate,
                 HUpdate = HUpdate,
                 xmin = Vector(-100.0),
                 xmax = Vector( 100.0))
elif crksph:
    hydro = CRKSPH(dataBase = db,
                   filter = filter,
                   cfl = cfl,
                   order = correctionOrder,
                   compatibleEnergyEvolution = compatibleEnergy,
                   evolveTotalEnergy = evolveTotalEnergy,
                   XSPH = XSPH,
                   densityUpdate = densityUpdate,
                   HUpdate = HUpdate)
elif psph:
    hydro = PSPH(dataBase = db,
                 cfl = cfl,
                 W = WT,
                 filter = filter,
                 compatibleEnergyEvolution = compatibleEnergy,
                 evolveTotalEnergy = evolveTotalEnergy,
                 densityUpdate = densityUpdate,
                 HUpdate = HUpdate,
                 XSPH = XSPH)
elif fsisph: 
    hydro = FSISPH(dataBase = db,
                   W = WT,
                   cfl = cfl,
                   sumDensityNodeLists = [nodes1],
                   compatibleEnergyEvolution = compatibleEnergy,          
                   evolveTotalEnergy = evolveTotalEnergy,
                   epsTensile = epsilonTensile)
else:
    hydro = SPH(dataBase = db,
                #Q=MonaghanGingoldViscosity(Cl,Cq),
                W = WT, 
                cfl = cfl,
                compatibleEnergyEvolution = compatibleEnergy,
                evolveTotalEnergy = evolveTotalEnergy,
                gradhCorrection = gradhCorrection,
                XSPH = XSPH,
                correctVelocityGradient=False,
                densityUpdate = densityUpdate,
                HUpdate = HUpdate,
                epsTensile = epsilonTensile,
                nTensile = nTensile)
output("hydro")

#-------------------------------------------------------------------------------
# Construct the artificial viscosity.
#-------------------------------------------------------------------------------
q = hydro.Q
q.Cl = Cl
q.Cq = Cq
q.epsilon2 = epsilon2
q.limiter = Qlimiter
output("q")
output("q.Cl")
output("q.Cq")
output("q.epsilon2")
output("q.limiter")

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
integrator = IntegratorConstructor(db, [hydro])
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
output("integrator.verbose")

# Special stuff for implicit integrators
if isinstance(integrator, ImplicitIntegrator):
    if beta:
        integrator.beta = beta
    if convergenceTolerance:
        integrator.convergenceTolerance = convergenceTolerance
    if maxIterations:
        integrator.maxIterations = maxIterations
    if maxAllowedDtMultiplier:
        integrator.maxAllowedDtMultiplier = maxAllowedDtMultiplier
    output("integrator.beta")
    output("integrator.convergenceTolerance")
    output("integrator.maxIterations")
    output("integrator.maxAllowedDtMultiplier")
if isinstance(integrator, BackwardEulerIntegrator):
    if ftol:
        integrator.ftol = ftol
    output("integrator.ftol")

def printTotalEnergy(cycle,time,dt):
    Etot=0.0
    Etherm=0.0
    Ekinetic=0.0
    Vtot = 0.0
    Ptot = Vector.zero
    mass00=db.fluidMass
    vel00=db.fluidVelocity
    eps00=db.fluidSpecificThermalEnergy
    rho00 = db.fluidMassDensity
    nodeLists = db.nodeLists
    for nodelisti in range(db.numNodeLists):

        for i in range(nodeLists[nodelisti].numInternalNodes):
            #print([mass00(nodelisti,i),eps00(nodelisti,i),vel00(nodelisti,i)])
            Vtot += mass00(nodelisti,i)/rho00(nodelisti,i)
            #Ptot += mass00(nodelisti,i)*vel00(nodelisti,i)
            Etot += mass00(nodelisti,i)*(0.5*vel00(nodelisti,i).magnitude2()+eps00(nodelisti,i))
            Etherm += mass00(nodelisti,i)*(eps00(nodelisti,i))
            Ekinetic += mass00(nodelisti,i)*(0.5*vel00(nodelisti,i).magnitude2())
    #Ptot = mpi.allreduce(Ptot,mpi.SUM)
    Vtot = mpi.allreduce(Vtot,mpi.SUM)
    Etot = mpi.allreduce(Etot,mpi.SUM)
    Etherm = mpi.allreduce(Etherm,mpi.SUM)
    Ekinetic = mpi.allreduce(Ekinetic,mpi.SUM)
    print((" TOTAL   VOLUME   : %.15g" % Vtot))
    #print([" TOTAL   MOMENTUM : ",Ptot])
    print((" TOTAL   ENERGY   : %.15g" % Etot))
    print((" KINETIC ENERGY   : %.15g" % Ekinetic))
    print((" THERMAL ENERGY   : %.15g" % Etherm))

#-------------------------------------------------------------------------------
# Make the problem controller.
#-------------------------------------------------------------------------------
print("Making controller.")
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle,
                            periodicWork=[(printTotalEnergy,1)])
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
# Compute the analytic answer.
#-------------------------------------------------------------------------------
import AcousticWaveSolution
xlocal = [pos.x for pos in nodes1.positions().internalValues()]
xprof = np.array(mpi.reduce(xlocal, mpi.SUM))
dx = (x1 - x0)/nx1
h1 = nPerh*dx
answer = AcousticWaveSolution.AcousticWaveSolution(eos, 
                                                   cs=cs, 
                                                   rho0=rho1, 
                                                   x0=x0, 
                                                   x1=x1, 
                                                   A=A, 
                                                   k=twopi*kfreq, 
                                                   h0=h1)
#print "\n\nPERIOD=",1.0/(kfreq*cs)

### Compute the simulated specific entropy.
##rho = mpi.allreduce(nodes1.massDensity().internalValues(), mpi.SUM)
##P = mpi.allreduce(nodes1.pressure().internalValues(), mpi.SUM)
##A = [Pi/rhoi**gamma for (Pi, rhoi) in zip(P, rho)]

### The analytic solution for the simulated entropy.
##xans, vans, uans, rhoans, Pans, hans = answer.solution(control.time(), xprof)
##Aans = [Pi/rhoi**gamma for (Pi, rhoi) in zip(Pans,  rhoans)]

#-------------------------------------------------------------------------------
# Plot the final state.
#-------------------------------------------------------------------------------

if graphics:
    from SpheralMatplotlib import *

    rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotState(db)
    plotAnswer(answer, control.time(),
               rhoPlot = rhoPlot,
               velPlot = velPlot,
               epsPlot = epsPlot,
               PPlot = PPlot,
               HPlot = HPlot)
    pE = plotEHistory(control.conserve)

    #csPlot = plotFieldList(cs, winTitle="Sound speed")
    #csPlot.plot(xans, csAns, "k-",
    #            label = "Analytic")

    #APlot = newFigure()
    #APlot.plot(xprof, A, "ro", label="Simulation")
    #APlot.plot(xans, Aans, "k-", label="Analytic")
    #plt.title("A entropy")

    plots = [(rhoPlot, "Sod-planar-rho.png"),
             (velPlot, "Sod-planar-vel.png"),
             (epsPlot, "Sod-planar-eps.png"),
             (PPlot, "Sod-planar-P.png"),
             (HPlot, "Sod-planar-h.png")]
             #(csPlot, "Sod-planar-cs.png"),
             #(APlot, "Sod-planar-entropy.png")]

    # Make hardcopies of the plots.
    for p, filename in plots:
        savefig(p, os.path.join(dataDir, filename))
Eerror = (control.conserve.EHistory[-1] - control.conserve.EHistory[1])/control.conserve.EHistory[1]
print("Total energy error: %g" % Eerror)

#-------------------------------------------------------------------------------
# If requested, write out the state in a global ordering to a file.
#-------------------------------------------------------------------------------
if outputFile:
    outputFile = os.path.join(dataDir, outputFile)
    from SpheralTestUtilities import multiSort
    mprof = np.array(mpi.reduce(nodes1.mass().internalValues(), mpi.SUM))
    rhoprof = np.array(mpi.reduce(nodes1.massDensity().internalValues(), mpi.SUM))
    P = ScalarField("pressure", nodes1)
    nodes1.pressure(P)
    Pprof = np.array(mpi.reduce(P.internalValues(), mpi.SUM))
    vprof = np.array(mpi.reduce([v.x for v in nodes1.velocity().internalValues()], mpi.SUM))
    epsprof = np.array(mpi.reduce(nodes1.specificThermalEnergy().internalValues(), mpi.SUM))
    hprof = np.array(mpi.reduce([1.0/H.xx for H in nodes1.Hfield().internalValues()], mpi.SUM))
    xans, vans, uans, rhoans, Pans, hans = answer.solution(control.time(), xprof)

    labels = ["x", "m", "rho", "P", "v", "eps", "h", 
              "rhoans", "Pans", "vans", "epsans", "hans"]
    stuff = [xprof, mprof, rhoprof, Pprof, vprof, epsprof, hprof, 
             rhoans, Pans, vans, uans, hans]


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
        xmin, xmax = x0, x1
        pickleDumpDict =[]
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
            # if name == "Mass Density":
            #     pickleDumpL1 = L1

        if normOutputFile:
            f.write("\n")
            f.close()

        # keep running tally of density norms in single
        # file for refinement study
        # filename = hydroname+'densityNorm.pkl'
        # if os.path.exists(filename):

        #     f=open(filename,'r')
        #     data=pickle.load(f)
        #     data.append(pickleDumpL1)
        #     f.close()
        # else:
        #     data = [pickleDumpL1]

        # with open(filename, 'wb') as f:
        #     pickle.dump(data,f)


if compatibleEnergy and abs(Eerror) > 1e-5:
    raise ValueError("Energy error outside allowed bounds.")




    
    
