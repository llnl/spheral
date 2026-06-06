#-------------------------------------------------------------------------------
# A 1-D acoustic wave test.  Here we propogate a simple sound wave round and
# round in a periodic box.  This specific example is based on the test case
# described in D.J. Price's dissertation as an example of the effect of the
# grad h terms.
#-------------------------------------------------------------------------------
from math import *
from Spheral1d import *
from SpheralTestUtilities import *
import mpi
import numpy as np
#import matplotlib.pyplot as plt

#from CRKSPH_mod_package import *

title("Acoustic wave propagation test.")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(nx1 = 100,
            x0 = 0.0,
            x1 = 1.0,

            rho1 = 1.0,
            eps1 = 1.0,
            A = 1e-3,
            kfreq = 1.0,

            cs2 = 1.0,
            mu = 1.0,

            nPerh = 2.01,

            Qconstructor = MonaghanGingoldViscosity,
            #Qconstructor = TensorMonaghanGingoldViscosity,
            Cl = 0.0,
            Cq = 0.0,
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
            HEvolution = IdealH,
            densityUpdate = RigorousSumDensity,
            compatibleEnergy = True,
            gradhCorrection = True,
            linearConsistent = False,
            ComputeL1Norm = False,

            restoreCycle = None,
            restartStep = 10000,
            restartBaseName = "dumps-AcousticWave-1d",

            graphics = True,

            checkReversibility = False,
            )

#-------------------------------------------------------------------------------
# Material properties.
#-------------------------------------------------------------------------------
eos = IsothermalEquationOfStateMKS(cs2, mu)

##gamma = 5.0/3.0
##eps1 = cs2/(gamma*(gamma - 1.0))
##eos = GammaLawGasMKS(gamma, mu)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
WT = TableKernel(BSplineKernel(), 1000)
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

twopi = 2.0*pi
cs = sqrt(cs2)
pos = nodes1.positions()
vel = nodes1.velocity()
rho = nodes1.massDensity()
mass = nodes1.mass()
dx = (x1-x0)/nx1
for i in range(nodes1.numInternalNodes):
    rho[i] = rho1*(1.0 + A*sin(twopi*kfreq*(pos[i].x-x0)/(x1-x0)))
    vel[i].x = A*cs*sin(twopi*kfreq*(pos[i].x-x0)/(x1-x0))
    mass[i] = rho[i]*(x1-x0)/nx1
    '''
    #xi = i * (x1-x0) / nodes1.numInternalNodes + x0
    xi0= (i-0.5)* (x1-x0) / nodes1.numInternalNodes + x0
    xi1= (i+0.5)* (x1-x0) / nodes1.numInternalNodes + x0
    pos[i].x = pos[i].x + 0.5*A * cos(twopi*kfreq*(xi1-x0)/(x1-x0))-0.5*A * cos(twopi*kfreq*(xi0-x0)/(x1-x0))
    rho[i] = rho1*(1.0 + A*sin(twopi*kfreq*(pos[i].x-x0)/(x1-x0)))
    vel[i].x = A*cs*sin(twopi*kfreq*(pos[i].x-x0)/(x1-x0))
    '''
print("position 0 has {0} and position N has {1}".format((pos[0].x - x0),pos[nx1-1].x-x1))


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
q = Qconstructor(Cl, Cq)
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
    hydro = SVPHFacetedHydro(W = WT, 
                             Q = q,
                             cfl = cfl,
                             compatibleEnergyEvolution = compatibleEnergy,
                             XSVPH = XSPH,
                             linearConsistent = linearConsistent,
                             densityUpdate = densityUpdate,
                             HUpdate = HEvolution,
                             xmin = Vector(-100.0),
                             xmax = Vector( 100.0))
elif CRKSPH:
    hydro = CRKSPHHydro(W = WT, 
                        Q = q,
                        filter = filter,
                        cfl = cfl,
                        compatibleEnergyEvolution = compatibleEnergy,
                        XSPH = XSPH,
                        densityUpdate = densityUpdate,
                        HUpdate = HEvolution)
#CRKSPH_mod = CRKSPH_mod_package()

elif TSPH:
    hydro = TaylorSPHHydro(W = WT, 
                           Q = q,
                           cfl = cfl,
                           compatibleEnergyEvolution = compatibleEnergy,
                           XSPH = XSPH,
                           HUpdate = HEvolution)
else:
    hydro = SPHHydro(W = WT, 
                     Q = q,
                     cfl = cfl,
                     compatibleEnergyEvolution = compatibleEnergy,
                     gradhCorrection = gradhCorrection,
                     XSPH = XSPH,
                     densityUpdate = densityUpdate,
                     HUpdate = HEvolution,
                     epsTensile = epsilonTensile,
                     nTensile = nTensile)
output("hydro")

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
# Compute the analytic answer.
#-------------------------------------------------------------------------------
'''
omegat  = kfreq*sqrt(cs2)*goalTime
Pans    = []

for i in xrange(nodes1.numInternalNodes):
    xi = i * (x1-x0) / nodes1.numInternalNodes + x0
    Pans.append(A*sin(xi))

'''

import AcousticWaveSolutionMod
xlocal = [pos.x for pos in nodes1.positions().internalValues()]
xglobal = mpi.reduce(xlocal, mpi.SUM)
dx = (x1 - x0)/nx1
h1 = 1.0/(nPerh*dx)
cs = sqrt(cs2)
answer = AcousticWaveSolutionMod.AcousticWaveSolutionMod(eos, cs, rho1, x0, x1, A, twopi*kfreq, h1)
#print "\n\nPERIOD=",1.0/(kfreq*cs)

### Compute the simulated specific entropy.
##rho = mpi.allreduce(nodes1.massDensity().internalValues(), mpi.SUM)
##P = mpi.allreduce(nodes1.pressure().internalValues(), mpi.SUM)
##A = [Pi/rhoi**gamma for (Pi, rhoi) in zip(P, rho)]

### The analytic solution for the simulated entropy.
##xans, vans, uans, rhoans, Pans, hans = answer.solution(control.time(), xglobal)
##Aans = [Pi/rhoi**gamma for (Pi, rhoi) in zip(Pans,  rhoans)]

#-------------------------------------------------------------------------------
# Plot the final state.
#-------------------------------------------------------------------------------
if graphics:
    from SpheralMatplotlib import *
    state = State(db, integrator.physicsPackages())
    rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotState(state)
    #if mpi.rank == 0:
    if goalTime > 0:
        plotAnswer(answer, control.time(), rhoPlot, velPlot, epsPlot, PPlot, HPlot, xglobal)
    cs = state.scalarFields(HydroFieldNames.soundSpeed)
    csPlot = plotFieldList(cs, winTitle="Sound speed", colorNodeLists=False)

    EPlot = plotEHistory(control.conserve)
    massPlot = generateNewGnuPlot()
    massPlot.plot(nodes1.mass(),nodes1.positions())

    #PansPlot = generateNewGnuPlot()
    #PansPlot.plot(Pans)
    # Plot the correction terms.

    # Plot the grad h correction term (omega)

    if SVPH:
        volPlot = plotFieldList(hydro.volume(),
                                winTitle = "volume",
                                colorNodeLists = False)
    else:
        omegaPlot = plotFieldList(hydro.omegaGradh(),
                                  winTitle = "grad h correction",
                                  colorNodeLists = False)
    if ComputeL1Norm:
       xans, vans, uans, rhoans, Pans, hans = answer.solution(control.time(), xglobal)
       #rho = hydro.massDensity() 
       fieldList = state.scalarFields(HydroFieldNames.massDensity)
       #rho = field.internalValues()
       for field in fieldList:
          rho = [eval("%s" % "y") for y in field.internalValues()]
          #plt.figure()
          #plt.plot(xans,rhoans)
          #plt.scatter(xans,rho)
          #plt.xlim([np.min(xans),np.max(xans)])
          #plt.ylim([np.min(rhoans),np.max(rhoans)])
          #plt.show()
          diff=np.array(rho)-np.array(rhoans)
          L1Norm=(1.0/len(diff))*np.sum(np.abs(diff))
          print("\n\nL1Norm=",L1Norm, "\n\n")

Eerror = (control.conserve.EHistory[-1] - control.conserve.EHistory[0])/control.conserve.EHistory[0]
print("Total energy error: %g" % Eerror)
if abs(Eerror) > 1e-13:
    raise ValueError("Energy error outside allowed bounds.")
