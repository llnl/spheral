#-------------------------------------------------------------------------------
# The spherical Sedov test case (3D in spherical coordinates)
#-------------------------------------------------------------------------------
import os, sys, shutil
from SphericalSpheral import *
from findLastRestart import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from GenerateNodeDistribution3d import *

import mpi

title("Spherical Sedov problem run in spherical coordinates")
#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(nr = 100,
            nPerh = 4.01,
            KernelConstructor = WendlandC4Kernel3d,
            kernelOrder = 5,

            rho0 = 1.0,
            eps0 = 0.0,
            smallPressure = False,
            Espike = 1.0,
            smoothSpike = True,
            topHatSpike = False,
            smoothSpikeScale = 4.0,
            gamma = 5.0/3.0,
            mu = 1.0,
            rhomin = 1e-10,

            # kernel
            HUpdate = IdealH,
            hmin = 1e-15,
            hmax = 1.0,

            # hydro parameters
            solid = False,
            XSPH = False,
            evolveTotalEnergy = False,
            compatibleEnergy = True,
            gradhCorrection = True,
            correctVelocityGradient = True,
            densityUpdate = RigorousSumDensity, 
            Qself = 2.0,

            # Integration
            IntegratorConstructor = VerletIntegrator,
            cfl = 0.5,
            useVelocityMagnitudeForDt = False,
            steps = None,
            goalTime = None,
            goalRadius = 0.8,
            dt = 1e-8,
            dtMin = 1.0e-10,
            dtMax = None,
            dtGrowth = 2.0,
            maxSteps = None,
            statsStep = 1,

            # IO
            restoreCycle = -1,
            restartStep = 1000,

            graphics = True,
            clearDirectories = True,
            dataDirBase = "dumps-spherical-Sedov",
            outputFile = None,
            )

if smallPressure:
    P0 = 1.0e-6
    eps0 = P0/((gamma - 1.0)*rho0)
    print("WARNING: smallPressure specified, so setting eps0=%g" % eps0)

# Figure out what our goal time should be.
import SedovAnalyticSolution
h0 = 1.0/nr*nPerh
answer = SedovAnalyticSolution.SedovSolution(nDim = 3,
                                             gamma = gamma,
                                             rho0 = rho0,
                                             E0 = Espike,
                                             h0 = h0)
if goalTime is None:
    assert not goalRadius is None
    nu1 = 1.0/(answer.nu + 2.0)
    nu2 = 2.0*nu1
    goalTime = (goalRadius*(answer.alpha*rho0/Espike)**nu1)**(1.0/nu2)
vs, r2, v2, rho2, P2 = answer.shockState(goalTime)
print("Predicted shock position %g at goal time %g." % (r2, goalTime))

#-------------------------------------------------------------------------------
# Path names.
#-------------------------------------------------------------------------------
hydroname = ""
if solid:
    hydroname += "Solid"
hydroname += "SPH"

dataDir = os.path.join(dataDirBase,
                       "1d",
                       hydroname,
                       "nperh=%4.2f" % nPerh,
                       "XSPH=%s" % XSPH,
                       "densityUpdate=%s" % densityUpdate,
                       "compatibleEnergy=%s" % compatibleEnergy,
                       "nr=%i" % nr)
restartDir = os.path.join(dataDir, "restarts")
restartBaseName = os.path.join(restartDir, "Sedov-spherical-1d")

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
eos = GammaLawGasMKS(gamma, mu)

#-------------------------------------------------------------------------------
# Create our interpolation kernels -- one for normal hydro interactions, and
# one for use with the artificial viscosity
#-------------------------------------------------------------------------------
if KernelConstructor is NBSplineKernel3d:
  W = NBSplineKernel3d(kernelOrder)
else:
  W = KernelConstructor()
output("W")
kernelExtent = W.kernelExtent

#-------------------------------------------------------------------------------
# Create a NodeList and associated Neighbor object.
#-------------------------------------------------------------------------------
if solid:
    nodeListConstructor = makeSolidNodeList
else:
    nodeListConstructor = makeFluidNodeList
nodes1 = nodeListConstructor("nodes1", eos, 
                             hmin = hmin,
                             hmax = hmax,
                             kernelExtent = kernelExtent,
                             nPerh = nPerh,
                             rhoMin = rhomin)

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
from GenerateSphericalNodeDistribution1d import GenerateSphericalNodeDistribution1d
from PeanoHilbertDistributeNodes import distributeNodes1d
generator = GenerateSphericalNodeDistribution1d(nr = nr,
                                                rho = rho0,
                                                rmin = 0.0,
                                                rmax = 1.0,
                                                nNodePerh = nPerh)
distributeNodes1d((nodes1, generator))
output("mpi.reduce(nodes1.numInternalNodes, mpi.MIN)")
output("mpi.reduce(nodes1.numInternalNodes, mpi.MAX)")
output("mpi.reduce(nodes1.numInternalNodes, mpi.SUM)")

# Set the point source of energy.
Espike /= 4.0*pi          # Cause our node volume element has the r^2 but not 4 pi
pos = nodes1.positions()
vel = nodes1.velocity()
mass = nodes1.mass()
eps = nodes1.specificThermalEnergy()
H = nodes1.Hfield()
Esum = 0.0
if smoothSpike or topHatSpike:
    Wsum = 0.0
    for nodeID in range(nodes1.numInternalNodes):
        Hi = H[nodeID]
        etaij = (Hi*pos[nodeID]).magnitude()
        if smoothSpike:
            Wi = W.kernelValue(etaij/smoothSpikeScale, 1.0)
        else:
            if etaij < smoothSpikeScale*kernelExtent:
                Wi = 1.0
            else:
                Wi = 0.0
        Ei = Wi*Espike
        epsi = Ei/mass[nodeID]
        eps[nodeID] = epsi
        Wsum += Wi
    Wsum = mpi.allreduce(Wsum, mpi.SUM)
    assert Wsum > 0.0
    for nodeID in range(nodes1.numInternalNodes):
        eps[nodeID] /= Wsum
        Esum += eps[nodeID]*mass[nodeID]
        eps[nodeID] += eps0
else:
    i = -1
    rmin = 1e50
    for nodeID in range(nodes1.numInternalNodes):
        rij = pos[nodeID].magnitude()
        if rij < rmin:
            i = nodeID
            rmin = rij
        eps[nodeID] = eps0
    rminglobal = mpi.allreduce(rmin, mpi.MIN)
    if fuzzyEqual(rmin, rminglobal):
        assert i >= 0 and i < nodes1.numInternalNodes
        eps[i] += Espike/mass[i]
        Esum += Espike
Eglobal = mpi.allreduce(Esum, mpi.SUM)
print("Initialized a total energy of", Eglobal * 4.0*pi)
assert fuzzyEqual(Eglobal, Espike)

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
hydro = SPH(dataBase = db,
            W = W, 
            cfl = cfl,
            compatibleEnergyEvolution = compatibleEnergy,
            evolveTotalEnergy = evolveTotalEnergy,
            gradhCorrection = gradhCorrection,
            correctVelocityGradient = correctVelocityGradient,
            densityUpdate = densityUpdate,
            XSPH = XSPH,
            HUpdate = HUpdate)
hydro.Qself = Qself
output("hydro")
output("hydro.cfl")
output("hydro.compatibleEnergyEvolution")
output("hydro.XSPH")
output("hydro.densityUpdate")
output("hydro.HEvolution")
output("hydro.Qself")

packages = [hydro]

#-------------------------------------------------------------------------------
# Construct a time integrator, and add the one physics package.
#-------------------------------------------------------------------------------
integrator = IntegratorConstructor(db)
for p in packages:
    integrator.appendPhysicsPackage(p)
integrator.lastDt = dt
if dtMin:
    integrator.dtMin = dtMin
if dtMax:
    integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
output("integrator")
output("integrator.havePhysicsPackage(hydro)")
output("integrator.dtGrowth")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")

#-------------------------------------------------------------------------------
# Build the controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, hydro.kernel,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle,
                            SPH = (not ASPH))
output("control")

#-------------------------------------------------------------------------------
# Finally run the problem and plot the results.
#-------------------------------------------------------------------------------
if steps is None:
    control.advance(goalTime, maxSteps)
else:
    control.step(steps)

# Output the energy conservation.
print("Energy conservation: ", ((control.conserve.EHistory[-1] -
                                 control.conserve.EHistory[0])/
                                control.conserve.EHistory[0]))

#-------------------------------------------------------------------------------
# Generate some error metrics comparing to the analytic solution.
#-------------------------------------------------------------------------------
# Report the error norms.
rmin, rmax = 0.0, 0.95
r = mpi.allreduce([x.magnitude() for x in nodes1.positions().internalValues()], mpi.SUM)
xprof = mpi.allreduce([x.x for x in nodes1.positions().internalValues()], mpi.SUM)
yprof = mpi.allreduce([x.y for x in nodes1.positions().internalValues()], mpi.SUM)
zprof = mpi.allreduce([x.z for x in nodes1.positions().internalValues()], mpi.SUM)
rho = mpi.allreduce(list(nodes1.massDensity().internalValues()), mpi.SUM)
mass = mpi.allreduce(list(nodes1.mass().internalValues()), mpi.SUM)
v = mpi.allreduce([x.magnitude() for x in nodes1.velocity().internalValues()], mpi.SUM)
eps = mpi.allreduce(list(nodes1.specificThermalEnergy().internalValues()), mpi.SUM)
Pf = ScalarField("pressure", nodes1)
nodes1.pressure(Pf)
P = mpi.allreduce(list(Pf.internalValues()), mpi.SUM)
A = mpi.allreduce([Pi/(rhoi**gamma) for (Pi, rhoi) in zip(Pf.internalValues(), nodes1.massDensity().internalValues())], mpi.SUM)

rans, vans, epsans, rhoans, Pans, Aans, hans = answer.solution(control.time(), r)
from SpheralTestUtilities import multiSort
multiSort(r, rho, v, eps, P, A, rhoans, vans, epsans, Pans, hans)

if mpi.rank == 0:
    from SpheralTestUtilities import multiSort
    import Pnorm
    multiSort(r, rho, v, eps, P, A)
    print("\tQuantity \t\tL1 \t\t\tL2 \t\t\tLinf")
    for (name, data, ans) in [("Mass Density", rho, rhoans),
                              ("Pressure", P, Pans),
                              ("Velocity", v, vans),
                              ("Thermal E", eps, epsans),
                              ("Entropy", A, Aans)]:
        assert len(data) == len(ans)
        error = [data[i] - ans[i] for i in range(len(data))]
        Pn = Pnorm.Pnorm(error, r)
        L1 = Pn.gridpnorm(1, rmin, rmax)
        L2 = Pn.gridpnorm(2, rmin, rmax)
        Linf = Pn.gridpnorm("inf", rmin, rmax)
        print("\t%s \t\t%g \t\t%g \t\t%g" % (name, L1, L2, Linf))

#-------------------------------------------------------------------------------
# If requested, write out the state in a global ordering to a file.
#-------------------------------------------------------------------------------
if outputFile and mpi.rank == 0:
    outputFile = os.path.join(dataDir, outputFile)
    f = open(outputFile, "w")
    f.write(("# " + 16*"%15s " + "\n") % ("r", "x", "y", "z", "rho", "m", "P", "v", "eps", "A",
                                          "rhoans", "Pans", "vans", "epsans", "Aans", "hrans"))
    for (ri, xi, yi, zi, rhoi, mi, Pi, vi, epsi, Ai, 
         rhoansi, Pansi, vansi, epsansi, Aansi, hansi)  in zip(r, xprof, yprof, zprof, rho, mass, P, v, eps, A,
                                                               rhoans, Pans, vans, epsans, Aans, hans):
         f.write((16*"%16.12e " + "\n") % (ri, xi, yi, zi, rhoi, mi, Pi, vi, epsi, Ai,
                                           rhoansi, Pansi, vansi, epsansi, Aansi, hansi))
    f.close()

#-------------------------------------------------------------------------------
# Plot the final state.
#-------------------------------------------------------------------------------
if graphics:
    rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotRadialState(db)
    plotAnswer(answer, control.time(),
               rhoPlot, velPlot, epsPlot, PPlot, HPlot)
    plots = [(rhoPlot, "Sedov-spherical-rho.png"),
             (velPlot, "Sedov-spherical-vel.png"),
             (epsPlot, "Sedov-spherical-eps.png"),
             (PPlot, "Sedov-spherical-P.png"),
             (HPlot, "Sedov-spherical-h.png")]

    # Plot the specific entropy.
    AsimData = Gnuplot.Data(xprof, A,
                            with_ = "points",
                            title = "Simulation",
                            inline = True)
    AansData = Gnuplot.Data(xprof, Aans,
                            with_ = "lines",
                            title = "Solution",
                            inline = True)
    
    Aplot = generateNewGnuPlot()
    Aplot.plot(AsimData)
    Aplot.replot(AansData)
    Aplot.title("Specific entropy")
    Aplot.refresh()
    plots.append((Aplot, "Sedov-spherical-entropy.png"))

    # Make hardcopies of the plots.
    for p, filename in plots:
        p.hardcopy(os.path.join(dataDir, filename), terminal="png")

