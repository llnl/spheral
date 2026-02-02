#-------------------------------------------------------------------------------
# The cylindrical Sedov test case (1-D).
#-------------------------------------------------------------------------------
import os, sys, shutil
from Spheral1d import *
from SpheralTestUtilities import *
#from SpheralGnuPlotUtilities import *

import mpi

title("1-D integrated hydro test -- planar Sedov problem")
#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(# discretization & domain
            nRadial = 50,
            rmin = 0.0,
            rmax = 1.0,

            # material properties
            rho0 = 1.0,
            eps0 = 0.0,
            Espike = 1.0,
            smoothSpike = True,
            topHatSpike = False,
            smoothSpikeScale = 0.5,
            gamma = 5.0/3.0,
            mu = 1.0,
            smallPressure = False,
            rhomin = 1e-10,

            #Kernel 
            HUpdate = IdealH,
            nPerh = 1.51,
            order = 5, 
            hmin = 1e-15,
            hmax = 1.0,

            # hydro type
            crksph = False,
            psph = False,
            gsph = False,
            fsisph = False,

            # hydro options
            solid = False,
            asph = False,                       # Selects the H update algorithm -- can be used with CRK, PSPH, SPH, etc.
            XSPH = False,
            evolveTotalEnergy = False,          # Only for SPH variants -- evolve total rather than specific energy
            compatibleEnergy = True,
            gradhCorrection = True,
            correctVelocityGradient = True,
            densityUpdate = RigorousSumDensity, # VolumeScaledDensity,
            filter = 0.0,
            
            # crksph options
            correctionOrder = LinearOrder,

            # gsph options
            RiemannGradientType = RiemannGradient, # (RiemannGradient,SPHGradient,HydroAccelerationGradient,OnlyDvDxGradient,MixedMethodGradient)
            linearReconstruction = True,

            # artificial viscosity options
            boolReduceViscosity = False,
            nh = 5.0,
            aMin = 0.1,
            aMax = 2.0,
            boolCullenViscosity = False,
            alphMax = 2.0,
            alphMin = 0.02,
            betaC = 0.7,
            betaD = 0.05,
            betaE = 1.0,
            fKern = 1.0/3.0,
            boolHopkinsCorrection = True,
            
            # integration
            cfl = 0.5,
            useVelocityMagnitudeForDt = False,
            steps = None,
            goalTime = None,
            goalRadius = 0.8,
            dt = 1e-8,
            dtMin = 1.0e-8,
            dtMax = None,
            dtGrowth = 2.0,
            dtverbose = False,
            maxSteps = None,
            statsStep = 1,
            smoothIters = 0,

            # IO
            restoreCycle = -1,
            restartStep = 1000,

            graphics = True,
            clearDirectories = False,
            dataDirBase = "dumps-planar-Sedov",
            outputFile = None,
            )

if smallPressure:
    P0 = 1.0e-6
    eps0 = P0/((gamma - 1.0)*rho0)
    print("WARNING: smallPressure specified, so setting eps0=%g" % eps0)

assert not(boolReduceViscosity and boolCullenViscosity)
assert not(gsph and (boolReduceViscosity or boolCullenViscosity))
assert not(fsisph and not solid)

# Figure out what our goal time should be.
import SedovAnalyticSolution
h0 = 1.0/nRadial*nPerh
answer = SedovAnalyticSolution.SedovSolution(nDim = 1,
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

# Scale the spike energy according to the boundary conditions we're using.
Espike *= 0.5

#-------------------------------------------------------------------------------
# Path names.
#-------------------------------------------------------------------------------
if crksph:
    hydroname = "CRKSPH"
elif psph:
    hydroname = "PSPH"
elif gsph:
    hydroname = "GSPH"
elif fsisph:
    hydroname = "FSISPH"
else:
    hydroname = "SPH"
if asph:
    hydroname = "A" + hydroname

dataDir = os.path.join(dataDirBase,
                       hydroname,
                       "nperh=%4.2f" % nPerh,
                       "XSPH=%s" % XSPH,
                       "densityUpdate=%s" % densityUpdate,
                       "compatibleEnergy=%s" % compatibleEnergy,
                       "Cullen=%s" % boolCullenViscosity,
                       "nr=%i" % nRadial)
restartDir = os.path.join(dataDir, "restarts")
restartBaseName = os.path.join(restartDir, "Sedov-planar-1d-%i" % nRadial)

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
WT = TableKernel(NBSplineKernel(order), 1000)
kernelExtent = WT.kernelExtent
output("WT")

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
                             nPerh = nPerh,
                             kernelExtent = kernelExtent,
                             rhoMin = rhomin)

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
pos = nodes1.positions()
vel = nodes1.velocity()
mass = nodes1.mass()
eps = nodes1.specificThermalEnergy()
H = nodes1.Hfield()
from DistributeNodes import distributeNodesInRange1d
distributeNodesInRange1d([(nodes1, nRadial, rho0, (0.0, rmax))])
output("mpi.reduce(nodes1.numInternalNodes, mpi.MIN)")
output("mpi.reduce(nodes1.numInternalNodes, mpi.MAX)")
output("mpi.reduce(nodes1.numInternalNodes, mpi.SUM)")

# Set the point source of energy.
Esum = 0.0
if smoothSpike or topHatSpike:
    Wsum = 0.0
    for nodeID in range(nodes1.numInternalNodes):
        Hi = H[nodeID]
        etaij = (Hi*pos[nodeID]).magnitude()
        if smoothSpike:
            Wi = WT.kernelValue(etaij/smoothSpikeScale, 1.0)
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
print("Initialized a total energy of", Eglobal)
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
if crksph:
    hydro = CRKSPH(dataBase = db,
                   order = correctionOrder,
                   filter = filter,
                   cfl = cfl,
                   compatibleEnergyEvolution = compatibleEnergy,
                   XSPH = XSPH,
                   densityUpdate = densityUpdate,
                   HUpdate = HUpdate,
                   ASPH = asph)
elif fsisph:
    hydro = FSISPH(dataBase = db,
                   W = WT,
                   cfl = cfl,
                   sumDensityNodeLists=[nodes1],                       
                   densityStabilizationCoefficient = 0.1,
                   specificThermalEnergyDiffusionCoefficient = 0.1,
                   useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                   compatibleEnergyEvolution = compatibleEnergy,
                   evolveTotalEnergy = evolveTotalEnergy,
                   linearCorrectGradients = correctVelocityGradient,
                   HUpdate = HUpdate) 
elif gsph:
    limiter = VanLeerLimiter()
    waveSpeed = DavisWaveSpeed()
    solver = HLLC(limiter,waveSpeed,linearReconstruction)
    hydro = GSPH(dataBase = db,
                riemannSolver = solver,
                W = WT,
                cfl=cfl,
                specificThermalEnergyDiffusionCoefficient = 0.00,
                compatibleEnergyEvolution = compatibleEnergy,
                correctVelocityGradient= correctVelocityGradient,
                evolveTotalEnergy = evolveTotalEnergy,
                XSPH = XSPH,
                ASPH = asph,
                densityUpdate=densityUpdate,
                HUpdate = HUpdate)
elif psph:
    hydro = PSPH(dataBase = db,
                      W = WT,
                      filter = filter,
                      cfl = cfl,
                      useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                      compatibleEnergyEvolution = compatibleEnergy,
                      evolveTotalEnergy = evolveTotalEnergy,
                      correctVelocityGradient = correctVelocityGradient,
                      densityUpdate = densityUpdate,
                      HUpdate = HUpdate,
                      XSPH = XSPH,
                      ASPH = asph)
else:
    hydro = SPH(dataBase = db,
                W = WT, 
                cfl = cfl,
                compatibleEnergyEvolution = compatibleEnergy,
                evolveTotalEnergy = evolveTotalEnergy,
                gradhCorrection = gradhCorrection,
                correctVelocityGradient = correctVelocityGradient,
                densityUpdate = densityUpdate,
                XSPH = XSPH,
                HUpdate = HUpdate,
                ASPH = asph)

packages = [hydro]
output("hydro")
output("hydro.cfl")
output("hydro.compatibleEnergyEvolution")
output("hydro.HEvolution")

if not gsph:
    q = hydro.Q
    output("q")
    output("q.Cl")
    output("q.Cq")
    output("q.epsilon2")
    output("q.limiter")
    output("q.balsaraShearCorrection")
    output("q.linearInExpansion")
    output("q.quadraticInExpansion")

    #-------------------------------------------------------------------------------
    # Construct the MMRV physics object.
    #-------------------------------------------------------------------------------
    if boolReduceViscosity:
        evolveReducingViscosityMultiplier = MorrisMonaghanReducingViscosity(nh,aMin,aMax)
        packages.append(evolveReducingViscosityMultiplier)
    elif boolCullenViscosity:
        evolveCullenViscosityMultiplier = CullenDehnenViscosity(WT,alphMax,alphMin,betaC,betaD,betaE,fKern,boolHopkinsCorrection)
        packages.append(evolveCullenViscosityMultiplier)

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
xPlane0 = Plane(Vector(0.0, 0.0), Vector(1.0, 0.0))
xbc0 = ReflectingBoundary(xPlane0)

for p in packages:
    p.appendBoundary(xbc0)

#-------------------------------------------------------------------------------
# Construct a time integrator, and add the one physics package.
#-------------------------------------------------------------------------------
integrator = CheapSynchronousRK2Integrator(db)
for p in packages:
    integrator.appendPhysicsPackage(p)
integrator.lastDt = dt
if dtMin:
    integrator.dtMin = dtMin
if dtMax:
    integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
integrator.verbose = dtverbose
output("integrator")
output("integrator.havePhysicsPackage(hydro)")
output("integrator.dtGrowth")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")

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
# Finally run the problem and plot the results.
#-------------------------------------------------------------------------------
if steps is None:
    control.advance(goalTime, maxSteps)
    if restoreCycle != control.totalSteps:
        control.dropRestartFile()
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
rho = mpi.allreduce(list(nodes1.massDensity().internalValues()), mpi.SUM)
v = mpi.allreduce([x.magnitude() for x in nodes1.velocity().internalValues()], mpi.SUM)
eps = mpi.allreduce(list(nodes1.specificThermalEnergy().internalValues()), mpi.SUM)
Pf = ScalarField("pressure", nodes1)
nodes1.pressure(Pf)
P = mpi.allreduce(list(Pf.internalValues()), mpi.SUM)
A = mpi.allreduce([Pi/(rhoi**gamma) for (Pi, rhoi) in zip(Pf.internalValues(), nodes1.massDensity().internalValues())], mpi.SUM)

Hinverse = db.newFluidSymTensorFieldList()
db.fluidHinverse(Hinverse)
hr = mpi.allreduce([x.xx for x in Hinverse[0].internalValues()], mpi.SUM)

Aans = None
if mpi.rank == 0:
    from SpheralTestUtilities import multiSort
    import Pnorm
    multiSort(r, rho, v, eps, P, A, hr)
    rans, vans, epsans, rhoans, Pans, Aans, hans = answer.solution(control.time(), r)
    print("\tQuantity \t\tL1 \t\t\tL2 \t\t\tLinf")
    #f = open("MCTesting.txt", "a")
    #f.write(("CL=%g, Cq=%g \t") %(Cl, Cq))
    for (name, data, ans) in [("Mass Density", rho, rhoans),
                              ("Pressure", P, Pans),
                              ("Velocity", v, vans),
                              ("Thermal E", eps, epsans),
                              ("Entropy", A, Aans),
                              ("hr", hr, hans)]:
        assert len(data) == len(ans)
        error = [data[i] - ans[i] for i in range(len(data))]
        Pn = Pnorm.Pnorm(error, r)
        L1 = Pn.gridpnorm(1, rmin, rmax)
        L2 = Pn.gridpnorm(2, rmin, rmax)
        Linf = Pn.gridpnorm("inf", rmin, rmax)
        print("\t%s \t\t%g \t\t%g \t\t%g" % (name, L1, L2, Linf))
        #f.write(("\t\t%g") % (L1))
    #f.write("\n")
Aans = mpi.bcast(Aans, 0)

#-------------------------------------------------------------------------------
# If requested, write out the state in a global ordering to a file.
#-------------------------------------------------------------------------------
if outputFile and mpi.rank == 0:
    outputFile = os.path.join(dataDir, outputFile)
    f = open(outputFile, "w")
    f.write(("# " + 14*"%15s " + "\n") % ("r", "x", "rho", "P", "v", "eps", "A", "hr",
                                          "rhoans", "Pans", "vans", "epsans", "Aans", "hrans"))
    for (ri, xi, rhoi, Pi, vi, epsi, Ai, hri, 
         rhoansi, Pansi, vansi, epsansi, Aansi, hansi)  in zip(r, xprof, rho, P, v, eps, A, hr,
                                                               rhoans, Pans, vans, epsans, Aans, hans):
         f.write((14*"%16.12e " + "\n") % (ri, xi, rhoi, Pi, vi, epsi, Ai, hri, 
                                           rhoansi, Pansi, vansi, epsansi, Aansi, hansi))
    f.close()

#-------------------------------------------------------------------------------
# Plot the final state.
#-------------------------------------------------------------------------------
if graphics:
    rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotRadialState(db)
    plotAnswer(answer, control.time(),
               rhoPlot, velPlot, epsPlot, PPlot, HPlot)
    plots = [(rhoPlot, "Sedov-planar-rho.png"),
             (velPlot, "Sedov-planar-vel.png"),
             (epsPlot, "Sedov-planar-eps.png"),
             (PPlot, "Sedov-planar-P.png"),
             (HPlot, "Sedov-planar-h.png")]

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
    plots.append((Aplot, "Sedov-planar-entropy.png"))

    if boolCullenViscosity:
        cullAlphaPlot = plotFieldList(q.ClMultiplier(),
                                      winTitle = "Cullen alpha")
        plots += [(cullAlphaPlot, "Sedov-planar-Cullen-alpha.png")]

    # Make hardcopies of the plots.
    for p, filename in plots:
        p.hardcopy(os.path.join(dataDir, filename), terminal="png")

