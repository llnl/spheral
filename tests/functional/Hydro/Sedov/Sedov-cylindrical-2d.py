#-------------------------------------------------------------------------------
# The cylindrical Sedov test case (2-D).
#-------------------------------------------------------------------------------
import os, sys, shutil, mpi
from Spheral2d import *
from SpheralTestUtilities import *
#from SpheralGnuPlotUtilities import *
from GenerateNodeDistribution2d import *
from CubicNodeGenerator import GenerateSquareNodeDistribution

if mpi.procs > 1:
    from VoronoiDistributeNodes import distributeNodes2d
    #from PeanoHilbertDistributeNodes import distributeNodes2d
else:
    from DistributeNodes import distributeNodes2d


title("2-D integrated hydro test -- planar Sedov problem")
#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(seed = "lattice",

            thetaFactor = 0.5,
            azimuthalOffsetFraction = 0.0,
            nRadial = 64,
            nTheta = 64,
            rmin = 0.0,
            rmax = 1.0,

            rho0 = 1.0,
            eps0 = 0.0,
            smallPressure = False,
            Espike = 1.0,
            smoothSpike = True,
            topHatSpike = False,
            smoothSpikeScale = 0.5,
            gamma = 5.0/3.0,
            mu = 1.0,
            rhomin = 1e-10,

            # kernel
            kernelConstructor = NBSplineKernel,
            nPerh = 1.00,
            order = 3,
            HUpdate = IdealH,
            hmin = 1e-15,
            hmax = 1.0,

            # hydros
            crksph = False,
            psph = False,
            fsisph = False,
            gsph = False,
            mfm = False,
            mfv = False,

            # hydro options
            solid = False,
            asph = False,
            XSPH = False,    
            evolveTotalEnergy = False,    
            compatibleEnergy = True,
            gradhCorrection = True,
            correctVelocityGradient = False,
            densityUpdate = RigorousSumDensity, 
            filter = 0.0,

            # crksph options
            correctionOrder = LinearOrder,
            volumeType = RKSumVolume,

            # gsph options
            RiemannGradientType = SPHSameTimeGradient, # (RiemannGradient,SPHGradient,HydroAccelerationGradient,OnlyDvDxGradient,MixedMethodGradient)
            linearReconstruction = True,
            nodeMotionType = "eulerian",

            # Artifical Viscosity
            boolReduceViscosity = False,
            nh = 5.0,
            aMin = 0.1,
            aMax = 2.0,
            Qhmult = 1.0,
            boolCullenViscosity = False,
            alphMax = 2.0,
            alphMin = 0.02,
            betaC = 0.7,
            betaD = 0.05,
            betaE = 1.0,
            fKern = 1.0/3.0,
            boolHopkinsCorrection = True,
    
            # Integration
            IntegratorConstructor = CheapSynchronousRK2Integrator,
            cfl = 0.25,
            steps = None,
            goalTime = None,
            goalRadius = 0.8,
            dt = 1e-8,
            dtMin = 1.0e-8,
            dtMax = None,
            dtGrowth = 2.0,
            maxSteps = None,
            statsStep = 1,
            smoothIters = 0,
            useVelocityMagnitudeForDt = False,
            dtverbose = False,

            # IO
            vizCycle = None,
            vizDerivs = False,
            vizTime = 0.1,
            restoreCycle = -1,
            restartStep = 1000,

            graphics = False,
            useVoronoiOutput = False,
            clearDirectories = False,
            dataDirBase = "dumps-cylindrical-Sedov",
            outputFile = None,
            serialDump=True,
            )

if smallPressure:
    P0 = 1.0e-6
    eps0 = P0/((gamma - 1.0)*rho0)
    print("WARNING: smallPressure specified, so setting eps0=%g" % eps0)

assert not(boolReduceViscosity and boolCullenViscosity)
assert thetaFactor in (0.5, 1.0, 2.0)
assert not((gsph or mfm or mfv) and (boolReduceViscosity or boolCullenViscosity))
assert not(fsisph and not solid)
theta = thetaFactor * pi

# Figure out what our goal time should be.
import SedovAnalyticSolution
h0 = 1.0/nRadial*nPerh
answer = SedovAnalyticSolution.SedovSolution(nDim = 2,
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
if thetaFactor == 0.5:
    Espike *= 0.25
elif thetaFactor == 1.0:
    Espike *= 0.5

#-------------------------------------------------------------------------------
# Path names.
#-------------------------------------------------------------------------------
if crksph:
    hydroname = "CRKSPH"
elif psph:
    hydroname = "PSPH"
elif fsisph:
    hydroname = "FSISPH"
elif gsph:
    hydroname = "GSPH"
elif mfm:
    hydroname = "MFM"
elif mfv:
    hydroname = "MFV"
    if nodeMotionType == "eulerian":
        hydroname += "_{0}".format(nodeMotionType)
        nodeMotionType = NodeMotionType.Eulerian
    elif nodeMotionType == "lagrangian":
        hydroname += "_{0}".format(nodeMotionType)
        nodeMotionType = NodeMotionType.Lagrangian
    elif nodeMotionType == "fician":
        hydroname += "_{0}".format(nodeMotionType)
        nodeMotionType = NodeMotionType.Fician
    elif nodeMotionType == "xsph":
        hydroname += "_{0}".format(nodeMotionType)
        nodeMotionType = NodeMotionType.XSPH
    else:
        raise ValueError ("WHAT DID YOU DO!!!???")
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
                       "seed=%s" % seed,
                       "nr=%i_nt=%i" % (nRadial, nTheta))
restartDir = os.path.join(dataDir, "restarts")
vizDir = os.path.join(dataDir, "visit")
restartBaseName = os.path.join(restartDir, "Sedov-cylindrical-2d-%i" % nRadial)

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
if mpi.rank == 0:
    if clearDirectories and os.path.exists(dataDir):
        shutil.rmtree(dataDir)
    if not os.path.exists(restartDir):
        os.makedirs(restartDir)
    if not os.path.exists(vizDir):
        os.makedirs(vizDir)
mpi.barrier()

#-------------------------------------------------------------------------------
# Material properties.
#-------------------------------------------------------------------------------
eos = GammaLawGasMKS(gamma, mu)

#-------------------------------------------------------------------------------
# Create our interpolation kernels -- one for normal hydro interactions, and
# one for use with the artificial viscosity
#-------------------------------------------------------------------------------
if kernelConstructor == NBSplineKernel:
    WT = TableKernel(NBSplineKernel(order), 1000)
else:
    WT = TableKernel(kernelConstructor(), 1000)

    output("WT")
    kernelExtent = WT.kernelExtent
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
pos = nodes1.positions()
vel = nodes1.velocity()
mass = nodes1.mass()
eps = nodes1.specificThermalEnergy()
H = nodes1.Hfield()
if thetaFactor == 0.5:
    xmin = (0.0, 0.0)
    xmax = (1.0, 1.0)
elif thetaFactor == 1.0:
    xmin = (-1.0, 0.0)
    xmax = (1.0, 1.0)
else:
    xmin = (-1.0, -1.0)
    xmax = (1.0, 1.0)

if seed == "square":
    generator = GenerateNodeDistribution2d(nRadial, nTheta, rho0, "lattice",
                                           rmin = rmin,
                                           rmax = rmax,
                                           xmin = xmin,
                                           xmax = xmax,
                                           theta = theta,
                                           #azimuthalOffsetFraction = azimuthalOffsetFraction,
                                           nNodePerh = nPerh,
                                           SPH = (not ASPH))
else:
    generator = GenerateNodeDistribution2d(nRadial, nTheta, rho0, seed,
                                           rmin = rmin,
                                           rmax = rmax,
                                           xmin = xmin,
                                           xmax = xmax,
                                           theta = theta,
                                           azimuthalOffsetFraction = azimuthalOffsetFraction,
                                           nNodePerh = nPerh,
                                           SPH = (not ASPH))

distributeNodes2d((nodes1, generator))
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
                gradientType = RiemannGradientType,
                XSPH = XSPH,
                ASPH = asph,
                densityUpdate=densityUpdate,
                HUpdate = HUpdate)
elif mfm:
    limiter = VanLeerLimiter()
    waveSpeed = DavisWaveSpeed()
    solver = HLLC(limiter,waveSpeed,linearReconstruction)
    hydro = MFM(dataBase = db,
                riemannSolver = solver,
                W = WT,
                cfl=cfl,
                specificThermalEnergyDiffusionCoefficient = 0.00,
                compatibleEnergyEvolution = compatibleEnergy,
                correctVelocityGradient= correctVelocityGradient,
                evolveTotalEnergy = evolveTotalEnergy,
                gradientType = RiemannGradientType,
                XSPH = XSPH,
                ASPH = asph,
                densityUpdate=densityUpdate,
                HUpdate = HUpdate)
elif mfv:
    limiter = VanLeerLimiter()
    waveSpeed = DavisWaveSpeed()
    solver = HLLC(limiter,waveSpeed,linearReconstruction)
    hydro = MFV(dataBase = db,
                riemannSolver = solver,
                W = WT,
                cfl=cfl,
                specificThermalEnergyDiffusionCoefficient = 0.00,
                compatibleEnergyEvolution = compatibleEnergy,
                correctVelocityGradient= correctVelocityGradient,
                evolveTotalEnergy = evolveTotalEnergy,
                gradientType = RiemannGradientType,
                nodeMotionType = nodeMotionType,
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
output("hydro.densityUpdate")
#output("hydro.HEvolution")

if not (gsph or mfm or mfv):
    q = hydro.Q
    q.Cq = 2
    q.Cl = 2
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
yPlane0 = Plane(Vector(0.0, 0.0), Vector(0.0, 1.0))
xbc0 = ReflectingBoundary(xPlane0)
ybc0 = ReflectingBoundary(yPlane0)

for p in packages:
    if thetaFactor in (0.5, ):
        p.appendBoundary(xbc0)
    if thetaFactor in (0.5, 1.0):
        p.appendBoundary(ybc0)

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
integrator.verbose = dtverbose
integrator.allowDtCheck = True
output("integrator")
output("integrator.havePhysicsPackage(hydro)")
output("integrator.dtGrowth")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.allowDtCheck")

#-------------------------------------------------------------------------------
# Build the controller.
#-------------------------------------------------------------------------------
vizMethod = None
if useVoronoiOutput:
    import SpheralVoronoiSiloDump
    vizMethod = SpheralVoronoiSiloDump.dumpPhysicsState
control = SpheralController(integrator, WT,
                            volumeType = volumeType,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle,
                            vizMethod = vizMethod,
                            vizBaseName = "Sedov-cylindrical-2d-%ix%i" % (nRadial, nTheta),
                            vizDerivs=vizDerivs,
                            vizDir = vizDir,
                            vizStep = vizCycle,
                            vizTime = vizTime,
                            SPH = (not ASPH))
output("control")

#-------------------------------------------------------------------------------
# Finally run the problem and plot the results.
#-------------------------------------------------------------------------------
if steps is None:
    control.advance(goalTime, maxSteps)
    if restoreCycle != control.totalSteps:
        control.updateViz(control.totalSteps, integrator.currentTime, 0.0)
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
yprof = mpi.allreduce([x.y for x in nodes1.positions().internalValues()], mpi.SUM)
rho = mpi.allreduce(list(nodes1.massDensity().internalValues()), mpi.SUM)
mass = mpi.allreduce(list(nodes1.mass().internalValues()), mpi.SUM)
v = mpi.allreduce([x.magnitude() for x in nodes1.velocity().internalValues()], mpi.SUM)
eps = mpi.allreduce(list(nodes1.specificThermalEnergy().internalValues()), mpi.SUM)
Pf = ScalarField("pressure", nodes1)
nodes1.pressure(Pf)
P = mpi.allreduce(list(Pf.internalValues()), mpi.SUM)
A = mpi.allreduce([Pi/(rhoi**gamma) for (Pi, rhoi) in zip(Pf.internalValues(), nodes1.massDensity().internalValues())], mpi.SUM)

Hinverse = db.newFluidSymTensorFieldList()
db.fluidHinverse(Hinverse)
hrfl = db.newFluidScalarFieldList()
htfl = db.newFluidScalarFieldList()
for Hfield, hrfield, htfield in zip(Hinverse,
                                    hrfl,
                                    htfl):
    n = Hfield.numElements
    assert hrfield.numElements == n
    assert htfield.numElements == n
    positions = Hfield.nodeList().positions()
    for i in range(n):
        runit = positions[i].unitVector()
        tunit = Vector(-(positions[i].y), positions[i].x).unitVector()
        hrfield[i] = (Hfield[i]*runit).magnitude()
        htfield[i] = (Hfield[i]*tunit).magnitude()
hr = mpi.allreduce(list(hrfl[0].internalValues()), mpi.SUM)
ht = mpi.allreduce(list(htfl[0].internalValues()), mpi.SUM)

Aans = None
if mpi.rank == 0:
    from SpheralTestUtilities import multiSort
    import Pnorm
    multiSort(r, rho, v, eps, P, A, hr, ht)
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
    f.write(("# " + 17*"%16s " + "\n") % ("r", "x", "y", "rho", "m", "P", "v", "eps", "A", "hr", "ht",
                                          "rhoans", "Pans", "vans", "epsans", "Aans", "hrans"))
    for (ri, xi, yi, rhoi, mi, Pi, vi, epsi, Ai, hri, hti, 
         rhoansi, Pansi, vansi, epsansi, Aansi, hansi)  in zip(r, xprof, yprof, rho, mass, P, v, eps, A, hr, ht,
                                                               rhoans, Pans, vans, epsans, Aans, hans):
         f.write((17*"%17.12e " + "\n") % (ri, xi, yi, rhoi, mi, Pi, vi, epsi, Ai, hri, hti, 
                                           rhoansi, Pansi, vansi, epsansi, Aansi, hansi))
    f.close()

if serialDump:
    procs = mpi.procs
    rank = mpi.rank
    serialData = []
    i,j = 0,0
    nodeSet = []
    nodeSet.append(nodes1)
    for i in range(procs):
        for nodeL in nodeSet:
            if rank == i:
                for j in range(nodeL.numInternalNodes):
                    serialData.append([nodeL.positions()[j],3.0/(nodeL.Hfield()[j].Trace()),nodeL.mass()[j],nodeL.massDensity()[j],nodeL.specificThermalEnergy()[j]])
    serialData = mpi.reduce(serialData,mpi.SUM)
    if rank == 0:
        f = open(dataDir + "/serialDump.ascii",'w')
        for i in range(len(serialData)):
            f.write("{0} {1} {2} {3} {4} {5} {6} {7}\n".format(i,serialData[i][0][0],serialData[i][0][1],0.0,serialData[i][1],serialData[i][2],serialData[i][3],serialData[i][4]))
        f.close()

#-------------------------------------------------------------------------------
# Plot the final state.
#-------------------------------------------------------------------------------
if graphics:
    from SpheralMatplotlib import *
    rPlot = plotNodePositions2d(db, colorNodeLists=0, colorDomains=1)
    rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotRadialState(db)
    plotAnswer(answer, control.time(),
               rhoPlot, velPlot, epsPlot, PPlot, HPlot)
    plots = [(rPlot, "Sedov-cylindrical-positions.png"),
             (rhoPlot, "Sedov-cylindrical-rho.png"),
             (velPlot, "Sedov-cylindrical-vel.png"),
             (epsPlot, "Sedov-cylindrical-eps.png"),
             (PPlot, "Sedov-cylindrical-P.png"),
             (HPlot, "Sedov-cylindrical-h.png")]

    # Plot the specific entropy.
    # AsimData = Gnuplot.Data(xprof, A,
    #                         with_ = "points",
    #                         title = "Simulation",
    #                         inline = True)
    # AansData = Gnuplot.Data(xprof, Aans,
    #                         with_ = "lines",
    #                         title = "Solution",
    #                         inline = True)
    
    # Aplot = generateNewGnuPlot()
    # Aplot.plot(AsimData)
    # Aplot.replot(AansData)
    # Aplot.title("Specific entropy")
    # Aplot.refresh()
    # plots.append((Aplot, "Sedov-cylindrical-entropy.png"))

    # if boolCullenViscosity:
    #     cullAlphaPlot = plotFieldList(q.ClMultiplier(),
    #                                   xFunction = "%s.magnitude()",
    #                                   plotStyle = "points",
    #                                   winTitle = "Cullen alpha")
    #     plots += [(cullAlphaPlot, "Sedov-planar-Cullen-alpha.png")]

    # Make hardcopies of the plots.
    for p, filename in plots:
        p.figure.savefig(os.path.join(dataDir, filename))

