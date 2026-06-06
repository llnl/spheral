#ATS:test(SELF, "--CRKSPH True --cfl 0.25 --Cl 1.0 --Cq 1.0 --clearDirectories True --filter 0.0 --goalTime 3.0", label="Gresho CRK, filter=0.0", np=10)
#ATS:test(SELF, "--CRKSPH True --cfl 0.25 --Cl 1.0 --Cq 1.0 --clearDirectories True --filter 0.01 --goalTime 3.0", label="Gresho CRK, filter=0.01", np=10)
#ATS:test(SELF, "--CRKSPH True --cfl 0.25 --Cl 1.0 --Cq 1.0 --clearDirectories True --filter 0.1 --goalTime 3.0", label="Gresho CRK, filter=0.1", np=10)
#ATS:test(SELF, "--CRKSPH True --cfl 0.25 --Cl 1.0 --Cq 1.0 --clearDirectories True --filter 0.2 --goalTime 3.0", label="Gresho CRK, filter=0.2", np=10)
#-------------------------------------------------------------------------------
# The Gresho-Vortex Test
#-------------------------------------------------------------------------------
import shutil
from math import *
from Spheral2d import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from findLastRestart import *
from GenerateNodeDistribution2d import *
from CubicNodeGenerator import GenerateSquareNodeDistribution
from CentroidalVoronoiRelaxation import *

import mpi
import DistributeNodes

title("2-D integrated hydro test --  Gresho-Vortex Test")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(
    rho = 1.0,
    gamma = 5.0/3.0,

    # Translation
    velTx=0.0,
    velTy=0.0,

    # Geometry of Box
    x0 = -0.5,
    x1 =  0.5,
    y0 = -0.5,
    y1 =  0.5,
   
    #Center of Vortex
    xc=0.0,
    yc=0.0,

    # Resolution and node seeding.
    nx1 = 64,
    ny1 = 64,
    seed = "lattice",

    # kernel 
    KernelConstructor = NBSplineKernel,
    order = 5,
    nPerh = 1.51,
    HUpdate = IdealH,
    hmin = 1e-5,
    hmax = 0.5,
    hminratio = 0.1,

    # hydro selection
    svph = False,
    crksph = False,
    psph = False,
    fsisph = False,
    
    # Hydro options
    asph = False, 
    solid = False, 
    XSPH = False,
    epsilonTensile = 0.0,
    nTensile = 8,
    densityUpdate = RigorousSumDensity,
    gradhCorrection = True,
    correctVelocityGradient = True,
    evolveTotalEnergy = False, 
    compatibleEnergy = True,    

    # CRK options
    filter = 0.0, 
    correctionOrder = LinearOrder,
    volumeType = RKSumVolume,

    # SVPH options
    fcentroidal = 0.0,
    fcellPressure = 0.0,

    # artificial viscosity
    Qconstructor = LimitedMonaghanGingoldViscosity,
    Cl = 1.0, 
    Cq = 0.75,
    boolReduceViscosity = False,
    nhQ = 5.0,
    nhL = 10.0,
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
    linearConsistent = False,
    linearInExpansion = False,
    Qlimiter = False,
    balsaraCorrection = False,
    epsilon2 = 1e-2,

    # integrator 
    IntegratorConstructor = CheapSynchronousRK2Integrator,
    cfl = 0.5,
    goalTime = 1.0,
    steps = None,
    vizCycle = 20,
    vizTime = 0.1,
    dt = 0.0001,
    dtMin = 1.0e-5, 
    dtMax = 0.1,
    dtGrowth = 2.0,
    maxSteps = None,
    statsStep = 10,
    domainIndependent = False,
    rigorousBoundaries = False,
    dtverbose = False,
 
    # output
    useVoronoiOutput = False,
    clearDirectories = False,
    restoreCycle = None,
    restartStep = 200,
    dataDir = "dumps-greshovortex-xy",
    graphics = True,
    smooth = None,
    outputFile = None,
    )

assert not(boolReduceViscosity and boolCullenViscosity)
# Decide on our hydro algorithm.

# Decide on our hydro algorithm.
if svph:
    hydroname = "SVPH"
elif crksph:
    Qconstructor = LimitedMonaghanGingoldViscosity
    hydroname = "CRKSPH"
elif psph:
    hydroname = "PSPH"
elif fsisph:
    hydroname = "FSISPH"
else:
    hydroname = "SPH"
if asph:
    hydroname = "A"+hydroname
if solid:
    hydroname = "solid"+hydroname

#-------------------------------------------------------------------------------
# Build our directory paths.
#-------------------------------------------------------------------------------
densityUpdateLabel = {IntegrateDensity : "IntegrateDensity",
                      SumDensity : "SumDensity",
                      RigorousSumDensity : "RigorousSumDensity",
                      SumVoronoiCellDensity : "SumVoronoiCellDensity",
                      VoronoiCellDensity : "VoronoiCellDensity"}
baseDir = os.path.join(dataDir,
                       hydroname,
                       Qconstructor.__name__,
                       KernelConstructor.__name__,
                       "Cl=%g_Cq=%g" % (Cl, Cq),
                       densityUpdateLabel[densityUpdate],
                       "compatibleEnergy=%s" % compatibleEnergy,
                       "Cullen=%s" % boolCullenViscosity,
                       "nPerh=%3.1f" % nPerh,
                       "fcentroidal=%f" % max(fcentroidal, filter),
                       "fcellPressure=%f" % fcellPressure,
                       "%ix%i" % (nx1, ny1))
restartDir = os.path.join(baseDir, "restarts")
restartBaseName = os.path.join(restartDir, "greshovortex-xy-%ix%i" % (nx1, ny1))

vizDir = os.path.join(baseDir, "visit")
if vizTime is None and vizCycle is None:
    vizBaseName = None
else:
    vizBaseName = "greshovortex-xy-%ix%i" % (nx1, ny1)

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
import os, sys
if mpi.rank == 0:
    if clearDirectories and os.path.exists(baseDir):
        shutil.rmtree(baseDir)
    if not os.path.exists(restartDir):
        os.makedirs(restartDir)
    if not os.path.exists(vizDir):
        os.makedirs(vizDir)
mpi.barrier()

#-------------------------------------------------------------------------------
# If we're restarting, find the set of most recent restart files.
#-------------------------------------------------------------------------------
if restoreCycle is None:
    restoreCycle = findLastRestart(restartBaseName)

#-------------------------------------------------------------------------------
# Material properties.
#-------------------------------------------------------------------------------
mu = 1.0
eos = GammaLawGasMKS(gamma, mu)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
if KernelConstructor==NBSplineKernel:
  WT = TableKernel(NBSplineKernel(order), 1000)
else:
  WT = TableKernel(KernelConstructor(), 1000)
output("WT")
kernelExtent = WT.kernelExtent

#-------------------------------------------------------------------------------
# Make the NodeLists.
#-------------------------------------------------------------------------------
if solid:
    nodeListConstructor = makeSolidNodeList
else:
    nodeListConstructor = makeFluidNodeList
nodes = nodeListConstructor("fluid", eos,
                               hmin = hmin,
                               hmax = hmax,
                               hminratio = hminratio,
                               kernelExtent = kernelExtent,
                               nPerh = nPerh)
output("nodes.name")
output("    nodes.hmin")
output("    nodes.hmax")
output("    nodes.hminratio")
output("    nodes.nodesPerSmoothingScale")

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
if restoreCycle is None:
    rmin = 0.0
    rmax = sqrt(2.0)*(x1-x0)
    
    if(seed=="latticeCylindrical"):
        rmin = x1-8.0*nPerh/nx1
        rmax = x1-2.0*nPerh/nx1
    
    generator = GenerateNodeDistribution2d(nx1, ny1, rho,
                                           distributionType = seed,
                                           xmin = (x0, y0),
                                           xmax = (x1, y1),
                                           #rmin = 0.0,
                                           theta = 2.0*pi,
                                           #rmax = sqrt(2.0)*(x1 - x0),
                                           rmax = rmax,
                                           rmin = rmin,
                                           nNodePerh = nPerh,
                                           SPH = (not ASPH))

    if mpi.procs > 1:
        from VoronoiDistributeNodes import distributeNodes2d
    else:
        from DistributeNodes import distributeNodes2d

    distributeNodes2d((nodes, generator))
    print(nodes.name, ":")
    output("    mpi.reduce(nodes.numInternalNodes, mpi.MIN)")
    output("    mpi.reduce(nodes.numInternalNodes, mpi.MAX)")
    output("    mpi.reduce(nodes.numInternalNodes, mpi.SUM)")

    #Set IC
    vel = nodes.velocity()
    eps = nodes.specificThermalEnergy()
    pos = nodes.positions()
    for i in range(nodes.numInternalNodes):
        xi, yi = pos[i]
        r2=(xi-xc)*(xi-xc)+(yi-yc)*(yi-yc)
        ri=sqrt(r2)
        vphi=0.0
        sinPhi=(yi-yc)/ri
        cosPhi=(xi-xc)/ri
        Pi=0.0
        if ri < 0.2:
           Pi=5.0+12.5*r2
 	   vphi=5*ri
        elif ri < 0.4 and ri >= 0.2:
           Pi=9.0+12.5*r2-20.0*ri+4.0*log(5.0*ri)
	   vphi=2.0-5.0*ri
        else:
           Pi=3.0+4*log(2.0)
           #vphi is zero
        velx=velTx-vphi*sinPhi #translation velocity + azimuthal velocity 
        vely=velTy+vphi*cosPhi
        vel[i]=Vector(velx,vely)
        eps0 = Pi/((gamma - 1.0)*rho)
        eps[i]=eps0

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node lists
#-------------------------------------------------------------------------------
db = DataBase()
output("db")
db.appendNodeList(nodes)
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Construct the artificial viscosity.
#-------------------------------------------------------------------------------
q = Qconstructor(Cl, Cq, linearInExpansion)
q.epsilon2 = epsilon2
q.limiter = Qlimiter
q.balsaraShearCorrection = balsaraCorrection
output("q")
output("q.Cl")
output("q.Cq")
output("q.epsilon2")
output("q.limiter")
output("q.balsaraShearCorrection")
output("q.linearInExpansion")
output("q.quadraticInExpansion")

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
packages=[]
if svph:
    hydro = SVPH(dataBase=db,
                W = WT,
                Q = q,
                cfl = cfl,
                compatibleEnergyEvolution = compatibleEnergy,
                densityUpdate = densityUpdate,
                XSVPH = XSPH,
                linearConsistent = linearConsistent,
                generateVoid = False,
                HUpdate = HUpdate,
                fcentroidal = fcentroidal,
                fcellPressure = fcellPressure,
                xmin = Vector(x0 - (x1 - x0), y0 - (y1 - y0)),
                xmax = Vector(x1 + (x1 - x0), y3 + (y1 - y0)))
elif crksph:
    hydro = CRKSPH(dataBase=db,
                    Q = q,
                    filter = filter,
                    epsTensile = epsilonTensile,
                    nTensile = nTensile,
                    cfl = cfl,
                    compatibleEnergyEvolution = compatibleEnergy,
                    XSPH = XSPH,
                    densityUpdate = densityUpdate,
                    HUpdate = HUpdate)
elif fsisph: 
    hydro = FSISPH(dataBase = db,
                Q=q,
                W = WT,
                cfl = cfl,                      
                densityStabilizationCoefficient = 0.1,              
                specificThermalEnergyDiffusionCoefficient = 0.1, 
                linearCorrectGradients = correctVelocityGradient,
                compatibleEnergyEvolution = compatibleEnergy,
                HUpdate = HUpdate,
                ASPH = asph,
                epsTensile = epsilonTensile,
                nTensile = nTensile)
elif psph:
    hydro = PSPH(dataBase=db,
                W=WT,
                Q=q,
                filter=filter,
                cfl=cfl,
                compatibleEnergyEvolution=compatibleEnergy,
                evolveTotalEnergy=evolveTotalEnergy,
                correctVelocityGradient=correctVelocityGradient,
                densityUpdate=densityUpdate,
                HUpdate=HUpdate,
                XSPH=XPSH)
else:
    hydro = SPH(dataBase=db,
                W = WT,
                Q = q,
                cfl = cfl,
                compatibleEnergyEvolution = compatibleEnergy,
                gradhCorrection = gradhCorrection,
                correctVelocityGradient = correctVelocityGradient,
                XSPH = XSPH,
                densityUpdate = densityUpdate,
                HUpdate = HUpdate,
                epsTensile = epsilonTensile,
                nTensile = nTensile)
output("hydro")
output("hydro.cfl")
output("hydro.compatibleEnergyEvolution")
output("hydro.densityUpdate")
output("hydro.HEvolution")

packages += [hydro]

#-------------------------------------------------------------------------------
# Construct the MMRV physics object.
#-------------------------------------------------------------------------------
if boolReduceViscosity:
    evolveReducingViscosityMultiplier = MorrisMonaghanReducingViscosity(nhQ,nhL,aMin,aMax)
    packages.append(evolveReducingViscosityMultiplier)
elif boolCullenViscosity:
    evolveCullenViscosityMultiplier = CullenDehnenViscosity(WT,alphMax,alphMin,betaC,betaD,betaE,fKern,boolHopkinsCorrection)
    packages.append(evolveCullenViscosityMultiplier)

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
xPlane0 = Plane(Vector(x0, y0), Vector( 1.0,  0.0))
xPlane1 = Plane(Vector(x1, y0), Vector(-1.0,  0.0))
yPlane0 = Plane(Vector(x0, y0), Vector( 0.0,  1.0))
yPlane1 = Plane(Vector(x0, y1), Vector( 0.0, -1.0))

xbc = PeriodicBoundary(xPlane0, xPlane1)
ybc = PeriodicBoundary(yPlane0, yPlane1)

xbc0 = ReflectingBoundary(xPlane0)
xbc1 = ReflectingBoundary(xPlane1)
ybc0 = ReflectingBoundary(yPlane0)
ybc1 = ReflectingBoundary(yPlane1)

bcSet = [xbc, ybc]
#bcSet = [xbc0, xbc1, ybc0, ybc1]

for p in packages:
    for bc in bcSet:
        p.appendBoundary(bc)

#-------------------------------------------------------------------------------
# Construct a time integrator, and add the physics packages.
#-------------------------------------------------------------------------------
integrator = IntegratorConstructor(db)
for p in packages:
    integrator.appendPhysicsPackage(p)
integrator.cullGhostNodes = False
integrator.lastDt = dt
integrator.dtMin = dtMin
integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
integrator.domainDecompositionIndependent = domainIndependent
integrator.verbose = dtverbose
integrator.rigorousBoundaries = rigorousBoundaries
output("integrator")
output("integrator.havePhysicsPackage(hydro)")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.dtGrowth")
output("integrator.domainDecompositionIndependent")
output("integrator.rigorousBoundaries")
output("integrator.verbose")

#-------------------------------------------------------------------------------
# If requested, smooth the initial conditions.
#-------------------------------------------------------------------------------
if smooth:
    for iter in range(smooth):
        db.updateConnectivityMap(False)
        cm = db.connectivityMap()
        position_fl = db.fluidPosition
        weight_fl = db.fluidMass
        H_fl = db.fluidHfield
        m0_fl = db.newFluidScalarFieldList(0.0, "m0")
        m1_fl = db.newFluidVectorFieldList(Vector.zero, "m1")
        m2_fl = db.newFluidSymTensorFieldList(SymTensor.zero, "m2")
        A0_fl = db.newFluidScalarFieldList(0.0, "A0")
        A_fl = db.newFluidScalarFieldList(0.0, "A")
        B_fl = db.newFluidVectorFieldList(Vector.zero, "B")
        C_fl = db.newFluidVectorFieldList(Vector.zero, "C")
        D_fl = db.newFluidTensorFieldList(Tensor.zero, "D")
        gradA0_fl = db.newFluidVectorFieldList(Vector.zero, "gradA0")
        gradA_fl = db.newFluidVectorFieldList(Vector.zero, "gradA")
        gradB_fl = db.newFluidTensorFieldList(Tensor.zero, "gradB")
        computeCRKSPHCorrections(cm, WT, weight_fl, position_fl, H_fl, True,
                               m0_fl, m1_fl, m2_fl,
                               A0_fl, A_fl, B_fl, C_fl, D_fl, gradA0_fl, gradA_fl, gradB_fl)
        eps0 = db.fluidSpecificThermalEnergy
        vel0 = db.fluidVelocity
        eps1 = interpolateCRKSPH(eps0, position_fl, weight_fl, H_fl, True, 
                               A_fl, B_fl, cm, WT)
        vel1 = interpolateCRKSPH(vel0, position_fl, weight_fl, H_fl, True, 
                               A_fl, B_fl, cm, WT)
        eps0.assignFields(eps1)
        vel0.assignFields(vel1)

#-------------------------------------------------------------------------------
# Make the problem controller.
#-------------------------------------------------------------------------------
# if useVoronoiOutput:
#     import SpheralVoronoiSiloDump
#     vizMethod = SpheralVoronoiSiloDump.dumpPhysicsState
# else:
#     vizMethod = None # default
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle,
                            #vizMethod = vizMethod,
                            vizBaseName = vizBaseName,
                            vizDir = vizDir,
                            vizStep = vizCycle,
                            vizTime = vizTime,
                            skipInitialPeriodicWork = svph)
output("control")

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
if not steps is None:
    control.step(steps)

else:
    control.advance(goalTime, maxSteps)
    control.updateViz(control.totalSteps, integrator.currentTime, 0.0)
    control.dropRestartFile()

# Plot the final velocity profile
if graphics:
    pos = nodes.positions()
    vel = nodes.velocity()
    vaz = db.newFluidScalarFieldList(0.0, "azimuthal velocity")
    for i in range(nodes.numInternalNodes):
        rhat = (pos[i] - Vector(xc, yc)).unitVector()
        vaz[0][i] = (vel[i] - vel[i].dot(rhat)*rhat).magnitude()
    paz = plotFieldList(vaz, xFunction="(%%s - Vector2d(%g,%g)).magnitude()" % (xc, yc), plotStyle="points", lineTitle="Simulation", winTitle="Azimuthal velocity")
    pmag = plotFieldList(db.fluidVelocity, xFunction="(%%s - Vector2d(%g,%g)).magnitude()" % (xc, yc), yFunction="%s.magnitude()", plotStyle="points", winTitle="Velocity magnitude")

    # Plot the analytic answer.
    xans = [0.0, 0.1, 0.2, 0.4, 1.0]
    yans = [0.0, 0.5, 1.0, 0.0, 0.0]
    ansData = Gnuplot.Data(xans, yans, title="Analytic", with_="lines lt 1 lw 3")
    paz.replot(ansData)
    pmag.replot(ansData)
    plots = [(paz, "GreshoVortex-velazimuthal.png"),
             (pmag, "GreshoVortex-velmag.png")]

    # Make hardcopies of the plots.
    for p, filename in plots:
        p.hardcopy(os.path.join(baseDir, filename), terminal="png")

#-------------------------------------------------------------------------------
# If requested, write out the state in a global ordering to a file.
#-------------------------------------------------------------------------------
if outputFile:
    outputFile = os.path.join(baseDir, outputFile)
    from SpheralTestUtilities import multiSort
    P = ScalarField("pressure", nodes)
    nodes.pressure(P)
    xprof = mpi.reduce([x.x for x in nodes.positions().internalValues()], mpi.SUM)
    yprof = mpi.reduce([x.y for x in nodes.positions().internalValues()], mpi.SUM)
    rhoprof = mpi.reduce(nodes.massDensity().internalValues(), mpi.SUM)
    Pprof = mpi.reduce(P.internalValues(), mpi.SUM)
    vprof = mpi.reduce([v.magnitude() for v in nodes.velocity().internalValues()], mpi.SUM)
    velx = mpi.reduce([v.x for v in nodes.velocity().internalValues()], mpi.SUM)
    vely = mpi.reduce([v.y for v in nodes.velocity().internalValues()], mpi.SUM)
    epsprof = mpi.reduce(nodes.specificThermalEnergy().internalValues(), mpi.SUM)
    hprof = mpi.reduce([1.0/sqrt(H.Determinant()) for H in nodes.Hfield().internalValues()], mpi.SUM)
    mof = mortonOrderIndices(db)
    mo = mpi.reduce(mof[0].internalValues(), mpi.SUM)

    if mpi.rank == 0:
        rprof = [sqrt(xi*xi + yi*yi) for xi, yi in zip(xprof, yprof)]
        multiSort(rprof, mo, xprof, yprof, rhoprof, Pprof, vprof, epsprof, hprof,velx,vely)
        L1 = 0.0
        vazprof = []
        for i in range(len(xprof)):
           rhat = (Vector(xprof[i],yprof[i]) - Vector(xc, yc)).unitVector()
           vel_vec = Vector(velx[i],vely[i])
           vazprof.append((vel_vec - vel_vec.dot(rhat)*rhat).magnitude())
           vans = 0.0
           if rprof[i] < 0.2:
             vans = 5*rprof[i]
           elif rprof[i] < 0.4:
             vans = 2.0-5.0*rprof[i]
           else:
             vans = 0.0
           L1 = L1 + abs(vazprof[i]-vans)
        L1 = L1/len(xprof)
        with open("Converge.txt.%s" % nPerh, "a") as myfile:
          myfile.write("%s\t %s\n" % (nx1,L1))
        f = open(outputFile, "w")
        f.write(("# " + 16*"%15s " + "\n") % ("r", "x", "y", "rho", "P", "v", "eps", "h", "mortonOrder",
                                              "x_uu", "y_uu", "rho_uu", "P_uu", "v_uu", "eps_uu", "h_uu"))
        for (ri, xi, yi, rhoi, Pi, vi, vaz, epsi, hi, mi)  in zip(rprof, xprof, yprof, rhoprof, Pprof, vprof, vazprof, epsprof, hprof, mo):
            f.write((9*"%16.12e " + "%i " + 7*"%i " + "\n") % (ri, xi, yi, rhoi, Pi, vi, vaz, epsi, hi, mi,
                                                               unpackElementUL(packElementDouble(xi)),
                                                               unpackElementUL(packElementDouble(yi)),
                                                               unpackElementUL(packElementDouble(rhoi)),
                                                               unpackElementUL(packElementDouble(Pi)),
                                                               unpackElementUL(packElementDouble(vi)),
                                                               unpackElementUL(packElementDouble(epsi)),
                                                               unpackElementUL(packElementDouble(hi))))
        f.close()
