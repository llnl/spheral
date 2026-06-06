import os, sys
import shutil
from SolidSpheral1d import *
from SpheralTestUtilities import *
from RiemannSolution import *
from DistributeNodes import distributeNodesInRange1d

title("1-D hydro test -- a variety of Riemann shocktubes")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(

    # Choose the node seeding
    nx = 500,
    nodeMatching = "equal_mass",     # ("equal_mass", "equal_volume")

    # Initial conditions
    problem = "Sod",                 # ("", "Sod", "123", "leftwc", "2shock_collision", "Stationary_contact", "Slow_shock", "shock_contact_shock", "LeBlanc")
    rho1 = 1.0,
    rho2 = 0.25,
    P1 = 1.0,
    P2 = 0.1795,
    v1 = 0.0,
    v2 = 0.0,
    x0 = -0.5,
    x1 = 0.0,
    x2 = 0.5,
    hsmooth = 0.0,                   # if >0, smooths intial discontinuity

    # Gamma-law gas 
    gammaGas = 5.0/3.0,
    mu = 1.0,

    # 1 or 2 NodeLists
    numNodeLists = 1,

    # Choose the hydro variant
    svph = False,
    crksph = False,
    psph = False,
    solid = False,                   # If true, use the fluid limit of the solid hydro option

    # All sorts of hydro parameters
    compatibleEnergy = True,
    evolveTotalEnergy = False,       # Only for SPH variants -- evolve total rather than specific energy
    HopkinsConductivity = False,     # PSPH
    correctionOrder = LinearOrder,   # CRK
    volumeType = RKSumVolume,       # CRK
    linearConsistent = False,        # SVPH
    correctVelocityGradient = True,  # SPH
    epsilonTensile = 0.0,            # SPH
    nTensile = 8,                    # SPH
    HUpdate = IdealH,
    densityUpdate = RigorousSumDensity,
    gradhCorrection = True,
    cfl = 0.25,
    XSPH = False,
    rhoMin = 0.01,
    filter = 0.00,
    order = 5,                       # Order of the B-spline kernel

    # Artificial viscosity
    Cl = None,
    Cq = None,

    # Smoothing scale
    nPerh = 1.35,
    hmin = 1e-10,
    hmax = 10.0,

    # Time integration
    IntegratorConstructor = VerletIntegrator,
    dtverbose = False,
    steps = None,
    goalTime = 0.15,
    dt = 1e-6,
    dtMin = 1.0e-7,
    dtMax = 0.1,
    dtGrowth = 2.0,
    maxSteps = None,

    # Output
    statsStep = 10,
    clearDirectories = False,
    restoreCycle = -1,
    restartStep = 10000,
    dataDirBase = "dumps-",
    outputFile = None,
    checkRestart = False,

    graphics = True,
)

assert numNodeLists in (1, 2)

if problem:
    assert problem.lower() in Riemann_packaged_problems
    x0, x2, x1, gammaGas, goalTime, rho1, v1, P1, rho2, v2, P2 = Riemann_packaged_problems[problem.lower()]
else:
    problem = "user"

# Node resolution
nx = int((x2 - x0)*nx + 0.49)
assert nodeMatching in ("equal_mass", "equal_volume")
if nodeMatching == "equal_mass":
    w1 = (x1 - x0)*rho1
    w2 = (x2 - x1)*rho2
else:
    w1 = x1 - x0
    w2 = x2 - x1
nx1 = int(w1/(w1 + w2)*nx + 0.5)
nx2 = nx - nx1
assert min(nx1, nx2) > 0
assert nx1 + nx2 == nx

print("Problem parameters:   (x0, x1, x2) = (%g, %g, %g)" % (x0, x1, x2))
print("                    (rho1, v1, P1) = (%g, %g, %g)" % (rho1, v1, P1))
print("                    (rho2, v2, P2) = (%g, %g, %g)" % (rho2, v2, P2))
print("                        (nx1, nx2) = (%i, %i)" % (nx1, nx2))

# Directories and naming.
if svph:
    hydroname = "SVPH"
elif crksph:
    hydroname = "CRKSPH"
elif psph:
    hydroname = "PSPH"
else:
    hydroname = "SPH"
if solid:
    hydroname = "Solid" + hydroname

dataDir = os.path.join(dataDirBase + problem, 
                       "numNodeLists=%i" % numNodeLists,
                       hydroname,
                       "nPerh=%f" % nPerh,
                       "densityUpdate=%s" % densityUpdate,
                       "compatibleEnergy=%s" % compatibleEnergy,
                       "evolveTotalEnergy=%s" % evolveTotalEnergy,
                       "correctionOrder=%s" % correctionOrder,
                       "filter=%f" % filter,
                       "nx=%i-%i" % (nx1, nx2))
restartDir = os.path.join(dataDir, "restarts")
restartBaseName = os.path.join(restartDir, problem)

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
eos = GammaLawGasMKS(gammaGas, mu)
strength = NullStrength()

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
WT = TableKernel(NBSplineKernel(order), 1000)
kernelExtent = WT.kernelExtent
output("WT")

#-------------------------------------------------------------------------------
# Make the NodeLists.
#-------------------------------------------------------------------------------
if solid:
    makeNL = makeSolidNodeList
else:
    makeNL = makeFluidNodeList
nodes1 = makeNL("nodes1", eos, 
                hmin = hmin,
                hmax = hmax,
                nPerh = nPerh,
                kernelExtent = kernelExtent,
                rhoMin = rhoMin)
nodes2 = makeNL("nodes2", eos, 
                hmin = hmin,
                hmax = hmax,
                nPerh = nPerh,
                kernelExtent = kernelExtent,
                rhoMin = rhoMin)
nodeSet = [nodes1, nodes2]

#-------------------------------------------------------------------------------
# Functions to specify the initial density and specific thermal energy with
# optional smoothing.
#-------------------------------------------------------------------------------
dx1 = (x1 - x0)/nx1
dx2 = (x2 - x1)/nx2
hfold = hsmooth*max(dx1, dx2)
def rho_initial(xi):
    if hfold > 0.0:
        return rho1 + (rho2 - rho1)/(1.0 + exp(-(xi - x1)/hfold))
    elif xi <= x1:
        return rho1
    else:
        return rho2

def vel_initial(xi):
    if hfold > 0.0:
        return v1 + (v2 - v1)/(1.0 + exp(-(xi - x1)/hfold))
    elif xi <= x1:
        return v1
    else:
        return v2

def specificEnergy(xi, rhoi):
    if hfold > 0.0:
        Pi = P1 + (P2 - P1)/(1.0 + exp(-(xi - x1)/hfold))
    elif xi <= x1:
        Pi = P1
    else:
        Pi = P2
    return Pi/((gammaGas - 1.0)*rhoi)

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
if numNodeLists == 1:
    distributeNodesInRange1d([(nodes1, nx1, rho1, (x0, x1)),
                              (nodes1, nx2, rho2, (x1, x2))])
else:
    distributeNodesInRange1d([(nodes1, nx1, rho1, (x0, x1)),
                              (nodes2, nx2, rho2, (x1, x2))])
output("nodes1.numNodes")
output("nodes2.numNodes")

# Set node specific thermal energies
for nodes in nodeSet:
    pos = nodes.positions()
    eps = nodes.specificThermalEnergy()
    rho = nodes.massDensity()
    vel = nodes.velocity()
    for i in range(nodes.numInternalNodes):
        eps[i] = specificEnergy(pos[i].x, rho[i])
        vel[i].x = vel_initial(pos[i].x)

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase()
output("db")
for nodes in nodeSet:
    output("db.appendNodeList(nodes)")
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
                 XSVPH = XSPH,
                 linearConsistent = linearConsistent,
                 generateVoid = False,
                 densityUpdate = densityUpdate,
                 HUpdate = HUpdate,
                 xmin = Vector(-100.0),
                 xmax = Vector( 100.0))
elif crksph:
    hydro = CRKSPH(dataBase = db,
                   W = WT, 
                   filter = filter,
                   cfl = cfl,
                   correctionOrder = correctionOrder,
                   volumeType = volumeType,
                   compatibleEnergyEvolution = compatibleEnergy,
                   evolveTotalEnergy = evolveTotalEnergy,
                   XSPH = XSPH,
                   densityUpdate = densityUpdate,
                   HUpdate = HUpdate)
elif psph:
    hydro = PSPH(dataBase = db,
                 W = WT,
                 cfl = cfl,
                 compatibleEnergyEvolution = compatibleEnergy,
                 evolveTotalEnergy = evolveTotalEnergy,
                 densityUpdate = densityUpdate,
                 HUpdate = HUpdate,
                 XSPH = XSPH,
                 correctVelocityGradient = correctVelocityGradient,
                 HopkinsConductivity = HopkinsConductivity)
else:
    hydro = SPH(dataBase = db,
                W = WT,
                cfl = cfl,
                compatibleEnergyEvolution = compatibleEnergy,
                evolveTotalEnergy = evolveTotalEnergy,
                gradhCorrection = gradhCorrection,
                densityUpdate = densityUpdate,
                HUpdate = HUpdate,
                XSPH = XSPH,
                correctVelocityGradient = correctVelocityGradient,
                epsTensile = epsilonTensile,
                nTensile = nTensile)
output("hydro")

packages = [hydro]

#-------------------------------------------------------------------------------
# Tweak the artificial viscosity.
#-------------------------------------------------------------------------------
q = hydro.Q
if Cl:
    q.Clinear = Cl
if Cq:
    q.Cquadratic = Cq
output("q")
output("q.Cl")
output("q.Cq")
output("q.limiter")
output("q.epsilon2")
output("q.linearInExpansion")
output("q.quadraticInExpansion")

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
if v1 == 0.0 and v2 == 0.0:
    xPlane0 = Plane(Vector(x0), Vector( 1.0))
    xPlane1 = Plane(Vector(x2), Vector(-1.0))
    xbc0 = ReflectingBoundary(xPlane0)
    xbc1 = ReflectingBoundary(xPlane1)

    for p in packages:
        p.appendBoundary(xbc0)
        p.appendBoundary(xbc1)

#-------------------------------------------------------------------------------
# Construct an integrator, and add the one physics package.
#-------------------------------------------------------------------------------
integrator = IntegratorConstructor(db)
for p in packages:
    integrator.appendPhysicsPackage(p)
integrator.lastDt = dt
integrator.dtGrowth = dtGrowth
if dtMin:
    integrator.dtMin = dtMin
if dtMax:
    integrator.dtMax = dtMax
integrator.verbose = dtverbose
output("integrator")
output("integrator.havePhysicsPackage(hydro)")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")

#-------------------------------------------------------------------------------
# Make the problem controller.
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
if not steps is None:
    control.step(steps)

    # Are we doing the restart test?
    if checkRestart:
        state0 = State(db, integrator.physicsPackages())
        state0.copyState()
        control.loadRestartFile(control.totalSteps)
        state1 = State(db, integrator.physicsPackages())
        if not state1 == state0:
            raise ValueError("The restarted state does not match!")
        else:
            print("Restart check PASSED.")

else:
    control.advance(goalTime, maxSteps)
    control.dropRestartFile()

#-------------------------------------------------------------------------------
# Plot the final state.
#-------------------------------------------------------------------------------
dx1 = (x1 - x0)/nx1
dx2 = (x2 - x1)/nx2
h1 = 1.0/(nPerh*dx1)
h2 = 1.0/(nPerh*dx2)
answer = RiemannSolution(problem = problem,
                         x0 = x0,
                         x1 = x2,
                         xdiaph = x1,
                         gamma_gas = gammaGas,
                         dl = rho1,
                         vl = v1,
                         pl = P1,
                         hl = 1.0/h1,
                         dr = rho2,
                         vr = v2,
                         pr = P2,
                         hr = 1.0/h2)

#cs = db.newFluidScalarFieldList(0.0, "sound speed")
#db.fluidSoundSpeed(cs)
cs = hydro.soundSpeed()

# Make a flat list from a FieldList
def createList(x):
    result = []
    for i in range(len(x)):
        for j in range(x[i].numInternalElements):
            result.append(x(i,j))
    return mpi.allreduce(result, mpi.SUM)

# Compute the simulated specific entropy.
rho = createList(db.fluidMassDensity)
P = createList(hydro.pressure())
A = [Pi/rhoi**gammaGas for (Pi, rhoi) in zip(P, rho)]

# The analytic solution for the simulated entropy.
xprof = [x.x for x in createList(db.fluidPosition)]
xans, vans, uans, rhoans, Pans, Aans, hans = answer.solution(control.time(), xprof)
csAns = [sqrt(gammaGas*Pi/rhoi) for (Pi, rhoi) in zip(Pans,  rhoans)]

if graphics:
    from SpheralMatplotlib import *
    
    rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotState(db, plotStyle="r-o")
    APlot = newFigure()
    APlot.plot(xprof, A, "r-o")
    plt.title("P/rho^\gamma")
    plotAnswer(answer, control.time(),
               rhoPlot, velPlot, epsPlot, PPlot, APlot, HPlot)
    pE = plotEHistory(control.conserve)

    csPlot = plotFieldList(cs, plotStyle="r-o", lineTitle="Simulation", winTitle="Sound speed")
    csPlot.plot(xans, csAns, "k-", label="Analytic")
    csPlot.axes.legend()

    suffix = "-%s-%s-compatibleEnergy=%s-densityUpdate=%s.pdf" % (problem, hydroname, compatibleEnergy, densityUpdate)
    plots = [(rhoPlot, "Sod-planar-rho" + suffix),
             (velPlot, "Sod-planar-vel" + suffix),
             (epsPlot, "Sod-planar-eps" + suffix),
             (PPlot, "Sod-planar-P" + suffix),
             (HPlot, "Sod-planar-h" + suffix),
             (csPlot, "Sod-planar-cs" + suffix),
             (APlot, "Sod-planar-entropy" + suffix)]
    
    if crksph:
        volPlot = plotFieldList(hydro.volume(), 
                                winTitle = "volume",
                                plotStyle = "r-o",
                                colorNodeLists = False, plotGhosts = False)
        aplot = plotFieldList(hydro.A(),
                              winTitle = "A",
                              plotStyle = "r-o",
                              colorNodeLists = False)
        bplot = plotFieldList(hydro.B(),
                              yFunction = "%s.x",
                              winTitle = "B",
                              plotStyle = "r-o",
                              colorNodeLists = False)
        splot = plotFieldList(hydro.surfacePoint(),
                              winTitle = "surface point",
                              plotStyle = "r-o",
                              colorNodeLists = False)
        voidplot = plotFieldList(hydro.voidPoint(),
                                 winTitle = "void point",
                                 plotGhosts = True,
                                 plotStyle = "r-o",
                                 colorNodeLists = False)
        plots += [(volPlot, "Sod-planar-vol" + suffix),
                   (aplot, "Sod-planar-ACRK" + suffix),
                   (bplot, "Sod-planar-BCRK" + suffix),
                   (splot, "Sod-planar-surfacePoint" + suffix),
                   (voidplot, "Sod-planar-voidPoint" + suffix)]
    
    viscPlot = plotFieldList(hydro.maxViscousPressure(),
                             winTitle = "$\max( \\rho^2 \Pi_{ij})$",
                             plotStyle = "r-o",
                             colorNodeLists = False)
    plots.append((viscPlot, "Sod-planar-viscosity" + suffix))
    
    # Make hardcopies of the plots.
    for p, filename in plots:
        p.set_xlim(0, 1)
        p.figure.savefig(os.path.join(dataDir, filename))
    pE.figure.savefig(os.path.join(dataDir, "Sod-planar-E" + suffix))

print("Energy conservation: original=%g, final=%g, error=%g" % (control.conserve.EHistory[0],
                                                                control.conserve.EHistory[-1],
                                                                (control.conserve.EHistory[-1] - control.conserve.EHistory[0])/control.conserve.EHistory[0]))

#-------------------------------------------------------------------------------
# If requested, write out the state in a global ordering to a file.
#-------------------------------------------------------------------------------
from SpheralTestUtilities import multiSort
mof = mortonOrderIndices(db)
mo = createList(mof)
rhoprof = createList(db.fluidMassDensity)
Pprof = createList(hydro.pressure())
vprof = [v.x for v in createList(db.fluidVelocity)]
epsprof = createList(db.fluidSpecificThermalEnergy)
hprof = [1.0/Hi.xx for Hi in createList(db.fluidHfield)]

rmin = x0
rmax = x2
if mpi.rank == 0:
    multiSort(mo, xprof, rhoprof, Pprof, vprof, epsprof, hprof)
    if outputFile:
        outputFile = os.path.join(dataDir, outputFile)
        f = open(outputFile, "w")
        f.write(("#  " + 19*"'%s' " + "\n") % ("x", "rho", "P", "v", "eps", "A", "h", "mo",
                                               "rhoans", "Pans", "vans", "Aans", "hans",
                                               "x_UU", "rho_UU", "P_UU", "v_UU", "eps_UU", "h_UU"))
        for (xi, rhoi, Pi, vi, epsi, Ai, hi, mi,
             rhoansi, Pansi, vansi, Aansi, hansi) in zip(xprof, rhoprof, Pprof, vprof, epsprof, A, hprof, mo,
                                                         rhoans, Pans, vans, Aans, hans):
            f.write((7*"%16.12e " + "%i " + 5*"%16.12e " + 6*"%i " + '\n') % 
                    (xi, rhoi, Pi, vi, epsi, Ai, hi, mi,
                     rhoansi, Pansi, vansi, Aansi, hansi,
                     unpackElementUL(packElementDouble(xi)),
                     unpackElementUL(packElementDouble(rhoi)),
                     unpackElementUL(packElementDouble(Pi)),
                     unpackElementUL(packElementDouble(vi)),
                     unpackElementUL(packElementDouble(epsi)),
                     unpackElementUL(packElementDouble(hi))))
        f.close()

    import Pnorm
    print("\tQuantity \t\tL1 \t\t\tL2 \t\t\tLinf")
    failure = False
    hD = []
    for (name, data, ans) in [("Mass Density", rhoprof, rhoans),
                              ("Pressure", Pprof, Pans),
                              ("Velocity", vprof, vans),
                              ("Thermal E", epsprof, uans),
                              ("Entropy", A, Aans),
                              ("h       ", hprof, hans)]:
        assert len(data) == len(ans)
        error = [data[i] - ans[i] for i in range(len(data))]
        Pn = Pnorm.Pnorm(error, xprof)
        L1 = Pn.gridpnorm(1, rmin, rmax)
        L2 = Pn.gridpnorm(2, rmin, rmax)
        Linf = Pn.gridpnorm("inf", rmin, rmax)
        print("\t%s \t\t%g \t\t%g \t\t%g" % (name, L1, L2, Linf))
        #f.write(("\t\t%g") % (L1))
        hD.append([L1,L2,Linf])
    #f.write("\n")

    print("%d\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t" % (nx1+nx2,hD[0][0],hD[1][0],hD[2][0],hD[3][0],
                                                                                hD[0][1],hD[1][1],hD[2][1],hD[3][1],
                                                                                hD[0][2],hD[1][2],hD[2][2],hD[3][2]))
