#
# SPH
#
#ATS:sph1 = test(        SELF, "--crksph False --cfl 0.25 --graphics None --clearDirectories True  --restartStep 20 --steps 40", label="Spherical Sod problem with SPH -- 1-D (serial)")
#ATS:sph2 = testif(sph1, SELF, "--crksph False --cfl 0.25 --graphics None --clearDirectories False --restartStep 20 --steps 20 --restoreCycle 20 --checkRestart True", label="Spherical Sod problem with SPH -- 1-D (serial) RESTART CHECK")
#
# CRK
#
#ATS:crk1 = test(        SELF, "--crksph True --cfl 0.25 --graphics None --clearDirectories True  --restartStep 20 --steps 40", label="Spherical Sod problem with CRK -- 1-D (serial)")
#ATS:crk2 = testif(crk1, SELF, "--crksph True --cfl 0.25 --graphics None --clearDirectories False --restartStep 20 --steps 20 --restoreCycle 20 --checkRestart True", label="Spherical Sod problem with CRK -- 1-D (serial) RESTART CHECK")
#
# Solid FSISPH
#
#ATS:fsisph1 = test(           SELF, "--crksph False --fsisph True --solid True --cfl 0.25 --graphics None --clearDirectories True  --restartStep 20 --steps 40", label="Spherical Sod problem with FSISPH -- 1-D (serial)", fsisph=True)
#ATS:fsisph2 = testif(fsisph1, SELF, "--crksph False --fsisph True --solid True --cfl 0.25 --graphics None --clearDirectories False --restartStep 20 --steps 20 --restoreCycle 20 --checkRestart True", label="Spherical Sod problem with FSISPH -- 1-D (serial) RESTART CHECK", fsisph=True)
#
# GSPH
#
#ATS:gsph1 = test(         SELF, "--gsph True --cfl 0.25 --graphics None --clearDirectories True  --restartStep 20 --steps 40", label="Spherical Sod problem with GSPH -- 1-D (serial)", gsph=True)
#ATS:gsph2 = testif(gsph1, SELF, "--gsph True --cfl 0.25 --graphics None --clearDirectories False --restartStep 20 --steps 20 --restoreCycle 20 --checkRestart True", label="Spherical Sod problem with GSPH -- 1-D (serial) RESTART CHECK", gsph=True)
#
import os, sys
import shutil
from SphericalSpheral import *
from SpheralTestUtilities import *
from SodAnalyticSolution import SodSolutionGasGas as SodSolution

title("1-D integrated hydro test -- spherical Sod problem")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(nr1 = 800,
            nr2 = 100,
            rho1 = 1.0,
            rho2 = 0.125,
            P1 = 1.0,
            P2 = 0.1,
            gamma1 = 1.4,
            gamma2 = 1.4,

            numNodeLists = 1,

            x0 = 0.0,
            x1 = 1.0,
            x2 = 2.0,

            hsmooth = 0.0,             # Optionally smooth initial discontinuity, expressed as particle spacings
            sumInitialDensity = False, # Optionally sum the initial density before setting the pressure and such

            nPerh = 4.00,

            mu = 1.0,
            
            crksph = False,
            psph = False,
            fsisph = False,
            gsph = False,

            evolveTotalEnergy = False,  # Only for SPH variants -- evolve total rather than specific energy
            solid = False,    # If true, use the fluid limit of the solid hydro option
            boolReduceViscosity = False,
            nh = 5.0,
            aMin = 0.1,
            aMax = 2.0,
            Cl = None,
            Cq = None,
            Qlimiter = None,
            epsilon2 = None,
            etaCritFrac = None,
            etaFoldFrac = None,
            linearInExpansion = None,
            quadraticInExpansion = None,
            boolCullenViscosity = False,
            alphMax = 2.0,
            alphMin = 0.02,
            betaC = 0.7,
            betaD = 0.05,
            betaE = 1.0,
            fKern = 1.0/3.0,
            hmin = 1e-10,
            hmax = 1.0e10,
            cfl = 0.25,
            XSPH = False,
            epsilonTensile = 0.0,
            nTensile = 8,
            rhoMin = 0.01,
            hourglass = None,
            hourglassOrder = 1,
            
            IntegratorConstructor = SynchronousRK2Integrator,
            dtverbose = False,
            steps = None,
            goalTime = 0.5,
            dt = 1e-6,
            dtMin = 1.0e-6,
            dtMax = 0.1,
            dtGrowth = 2.0,
            rigorousBoundaries = False,
            maxSteps = None,
            statsStep = 10,
            HUpdate = IdealH,
            correctionOrder = LinearOrder,
            volumeType = RKSumVolume,
            limitMultimaterialTopology = True,
            densityUpdate = RigorousSumDensity,
            compatibleEnergy = True,
            correctVelocityGradient = True,
            gradhCorrection = False,
            linearConsistent = False,

            useRefinement = False,

            clearDirectories = True,
            restoreCycle = -1,
            restartStep = 10000,
            dataDirBase = "dumps-Sod-spherical",
            restartBaseName = "Sod-spherical-1d-restart",
            outputFile = None,
            checkRestart = False,

            graphics = True,
            )

assert not(boolReduceViscosity and boolCullenViscosity)
assert not (GSPH and (boolReduceViscosity or boolCullenViscosity))
assert not (fsisph and not solid)
assert numNodeLists in (1, 2)

if crksph:
    hydroname = os.path.join("CRKSPH",
                             str(volumeType),
                             str(correctionOrder))
elif psph:
    hydroname = "PSPH"
elif fsisph:
    hydroname = "FSISPH"
elif gsph:
    hydroname = "GSPH"
else:
    hydroname = "SPH"
if solid:
    hydroname = "Solid" + hydroname

dataDir = os.path.join(dataDirBase, 
                       "numNodeLists=%i" % numNodeLists,
                       hydroname,
                       "nPerh=%f" % nPerh,
                       "compatibleEnergy=%s" % compatibleEnergy,
                       "%i" % (nr1 + nr2))
restartDir = os.path.join(dataDir, "restarts")
restartBaseName = os.path.join(restartDir, "Sod-spherical-1d-%i" % (nr1 + nr2))

# if two gammas were given but one node list we need to sort things out
if numNodeLists == 1 and abs(gamma1-gamma2) > 0.001:
    gamma2 = gamma1
    print(" ")
    print("============================================================")
    print("warning: because we're running 1 nodelist, gamma2 had to be ")
    print("         set equal to gamma1")
    print("============================================================")
    print(" ")
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
eos1 = GammaLawGasMKS(gamma1, mu)
eos2 = GammaLawGasMKS(gamma2, mu)
strength = NullStrength()

#-------------------------------------------------------------------------------
# Interpolation ernels.
#-------------------------------------------------------------------------------
W = WendlandC4Kernel3d()
kernelExtent = W.kernelExtent
output("W")

#-------------------------------------------------------------------------------
# Make the NodeLists.
#-------------------------------------------------------------------------------
if solid:
    makeNL = makeSolidNodeList
else:
    makeNL = makeFluidNodeList
nodes1 = makeNL("nodes1", eos1, 
                hmin = hmin,
                hmax = hmax,
                nPerh = nPerh,
                kernelExtent = kernelExtent,
                rhoMin = rhoMin)
nodes2 = makeNL("nodes2", eos2, 
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
dx1 = (x1 - x0)/nr1
dx2 = (x2 - x1)/nr2
hfold = hsmooth*max(dx1, dx2)
def rho_initial(xi):
    if hfold > 0.0:
        return rho1 + (rho2 - rho1)/(1.0 + exp(-(xi - x1)/hfold))
    elif xi <= x1:
        return rho1
    else:
        return rho2

def specificEnergy(xi, rhoi, gammai):
    if hfold > 0.0:
        Pi = P1 + (P2 - P1)/(1.0 + exp((-xi - x1)/hfold))
    elif xi <= x1:
        Pi = P1
    else:
        Pi = P2
    return Pi/((gammai - 1.0)*rhoi)

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
from GenerateSphericalNodeProfile1d import GenerateSphericalNodeProfile1d
from GenerateSphericalNodeDistribution1d import GenerateSphericalNodeDistribution1d
from CompositeNodeDistribution import CompositeNodeDistribution
from PeanoHilbertDistributeNodes import distributeNodes1d
if numNodeLists == 1:
    # gen = GenerateSphericalNodeProfile1d(nr = nr1 + nr2,
    #                                      rho = rho_initial,
    #                                      rmin = x0,
    #                                      rmax = x2,
    #                                      nNodePerh = nPerh)
    gen1 = GenerateSphericalNodeDistribution1d(nr = nr1,
                                               rho = rho1,
                                               rmin = x0,
                                               rmax = x1,
                                               nNodePerh = nPerh)
    gen2 = GenerateSphericalNodeDistribution1d(nr = nr2,
                                               rho = rho2,
                                               rmin = x1,
                                               rmax = x2,
                                               nNodePerh = nPerh)
    gen = CompositeNodeDistribution(gen1, gen2)
    distributeNodes1d((nodes1, gen))
else:
    if hsmooth > 0:
        gen1 = GenerateSphericalNodeProfile1d(nr = nr1,
                                              rho = rho_initial,
                                              rmin = x0,
                                              rmax = x1,
                                              nNodePerh = nPerh)
        gen2 = GenerateSphericalNodeProfile1d(nr = nr2,
                                              rho = rho_initial,
                                              rmin = x1,
                                              rmax = x2,
                                              nNodePerh = nPerh)
    else:
        gen1 = GenerateSphericalNodeProfile1d(nr = nr1,
                                              rho = rho1,
                                              rmin = x0,
                                              rmax = x1,
                                              nNodePerh = nPerh)
        gen2 = GenerateSphericalNodeProfile1d(nr = nr2,
                                              rho = rho2,
                                              rmin = x1,
                                              rmax = x2,
                                              nNodePerh = nPerh)
    distributeNodes1d((nodes1, gen1),
                      (nodes2, gen2))
output("nodes1.numNodes")
output("nodes2.numNodes")

# Set node specific thermal energies
for nodes in nodeSet:

    pos = nodes.positions()
    eps = nodes.specificThermalEnergy()
    rho = nodes.massDensity()
    gammai = nodes.eos.gamma

    for i in range(nodes.numInternalNodes):
        eps[i] = specificEnergy(pos[i].x, rho[i], gammai)

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
if crksph:
    hydro = CRKSPH(dataBase = db,
                   order = correctionOrder,
                   cfl = cfl,
                   compatibleEnergyEvolution = compatibleEnergy,
                   evolveTotalEnergy = evolveTotalEnergy,
                   XSPH = XSPH,
                   densityUpdate = densityUpdate,
                   HUpdate = HUpdate)
elif psph:
    hydro = PSPH(dataBase = db,
                 W = W,
                 cfl = cfl,
                 compatibleEnergyEvolution = compatibleEnergy,
                 evolveTotalEnergy = evolveTotalEnergy,
                 densityUpdate = densityUpdate,
                 HUpdate = HUpdate,
                 XSPH = XSPH,
                 correctVelocityGradient = correctVelocityGradient)
elif fsisph:
    sumDensityNodeLists = [nodes1]
    if numNodeLists == 2:
        sumDensityNodeLists += [nodes2]
    hydro = FSISPH(dataBase = db,
                   W = W,
                   cfl = cfl,
                   sumDensityNodeLists=sumDensityNodeLists,                       
                   densityStabilizationCoefficient = 0.00,
                   specificThermalEnergyDiffusionCoefficient = 0.00,
                   interfaceMethod = HLLCInterface,
                   compatibleEnergyEvolution = compatibleEnergy,
                   evolveTotalEnergy = evolveTotalEnergy,
                   correctVelocityGradient = correctVelocityGradient,
                   HUpdate = HUpdate)
elif gsph:
    limiter = VanLeerLimiter()
    waveSpeed = DavisWaveSpeed()
    solver = HLLC(limiter,waveSpeed,True,RiemannGradient)
    hydro = GSPH(dataBase = db,
                riemannSolver = solver,
                W = W,
                cfl=cfl,
                compatibleEnergyEvolution = compatibleEnergy,
                correctVelocityGradient=correctVelocityGradient,
                evolveTotalEnergy = evolveTotalEnergy,
                XSPH = XSPH,
                densityUpdate=densityUpdate,
                HUpdate = IdealH,
                epsTensile = epsilonTensile,
                nTensile = nTensile)
else:
    hydro = SPH(dataBase = db,
                W = W,
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
output("hydro.compatibleEnergyEvolution")
packages = [hydro]

#-------------------------------------------------------------------------------
# Tweak the artificial viscosity.
#-------------------------------------------------------------------------------
if not gsph:
    q = hydro.Q
    if not Cl is None:
        q.Cl = Cl
    if not Cq is None:
        q.Cq = Cq
    if not linearInExpansion is None:
        q.linearInExpansion = linearInExpansion
    if not quadraticInExpansion is None:
        q.quadraticInExpansion = quadraticInExpansion
    if not Qlimiter is None:
        q.limiter = Qlimiter
    if not epsilon2 is None:
        q.epsilon2 = epsilon2
    if not etaCritFrac is None:
        q.etaCritFrac = etaCritFrac
    if not etaFoldFrac is None:
        q.etaFoldFrac = etaFoldFrac
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
rPlane1 = Plane(Vector(x2), Vector(-1.0))
rbc1 = ReflectingBoundary(rPlane1)

for p in packages:
    p.appendBoundary(rbc1)

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
integrator.rigorousBoundaries = rigorousBoundaries
integrator.verbose = dtverbose
output("integrator")
output("integrator.havePhysicsPackage(hydro)")
if hourglass:
    output("integrator.havePhysicsPackage(hg)")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.rigorousBoundaries")

#-------------------------------------------------------------------------------
# Make the problem controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator,
                            kernel = hydro.kernel,
                            volumeType = volumeType,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle)
output("control")

#-------------------------------------------------------------------------------
# If requested, reset the intial energies and densities.
#-------------------------------------------------------------------------------
if sumInitialDensity and control.totalSteps == 0:
    packages = integrator.physicsPackages()
    state = State(db, packages)
    derivs = StateDerivatives(db, packages)
    integrator.setGhostNodes()
    integrator.preStepInitialize(state, derivs)
    cm = db.connectivityMap()
    pos = db.fluidPosition
    mass = db.fluidMass
    H = db.fluidHfield
    rho = db.fluidMassDensity
    computeSPHSumMassDensity(cm, W, True, pos, mass, H, rho)
    for nodes in nodeSet:
        pos = nodes.positions()
        eps = nodes.specificThermalEnergy()
        rho = nodes.massDensity()
        gammai = nodes.eos.gamma
        for i in range(nodes.numInternalNodes):
            eps[i] = specificEnergy(pos[i].x, rho[i],gammai)

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
if not steps is None:
    if checkRestart:
        control.setRestartBaseName(restartBaseName + "_CHECK")
    control.step(steps)
    if checkRestart:
        control.setRestartBaseName(restartBaseName)

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
dx1 = (x1 - x0)/nr1
dx2 = (x2 - x1)/nr2
h1 = 1.0/(nPerh*dx1)
h2 = 1.0/(nPerh*dx2)
answer = SodSolution(nPoints=nr1 + nr2,
                     gamma1 = gamma1,
                     gamma2 = gamma2,
                     rho1 = rho1,
                     P1 = P1,
                     rho2 = rho2,
                     P2 = P2,
                     x0 = x0,
                     x1 = x1,
                     x2 = x2,
                     h1 = 1.0/h1,
                     h2 = 1.0/h2)

#cs = db.newFluidScalarFieldList(0.0, "sound speed")
#db.fluidSoundSpeed(cs)
cs = hydro.soundSpeed

# Make a flat list from a FieldList
def createList(x):
    result = []
    for i in range(len(x)):
        for j in range(x[i].numInternalElements):
            result.append(x(i,j))
    return mpi.allreduce(result, mpi.SUM)

# Compute the simulated specific entropy.
gamma = db.newFluidScalarFieldList(0.0, "specific heat ratio")
db.fluidGamma(gamma)
gammaList = createList(gamma)

rho = createList(db.fluidMassDensity)
P = createList(hydro.pressure)
A = [Pi/rhoi**gammai for (Pi, rhoi,gammai) in zip(P, rho,gammaList)]

if graphics:
    from SpheralMatplotlib import *

    rhoPlot, velPlot, epsPlot, PPlot, HPlot = plotState(db)
    pE = plotEHistory(control.conserve)
    csPlot = plotFieldList(cs, winTitle="Sound speed")
    plots = [(rhoPlot, "Sod-spherical-rho.png"),
             (velPlot, "Sod-spherical-vel.png"),
             (epsPlot, "Sod-spherical-eps.png"),
             (PPlot, "Sod-spherical-P.png"),
             (HPlot, "Sod-spherical-h.png"),
             (csPlot, "Sod-spherical-cs.png")]
    
    if crksph:
        volPlot = plotFieldList(control.RKCorrections.volume, 
                                winTitle = "volume",
                                colorNodeLists = False, plotGhosts = False)
        splot = plotFieldList(control.RKCorrections.surfacePoint,
                              winTitle = "surface point",
                              colorNodeLists = False)
        plots += [(volPlot, "Sod-spherical-vol.png"),
                   (splot, "Sod-spherical-surfacePoint.png")]
    
    if not gsph:
        viscPlot = plotFieldList(hydro.maxViscousPressure,
                             winTitle = "max($\\rho^2 \pi_{ij}$)",
                             colorNodeLists = False)
        plots.append((viscPlot, "Sod-spherical-viscosity.png"))
    
    if boolCullenViscosity:
        cullAlphaPlot = plotFieldList(q.ClMultiplier,
                                      winTitle = "Cullen alpha")
        cullDalphaPlot = plotFieldList(evolveCullenViscosityMultiplier.DalphaDt,
                                       winTitle = "Cullen DalphaDt")
        plots += [(cullAlphaPlot, "Sod-spherical-Cullen-alpha.png"),
                  (cullDalphaPlot, "Sod-spherical-Cullen-DalphaDt.png")]

    if boolReduceViscosity:
        alphaPlot = plotFieldList(q.ClMultiplier,
                                  winTitle = "rvAlpha",
                                  colorNodeLists = False)

    # Make hardcopies of the plots.
    for p, filename in plots:
        savefig(p, os.path.join(dataDir, filename))

print("Energy conservation: original=%g, final=%g, error=%g" % (control.conserve.EHistory[0],
                                                                control.conserve.EHistory[-1],
                                                                (control.conserve.EHistory[-1] - control.conserve.EHistory[0])/control.conserve.EHistory[0]))

# #-------------------------------------------------------------------------------
# # If requested, write out the state in a global ordering to a file.
# #-------------------------------------------------------------------------------
# from SpheralTestUtilities import multiSort
# mof = mortonOrderIndices(db)
# mo = createList(mof)
# rhoprof = createList(db.fluidMassDensity)
# Pprof = createList(hydro.pressure)
# vprof = [v.x for v in createList(db.fluidVelocity)]
# epsprof = createList(db.fluidSpecificThermalEnergy)
# hprof = [1.0/Hi.xx for Hi in createList(db.fluidHfield)]

# rmin = x0
# rmax = x2
# if mpi.rank == 0:
#     multiSort(mo, xprof, rhoprof, Pprof, vprof, epsprof, hprof)
#     if outputFile:
#         outputFile = os.path.join(dataDir, outputFile)
#         f = open(outputFile, "w")
#         f.write(("#  " + 19*"'%s' " + "\n") % ("x", "rho", "P", "v", "eps", "A", "h", "mo",
#                                                "rhoans", "Pans", "vans", "Aans", "hans",
#                                                "x_UU", "rho_UU", "P_UU", "v_UU", "eps_UU", "h_UU"))
#         for (xi, rhoi, Pi, vi, epsi, Ai, hi, mi,
#              rhoansi, Pansi, vansi, Aansi, hansi) in zip(xprof, rhoprof, Pprof, vprof, epsprof, A, hprof, mo,
#                                                          rhoans, Pans, vans, Aans, hans):
#             f.write((7*"%16.12e " + "%i " + 5*"%16.12e " + 6*"%i " + '\n') % 
#                     (xi, rhoi, Pi, vi, epsi, Ai, hi, mi,
#                      rhoansi, Pansi, vansi, Aansi, hansi,
#                      unpackElementUL(packElementDouble(xi)),
#                      unpackElementUL(packElementDouble(rhoi)),
#                      unpackElementUL(packElementDouble(Pi)),
#                      unpackElementUL(packElementDouble(vi)),
#                      unpackElementUL(packElementDouble(epsi)),
#                      unpackElementUL(packElementDouble(hi))))
#         f.close()

#     import Pnorm
#     print "\tQuantity \t\tL1 \t\t\tL2 \t\t\tLinf"
#     failure = False
#     hD = []
#     for (name, data, ans) in [("Mass Density", rhoprof, rhoans),
#                               ("Pressure", Pprof, Pans),
#                               ("Velocity", vprof, vans),
#                               ("Thermal E", epsprof, uans),
#                               ("Entropy", A, Aans),
#                               ("h       ", hprof, hans)]:
#         assert len(data) == len(ans)
#         error = [data[i] - ans[i] for i in xrange(len(data))]
#         Pn = Pnorm.Pnorm(error, xprof)
#         L1 = Pn.gridpnorm(1, rmin, rmax)
#         L2 = Pn.gridpnorm(2, rmin, rmax)
#         Linf = Pn.gridpnorm("inf", rmin, rmax)
#         print "\t%s \t\t%g \t\t%g \t\t%g" % (name, L1, L2, Linf)
#         #f.write(("\t\t%g") % (L1))
#         hD.append([L1,L2,Linf])
#     #f.write("\n")

#     print "%d\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t" % (nr1+nr2,hD[0][0],hD[1][0],hD[2][0],hD[3][0],
#                                                                                 hD[0][1],hD[1][1],hD[2][1],hD[3][1],
#                                                                                 hD[0][2],hD[1][2],hD[2][2],hD[3][2])
