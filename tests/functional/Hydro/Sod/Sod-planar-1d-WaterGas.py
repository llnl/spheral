
# Solid FSISPH
#
#ATS:fsisph1 = test(           SELF, "--fsisph True --solid True --nx1 500 --nx2 30 --cfl 0.45 --graphics None --clearDirectories True  --restartStep 20 --steps 40", label="Planar Water-Gas Sod problem with FSISPH -- 1-D (serial)", fsisph=True)
#ATS:fsisph2 = testif(fsisph1, SELF, "--fsisph True --solid True --nx1 500 --nx2 30 --cfl 0.45 --graphics None --clearDirectories False --restartStep 20 --steps 20 --restoreCycle 20 --checkRestart True", label="Planar Water-Gas Sod problem with FSISPH -- 1-D (serial) RESTART CHECK", fsisph=True)
#
# GSPH
#
#ATS:gsph1 = test(         SELF, "--gsph True --nx1 500 --nx2 30 --cfl 0.45 --graphics None --clearDirectories True  --restartStep 20 --steps 40", label="Planar Water-Gas Sod problem with GSPH -- 1-D (serial)", gsph=True)
#ATS:gsph2 = testif(gsph1, SELF, "--gsph True --nx1 500 --nx2 30 --cfl 0.45 --graphics None --clearDirectories False --restartStep 20 --steps 20 --restoreCycle 20 --checkRestart True", label="Planar Water-Gas Sod problem with GSPH -- 1-D (serial) RESTART CHECK", gsph=True)
#

import os, sys
import shutil
from SolidSpheral1d import *
from SpheralTestUtilities import *
from SodAnalyticSolution import SodSolutionStiffGasStiffGas as SodSolution

from GenerateNodeProfile import GenerateNodeProfile1d
from VoronoiDistributeNodes import distributeNodes1d

title("1-D integrated hydro test -- planar Sod problem Water-Gas")

#-------------------------------------------------------------------------------
# Generic problem parameters -- defaults are from the gas-water shock test of 
# Liang & Chen, "Flow visualization of shock/water column interactions" Shocks
# Waves, 2008.
#-------------------------------------------------------------------------------
            # materials
commandLine(nx1 = 1000,     # number of nodes
            nx2 = 60,
            rho1 = 1000.0,  # density
            rho2 = 50.0,
            P1 = 1.0e9,     # Pressure
            P2 = 1.0e5,
            Po1 = 6.6e8,    # stiffening pressure
            Po2 = 0.0,
            gamma1 = 4.4,   # specific heat ratio
            gamma2 = 1.4,
            Cv1 = 4184.0,   # specific heat (irrevlant on this prob need for initialization)
            Cv2 = 50.0,     
            
            # domain
            x0 = -0.5,    # bound1 
            x1 = 0.2,     # interface    
            x2 = 0.5,     # bound2

            # kernel things
            KernelConstructor = NBSplineKernel,  # (NBSplineKernel, WendlandC2Kernel, WendlandC4Kernel)
            order = 5,                           # spline order for NBSplineKernel
            nPerh = 1.35,                        # neighbors per smoothing scale
            HUpdate = IdealH,                    # (IdealH, IntegrateH) how do we update smoothing scale?
            iterateInitialH = True,              # do we want to find an ideal H to start 
            gradhCorrection = True,              # correct for adaptive-h (not implemented in FSISPH)
            hmin = 1e-10,                    
            hmax = 1.0,

            # hydro type (if all false default to SPH)
            crksph = False,   # based on conservative formulation w/ repoducing kernels
            psph = False,     # pressure based sph
            fsisph = False,   # multimaterial patching method
            gsph = False,     # Convolution-free godunov-SPH
            mfm = False,      # moving finite mass -- hopkins 2015

            # FSI parameters
            fsiRhoStabilizeCoeff = 0.1,         # diffusion operating through the vel-gradient
            fsiEpsDiffuseCoeff = 0.1,           # diffusion coeff for specific thermal energy
            fsiXSPHCoeff = 0.0,                 # ramps xsph up
            fsiInterfaceMethod = HLLCInterface, # (HLLCInterface, ModulusInterface, NoInterface)

            # CRK parameters
            correctionOrder = LinearOrder,
            volumeType = RKSumVolume,

            # GSPH parameters
            gsphReconstructionGradient = RiemannGradient,
            linearReconstruction = True,

            # general hydro parameters
            densityUpdate = IntegrateDensity,     # (RigorousSumDensity, IntegrateDensity, CorrectedSumDensity)
            correctVelocityGradient = True,       # linear correction to velocity gradient
            compatibleEnergy = True,              # conservative specific energy updat method (reccomended)
            evolveTotalEnergy = False,            # evolve total rather than specific energy
            solid = False,                        # turns on solid nodelist & hydro objects
            XSPH = False,                         # position updated w/ smoothed velocity
            epsilonTensile = 0.0,
            nTensile = 8,
            rhoMin = 0.01,
            filter = 0.00,

            # artificial viscosity
            Cl = None,
            Cq = None,
            boolReduceViscosity = False,
            nh = 5.0,
            aMin = 0.1,
            aMax = 2.0,
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
            boolHopkinsCorrection = True,
            
            # integrator settings
            IntegratorConstructor = CheapSynchronousRK2Integrator,            
            cfl = 0.25,
            goalTime = 240e-6,
            dt = 1.0e-12,
            dtMin = 1.0e-12,
            dtMax = 0.1,
            dtGrowth = 2.0,
            steps = None,
            maxSteps = None,
            statsStep = 10,
            rigorousBoundaries = False,
            dtverbose = False,
            
            # IO settings
            clearDirectories = False,
            restoreCycle = -1,
            restartStep = 10000,
            dataDirBase = "dumps-Sod-planar",
            restartBaseName = "Sod-planar-1d-restart",
            outputFile = None,
            checkRestart = False,
            graphics = True,
            )

assert not(boolReduceViscosity and boolCullenViscosity)
assert not (fsisph and not solid)              # only implemented for solid
assert (mfm + fsisph + crksph + psph + gsph <= 1)    # only one hydro selection

if crksph:
    hydroname = os.path.join("CRKSPH",
                             str(volumeType),
                             str(correctionOrder))
elif psph:
    hydroname = "PSPH"
elif fsisph:
    hydroname = "FSISPH"
elif gsph:
    hydroname = os.path.join("GSPH",str(gsphReconstructionGradient))
elif mfm:
    hydroname = os.path.join("MFM",str(gsphReconstructionGradient))
else:
    hydroname = "SPH"
if solid:
    hydroname = "Solid" + hydroname

dataDir = os.path.join(dataDirBase, 
                       hydroname,
                       "nPerh=%f" % nPerh,
                       "compatibleEnergy=%s" % compatibleEnergy,
                       "Cullen=%s" % boolCullenViscosity,
                       "%i" % (nx1 + nx2))
restartDir = os.path.join(dataDir, "restarts")
restartBaseName = os.path.join(restartDir, "Sod-planar-1d-WG-%i" % (nx1 + nx2))

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
units = MKS()
eos1 = StiffenedGas(gamma=gamma1, 
                    P0=Po1, 
                    Cv=Cv1,
                    constants=units)
eos2 = StiffenedGas(gamma=gamma2, 
                    P0=Po2, 
                    Cv=Cv2,
                    constants=units)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
if KernelConstructor==NBSplineKernel:
    assert order in (3,5,7)
    WT = TableKernel(KernelConstructor(order), 1000)
else:
    WT = TableKernel(KernelConstructor(), 1000) 
kernelExtent = WT.kernelExtent
output("WT")

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
dx1 = (x1 - x0)/nx1
dx2 = (x2 - x1)/nx2

def rho_initial(xi):
    if xi < x1:
        return rho1
    else:
        return rho2

def specificEnergy(xi, rhoi):
    if xi < x1:
        Pi = P1
        Poi = Po1
        gammai = gamma1
    else:
        Pi = P2
        Poi = Po2
        gammai = gamma2
    return (Pi+gammai*Poi)/((gammai - 1.0)*rhoi)

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
gen1 = GenerateNodeProfile1d(nx = nx1,
                             rho = rho1,
                             xmin = x0,
                             xmax = x1,
                             nNodePerh = nPerh)
gen2 = GenerateNodeProfile1d(nx = nx2,
                             rho = rho2,
                             xmin = x1,
                             xmax = x2,
                             nNodePerh = nPerh)
distributeNodes1d((nodes1, gen1),(nodes2, gen2))
output("nodes1.numNodes")
output("nodes2.numNodes")

# Set node specific thermal energies
for nodes in nodeSet:
    pos = nodes.positions()
    eps = nodes.specificThermalEnergy()
    rho = nodes.massDensity()
    for i in range(nodes.numInternalNodes):
        eps[i] = specificEnergy(pos[i].x, rho[i])

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
                   filter = filter,
                   cfl = cfl,
                   compatibleEnergyEvolution = compatibleEnergy,
                   evolveTotalEnergy = evolveTotalEnergy,
                   XSPH = XSPH,
                   densityUpdate = RigorousSumDensity,#densityUpdate,
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
                 correctVelocityGradient = correctVelocityGradient)
elif fsisph:
    hydro = FSISPH(dataBase = db,
                   W = WT,
                   cfl = cfl,
                   sumDensityNodeLists=[],                       
                   densityStabilizationCoefficient = fsiRhoStabilizeCoeff,
                   specificThermalEnergyDiffusionCoefficient = fsiEpsDiffuseCoeff,
                   xsphCoefficient = fsiXSPHCoeff,
                   interfaceMethod = fsiInterfaceMethod,
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
                compatibleEnergyEvolution = compatibleEnergy,
                correctVelocityGradient=correctVelocityGradient,
                evolveTotalEnergy = evolveTotalEnergy,
                XSPH = XSPH,
                gradientType = gsphReconstructionGradient,
                densityUpdate=densityUpdate,
                HUpdate = IdealH,
                epsTensile = epsilonTensile,
                nTensile = nTensile)
elif mfm:
    limiter = VanLeerLimiter()
    waveSpeed = DavisWaveSpeed()
    solver = HLLC(limiter,waveSpeed,linearReconstruction)
    hydro = MFM(dataBase = db,
                riemannSolver = solver,
                W = WT,
                cfl=cfl,
                compatibleEnergyEvolution = compatibleEnergy,
                correctVelocityGradient=correctVelocityGradient,
                evolveTotalEnergy = evolveTotalEnergy,
                XSPH = XSPH,
                gradientType = gsphReconstructionGradient,
                densityUpdate=densityUpdate,
                HUpdate = IdealH,
                epsTensile = epsilonTensile,
                nTensile = nTensile)
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
output("hydro.compatibleEnergyEvolution")
packages = [hydro]

#-------------------------------------------------------------------------------
# Tweak the artificial viscosity.
#-------------------------------------------------------------------------------
if not (gsph or mfm):
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
    # Construct the MMRV physics object.
    #-------------------------------------------------------------------------------
    if boolReduceViscosity:
        evolveReducingViscosityMultiplier = MorrisMonaghanReducingViscosity(nh,nh,aMin,aMax)
        packages.append(evolveReducingViscosityMultiplier)
    elif boolCullenViscosity:
        evolveCullenViscosityMultiplier = CullenDehnenViscosity(WT,alphMax,alphMin,betaC,betaD,betaE,fKern,boolHopkinsCorrection)
        packages.append(evolveCullenViscosityMultiplier)



#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
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
integrator.rigorousBoundaries = rigorousBoundaries
integrator.verbose = dtverbose
output("integrator")
output("integrator.havePhysicsPackage(hydro)")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.rigorousBoundaries")

#-------------------------------------------------------------------------------
# Make the problem controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator,
                            iterateInitialH=iterateInitialH,
                            kernel = WT,
                            volumeType = volumeType,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle)
output("control")


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
dx1 = (x1 - x0)/nx1
dx2 = (x2 - x1)/nx2
h1 = 1.0/(nPerh*dx1)
h2 = 1.0/(nPerh*dx2)
answer = SodSolution(nPoints=nx1 + nx2,
                     gamma1 = gamma1,
                     gamma2 = gamma2,
                     rho1 = rho1,
                     P1 = P1,
                     Po1 = Po1,
                     rho2 = rho2,
                     P2 = P2,
                     Po2 =Po2,
                     x0 = x0,
                     x1 = x1,
                     x2 = x2,
                     h1 = 1.0/h1,
                     h2 = 1.0/h2)



# Make a flat list from a FieldList
def createList(x):
    result = []
    for i in range(len(x)):
        for j in range(x[i].numInternalElements):
            result.append(x(i,j))
    return mpi.allreduce(result, mpi.SUM)

# Compute the simulated specific entropy.
#gamma = db.newFluidScalarFieldList(0.0, "specific heat ratio")
#db.fluidGamma(gamma)
#gammaList = createList(gamma)

#cs = hydro.soundSpeed

#rho = createList(db.fluidMassDensity)
#P = createList(hydro.pressure)
#A = [Pi/rhoi**gammai for (Pi, rhoi,gammai) in zip(P, rho,gammaList)]

# The analytic solution for the simulated entropy.
xprof = [x.x for x in createList(db.fluidPosition)]
xans, vans, uans, rhoans, Pans, hans, gammaans = answer.solution(control.time(), xprof)
#Aans = [Pi/rhoi**gammai for (Pi, rhoi, gammai) in zip(Pans,  rhoans, gammaans)]
#csAns = [sqrt(gammai*Pi/rhoi) for (Pi, rhoi, gammai) in zip(Pans,  rhoans, gammaans)]

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
    
    if crksph:
        volPlot = plotFieldList(control.RKCorrections.volume, 
                                winTitle = "volume",
                                colorNodeLists = False, plotGhosts = False)
        splot = plotFieldList(control.RKCorrections.surfacePoint,
                              winTitle = "surface point",
                              colorNodeLists = False)
        plots += [(volPlot, "Sod-planar-vol.png"),
                   (splot, "Sod-planar-surfacePoint.png")]
    
    viscPlot = plotFieldList(hydro.maxViscousPressure,
                             winTitle = "max($\\rho^2 \pi_{ij}$)",
                             colorNodeLists = False)
    plots.append((viscPlot, "Sod-planar-viscosity.png"))
    
    if boolCullenViscosity:
        cullAlphaPlot = plotFieldList(q.ClMultiplier,
                                      winTitle = "Cullen alpha")
        cullDalphaPlot = plotFieldList(evolveCullenViscosityMultiplier.DalphaDt,
                                       winTitle = "Cullen DalphaDt")
        plots += [(cullAlphaPlot, "Sod-planar-Cullen-alpha.png"),
                  (cullDalphaPlot, "Sod-planar-Cullen-DalphaDt.png")]

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

#-------------------------------------------------------------------------------
# If requested, write out the state in a global ordering to a file.
#-------------------------------------------------------------------------------
from SpheralTestUtilities import multiSort
mof = mortonOrderIndices(db)
mo = createList(mof)
rhoprof = createList(db.fluidMassDensity)
Pprof = createList(hydro.pressure)
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
        f.write(("#  " + 17*"'%s' " + "\n") % ("x", "rho", "P", "v", "eps", "h", "mo",
                                               "rhoans", "Pans", "vans", "hans",
                                               "x_UU", "rho_UU", "P_UU", "v_UU", "eps_UU", "h_UU"))
        for (xi, rhoi, Pi, vi, epsi, hi, mi,
             rhoansi, Pansi, vansi, hansi) in zip(xprof, rhoprof, Pprof, vprof, epsprof, hprof, mo,
                                                         rhoans, Pans, vans, hans):
            f.write((6*"%16.12e " + "%i " + 4*"%16.12e " + 6*"%i " + '\n') % 
                    (xi, rhoi, Pi, vi, epsi, hi, mi,
                     rhoansi, Pansi, vansi,  hansi,
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
