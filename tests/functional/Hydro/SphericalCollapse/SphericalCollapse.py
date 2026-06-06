#-------------------------------------------------------------------------------
# The spherical Sedov test case (3-D).
#-------------------------------------------------------------------------------
import os, sys, shutil
from Spheral3d import *
from findLastRestart import *
from SpheralTestUtilities import *
from SpheralGnuPlotUtilities import *
from GenerateNodeDistribution3d import *

import mpi

title("3-D integrated hydro test -- spherical collapse to disk")
#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(seed = "lattice",

            nx = 50,
            ny = 50,
            nz = 50,
            ns = 100000,
            nPerh = 1.51,
            KernelConstructor = BSplineKernel,
            order = 5,

            M0 = 1.0,
            G0 = 1.0,
            Rc = 0.01,
            R0 = Vector(0.0,0.0,0.0),

            rmin = 0.0,
            rmax = 1.0,

            fvelSupport = 1.0,

            rho0 = 1.0e-4,
            eps0 = 0.5,
            smallPressure = False,
            Espike = 1.0,
            smoothSpike = True,
            topHatSpike = False,
            smoothSpikeScale = 0.5,
            gamma = 5.0/3.0,
            mu = 1.0,

            Cl = 1.0,
            Cq = 0.75,
            epsilon2 = 1e-2,
            Qlimiter = False,
            balsaraCorrection = False,
            linearInExpansion = False,

            CRKSPH = False,
            PSPH = False,
            ASPH = False,     # Only for H evolution, not hydro algorithm
            Qconstructor = MonaghanGingoldViscosity,
            correctionOrder = LinearOrder,
            densityUpdate = RigorousSumDensity, # VolumeScaledDensity,
            HUpdate = IdealH,
            filter = 0.0,
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
            HopkinsConductivity = False,     # For PSPH
            evolveTotalEnergy = False,       # Only for SPH variants -- evolve total rather than specific energy

            hmin = 1e-15,
            hmax = 1.0,
            cfl = 0.5,
            deltaPhi = 0.01,
            useVelocityMagnitudeForDt = True,
            XSPH = False,
            rhomin = 1e-10,

            IntegratorConstructor = CheapSynchronousRK2Integrator,
            steps = None,
            goalTime = 10000.0,
            goalRadius = 0.8,
            dt = 1e-8,
            dtMin = 1.0e-8,
            dtMax = None,
            dtGrowth = 2.0,
            vizCycle = None,
            vizTime = 0.1,
            maxSteps = None,
            statsStep = 1,
            smoothIters = 0,
            HEvolution = IdealH,
            compatibleEnergy = True,
            gradhCorrection = True,
            correctVelocityGradient = True,

            restoreCycle = None,
            restartStep = 1000,

            graphics = True,
            clearDirectories = False,
            dataRoot = "dumps-spherical-collapse",
            outputFile = None,
            )


assert not(boolReduceViscosity and boolCullenViscosity)



#-------------------------------------------------------------------------------
# Set the hydro choice.
#-------------------------------------------------------------------------------
if CRKSPH:
    Qconstructor = LimitedMonaghanGingoldViscosity
    if ASPH:
        HydroConstructor = ACRKSPHHydro
    else:
        HydroConstructor = CRKSPHHydro
elif PSPH:
    if ASPH:
        HydroConstructor = APSPHHydro
    else:
        HydroConstructor = PSPHHydro
else:
    if ASPH:
        HydroConstructor = ASPHHydro
    else:
        HydroConstructor = SPHHydro

#-------------------------------------------------------------------------------
# Path names.
#-------------------------------------------------------------------------------
dataDir = os.path.join(dataRoot,
                       HydroConstructor.__name__,
                       Qconstructor.__name__,
                       "nperh=%4.2f" % nPerh,
                       "XSPH=%s" % XSPH,
                       "densityUpdate=%s" % densityUpdate,
                       "compatibleEnergy=%s" % compatibleEnergy,
                       "Cullen=%s" % boolCullenViscosity,
                       "seed=%s" % seed,
                       "nx=%i_ny=%i_nz=%i" % (nx, ny, nz))
restartDir = os.path.join(dataDir, "restarts")
vizDir = os.path.join(dataDir, "visit")
restartBaseName = os.path.join(restartDir, "Sedov-spherical-3d")

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
import os, sys
if mpi.rank == 0:
    if clearDirectories and os.path.exists(dataDir):
        shutil.rmtree(dataDir)
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
eos = GammaLawGasMKS(gamma, mu)

#-------------------------------------------------------------------------------
# Create our interpolation kernels -- one for normal hydro interactions, and
# one for use with the artificial viscosity
#-------------------------------------------------------------------------------
if KernelConstructor==NBSplineKernel:
  WT = TableKernel(NBSplineKernel(order), 1000)
else:
  WT = TableKernel(KernelConstructor(), 1000)
output("WT")
kernelExtent = WT.kernelExtent

#-------------------------------------------------------------------------------
# Create a NodeList and associated Neighbor object.
#-------------------------------------------------------------------------------
nodes1 = makeFluidNodeList("nodes1", eos, 
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
if restoreCycle is None:
    if seed == "lattice":
        generator = GenerateNodeDistribution3d(nx, ny, nz,
                                               rho0, seed,
                                               xmin = (-1.0, -1.0, -1.0),
                                               xmax = (1.0, 1.0, 1.0),
                                               rmin = 0.0,
                                               rmax = 1.0,
                                               nNodePerh = nPerh,
                                               SPH = (not ASPH))
    else:
        generator = GenerateIcosahedronMatchingProfile3d(n=ns,
                                                         densityProfileMethod=rho0,
                                                         rmin=rmin,rmax=rmax,
                                                         nNodePerh=nPerh)

    if mpi.procs > 1:
        from VoronoiDistributeNodes import distributeNodes3d
    else:
        from DistributeNodes import distributeNodes3d

    distributeNodes3d((nodes1, generator))
    output("mpi.reduce(nodes1.numInternalNodes, mpi.MIN)")
    output("mpi.reduce(nodes1.numInternalNodes, mpi.MAX)")
    output("mpi.reduce(nodes1.numInternalNodes, mpi.SUM)")

    # Set the velocities and energies
    v = nodes1.velocity()
    r = nodes1.positions()
    m = nodes1.mass()
    e = nodes1.specificThermalEnergy()
    omega = sqrt(G0*M0/rmax**3)
    for i in range(nodes1.numInternalNodes):
        x = r[i].x
        y = r[i].y
        z = r[i].z
        vx = -omega*y
        vy = omega*x
        v[i].x = fvelSupport*vx
        v[i].y = fvelSupport*vy
        rr = r[i].magnitude()
        rz = sqrt(x*x+y*y)
        theta = atan2(rz,z)
        if (rr > 0.0):
            P = rho0*G0*M0*(1.0/rr-1.0/rmax)*(1.0-fvelSupport*sin(theta))
            eps = P/rho0*(1.0/(gamma-1.0))
            e[i] = eps
        
    # Set the point source of energy.


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
# Create the gravity physics object.
#-------------------------------------------------------------------------------
gravity = PointPotential(G0, M0, Rc, R0)
gravity.deltaPotentialFraction = deltaPhi
output("gravity.G")
output("gravity.mass")
output("gravity.coreRadius")
output("gravity.origin")
output("gravity.deltaPotentialFraction")

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
if CRKSPH:
    hydro = HydroConstructor(W = WT,
                             Q = q,
                             filter = filter,
                             cfl = cfl,
                             compatibleEnergyEvolution = compatibleEnergy,
                             XSPH = XSPH,
                             correctionOrder = correctionOrder,
                             densityUpdate = densityUpdate,
                             HUpdate = HUpdate)
elif PSPH:
    hydro = HydroConstructor(W = WT,
                             Q = q,
                             filter = filter,
                             cfl = cfl,
                             useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                             compatibleEnergyEvolution = compatibleEnergy,
                             evolveTotalEnergy = evolveTotalEnergy,
                             HopkinsConductivity = HopkinsConductivity,
                             densityUpdate = densityUpdate,
                             HUpdate = HUpdate,
                             XSPH = XSPH)
else:
    hydro = HydroConstructor(W = WT, 
                             Q = q,
                             cfl = cfl,
                             compatibleEnergyEvolution = compatibleEnergy,
                             evolveTotalEnergy = evolveTotalEnergy,
                             gradhCorrection = gradhCorrection,
                             correctVelocityGradient = correctVelocityGradient,
                             densityUpdate = densityUpdate,
                             XSPH = XSPH,
                             HUpdate = HEvolution)
output("hydro")
output("hydro.kernel()")
output("hydro.PiKernel()")
output("hydro.cfl")
output("hydro.compatibleEnergyEvolution")
output("hydro.XSPH")
output("hydro.densityUpdate")
output("hydro.HEvolution")

packages = [hydro]

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
# Construct a time integrator, and add the one physics package.
#-------------------------------------------------------------------------------
packages.append(gravity)
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
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            restoreCycle = restoreCycle,
                            vizBaseName = "SphericalCollapse",
                            vizDir = vizDir,
                            vizStep = vizCycle,
                            vizTime = vizTime,
                            SPH = (not ASPH))
output("control")

#-------------------------------------------------------------------------------
# Finally run the problem and plot the results.
#-------------------------------------------------------------------------------

v = nodes1.velocity()
r = nodes1.positions()
m = nodes1.mass()
totalL = 0
for i in range(nodes1.numInternalNodes):
    x = r[i].x
    y = r[i].y
    z = r[i].z
    vx = v[i].x
    vy = v[i].y
    vz = v[i].z
    totalL += m[i]*sqrt((y*vz-z*vy)**2+(z*vx-x*vz)**2+(x*vy-y*vx)**2)

print("Total L=%e" % totalL)

if steps is None:
    control.advance(goalTime, maxSteps)
    if restoreCycle != control.totalSteps:
        control.updateViz(control.totalSteps, integrator.currentTime, 0.0)
        control.dropRestartFile()
else:
    control.step(steps)

control.conserve.writeHistory("collapseHistory-CRK-%s" % CRKSPH)

# Output the energy conservation.
print("Energy conservation: ", ((control.conserve.EHistory[-1] -
                                 control.conserve.EHistory[0])/
                                control.conserve.EHistory[0]))

totalL = 0
for i in range(nodes1.numInternalNodes):
    x = r[i].x
    y = r[i].y
    z = r[i].z
    vx = v[i].x
    vy = v[i].y
    vz = v[i].z
    totalL += m[i]*sqrt((y*vz-z*vy)**2+(z*vx-x*vz)**2+(x*vy-y*vx)**2)

print("Total L=%e" % totalL)

