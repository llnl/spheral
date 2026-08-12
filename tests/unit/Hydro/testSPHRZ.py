#-------------------------------------------------------------------------------
# Unit test of derivatives for analytic initial conditions
#-------------------------------------------------------------------------------
import os, shutil, mpi
import numpy as np
from SpheralRZ import *
from SpheralTestUtilities import *

from GenerateNodeDistribution2d import *
from DistributeNodes import distributeNodes2d
from SpheralMatplotlib import *

title("RZ SPH hydro unit test")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(geometry = "planar",     # one of (planar, cylindrical, spherical)
            problem = "linear",
            KernelConstructor = WendlandC4Kernel,

            n1 = 100,
            n2 = 20,

            seed = "lattice",
            nPerh = 4.01,

            gamma = 5.0/3.0,
            mu = 1.0,

            densityUpdate = RigorousSumDensity,

            # Problem parameters
            rho0 = 1.0,
            P0 = 5.0,
            Pslope = -1.0,
            )

#-------------------------------------------------------------------------------
# Material properties.
#-------------------------------------------------------------------------------
eos = GammaLawGasMKS(gamma, mu, minimumPressure=0.0)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
WT = TableKernel(KernelConstructor(), 200)
output("WT")

#-------------------------------------------------------------------------------
# Make the NodeList.
#-------------------------------------------------------------------------------
nodes1 = makeFluidNodeList("nodes1", eos, 
                           nPerh = nPerh,
                           xmin = Vector(-10.0, -10.0),
                           xmax = Vector( 10.0,  10.0),
                           kernelExtent = WT.kernelExtent)
    
output("nodes1")
output("nodes1.hmin")
output("nodes1.hmax")
output("nodes1.nodesPerSmoothingScale")

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
if geometry == "planar":
    nz = n1
    nr = n2
    z0, z1 = 0.0, 1.0
    r0, r1 = 0.0, 0.2
    rmin, rmax = None, None
elif geometry == "cylindrical":
    nz = n2
    nr = n1
    z0, z1 = 0.0, 0.2
    r0, r1 = 0.0, 1.0
    rmin, rmax = None, None
else:
    assert geometry == "spherical"
    nz = n1
    nr = n1
    rmin, rmax = 0.0, 1.0
    z0, z1 = 0.0, 1.0
    r0, r1 = 0.0, 1.0

generator = GenerateNodeDistribution2d(nz, nr, rho0, seed,
                                       xmin = (z0, r0),
                                       xmax = (z1, r1),
                                       rmin = rmin,
                                       rmax = rmax,
                                       nNodePerh = nPerh)

distributeNodes2d((nodes1, generator))
output("mpi.reduce(nodes1.numInternalNodes, mpi.MIN)")
output("mpi.reduce(nodes1.numInternalNodes, mpi.MAX)")
output("mpi.reduce(nodes1.numInternalNodes, mpi.SUM)")

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
            W = WT,
            densityUpdate = densityUpdate)
output("hydro")
output("hydro.cfl")
output("hydro.compatibleEnergyEvolution")
output("hydro.densityUpdate")
output("hydro.XSPH")
output("hydro._smoothingScaleMethod.HEvolution")

packages = [hydro]

#-------------------------------------------------------------------------------
# Create boundary conditions.
#-------------------------------------------------------------------------------
if geometry == "planar":
    bcs = [ReflectingBoundary(Plane(Vector(z0, r1), Vector( 0.0, -1.0))),
           ReflectingBoundary(Plane(Vector(z1, r0), Vector(-1.0,  0.0)))]
    if z0 == 0.0:
        bcs.append(ReflectingBoundary(Plane(Vector(z0, r0), Vector( 1.0,  0.0))))
    if r0 != 0.0:
        bcs.append(ReflectingBoundary(Plane(Vector(z0, r0), Vector( 0.0, 1.0))))
elif geometry == "cylindrical":
    bcs = [ReflectingBoundary(Plane(Vector(z0, r0), Vector( 1.0,  0.0))),
           ReflectingBoundary(Plane(Vector(z1, r0), Vector(-1.0,  0.0))),
           ReflectingBoundary(Plane(Vector(z0, r1), Vector( 0.0, -1.0)))]
else:
    assert geometry == "spherical"
    bcs = [ReflectingBoundary(Plane(Vector(z0, r0), Vector( 1.0,  0.0)))]

for bc in bcs:
    for p in packages:
        p.appendBoundary(bc)

#-------------------------------------------------------------------------------
# Construct an integrator.
#-------------------------------------------------------------------------------
integrator = ForwardEulerIntegrator(db, packages)
output("integrator")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.dtGrowth")
output("integrator.domainDecompositionIndependent")
output("integrator.verbose")
output("integrator.allowDtCheck")

#-------------------------------------------------------------------------------
# Make the problem controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT)
output("control")

#-------------------------------------------------------------------------------
# Solution object
#-------------------------------------------------------------------------------
class Answer:
    def pressure(self, x):
        return P0 + Pslope*x
    def rho(self, x):
        if x is np.array:
            return np.array([rho0]*len(x))
        else:
            return rho0
    def eps(self, x):
        return self.pressure(x)/((gamma - 1.0)*rho0)
    def vel(self, x):
        if x is np.array:
            return np.array([0.0]*len(x))
        else:
            return 0.0
    def answer(self, t, x):
        return x, self.vel(x), self.eps(x), np.array([rho0]*len(x)), self.pressure(x), np.array([h0]*len(x))
            
answer = Answer()

#-------------------------------------------------------------------------------
# Generate initial conditions
#-------------------------------------------------------------------------------
pos = nodes1.positions()
eps = nodes1.specificThermalEnergy()
if geometry == "planar":
    x = np.array([x.x for x in pos])
    eps.assign(list(answer.eps(x)))
elif geometry == "cylindrical":
    x = np.array([x.y for x in pos])
    eps.assign(list(answer.eps(x)))

#-------------------------------------------------------------------------------
# Generate the initial derivatives
#-------------------------------------------------------------------------------
packages = integrator.physicsPackages()
state = State(db, packages)
derivs = StateDerivatives(db, packages)
hydro.initializeProblemStartupDependencies(db, state, derivs)
integrator.setGhostNodes()
integrator.applyGhostBoundaries(state, derivs)
integrator.preStepInitialize(state, derivs)
dt = integrator.selectDt(integrator.dtMin, integrator.dtMax, state, derivs)
integrator.initializeDerivatives(0.0, dt, state, derivs)
derivs.Zero()
integrator.evaluateDerivatives(0.0, dt, db, state, derivs)
integrator.finalizeDerivatives(0.0, dt, db, state, derivs)

# Plot stuff
if geometry == "planar":
    xfunc = "%s.x"
    yfunc = "%s.x"
elif geometry == "cylindrical":
    xfunc = "%s.y"
    yfunc = "%s.y"

rhoPlot, velPlot, epsPlot, Pplot, Hplot = plotState(state,
                                                    xFunction = xfunc,
                                                    vecyFunction = yfunc)
DvDt = derivs.vectorFields(HydroFieldNames.hydroAcceleration)
zaccPlot = plotFieldList(DvDt,
                         xFunction = xfunc,
                         yFunction = "%s.x",
                         winTitle = "Acceleration (z)")
raccPlot = plotFieldList(DvDt,
                         xFunction = xfunc,
                         yFunction = "%s.y",
                         winTitle = "Acceleration (r)")
rhoRZ = state.scalarFields(HydroFieldNames.massDensityRZ)
rhoRZplot = plotFieldList(rhoRZ,
                          xFunction = xfunc,
                          winTitle = "Mass density (RZ)")
