#-------------------------------------------------------------------------------
# testJaumannRate
#
# Impose a pure-rotational velocity field and test that the Jaumann rate
# formulation for the deviatoric stress evoultion is correct
#
# Based on an example originally submitted in a Spheral bug report:
# https://github.com/llnl/spheral/issues/534
#-------------------------------------------------------------------------------

#ATS:t0 = test(SELF, "--hydroType SPH",             label="Jaumann rate rotation unit test SPH (serial)")
#ATS:t1 = test(SELF, "--hydroType SPH",             label="Jaumann rate rotation unit test SPH (4 proc)", np=4)
#ATS:t2 = test(SELF, "--hydroType SPH --raja True", label="Jaumann rate rotation unit test SPH and RAJA (serial)")
#ATS:t3 = test(SELF, "--hydroType CRKSPH",          label="Jaumann rate rotation unit test CRKSPH (serial)")
#ATS:t4 = test(SELF, "--hydroType CRKSPH",          label="Jaumann rate rotation unit test CRKSPH (4 proc)", np=4)
#ATS:t5 = test(SELF, "--hydroType FSISPH",          label="Jaumann rate rotation unit test FSISPH (serial)")
#ATS:t6 = test(SELF, "--hydroType FSISPH",          label="Jaumann rate rotation unit test FSISPH (4 proc)", np=4)

from Spheral2d import *
from SpheralController import SpheralController
from SpheralTestUtilities import *

from GenerateNodeDistribution2d import *
from DistributeNodes import distributeNodes2d

import os, mpi
import numpy as np

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(hydroType = "SPH",
            correctionOrder = LinearOrder,
            raja = False,
            omega = 3.0,
            S0 = 1.0,
            nx = 40,
            x0 = -1.0,
            x1 =  1.0,
            y0 = -1.0,
            y1 =  1.0,
            rho0 = 1.0,
            nPerh = 4.01,
            tol = 1.0e-8,
            vizName = None,
            vizDir = None)

hydroType = hydroType.upper()
if vizName:
    vizName = vizName + "_" + hydroType
    if raja:
        vizName += "_RAJA"

#-------------------------------------------------------------------------------
# If needed, prepare viz directory
#-------------------------------------------------------------------------------
if mpi.rank == 0 and vizDir:
    if not os.path.exists(vizDir):
        os.makedirs(vizDir)
mpi.barrier()

#-------------------------------------------------------------------------------
# Material properties
#-------------------------------------------------------------------------------
units = CGuS()
eos = LinearPolynomialEquationOfState(rho0, 0.01, 100.0,
                                      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 55.0, units)
strength = ConstantStrength(10.0, 1.0e30)

#-------------------------------------------------------------------------------
# Interpolation kernel
#-------------------------------------------------------------------------------
WT = TableKernel(WendlandC4Kernel(), 100)

#-------------------------------------------------------------------------------
# Make the NodeList
#-------------------------------------------------------------------------------
nodes = makeSolidNodeList("patch", eos, strength,
                          nPerh = nPerh,
                          hmin = 1e-5,
                          hmax = 1.0,
                          rhoMin = 0.01,
                          rhoMax = 100.0)

#-------------------------------------------------------------------------------
# Create nodes and initial conditions
#-------------------------------------------------------------------------------
gen = GenerateNodeDistribution2d(nx, nx, rho0, "lattice",
                                 xmin = (x0, y0),
                                 xmax = (x1, y1),
                                 nNodePerh = nPerh)

distributeNodes2d((nodes, gen))
output("mpi.reduce(nodes.numInternalNodes, mpi.MIN)")
output("mpi.reduce(nodes.numInternalNodes, mpi.MAX)")
output("mpi.reduce(nodes.numInternalNodes, mpi.SUM)")

db = DataBase()
db.appendNodeList(nodes)

pos = nodes.positions()
vel = nodes.velocity()
S = nodes.deviatoricStress()
for i in range(nodes.numInternalNodes):
    posi = pos[i]
    vel[i] = Vector(-omega * posi.y, omega * posi.x)
    S[i] = SymTensor(S0, 0.0, 0.0, -S0)

#-------------------------------------------------------------------------------
# Hydro package
#-------------------------------------------------------------------------------
if hydroType == "SPH":
    hydro = SPH(dataBase = db,
                W = WT,
                RAJA = raja)
elif hydroType == "CRKSPH":
    hydro = CRKSPH(dataBase = db,
                   W = WT,
                   order = correctionOrder)
elif hydroType == "FSISPH":
    hydro = FSISPH(dataBase = db,
                   W = WT,
                   interfaceMethod = HLLCInterface,
                   sumDensityNodeLists=[nodes],                       
                   densityStabilizationCoefficient = 0.00)
else:
    raise ValueError("Unknown hydro type passed: " + hydroType)

#-------------------------------------------------------------------------------
# Time integrator
#-------------------------------------------------------------------------------
integrator = CheapSynchronousRK2Integrator(db)
integrator.appendPhysicsPackage(hydro)
integrator.lastDt = 1.0e-7

#-------------------------------------------------------------------------------
# Build the controller and take a tiny step to test the deviatoric stress evolution
#-------------------------------------------------------------------------------
control = SpheralController(integrator,
                            WT,
                            vizDir = vizDir,
                            vizBaseName = vizName,
                            vizStep = 1,
                            vizDerivs = True)

before = np.array([Si.xy for Si in nodes.deviatoricStress().internalValues()])
control.advance(control.time() + 1.0e-7)
after = np.array([Si.xy for Si in nodes.deviatoricStress().internalValues()])

rate = (after - before) / control.time()
ans = 2.0 * omega * S0
print("\n  dS_xy/dt measured       = ", rate)
print("  co-rotation +2*omega*S0 = ", ans)
print("  ratio to co-rotation    = ", (rate / ans), "\n")

if np.any(abs(rate/ans - 1.0) > tol):
    raise ValueError("Ratio of (dS_xy/dt)/ans out of bounds")
else:
    print("PASSED")
