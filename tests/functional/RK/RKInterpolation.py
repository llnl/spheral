#ATS:test(SELF, "--dimension 1 --correctionOrder ZerothOrder --funcType constant", label="RK interpolation - 1D zeroth")
#ATS:test(SELF, "--dimension 2 --correctionOrder ZerothOrder --funcType constant --numToCheck 10", label="RK interpolation - 2D zeroth")
#ATS:test(SELF, "--dimension 3 --correctionOrder ZerothOrder --funcType constant --numToCheck 10", label="RK interpolation - 3D zeroth")

#ATS:test(SELF, "--dimension 2 --correctionOrder LinearOrder --funcType linear", label="RK interpolation - 2D linear")
#ATS:test(SELF, "--dimension 3 --correctionOrder LinearOrder --funcType linear --numToCheck 10", label="RK interpolation - 3D linear")
#ATS:test(SELF, "--dimension 1 --correctionOrder LinearOrder --funcType linear --numToCheck 10", label="RK interpolation - 1D linear")

#ATS:test(SELF, "--dimension 1 --correctionOrder QuadraticOrder --funcType quadratic --testHessian True", label="RK interpolation - 1D quadratic")
#ATS:test(SELF, "--dimension 2 --correctionOrder QuadraticOrder --funcType quadratic --testHessian True --numToCheck 10", label="RK interpolation - 2D quadratic")
#ATS:test(SELF, "--dimension 3 --correctionOrder QuadraticOrder --funcType quadratic --testHessian True --numToCheck 10", label="RK interpolation - 3D quadratic")

#ATS:test(SELF, "--dimension 1 --correctionOrder CubicOrder --funcType cubic --testHessian True", label="RK interpolation - 1D cubic")
#ATS:test(SELF, "--dimension 2 --correctionOrder CubicOrder --funcType cubic --testHessian True --numToCheck 10", label="RK interpolation - 2D cubic")
#ATS:test(SELF, "--dimension 3 --correctionOrder CubicOrder --funcType cubic --testHessian True --numToCheck 10", label="RK interpolation - 3D cubic")

#ATS:test(SELF, "--dimension 1 --correctionOrder QuarticOrder --funcType quartic --nPerh 5.01 --testHessian True", label="RK interpolation - 1D quartic")
#ATS:test(SELF, "--dimension 2 --correctionOrder QuarticOrder --funcType quartic --nPerh 5.01 --testHessian True --numToCheck 10", level=100, label="RK interpolation - 2D quartic")
#ATS:test(SELF, "--dimension 3 --correctionOrder QuarticOrder --funcType quartic --nPerh 5.01 --testHessian True --numToCheck 10", level=100, label="RK interpolation - 3D quartic")

#ATS:test(SELF, "--dimension 1 --correctionOrder QuinticOrder --funcType quintic --nPerh 6.01 --tolerance 1.e-11 --testHessian True", label="RK interpolation - 1D quintic")
#ATS:test(SELF, "--dimension 2 --correctionOrder QuinticOrder --funcType quintic --nPerh 6.01 --tolerance 1.e-11 --testHessian True --numToCheck 10", level=100, label="RK interpolation - 2D quintic")
#ATS:test(SELF, "--dimension 3 --correctionOrder QuinticOrder --funcType quintic --nPerh 6.01 --tolerance 1.e-11 --testHessian True --numToCheck 10", level=100, label="RK interpolation - 3D quintic")

#ATS:test(SELF, "--dimension 1 --correctionOrder SexticOrder --funcType sextic --nPerh 7.01 --tolerance 1.e-10 --testHessian True", label="RK interpolation - 1D sextic")
#ATS:test(SELF, "--dimension 2 --correctionOrder SexticOrder --funcType sextic --nPerh 7.01 --tolerance 1.e-10 --testHessian True --numToCheck 10", level=100, label="RK interpolation - 2D sextic")
#ATS:test(SELF, "--dimension 3 --correctionOrder SexticOrder --funcType sextic --nPerh 7.01 --tolerance 1.e-10 --testHessian True --numToCheck 10", level=100, label="RK interpolation - 3D sextic")

#ATS:test(SELF, "--dimension 1 --correctionOrder SepticOrder --funcType septic --nx 12 --nPerh 8.01 --tolerance 2.e-8 --testHessian True", label="RK interpolation - 1D septic")
#ATS:test(SELF, "--dimension 2 --correctionOrder SepticOrder --funcType septic --nPerh 8.01 --tolerance 2.e-8 --testHessian True --numToCheck 10", level=100, label="RK interpolation - 2D septic")
#ATS:test(SELF, "--dimension 3 --correctionOrder SepticOrder --funcType septic --nPerh 8.01 --tolerance 2.e-8 --testHessian True --numToCheck 10", level=100, label="RK interpolation - 3D septic")

#-------------------------------------------------------------------------------
# Test of interpolation for reproducing kernels
#-------------------------------------------------------------------------------
from Spheral import *
from SpheralTestUtilities import *
import os, shutil
import numpy as np
import time
from matplotlib import pyplot as plt
title("Reproducing kernel test")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(
    # Spatial parameters
    dimension = 1,
    nx = 10,
    ny = 10,
    nz = 10,
    x0 = -2.0,
    x1 = 2.0,
    y0 = -2.0,
    y1 = 2.0,
    z0 = -2.0,
    z1 = 2.0,

    # RK options
    correctionOrder = LinearOrder,
    volumeType = RKMassOverDensity,

    # Testing options
    randomizeNodes = True,
    ranfrac = 0.2,
    testHessian = False,
    useOldKernel = False, # Test using old kernel
    useBaseKernel = False, # Test using standard SPH kernel
    printErrors = False,
    quitAfterTiming = False,
    numToCheck = -1, # -1 for all nodes, positive int for a few random nodes
    
    # Manufactured parameters
    funcType = "linear",
    
    # Interpolation kernel choice
    nPerh = 4.01,
    hminmult = 1.e-3,
    hmaxmult = 1.e3,

    # Options for checking
    checkConditions = True,
    
    # Material parameters
    rho0 = 2.5e-7,
    gamma = 5.0/3.0,
    mu = 1.0,

    # Data dir
    dataDirBase = "dumps-RKInterpolation",

    # Error
    tolerance = 1.e-12,
    
    # Plotting
    plot = False,
)

if useOldKernel and useBaseKernel:
    raise ValueError("cannot use both old and base kernel")
if useOldKernel and correctionOrder > QuadraticOrder:
    raise ValueError("correction order must be quadratic to use old kernel")

if nPerh < int(correctionOrder):
    print("nPerh is not large enough for correction order: {} < {}".format(nPerh, int(correctionOrder)))
    
if mpi.procs > 1:
    raise ValueError("parallel node generation not working")
    
#-------------------------------------------------------------------------------
# Choose correct dimension aliases
#-------------------------------------------------------------------------------
exec("from Spheral%id import *" % dimension, globals())

#-------------------------------------------------------------------------------
# Set up data
#-------------------------------------------------------------------------------
# Set limits
lims = [[x0, x1], [y0, y1], [z0, z1]]

# Get point spacing
units = MKS()
delta = [(x1-x0)/nx]
if dimension > 1:
    delta.append((y1 -y0)/ny)
if dimension > 2:
    delta.append((z1 - z0)/nz)
deltaMin = min(delta)
deltaMax = max(delta)
hmin = deltaMin * nPerh * hminmult
hmax = deltaMax * nPerh * hmaxmult

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
dataDir = os.path.join(dataDirBase,
                       "correctionOrder={}".format(correctionOrder),
                       "funcType={}".format(funcType),
                       "dim={}".format(dimension),
                       "nx={}".format(nx))
restartDir = os.path.join(dataDir, "restarts")
vizDir = os.path.join(dataDir, "viz")
restartBaseName = os.path.join(restartDir, "RKInterpolation")

import os, sys
if mpi.rank == 0:
    if not os.path.exists(dataDir):
        os.makedirs(dataDir)
    if not os.path.exists(vizDir):
        os.makedirs(vizDir)
    if not os.path.exists(restartDir):
        os.makedirs(restartDir)
        mpi.barrier()

#-------------------------------------------------------------------------------
# Material properties.
#-------------------------------------------------------------------------------
eos = GammaLawGas(gamma, mu, units)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
WT = TableKernel(WendlandC4Kernel(), 1000)
kernelExtent = WT.kernelExtent
output("WT")

#-------------------------------------------------------------------------------
# Make the NodeList.
#-------------------------------------------------------------------------------
nodes = makeFluidNodeList("nodes", eos, 
                          hmin = hmin,
                          hmax = hmax,
                          nPerh = nPerh,
                          kernelExtent = kernelExtent)
output("nodes")
output("nodes.hmin")
output("nodes.hmax")
output("nodes.nodesPerSmoothingScale")

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
dataBase = DataBase()
dataBase.appendNodeList(nodes)
output("dataBase")
output("dataBase.numNodeLists")
output("dataBase.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Seed the nodes
#-------------------------------------------------------------------------------
if dimension == 1:
    from DistributeNodes import distributeNodesInRange1d
    distributeNodesInRange1d([(nodes, nx, rho0, (x0, x1))],
                             nPerh = nPerh)
elif dimension == 2:
    from GenerateNodeDistribution2d import *
    generator = GenerateNodeDistribution2d(distributionType="lattice",
                                           nRadial = nx,
                                           nTheta = ny,
                                           xmin = (x0, y0),
                                           xmax = (x1, y1),
                                           rho = rho0,
                                           nNodePerh = nPerh)
    if mpi.procs > 1:
        from VoronoiDistributeNodes import distributeNodes2d
    else:
        from DistributeNodes import distributeNodes2d
        distributeNodes2d((nodes, generator))
else:
    from GenerateNodeDistribution3d import *
    generator = GenerateNodeDistribution3d(distributionType="lattice",
                                           n1 = nx,
                                           n2 = ny,
                                           n3 = nz,
                                           xmin = (x0, y0, z0),
                                           xmax = (x1, y1, z1),
                                           rho=rho0,
                                           nNodePerh = nPerh)
    if mpi.procs > 1:
        from VoronoiDistributeNodes import distributeNodes3d
    else:
        from DistributeNodes import distributeNodes3d
        distributeNodes3d((nodes, generator))
        
output("nodes.numNodes")
numLocal = nodes.numInternalNodes
output("numLocal")

#-------------------------------------------------------------------------------
# Randomize nodes
#-------------------------------------------------------------------------------
import random
seed = 459297849234
random.seed(seed)

if randomizeNodes:
    dx = (x1 - x0)/nx
    dy = (y1 - y0)/ny
    dz = (z1 - z0)/nz
    pos = nodes.positions()
    for i in range(nodes.numInternalNodes):
        if dimension == 1:
            pos[i].x += ranfrac * dx * random.uniform(-1.0, 1.0)
        elif dimension == 2:
            pos[i].x += ranfrac * dx * random.uniform(-1.0, 1.0)
            pos[i].y += ranfrac * dy * random.uniform(-1.0, 1.0)
        elif dimension == 3:
            pos[i].x += ranfrac * dx * random.uniform(-1.0, 1.0)
            pos[i].y += ranfrac * dy * random.uniform(-1.0, 1.0)
            pos[i].z += ranfrac * dz * random.uniform(-1.0, 1.0)
            
#-------------------------------------------------------------------------------
# Iterate h
#-------------------------------------------------------------------------------
method = SPHSmoothingScale(IdealH, WT)
iterateIdealH(dataBase,
              vector_of_Physics([method]),
              vector_of_Boundary(),
              100, # max h iterations
              1.e-4) # h tolerance
dataBase.updateConnectivityMap(True)

#-------------------------------------------------------------------------------
# Get interpolant
#-------------------------------------------------------------------------------
if funcType == "constant":
    def func(x):
        return 2.0
    def dfunc(x):
        return [0.0 for v in x]
    def ddfunc(x):
        return [[0.0 for v in x] for v in x]
elif funcType == "linear":
    if dimension == 1:
        def func(x):
            return -10 - 5*x[0]
        def dfunc(x):
            return [-5]
        def ddfunc(x):
            return [[0]]
    elif dimension == 2:
        def func(x):
            return 5 + 8*x[0] - 4*x[1]
        def dfunc(x):
            return [8, -4]
        def ddfunc(x):
            return [[0, 0], [0, 0]]
    else:
        def func(x):
            return -5 + 2*x[0] - 4*x[1] - 9*x[2]
        def dfunc(x):
            return [2, -4, -9]
        def ddfunc(x):
            return [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
elif funcType == "quadratic":
    if dimension == 1:
        def func(x):
            return 9 - 9*x[0] - 9*np.power(x[0],2)
        def dfunc(x):
            return [-9 - 18*x[0]]
        def ddfunc(x):
            return [[-18]]
    elif dimension == 2:
        def func(x):
            return 4 - 7*x[0] - np.power(x[0],2) - 6*x[1] + 9*x[0]*x[1] + 10*np.power(x[1],2)
        def dfunc(x):
            return [-7 - 2*x[0] + 9*x[1], -6 + 9*x[0] + 20*x[1]]
        def ddfunc(x):
            return [[-2, 9], [9, 20]]
    else:
        def func(x):
            return 10 - 6*x[0] - np.power(x[0],2) - 6*x[1] + 4*x[0]*x[1] + np.power(x[1],2) + 8*x[2] + 6*x[0]*x[2] - 7*x[1]*x[2] + 3*np.power(x[2],2)
        def dfunc(x):
            return [-6 - 2*x[0] + 4*x[1] + 6*x[2], -6 + 4*x[0] + 2*x[1] - 7*x[2], 8 + 6*x[0] - 7*x[1] + 6*x[2]]
        def ddfunc(x):
            return [[-2, 4, 6], [4, 2, -7], [6, -7, 6]]
elif funcType == "cubic":
    if dimension == 1:
        def func(x):
            return 9 + 5*x[0] + 8*np.power(x[0],2) - 9*np.power(x[0],3)
        def dfunc(x):
            return [5 + 16*x[0] - 27*np.power(x[0],2)]
        def ddfunc(x):
            return [[16 - 54*x[0]]]
    elif dimension == 2:
        def func(x):
            return 1 + 10*x[0] - 4*np.power(x[0],2) + 8*np.power(x[0],3) + 10*x[1] - 2*x[0]*x[1] + 5*np.power(x[0],2)*x[1] + 6*np.power(x[1],2) + 8*x[0]*np.power(x[1],2) - 8*np.power(x[1],3)
        def dfunc(x):
            return [10 - 8*x[0] + 24*np.power(x[0],2) - 2*x[1] + 10*x[0]*x[1] + 8*np.power(x[1],2), 10 - 2*x[0] + 5*np.power(x[0],2) + 12*x[1] + 16*x[0]*x[1] - 24*np.power(x[1],2)]
        def ddfunc(x):
            return [[-8 + 48*x[0] + 10*x[1], -2 + 10*x[0] + 16*x[1]], [-2 + 10*x[0] + 16*x[1], 12 + 16*x[0] - 48*x[1]]]
    else:
        def func(x):
            return 6 + 8*x[0] + 5*np.power(x[0],2) + 2*np.power(x[0],3) - 5*x[1] - 8*x[0]*x[1] - np.power(x[0],2)*x[1] - 10*np.power(x[1],2) - 9*x[0]*np.power(x[1],2) - 4*np.power(x[1],3) - 8*x[2] + 9*x[0]*x[2] - 10*np.power(x[0],2)*x[2] - 7*x[1]*x[2] - 7*x[0]*x[1]*x[2] - 5*np.power(x[1],2)*x[2] + 2*np.power(x[2],2) + 5*x[0]*np.power(x[2],2) - 4*x[1]*np.power(x[2],2) + 10*np.power(x[2],3)
        def dfunc(x):
            return [8 + 10*x[0] + 6*np.power(x[0],2) - 8*x[1] - 2*x[0]*x[1] - 9*np.power(x[1],2) + 9*x[2] - 20*x[0]*x[2] - 7*x[1]*x[2] + 5*np.power(x[2],2), -5 - 8*x[0] - np.power(x[0],2) - 20*x[1] - 18*x[0]*x[1] - 12*np.power(x[1],2) - 7*x[2] - 7*x[0]*x[2] - 10*x[1]*x[2] - 4*np.power(x[2],2), -8 + 9*x[0] - 10*np.power(x[0],2) - 7*x[1] - 7*x[0]*x[1] - 5*np.power(x[1],2) + 4*x[2] + 10*x[0]*x[2] - 8*x[1]*x[2] + 30*np.power(x[2],2)]
        def ddfunc(x):
            return [[10 + 12*x[0] - 2*x[1] - 20*x[2], -8 - 2*x[0] - 18*x[1] - 7*x[2], 9 - 20*x[0] - 7*x[1] + 10*x[2]], [-8 - 2*x[0] - 18*x[1] - 7*x[2], -20 - 18*x[0] - 24*x[1] - 10*x[2], -7 - 7*x[0] - 10*x[1] - 8*x[2]], [9 - 20*x[0] - 7*x[1] + 10*x[2], -7 - 7*x[0] - 10*x[1] - 8*x[2], 4 + 10*x[0] - 8*x[1] + 60*x[2]]]
elif funcType == "quartic":
    if dimension == 1:
        def func(x):
            return 10 + 7*x[0] + 3*np.power(x[0],2) - 5*np.power(x[0],3) - 9*np.power(x[0],4)
        def dfunc(x):
            return [7 + 6*x[0] - 15*np.power(x[0],2) - 36*np.power(x[0],3)]
        def ddfunc(x):
            return [[6 - 30*x[0] - 108*np.power(x[0],2)]]
    elif dimension == 2:
        def func(x):
            return 4 + 9*x[0] - 2*np.power(x[0],2) + 3*np.power(x[0],3) + 2*np.power(x[0],4) - 9*x[1] + 7*x[0]*x[1] - 2*np.power(x[0],2)*x[1] - 3*np.power(x[0],3)*x[1] - 7*np.power(x[1],2) - 2*x[0]*np.power(x[1],2) - 2*np.power(x[0],2)*np.power(x[1],2) - 2*np.power(x[1],3) + 4*x[0]*np.power(x[1],3) + 4*np.power(x[1],4)
        def dfunc(x):
            return [9 - 4*x[0] + 9*np.power(x[0],2) + 8*np.power(x[0],3) + 7*x[1] - 4*x[0]*x[1] - 9*np.power(x[0],2)*x[1] - 2*np.power(x[1],2) - 4*x[0]*np.power(x[1],2) + 4*np.power(x[1],3), -9 + 7*x[0] - 2*np.power(x[0],2) - 3*np.power(x[0],3) - 14*x[1] - 4*x[0]*x[1] - 4*np.power(x[0],2)*x[1] - 6*np.power(x[1],2) + 12*x[0]*np.power(x[1],2) + 16*np.power(x[1],3)]
        def ddfunc(x):
            return [[-4 + 18*x[0] + 24*np.power(x[0],2) - 4*x[1] - 18*x[0]*x[1] - 4*np.power(x[1],2), 7 - 4*x[0] - 9*np.power(x[0],2) - 4*x[1] - 8*x[0]*x[1] + 12*np.power(x[1],2)], [7 - 4*x[0] - 9*np.power(x[0],2) - 4*x[1] - 8*x[0]*x[1] + 12*np.power(x[1],2), -14 - 4*x[0] - 4*np.power(x[0],2) - 12*x[1] + 24*x[0]*x[1] + 48*np.power(x[1],2)]]
    else:
        def func(x):
            return 9 - x[0] - np.power(x[0],2) + 8*np.power(x[0],3) - 7*np.power(x[0],4) - 7*x[1] - 3*x[0]*x[1] - 3*np.power(x[0],2)*x[1] + 6*np.power(x[0],3)*x[1] + 5*np.power(x[1],2) - 6*x[0]*np.power(x[1],2) - 7*np.power(x[0],2)*np.power(x[1],2) + np.power(x[1],3) - 2*x[0]*np.power(x[1],3) + 4*np.power(x[1],4) - x[2] + 8*x[0]*x[2] - 8*np.power(x[0],2)*x[2] + 4*np.power(x[0],3)*x[2] - 6*x[1]*x[2] + 2*x[0]*x[1]*x[2] - 2*np.power(x[0],2)*x[1]*x[2] + 7*np.power(x[1],2)*x[2] - 5*x[0]*np.power(x[1],2)*x[2] - 9*np.power(x[1],3)*x[2] - np.power(x[2],2) - 4*x[0]*np.power(x[2],2) - 2*np.power(x[0],2)*np.power(x[2],2) + 2*x[1]*np.power(x[2],2) + 2*x[0]*x[1]*np.power(x[2],2) + 8*np.power(x[1],2)*np.power(x[2],2) - 7*np.power(x[2],3) - 10*x[0]*np.power(x[2],3) + 2*x[1]*np.power(x[2],3) + 7*np.power(x[2],4)
        def dfunc(x):
            return [-1 - 2*x[0] + 24*np.power(x[0],2) - 28*np.power(x[0],3) - 3*x[1] - 6*x[0]*x[1] + 18*np.power(x[0],2)*x[1] - 6*np.power(x[1],2) - 14*x[0]*np.power(x[1],2) - 2*np.power(x[1],3) + 8*x[2] - 16*x[0]*x[2] + 12*np.power(x[0],2)*x[2] + 2*x[1]*x[2] - 4*x[0]*x[1]*x[2] - 5*np.power(x[1],2)*x[2] - 4*np.power(x[2],2) - 4*x[0]*np.power(x[2],2) + 2*x[1]*np.power(x[2],2) - 10*np.power(x[2],3), -7 - 3*x[0] - 3*np.power(x[0],2) + 6*np.power(x[0],3) + 10*x[1] - 12*x[0]*x[1] - 14*np.power(x[0],2)*x[1] + 3*np.power(x[1],2) - 6*x[0]*np.power(x[1],2) + 16*np.power(x[1],3) - 6*x[2] + 2*x[0]*x[2] - 2*np.power(x[0],2)*x[2] + 14*x[1]*x[2] - 10*x[0]*x[1]*x[2] - 27*np.power(x[1],2)*x[2] + 2*np.power(x[2],2) + 2*x[0]*np.power(x[2],2) + 16*x[1]*np.power(x[2],2) + 2*np.power(x[2],3), -1 + 8*x[0] - 8*np.power(x[0],2) + 4*np.power(x[0],3) - 6*x[1] + 2*x[0]*x[1] - 2*np.power(x[0],2)*x[1] + 7*np.power(x[1],2) - 5*x[0]*np.power(x[1],2) - 9*np.power(x[1],3) - 2*x[2] - 8*x[0]*x[2] - 4*np.power(x[0],2)*x[2] + 4*x[1]*x[2] + 4*x[0]*x[1]*x[2] + 16*np.power(x[1],2)*x[2] - 21*np.power(x[2],2) - 30*x[0]*np.power(x[2],2) + 6*x[1]*np.power(x[2],2) + 28*np.power(x[2],3)]
        def ddfunc(x):
            return [[-2 + 48*x[0] - 84*np.power(x[0],2) - 6*x[1] + 36*x[0]*x[1] - 14*np.power(x[1],2) - 16*x[2] + 24*x[0]*x[2] - 4*x[1]*x[2] - 4*np.power(x[2],2), -3 - 6*x[0] + 18*np.power(x[0],2) - 12*x[1] - 28*x[0]*x[1] - 6*np.power(x[1],2) + 2*x[2] - 4*x[0]*x[2] - 10*x[1]*x[2] + 2*np.power(x[2],2), 8 - 16*x[0] + 12*np.power(x[0],2) + 2*x[1] - 4*x[0]*x[1] - 5*np.power(x[1],2) - 8*x[2] - 8*x[0]*x[2] + 4*x[1]*x[2] - 30*np.power(x[2],2)], [-3 - 6*x[0] + 18*np.power(x[0],2) - 12*x[1] - 28*x[0]*x[1] - 6*np.power(x[1],2) + 2*x[2] - 4*x[0]*x[2] - 10*x[1]*x[2] + 2*np.power(x[2],2), 10 - 12*x[0] - 14*np.power(x[0],2) + 6*x[1] - 12*x[0]*x[1] + 48*np.power(x[1],2) + 14*x[2] - 10*x[0]*x[2] - 54*x[1]*x[2] + 16*np.power(x[2],2), -6 + 2*x[0] - 2*np.power(x[0],2) + 14*x[1] - 10*x[0]*x[1] - 27*np.power(x[1],2) + 4*x[2] + 4*x[0]*x[2] + 32*x[1]*x[2] + 6*np.power(x[2],2)], [8 - 16*x[0] + 12*np.power(x[0],2) + 2*x[1] - 4*x[0]*x[1] - 5*np.power(x[1],2) - 8*x[2] - 8*x[0]*x[2] + 4*x[1]*x[2] - 30*np.power(x[2],2), -6 + 2*x[0] - 2*np.power(x[0],2) + 14*x[1] - 10*x[0]*x[1] - 27*np.power(x[1],2) + 4*x[2] + 4*x[0]*x[2] + 32*x[1]*x[2] + 6*np.power(x[2],2), -2 - 8*x[0] - 4*np.power(x[0],2) + 4*x[1] + 4*x[0]*x[1] + 16*np.power(x[1],2) - 42*x[2] - 60*x[0]*x[2] + 12*x[1]*x[2] + 84*np.power(x[2],2)]]
elif funcType == "quintic":
    if dimension == 1:
        def func(x):
            return -3 - 10*x[0] + np.power(x[0],2) + 8*np.power(x[0],3) - 6*np.power(x[0],4) + 5*np.power(x[0],5)
        def dfunc(x):
            return [-10 + 2*x[0] + 24*np.power(x[0],2) - 24*np.power(x[0],3) + 25*np.power(x[0],4)]
        def ddfunc(x):
            return [[2 + 48*x[0] - 72*np.power(x[0],2) + 100*np.power(x[0],3)]]
    elif dimension == 2:
        def func(x):
            return 4 + 6*x[0] - 7*np.power(x[0],2) - 2*np.power(x[0],3) + np.power(x[0],4) + 4*np.power(x[0],5) + 2*x[1] - 7*x[0]*x[1] + 5*np.power(x[0],2)*x[1] - np.power(x[0],3)*x[1] - 6*np.power(x[0],4)*x[1] + 10*np.power(x[1],2) - 6*x[0]*np.power(x[1],2) - 4*np.power(x[0],2)*np.power(x[1],2) + np.power(x[0],3)*np.power(x[1],2) - 8*np.power(x[1],3) - 2*x[0]*np.power(x[1],3) - 10*np.power(x[0],2)*np.power(x[1],3) + 2*np.power(x[1],4) - 7*x[0]*np.power(x[1],4) - 8*np.power(x[1],5)
        def dfunc(x):
            return [6 - 14*x[0] - 6*np.power(x[0],2) + 4*np.power(x[0],3) + 20*np.power(x[0],4) - 7*x[1] + 10*x[0]*x[1] - 3*np.power(x[0],2)*x[1] - 24*np.power(x[0],3)*x[1] - 6*np.power(x[1],2) - 8*x[0]*np.power(x[1],2) + 3*np.power(x[0],2)*np.power(x[1],2) - 2*np.power(x[1],3) - 20*x[0]*np.power(x[1],3) - 7*np.power(x[1],4), 2 - 7*x[0] + 5*np.power(x[0],2) - np.power(x[0],3) - 6*np.power(x[0],4) + 20*x[1] - 12*x[0]*x[1] - 8*np.power(x[0],2)*x[1] + 2*np.power(x[0],3)*x[1] - 24*np.power(x[1],2) - 6*x[0]*np.power(x[1],2) - 30*np.power(x[0],2)*np.power(x[1],2) + 8*np.power(x[1],3) - 28*x[0]*np.power(x[1],3) - 40*np.power(x[1],4)]
        def ddfunc(x):
            return [[-14 - 12*x[0] + 12*np.power(x[0],2) + 80*np.power(x[0],3) + 10*x[1] - 6*x[0]*x[1] - 72*np.power(x[0],2)*x[1] - 8*np.power(x[1],2) + 6*x[0]*np.power(x[1],2) - 20*np.power(x[1],3), -7 + 10*x[0] - 3*np.power(x[0],2) - 24*np.power(x[0],3) - 12*x[1] - 16*x[0]*x[1] + 6*np.power(x[0],2)*x[1] - 6*np.power(x[1],2) - 60*x[0]*np.power(x[1],2) - 28*np.power(x[1],3)], [-7 + 10*x[0] - 3*np.power(x[0],2) - 24*np.power(x[0],3) - 12*x[1] - 16*x[0]*x[1] + 6*np.power(x[0],2)*x[1] - 6*np.power(x[1],2) - 60*x[0]*np.power(x[1],2) - 28*np.power(x[1],3), 20 - 12*x[0] - 8*np.power(x[0],2) + 2*np.power(x[0],3) - 48*x[1] - 12*x[0]*x[1] - 60*np.power(x[0],2)*x[1] + 24*np.power(x[1],2) - 84*x[0]*np.power(x[1],2) - 160*np.power(x[1],3)]]
    else:
        def func(x):
            return -8 + 10*x[0] + 4*np.power(x[0],2) + 4*np.power(x[0],3) + 6*np.power(x[0],4) + 9*np.power(x[0],5) + 9*x[1] - 10*x[0]*x[1] + 2*np.power(x[0],2)*x[1] + np.power(x[0],3)*x[1] + 6*np.power(x[0],4)*x[1] - 8*np.power(x[1],2) - 4*x[0]*np.power(x[1],2) + np.power(x[0],2)*np.power(x[1],2) + 2*np.power(x[0],3)*np.power(x[1],2) + 10*np.power(x[1],3) - 10*x[0]*np.power(x[1],3) - 7*np.power(x[0],2)*np.power(x[1],3) - 10*np.power(x[1],4) + 10*x[0]*np.power(x[1],4) - 2*np.power(x[1],5) + 10*x[2] - 4*x[0]*x[2] - 5*np.power(x[0],2)*x[2] - 5*np.power(x[0],3)*x[2] + 2*np.power(x[0],4)*x[2] + 8*x[1]*x[2] - 6*x[0]*x[1]*x[2] + np.power(x[0],2)*x[1]*x[2] + 5*np.power(x[0],3)*x[1]*x[2] - 7*np.power(x[1],2)*x[2] - x[0]*np.power(x[1],2)*x[2] - 8*np.power(x[0],2)*np.power(x[1],2)*x[2] - 5*np.power(x[1],3)*x[2] + 8*x[0]*np.power(x[1],3)*x[2] + 10*np.power(x[1],4)*x[2] - 10*np.power(x[2],2) - 3*x[0]*np.power(x[2],2) + 6*np.power(x[0],2)*np.power(x[2],2) - 7*np.power(x[0],3)*np.power(x[2],2) - 2*x[1]*np.power(x[2],2) - 6*x[0]*x[1]*np.power(x[2],2) + 8*np.power(x[0],2)*x[1]*np.power(x[2],2) + 6*np.power(x[1],2)*np.power(x[2],2) - 7*x[0]*np.power(x[1],2)*np.power(x[2],2) + 2*np.power(x[1],3)*np.power(x[2],2) - 5*np.power(x[2],3) - 4*x[0]*np.power(x[2],3) - 7*np.power(x[0],2)*np.power(x[2],3) - 4*x[1]*np.power(x[2],3) + x[0]*x[1]*np.power(x[2],3) - 6*np.power(x[1],2)*np.power(x[2],3) - 3*np.power(x[2],4) - 10*x[0]*np.power(x[2],4) + 7*x[1]*np.power(x[2],4) - 8*np.power(x[2],5)
        def dfunc(x):
            return [10 + 8*x[0] + 12*np.power(x[0],2) + 24*np.power(x[0],3) + 45*np.power(x[0],4) - 10*x[1] + 4*x[0]*x[1] + 3*np.power(x[0],2)*x[1] + 24*np.power(x[0],3)*x[1] - 4*np.power(x[1],2) + 2*x[0]*np.power(x[1],2) + 6*np.power(x[0],2)*np.power(x[1],2) - 10*np.power(x[1],3) - 14*x[0]*np.power(x[1],3) + 10*np.power(x[1],4) - 4*x[2] - 10*x[0]*x[2] - 15*np.power(x[0],2)*x[2] + 8*np.power(x[0],3)*x[2] - 6*x[1]*x[2] + 2*x[0]*x[1]*x[2] + 15*np.power(x[0],2)*x[1]*x[2] - np.power(x[1],2)*x[2] - 16*x[0]*np.power(x[1],2)*x[2] + 8*np.power(x[1],3)*x[2] - 3*np.power(x[2],2) + 12*x[0]*np.power(x[2],2) - 21*np.power(x[0],2)*np.power(x[2],2) - 6*x[1]*np.power(x[2],2) + 16*x[0]*x[1]*np.power(x[2],2) - 7*np.power(x[1],2)*np.power(x[2],2) - 4*np.power(x[2],3) - 14*x[0]*np.power(x[2],3) + x[1]*np.power(x[2],3) - 10*np.power(x[2],4), 9 - 10*x[0] + 2*np.power(x[0],2) + np.power(x[0],3) + 6*np.power(x[0],4) - 16*x[1] - 8*x[0]*x[1] + 2*np.power(x[0],2)*x[1] + 4*np.power(x[0],3)*x[1] + 30*np.power(x[1],2) - 30*x[0]*np.power(x[1],2) - 21*np.power(x[0],2)*np.power(x[1],2) - 40*np.power(x[1],3) + 40*x[0]*np.power(x[1],3) - 10*np.power(x[1],4) + 8*x[2] - 6*x[0]*x[2] + np.power(x[0],2)*x[2] + 5*np.power(x[0],3)*x[2] - 14*x[1]*x[2] - 2*x[0]*x[1]*x[2] - 16*np.power(x[0],2)*x[1]*x[2] - 15*np.power(x[1],2)*x[2] + 24*x[0]*np.power(x[1],2)*x[2] + 40*np.power(x[1],3)*x[2] - 2*np.power(x[2],2) - 6*x[0]*np.power(x[2],2) + 8*np.power(x[0],2)*np.power(x[2],2) + 12*x[1]*np.power(x[2],2) - 14*x[0]*x[1]*np.power(x[2],2) + 6*np.power(x[1],2)*np.power(x[2],2) - 4*np.power(x[2],3) + x[0]*np.power(x[2],3) - 12*x[1]*np.power(x[2],3) + 7*np.power(x[2],4), 10 - 4*x[0] - 5*np.power(x[0],2) - 5*np.power(x[0],3) + 2*np.power(x[0],4) + 8*x[1] - 6*x[0]*x[1] + np.power(x[0],2)*x[1] + 5*np.power(x[0],3)*x[1] - 7*np.power(x[1],2) - x[0]*np.power(x[1],2) - 8*np.power(x[0],2)*np.power(x[1],2) - 5*np.power(x[1],3) + 8*x[0]*np.power(x[1],3) + 10*np.power(x[1],4) - 20*x[2] - 6*x[0]*x[2] + 12*np.power(x[0],2)*x[2] - 14*np.power(x[0],3)*x[2] - 4*x[1]*x[2] - 12*x[0]*x[1]*x[2] + 16*np.power(x[0],2)*x[1]*x[2] + 12*np.power(x[1],2)*x[2] - 14*x[0]*np.power(x[1],2)*x[2] + 4*np.power(x[1],3)*x[2] - 15*np.power(x[2],2) - 12*x[0]*np.power(x[2],2) - 21*np.power(x[0],2)*np.power(x[2],2) - 12*x[1]*np.power(x[2],2) + 3*x[0]*x[1]*np.power(x[2],2) - 18*np.power(x[1],2)*np.power(x[2],2) - 12*np.power(x[2],3) - 40*x[0]*np.power(x[2],3) + 28*x[1]*np.power(x[2],3) - 40*np.power(x[2],4)]
        def ddfunc(x):
            return [[8 + 24*x[0] + 72*np.power(x[0],2) + 180*np.power(x[0],3) + 4*x[1] + 6*x[0]*x[1] + 72*np.power(x[0],2)*x[1] + 2*np.power(x[1],2) + 12*x[0]*np.power(x[1],2) - 14*np.power(x[1],3) - 10*x[2] - 30*x[0]*x[2] + 24*np.power(x[0],2)*x[2] + 2*x[1]*x[2] + 30*x[0]*x[1]*x[2] - 16*np.power(x[1],2)*x[2] + 12*np.power(x[2],2) - 42*x[0]*np.power(x[2],2) + 16*x[1]*np.power(x[2],2) - 14*np.power(x[2],3), -10 + 4*x[0] + 3*np.power(x[0],2) + 24*np.power(x[0],3) - 8*x[1] + 4*x[0]*x[1] + 12*np.power(x[0],2)*x[1] - 30*np.power(x[1],2) - 42*x[0]*np.power(x[1],2) + 40*np.power(x[1],3) - 6*x[2] + 2*x[0]*x[2] + 15*np.power(x[0],2)*x[2] - 2*x[1]*x[2] - 32*x[0]*x[1]*x[2] + 24*np.power(x[1],2)*x[2] - 6*np.power(x[2],2) + 16*x[0]*np.power(x[2],2) - 14*x[1]*np.power(x[2],2) + np.power(x[2],3), -4 - 10*x[0] - 15*np.power(x[0],2) + 8*np.power(x[0],3) - 6*x[1] + 2*x[0]*x[1] + 15*np.power(x[0],2)*x[1] - np.power(x[1],2) - 16*x[0]*np.power(x[1],2) + 8*np.power(x[1],3) - 6*x[2] + 24*x[0]*x[2] - 42*np.power(x[0],2)*x[2] - 12*x[1]*x[2] + 32*x[0]*x[1]*x[2] - 14*np.power(x[1],2)*x[2] - 12*np.power(x[2],2) - 42*x[0]*np.power(x[2],2) + 3*x[1]*np.power(x[2],2) - 40*np.power(x[2],3)], [-10 + 4*x[0] + 3*np.power(x[0],2) + 24*np.power(x[0],3) - 8*x[1] + 4*x[0]*x[1] + 12*np.power(x[0],2)*x[1] - 30*np.power(x[1],2) - 42*x[0]*np.power(x[1],2) + 40*np.power(x[1],3) - 6*x[2] + 2*x[0]*x[2] + 15*np.power(x[0],2)*x[2] - 2*x[1]*x[2] - 32*x[0]*x[1]*x[2] + 24*np.power(x[1],2)*x[2] - 6*np.power(x[2],2) + 16*x[0]*np.power(x[2],2) - 14*x[1]*np.power(x[2],2) + np.power(x[2],3), -16 - 8*x[0] + 2*np.power(x[0],2) + 4*np.power(x[0],3) + 60*x[1] - 60*x[0]*x[1] - 42*np.power(x[0],2)*x[1] - 120*np.power(x[1],2) + 120*x[0]*np.power(x[1],2) - 40*np.power(x[1],3) - 14*x[2] - 2*x[0]*x[2] - 16*np.power(x[0],2)*x[2] - 30*x[1]*x[2] + 48*x[0]*x[1]*x[2] + 120*np.power(x[1],2)*x[2] + 12*np.power(x[2],2) - 14*x[0]*np.power(x[2],2) + 12*x[1]*np.power(x[2],2) - 12*np.power(x[2],3), 8 - 6*x[0] + np.power(x[0],2) + 5*np.power(x[0],3) - 14*x[1] - 2*x[0]*x[1] - 16*np.power(x[0],2)*x[1] - 15*np.power(x[1],2) + 24*x[0]*np.power(x[1],2) + 40*np.power(x[1],3) - 4*x[2] - 12*x[0]*x[2] + 16*np.power(x[0],2)*x[2] + 24*x[1]*x[2] - 28*x[0]*x[1]*x[2] + 12*np.power(x[1],2)*x[2] - 12*np.power(x[2],2) + 3*x[0]*np.power(x[2],2) - 36*x[1]*np.power(x[2],2) + 28*np.power(x[2],3)], [-4 - 10*x[0] - 15*np.power(x[0],2) + 8*np.power(x[0],3) - 6*x[1] + 2*x[0]*x[1] + 15*np.power(x[0],2)*x[1] - np.power(x[1],2) - 16*x[0]*np.power(x[1],2) + 8*np.power(x[1],3) - 6*x[2] + 24*x[0]*x[2] - 42*np.power(x[0],2)*x[2] - 12*x[1]*x[2] + 32*x[0]*x[1]*x[2] - 14*np.power(x[1],2)*x[2] - 12*np.power(x[2],2) - 42*x[0]*np.power(x[2],2) + 3*x[1]*np.power(x[2],2) - 40*np.power(x[2],3), 8 - 6*x[0] + np.power(x[0],2) + 5*np.power(x[0],3) - 14*x[1] - 2*x[0]*x[1] - 16*np.power(x[0],2)*x[1] - 15*np.power(x[1],2) + 24*x[0]*np.power(x[1],2) + 40*np.power(x[1],3) - 4*x[2] - 12*x[0]*x[2] + 16*np.power(x[0],2)*x[2] + 24*x[1]*x[2] - 28*x[0]*x[1]*x[2] + 12*np.power(x[1],2)*x[2] - 12*np.power(x[2],2) + 3*x[0]*np.power(x[2],2) - 36*x[1]*np.power(x[2],2) + 28*np.power(x[2],3), -20 - 6*x[0] + 12*np.power(x[0],2) - 14*np.power(x[0],3) - 4*x[1] - 12*x[0]*x[1] + 16*np.power(x[0],2)*x[1] + 12*np.power(x[1],2) - 14*x[0]*np.power(x[1],2) + 4*np.power(x[1],3) - 30*x[2] - 24*x[0]*x[2] - 42*np.power(x[0],2)*x[2] - 24*x[1]*x[2] + 6*x[0]*x[1]*x[2] - 36*np.power(x[1],2)*x[2] - 36*np.power(x[2],2) - 120*x[0]*np.power(x[2],2) + 84*x[1]*np.power(x[2],2) - 160*np.power(x[2],3)]]
elif funcType == "sextic":
    if dimension == 1:
        def func(x):
            return 5 - 3*x[0] - 8*(x[0]*x[0]) - 9*(x[0]*x[0]*x[0]) - 7*(x[0]*x[0]*x[0]*x[0]) + 2*(x[0]*x[0]*x[0]*x[0]*x[0]) + 9*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0])
        def dfunc(x):
            return [-3 - 16*x[0] - 27*(x[0]*x[0]) - 28*(x[0]*x[0]*x[0]) + 10*(x[0]*x[0]*x[0]*x[0]) + 54*(x[0]*x[0]*x[0]*x[0]*x[0])]
        def ddfunc(x):
            return [[-16 - 54*x[0] - 84*(x[0]*x[0]) + 40*(x[0]*x[0]*x[0]) + 270*(x[0]*x[0]*x[0]*x[0])]]
    elif dimension == 2:
        def func(x):
            return -1 + x[0] - x[0]*x[0] - 6*(x[0]*x[0]*x[0]) + 3*(x[0]*x[0]*x[0]*x[0]) - 6*(x[0]*x[0]*x[0]*x[0]*x[0]) + 9*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0]) - 3*x[1] + 9*x[0]*x[1] - 4*(x[0]*x[0])*x[1] - 6*(x[0]*x[0]*x[0])*x[1] - 4*(x[0]*x[0]*x[0]*x[0])*x[1] - 9*(x[0]*x[0]*x[0]*x[0]*x[0])*x[1] + 9*(x[1]*x[1]) + 6*x[0]*(x[1]*x[1]) + 6*(x[0]*x[0])*(x[1]*x[1]) + 5*(x[0]*x[0]*x[0])*(x[1]*x[1]) + 10*(x[0]*x[0]*x[0]*x[0])*(x[1]*x[1]) + 9*(x[1]*x[1]*x[1]) - 7*x[0]*(x[1]*x[1]*x[1]) + 4*(x[0]*x[0])*(x[1]*x[1]*x[1]) - 3*(x[0]*x[0]*x[0])*(x[1]*x[1]*x[1]) + 9*(x[1]*x[1]*x[1]*x[1]) + 4*x[0]*(x[1]*x[1]*x[1]*x[1]) + 7*(x[0]*x[0])*(x[1]*x[1]*x[1]*x[1]) - 2*(x[1]*x[1]*x[1]*x[1]*x[1]) - 4*x[0]*(x[1]*x[1]*x[1]*x[1]*x[1]) + x[1]*x[1]*x[1]*x[1]*x[1]*x[1]
        def dfunc(x):
            return [1 - 2*x[0] - 18*(x[0]*x[0]) + 12*(x[0]*x[0]*x[0]) - 30*(x[0]*x[0]*x[0]*x[0]) + 54*(x[0]*x[0]*x[0]*x[0]*x[0]) + 9*x[1] - 8*x[0]*x[1] - 18*(x[0]*x[0])*x[1] - 16*(x[0]*x[0]*x[0])*x[1] - 45*(x[0]*x[0]*x[0]*x[0])*x[1] + 6*(x[1]*x[1]) + 12*x[0]*(x[1]*x[1]) + 15*(x[0]*x[0])*(x[1]*x[1]) + 40*(x[0]*x[0]*x[0])*(x[1]*x[1]) - 7*(x[1]*x[1]*x[1]) + 8*x[0]*(x[1]*x[1]*x[1]) - 9*(x[0]*x[0])*(x[1]*x[1]*x[1]) + 4*(x[1]*x[1]*x[1]*x[1]) + 14*x[0]*(x[1]*x[1]*x[1]*x[1]) - 4*(x[1]*x[1]*x[1]*x[1]*x[1]), -3 + 9*x[0] - 4*(x[0]*x[0]) - 6*(x[0]*x[0]*x[0]) - 4*(x[0]*x[0]*x[0]*x[0]) - 9*(x[0]*x[0]*x[0]*x[0]*x[0]) + 18*x[1] + 12*x[0]*x[1] + 12*(x[0]*x[0])*x[1] + 10*(x[0]*x[0]*x[0])*x[1] + 20*(x[0]*x[0]*x[0]*x[0])*x[1] + 27*(x[1]*x[1]) - 21*x[0]*(x[1]*x[1]) + 12*(x[0]*x[0])*(x[1]*x[1]) - 9*(x[0]*x[0]*x[0])*(x[1]*x[1]) + 36*(x[1]*x[1]*x[1]) + 16*x[0]*(x[1]*x[1]*x[1]) + 28*(x[0]*x[0])*(x[1]*x[1]*x[1]) - 10*(x[1]*x[1]*x[1]*x[1]) - 20*x[0]*(x[1]*x[1]*x[1]*x[1]) + 6*(x[1]*x[1]*x[1]*x[1]*x[1])]
        def ddfunc(x):
            return [[-2 - 36*x[0] + 36*(x[0]*x[0]) - 120*(x[0]*x[0]*x[0]) + 270*(x[0]*x[0]*x[0]*x[0]) - 8*x[1] - 36*x[0]*x[1] - 48*(x[0]*x[0])*x[1] - 180*(x[0]*x[0]*x[0])*x[1] + 12*(x[1]*x[1]) + 30*x[0]*(x[1]*x[1]) + 120*(x[0]*x[0])*(x[1]*x[1]) + 8*(x[1]*x[1]*x[1]) - 18*x[0]*(x[1]*x[1]*x[1]) + 14*(x[1]*x[1]*x[1]*x[1]), 9 - 8*x[0] - 18*(x[0]*x[0]) - 16*(x[0]*x[0]*x[0]) - 45*(x[0]*x[0]*x[0]*x[0]) + 12*x[1] + 24*x[0]*x[1] + 30*(x[0]*x[0])*x[1] + 80*(x[0]*x[0]*x[0])*x[1] - 21*(x[1]*x[1]) + 24*x[0]*(x[1]*x[1]) - 27*(x[0]*x[0])*(x[1]*x[1]) + 16*(x[1]*x[1]*x[1]) + 56*x[0]*(x[1]*x[1]*x[1]) - 20*(x[1]*x[1]*x[1]*x[1])], [9 - 8*x[0] - 18*(x[0]*x[0]) - 16*(x[0]*x[0]*x[0]) - 45*(x[0]*x[0]*x[0]*x[0]) + 12*x[1] + 24*x[0]*x[1] + 30*(x[0]*x[0])*x[1] + 80*(x[0]*x[0]*x[0])*x[1] - 21*(x[1]*x[1]) + 24*x[0]*(x[1]*x[1]) - 27*(x[0]*x[0])*(x[1]*x[1]) + 16*(x[1]*x[1]*x[1]) + 56*x[0]*(x[1]*x[1]*x[1]) - 20*(x[1]*x[1]*x[1]*x[1]), 18 + 12*x[0] + 12*(x[0]*x[0]) + 10*(x[0]*x[0]*x[0]) + 20*(x[0]*x[0]*x[0]*x[0]) + 54*x[1] - 42*x[0]*x[1] + 24*(x[0]*x[0])*x[1] - 18*(x[0]*x[0]*x[0])*x[1] + 108*(x[1]*x[1]) + 48*x[0]*(x[1]*x[1]) + 84*(x[0]*x[0])*(x[1]*x[1]) - 40*(x[1]*x[1]*x[1]) - 80*x[0]*(x[1]*x[1]*x[1]) + 30*(x[1]*x[1]*x[1]*x[1])]]
    else:
        def func(x):
            return 7 - 7*x[0] - 2*(x[0]*x[0]) - 3*(x[0]*x[0]*x[0]) - x[0]*x[0]*x[0]*x[0] - 10*(x[0]*x[0]*x[0]*x[0]*x[0]) - 7*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0]) + 7*x[1] - 5*x[0]*x[1] - 4*(x[0]*x[0])*x[1] - x[0]*x[0]*x[0]*x[1] + 4*(x[0]*x[0]*x[0]*x[0])*x[1] - 9*(x[0]*x[0]*x[0]*x[0]*x[0])*x[1] + 4*(x[1]*x[1]) - 7*x[0]*(x[1]*x[1]) + 9*(x[0]*x[0])*(x[1]*x[1]) + x[0]*x[0]*x[0]*(x[1]*x[1]) + 7*(x[0]*x[0]*x[0]*x[0])*(x[1]*x[1]) - 10*(x[1]*x[1]*x[1]) - 9*x[0]*(x[1]*x[1]*x[1]) - 4*(x[0]*x[0])*(x[1]*x[1]*x[1]) - 4*(x[0]*x[0]*x[0])*(x[1]*x[1]*x[1]) + 7*(x[1]*x[1]*x[1]*x[1]) - 4*x[0]*(x[1]*x[1]*x[1]*x[1]) - 9*(x[0]*x[0])*(x[1]*x[1]*x[1]*x[1]) + 6*(x[1]*x[1]*x[1]*x[1]*x[1]) - 8*x[0]*(x[1]*x[1]*x[1]*x[1]*x[1]) + 7*(x[1]*x[1]*x[1]*x[1]*x[1]*x[1]) + 3*x[2] + 8*x[0]*x[2] + x[0]*x[0]*x[2] + 2*(x[0]*x[0]*x[0])*x[2] - 3*(x[0]*x[0]*x[0]*x[0])*x[2] + 10*(x[0]*x[0]*x[0]*x[0]*x[0])*x[2] - 3*x[1]*x[2] - 2*x[0]*x[1]*x[2] + 3*(x[0]*x[0])*x[1]*x[2] + 9*(x[0]*x[0]*x[0])*x[1]*x[2] + 5*(x[0]*x[0]*x[0]*x[0])*x[1]*x[2] + 7*(x[1]*x[1])*x[2] + 4*x[0]*(x[1]*x[1])*x[2] - 3*(x[0]*x[0])*(x[1]*x[1])*x[2] + 6*(x[0]*x[0]*x[0])*(x[1]*x[1])*x[2] + 10*(x[1]*x[1]*x[1])*x[2] - 6*x[0]*(x[1]*x[1]*x[1])*x[2] + x[0]*x[0]*(x[1]*x[1]*x[1])*x[2] - 8*(x[1]*x[1]*x[1]*x[1])*x[2] + 10*x[0]*(x[1]*x[1]*x[1]*x[1])*x[2] + 10*(x[1]*x[1]*x[1]*x[1]*x[1])*x[2] + x[2]*x[2] - 9*x[0]*(x[2]*x[2]) - x[0]*x[0]*(x[2]*x[2]) - 7*(x[0]*x[0]*x[0])*(x[2]*x[2]) + 8*(x[0]*x[0]*x[0]*x[0])*(x[2]*x[2]) - 2*x[1]*(x[2]*x[2]) + 7*x[0]*x[1]*(x[2]*x[2]) - 2*(x[0]*x[0])*x[1]*(x[2]*x[2]) + 5*(x[0]*x[0]*x[0])*x[1]*(x[2]*x[2]) - x[1]*x[1]*(x[2]*x[2]) - 7*x[0]*(x[1]*x[1])*(x[2]*x[2]) - 9*(x[0]*x[0])*(x[1]*x[1])*(x[2]*x[2]) - x[1]*x[1]*x[1]*(x[2]*x[2]) - 10*x[0]*(x[1]*x[1]*x[1])*(x[2]*x[2]) + 3*(x[1]*x[1]*x[1]*x[1])*(x[2]*x[2]) - x[2]*x[2]*x[2] + 9*x[0]*(x[2]*x[2]*x[2]) - 10*(x[0]*x[0])*(x[2]*x[2]*x[2]) + 8*(x[0]*x[0]*x[0])*(x[2]*x[2]*x[2]) - x[1]*(x[2]*x[2]*x[2]) + 3*x[0]*x[1]*(x[2]*x[2]*x[2]) - 8*(x[0]*x[0])*x[1]*(x[2]*x[2]*x[2]) - 3*(x[1]*x[1])*(x[2]*x[2]*x[2]) - 7*x[0]*(x[1]*x[1])*(x[2]*x[2]*x[2]) + 8*(x[1]*x[1]*x[1])*(x[2]*x[2]*x[2]) + 4*(x[2]*x[2]*x[2]*x[2]) - 3*x[0]*(x[2]*x[2]*x[2]*x[2]) - 7*(x[0]*x[0])*(x[2]*x[2]*x[2]*x[2]) + 6*x[1]*(x[2]*x[2]*x[2]*x[2]) - 4*x[0]*x[1]*(x[2]*x[2]*x[2]*x[2]) - 4*(x[1]*x[1])*(x[2]*x[2]*x[2]*x[2]) + 9*(x[2]*x[2]*x[2]*x[2]*x[2]) - 10*x[0]*(x[2]*x[2]*x[2]*x[2]*x[2]) - x[1]*(x[2]*x[2]*x[2]*x[2]*x[2]) - 5*(x[2]*x[2]*x[2]*x[2]*x[2]*x[2])
        def dfunc(x):
            return [-7 - 4*x[0] - 9*(x[0]*x[0]) - 4*(x[0]*x[0]*x[0]) - 50*(x[0]*x[0]*x[0]*x[0]) - 42*(x[0]*x[0]*x[0]*x[0]*x[0]) - 5*x[1] - 8*x[0]*x[1] - 3*(x[0]*x[0])*x[1] + 16*(x[0]*x[0]*x[0])*x[1] - 45*(x[0]*x[0]*x[0]*x[0])*x[1] - 7*(x[1]*x[1]) + 18*x[0]*(x[1]*x[1]) + 3*(x[0]*x[0])*(x[1]*x[1]) + 28*(x[0]*x[0]*x[0])*(x[1]*x[1]) - 9*(x[1]*x[1]*x[1]) - 8*x[0]*(x[1]*x[1]*x[1]) - 12*(x[0]*x[0])*(x[1]*x[1]*x[1]) - 4*(x[1]*x[1]*x[1]*x[1]) - 18*x[0]*(x[1]*x[1]*x[1]*x[1]) - 8*(x[1]*x[1]*x[1]*x[1]*x[1]) + 8*x[2] + 2*x[0]*x[2] + 6*(x[0]*x[0])*x[2] - 12*(x[0]*x[0]*x[0])*x[2] + 50*(x[0]*x[0]*x[0]*x[0])*x[2] - 2*x[1]*x[2] + 6*x[0]*x[1]*x[2] + 27*(x[0]*x[0])*x[1]*x[2] + 20*(x[0]*x[0]*x[0])*x[1]*x[2] + 4*(x[1]*x[1])*x[2] - 6*x[0]*(x[1]*x[1])*x[2] + 18*(x[0]*x[0])*(x[1]*x[1])*x[2] - 6*(x[1]*x[1]*x[1])*x[2] + 2*x[0]*(x[1]*x[1]*x[1])*x[2] + 10*(x[1]*x[1]*x[1]*x[1])*x[2] - 9*(x[2]*x[2]) - 2*x[0]*(x[2]*x[2]) - 21*(x[0]*x[0])*(x[2]*x[2]) + 32*(x[0]*x[0]*x[0])*(x[2]*x[2]) + 7*x[1]*(x[2]*x[2]) - 4*x[0]*x[1]*(x[2]*x[2]) + 15*(x[0]*x[0])*x[1]*(x[2]*x[2]) - 7*(x[1]*x[1])*(x[2]*x[2]) - 18*x[0]*(x[1]*x[1])*(x[2]*x[2]) - 10*(x[1]*x[1]*x[1])*(x[2]*x[2]) + 9*(x[2]*x[2]*x[2]) - 20*x[0]*(x[2]*x[2]*x[2]) + 24*(x[0]*x[0])*(x[2]*x[2]*x[2]) + 3*x[1]*(x[2]*x[2]*x[2]) - 16*x[0]*x[1]*(x[2]*x[2]*x[2]) - 7*(x[1]*x[1])*(x[2]*x[2]*x[2]) - 3*(x[2]*x[2]*x[2]*x[2]) - 14*x[0]*(x[2]*x[2]*x[2]*x[2]) - 4*x[1]*(x[2]*x[2]*x[2]*x[2]) - 10*(x[2]*x[2]*x[2]*x[2]*x[2]), 7 - 5*x[0] - 4*(x[0]*x[0]) - x[0]*x[0]*x[0] + 4*(x[0]*x[0]*x[0]*x[0]) - 9*(x[0]*x[0]*x[0]*x[0]*x[0]) + 8*x[1] - 14*x[0]*x[1] + 18*(x[0]*x[0])*x[1] + 2*(x[0]*x[0]*x[0])*x[1] + 14*(x[0]*x[0]*x[0]*x[0])*x[1] - 30*(x[1]*x[1]) - 27*x[0]*(x[1]*x[1]) - 12*(x[0]*x[0])*(x[1]*x[1]) - 12*(x[0]*x[0]*x[0])*(x[1]*x[1]) + 28*(x[1]*x[1]*x[1]) - 16*x[0]*(x[1]*x[1]*x[1]) - 36*(x[0]*x[0])*(x[1]*x[1]*x[1]) + 30*(x[1]*x[1]*x[1]*x[1]) - 40*x[0]*(x[1]*x[1]*x[1]*x[1]) + 42*(x[1]*x[1]*x[1]*x[1]*x[1]) - 3*x[2] - 2*x[0]*x[2] + 3*(x[0]*x[0])*x[2] + 9*(x[0]*x[0]*x[0])*x[2] + 5*(x[0]*x[0]*x[0]*x[0])*x[2] + 14*x[1]*x[2] + 8*x[0]*x[1]*x[2] - 6*(x[0]*x[0])*x[1]*x[2] + 12*(x[0]*x[0]*x[0])*x[1]*x[2] + 30*(x[1]*x[1])*x[2] - 18*x[0]*(x[1]*x[1])*x[2] + 3*(x[0]*x[0])*(x[1]*x[1])*x[2] - 32*(x[1]*x[1]*x[1])*x[2] + 40*x[0]*(x[1]*x[1]*x[1])*x[2] + 50*(x[1]*x[1]*x[1]*x[1])*x[2] - 2*(x[2]*x[2]) + 7*x[0]*(x[2]*x[2]) - 2*(x[0]*x[0])*(x[2]*x[2]) + 5*(x[0]*x[0]*x[0])*(x[2]*x[2]) - 2*x[1]*(x[2]*x[2]) - 14*x[0]*x[1]*(x[2]*x[2]) - 18*(x[0]*x[0])*x[1]*(x[2]*x[2]) - 3*(x[1]*x[1])*(x[2]*x[2]) - 30*x[0]*(x[1]*x[1])*(x[2]*x[2]) + 12*(x[1]*x[1]*x[1])*(x[2]*x[2]) - x[2]*x[2]*x[2] + 3*x[0]*(x[2]*x[2]*x[2]) - 8*(x[0]*x[0])*(x[2]*x[2]*x[2]) - 6*x[1]*(x[2]*x[2]*x[2]) - 14*x[0]*x[1]*(x[2]*x[2]*x[2]) + 24*(x[1]*x[1])*(x[2]*x[2]*x[2]) + 6*(x[2]*x[2]*x[2]*x[2]) - 4*x[0]*(x[2]*x[2]*x[2]*x[2]) - 8*x[1]*(x[2]*x[2]*x[2]*x[2]) - x[2]*x[2]*x[2]*x[2]*x[2], 3 + 8*x[0] + x[0]*x[0] + 2*(x[0]*x[0]*x[0]) - 3*(x[0]*x[0]*x[0]*x[0]) + 10*(x[0]*x[0]*x[0]*x[0]*x[0]) - 3*x[1] - 2*x[0]*x[1] + 3*(x[0]*x[0])*x[1] + 9*(x[0]*x[0]*x[0])*x[1] + 5*(x[0]*x[0]*x[0]*x[0])*x[1] + 7*(x[1]*x[1]) + 4*x[0]*(x[1]*x[1]) - 3*(x[0]*x[0])*(x[1]*x[1]) + 6*(x[0]*x[0]*x[0])*(x[1]*x[1]) + 10*(x[1]*x[1]*x[1]) - 6*x[0]*(x[1]*x[1]*x[1]) + x[0]*x[0]*(x[1]*x[1]*x[1]) - 8*(x[1]*x[1]*x[1]*x[1]) + 10*x[0]*(x[1]*x[1]*x[1]*x[1]) + 10*(x[1]*x[1]*x[1]*x[1]*x[1]) + 2*x[2] - 18*x[0]*x[2] - 2*(x[0]*x[0])*x[2] - 14*(x[0]*x[0]*x[0])*x[2] + 16*(x[0]*x[0]*x[0]*x[0])*x[2] - 4*x[1]*x[2] + 14*x[0]*x[1]*x[2] - 4*(x[0]*x[0])*x[1]*x[2] + 10*(x[0]*x[0]*x[0])*x[1]*x[2] - 2*(x[1]*x[1])*x[2] - 14*x[0]*(x[1]*x[1])*x[2] - 18*(x[0]*x[0])*(x[1]*x[1])*x[2] - 2*(x[1]*x[1]*x[1])*x[2] - 20*x[0]*(x[1]*x[1]*x[1])*x[2] + 6*(x[1]*x[1]*x[1]*x[1])*x[2] - 3*(x[2]*x[2]) + 27*x[0]*(x[2]*x[2]) - 30*(x[0]*x[0])*(x[2]*x[2]) + 24*(x[0]*x[0]*x[0])*(x[2]*x[2]) - 3*x[1]*(x[2]*x[2]) + 9*x[0]*x[1]*(x[2]*x[2]) - 24*(x[0]*x[0])*x[1]*(x[2]*x[2]) - 9*(x[1]*x[1])*(x[2]*x[2]) - 21*x[0]*(x[1]*x[1])*(x[2]*x[2]) + 24*(x[1]*x[1]*x[1])*(x[2]*x[2]) + 16*(x[2]*x[2]*x[2]) - 12*x[0]*(x[2]*x[2]*x[2]) - 28*(x[0]*x[0])*(x[2]*x[2]*x[2]) + 24*x[1]*(x[2]*x[2]*x[2]) - 16*x[0]*x[1]*(x[2]*x[2]*x[2]) - 16*(x[1]*x[1])*(x[2]*x[2]*x[2]) + 45*(x[2]*x[2]*x[2]*x[2]) - 50*x[0]*(x[2]*x[2]*x[2]*x[2]) - 5*x[1]*(x[2]*x[2]*x[2]*x[2]) - 30*(x[2]*x[2]*x[2]*x[2]*x[2])]
        def ddfunc(x):
            return [[-4 - 18*x[0] - 12*(x[0]*x[0]) - 200*(x[0]*x[0]*x[0]) - 210*(x[0]*x[0]*x[0]*x[0]) - 8*x[1] - 6*x[0]*x[1] + 48*(x[0]*x[0])*x[1] - 180*(x[0]*x[0]*x[0])*x[1] + 18*(x[1]*x[1]) + 6*x[0]*(x[1]*x[1]) + 84*(x[0]*x[0])*(x[1]*x[1]) - 8*(x[1]*x[1]*x[1]) - 24*x[0]*(x[1]*x[1]*x[1]) - 18*(x[1]*x[1]*x[1]*x[1]) + 2*x[2] + 12*x[0]*x[2] - 36*(x[0]*x[0])*x[2] + 200*(x[0]*x[0]*x[0])*x[2] + 6*x[1]*x[2] + 54*x[0]*x[1]*x[2] + 60*(x[0]*x[0])*x[1]*x[2] - 6*(x[1]*x[1])*x[2] + 36*x[0]*(x[1]*x[1])*x[2] + 2*(x[1]*x[1]*x[1])*x[2] - 2*(x[2]*x[2]) - 42*x[0]*(x[2]*x[2]) + 96*(x[0]*x[0])*(x[2]*x[2]) - 4*x[1]*(x[2]*x[2]) + 30*x[0]*x[1]*(x[2]*x[2]) - 18*(x[1]*x[1])*(x[2]*x[2]) - 20*(x[2]*x[2]*x[2]) + 48*x[0]*(x[2]*x[2]*x[2]) - 16*x[1]*(x[2]*x[2]*x[2]) - 14*(x[2]*x[2]*x[2]*x[2]), -5 - 8*x[0] - 3*(x[0]*x[0]) + 16*(x[0]*x[0]*x[0]) - 45*(x[0]*x[0]*x[0]*x[0]) - 14*x[1] + 36*x[0]*x[1] + 6*(x[0]*x[0])*x[1] + 56*(x[0]*x[0]*x[0])*x[1] - 27*(x[1]*x[1]) - 24*x[0]*(x[1]*x[1]) - 36*(x[0]*x[0])*(x[1]*x[1]) - 16*(x[1]*x[1]*x[1]) - 72*x[0]*(x[1]*x[1]*x[1]) - 40*(x[1]*x[1]*x[1]*x[1]) - 2*x[2] + 6*x[0]*x[2] + 27*(x[0]*x[0])*x[2] + 20*(x[0]*x[0]*x[0])*x[2] + 8*x[1]*x[2] - 12*x[0]*x[1]*x[2] + 36*(x[0]*x[0])*x[1]*x[2] - 18*(x[1]*x[1])*x[2] + 6*x[0]*(x[1]*x[1])*x[2] + 40*(x[1]*x[1]*x[1])*x[2] + 7*(x[2]*x[2]) - 4*x[0]*(x[2]*x[2]) + 15*(x[0]*x[0])*(x[2]*x[2]) - 14*x[1]*(x[2]*x[2]) - 36*x[0]*x[1]*(x[2]*x[2]) - 30*(x[1]*x[1])*(x[2]*x[2]) + 3*(x[2]*x[2]*x[2]) - 16*x[0]*(x[2]*x[2]*x[2]) - 14*x[1]*(x[2]*x[2]*x[2]) - 4*(x[2]*x[2]*x[2]*x[2]), 8 + 2*x[0] + 6*(x[0]*x[0]) - 12*(x[0]*x[0]*x[0]) + 50*(x[0]*x[0]*x[0]*x[0]) - 2*x[1] + 6*x[0]*x[1] + 27*(x[0]*x[0])*x[1] + 20*(x[0]*x[0]*x[0])*x[1] + 4*(x[1]*x[1]) - 6*x[0]*(x[1]*x[1]) + 18*(x[0]*x[0])*(x[1]*x[1]) - 6*(x[1]*x[1]*x[1]) + 2*x[0]*(x[1]*x[1]*x[1]) + 10*(x[1]*x[1]*x[1]*x[1]) - 18*x[2] - 4*x[0]*x[2] - 42*(x[0]*x[0])*x[2] + 64*(x[0]*x[0]*x[0])*x[2] + 14*x[1]*x[2] - 8*x[0]*x[1]*x[2] + 30*(x[0]*x[0])*x[1]*x[2] - 14*(x[1]*x[1])*x[2] - 36*x[0]*(x[1]*x[1])*x[2] - 20*(x[1]*x[1]*x[1])*x[2] + 27*(x[2]*x[2]) - 60*x[0]*(x[2]*x[2]) + 72*(x[0]*x[0])*(x[2]*x[2]) + 9*x[1]*(x[2]*x[2]) - 48*x[0]*x[1]*(x[2]*x[2]) - 21*(x[1]*x[1])*(x[2]*x[2]) - 12*(x[2]*x[2]*x[2]) - 56*x[0]*(x[2]*x[2]*x[2]) - 16*x[1]*(x[2]*x[2]*x[2]) - 50*(x[2]*x[2]*x[2]*x[2])], [-5 - 8*x[0] - 3*(x[0]*x[0]) + 16*(x[0]*x[0]*x[0]) - 45*(x[0]*x[0]*x[0]*x[0]) - 14*x[1] + 36*x[0]*x[1] + 6*(x[0]*x[0])*x[1] + 56*(x[0]*x[0]*x[0])*x[1] - 27*(x[1]*x[1]) - 24*x[0]*(x[1]*x[1]) - 36*(x[0]*x[0])*(x[1]*x[1]) - 16*(x[1]*x[1]*x[1]) - 72*x[0]*(x[1]*x[1]*x[1]) - 40*(x[1]*x[1]*x[1]*x[1]) - 2*x[2] + 6*x[0]*x[2] + 27*(x[0]*x[0])*x[2] + 20*(x[0]*x[0]*x[0])*x[2] + 8*x[1]*x[2] - 12*x[0]*x[1]*x[2] + 36*(x[0]*x[0])*x[1]*x[2] - 18*(x[1]*x[1])*x[2] + 6*x[0]*(x[1]*x[1])*x[2] + 40*(x[1]*x[1]*x[1])*x[2] + 7*(x[2]*x[2]) - 4*x[0]*(x[2]*x[2]) + 15*(x[0]*x[0])*(x[2]*x[2]) - 14*x[1]*(x[2]*x[2]) - 36*x[0]*x[1]*(x[2]*x[2]) - 30*(x[1]*x[1])*(x[2]*x[2]) + 3*(x[2]*x[2]*x[2]) - 16*x[0]*(x[2]*x[2]*x[2]) - 14*x[1]*(x[2]*x[2]*x[2]) - 4*(x[2]*x[2]*x[2]*x[2]), 8 - 14*x[0] + 18*(x[0]*x[0]) + 2*(x[0]*x[0]*x[0]) + 14*(x[0]*x[0]*x[0]*x[0]) - 60*x[1] - 54*x[0]*x[1] - 24*(x[0]*x[0])*x[1] - 24*(x[0]*x[0]*x[0])*x[1] + 84*(x[1]*x[1]) - 48*x[0]*(x[1]*x[1]) - 108*(x[0]*x[0])*(x[1]*x[1]) + 120*(x[1]*x[1]*x[1]) - 160*x[0]*(x[1]*x[1]*x[1]) + 210*(x[1]*x[1]*x[1]*x[1]) + 14*x[2] + 8*x[0]*x[2] - 6*(x[0]*x[0])*x[2] + 12*(x[0]*x[0]*x[0])*x[2] + 60*x[1]*x[2] - 36*x[0]*x[1]*x[2] + 6*(x[0]*x[0])*x[1]*x[2] - 96*(x[1]*x[1])*x[2] + 120*x[0]*(x[1]*x[1])*x[2] + 200*(x[1]*x[1]*x[1])*x[2] - 2*(x[2]*x[2]) - 14*x[0]*(x[2]*x[2]) - 18*(x[0]*x[0])*(x[2]*x[2]) - 6*x[1]*(x[2]*x[2]) - 60*x[0]*x[1]*(x[2]*x[2]) + 36*(x[1]*x[1])*(x[2]*x[2]) - 6*(x[2]*x[2]*x[2]) - 14*x[0]*(x[2]*x[2]*x[2]) + 48*x[1]*(x[2]*x[2]*x[2]) - 8*(x[2]*x[2]*x[2]*x[2]), -3 - 2*x[0] + 3*(x[0]*x[0]) + 9*(x[0]*x[0]*x[0]) + 5*(x[0]*x[0]*x[0]*x[0]) + 14*x[1] + 8*x[0]*x[1] - 6*(x[0]*x[0])*x[1] + 12*(x[0]*x[0]*x[0])*x[1] + 30*(x[1]*x[1]) - 18*x[0]*(x[1]*x[1]) + 3*(x[0]*x[0])*(x[1]*x[1]) - 32*(x[1]*x[1]*x[1]) + 40*x[0]*(x[1]*x[1]*x[1]) + 50*(x[1]*x[1]*x[1]*x[1]) - 4*x[2] + 14*x[0]*x[2] - 4*(x[0]*x[0])*x[2] + 10*(x[0]*x[0]*x[0])*x[2] - 4*x[1]*x[2] - 28*x[0]*x[1]*x[2] - 36*(x[0]*x[0])*x[1]*x[2] - 6*(x[1]*x[1])*x[2] - 60*x[0]*(x[1]*x[1])*x[2] + 24*(x[1]*x[1]*x[1])*x[2] - 3*(x[2]*x[2]) + 9*x[0]*(x[2]*x[2]) - 24*(x[0]*x[0])*(x[2]*x[2]) - 18*x[1]*(x[2]*x[2]) - 42*x[0]*x[1]*(x[2]*x[2]) + 72*(x[1]*x[1])*(x[2]*x[2]) + 24*(x[2]*x[2]*x[2]) - 16*x[0]*(x[2]*x[2]*x[2]) - 32*x[1]*(x[2]*x[2]*x[2]) - 5*(x[2]*x[2]*x[2]*x[2])], [8 + 2*x[0] + 6*(x[0]*x[0]) - 12*(x[0]*x[0]*x[0]) + 50*(x[0]*x[0]*x[0]*x[0]) - 2*x[1] + 6*x[0]*x[1] + 27*(x[0]*x[0])*x[1] + 20*(x[0]*x[0]*x[0])*x[1] + 4*(x[1]*x[1]) - 6*x[0]*(x[1]*x[1]) + 18*(x[0]*x[0])*(x[1]*x[1]) - 6*(x[1]*x[1]*x[1]) + 2*x[0]*(x[1]*x[1]*x[1]) + 10*(x[1]*x[1]*x[1]*x[1]) - 18*x[2] - 4*x[0]*x[2] - 42*(x[0]*x[0])*x[2] + 64*(x[0]*x[0]*x[0])*x[2] + 14*x[1]*x[2] - 8*x[0]*x[1]*x[2] + 30*(x[0]*x[0])*x[1]*x[2] - 14*(x[1]*x[1])*x[2] - 36*x[0]*(x[1]*x[1])*x[2] - 20*(x[1]*x[1]*x[1])*x[2] + 27*(x[2]*x[2]) - 60*x[0]*(x[2]*x[2]) + 72*(x[0]*x[0])*(x[2]*x[2]) + 9*x[1]*(x[2]*x[2]) - 48*x[0]*x[1]*(x[2]*x[2]) - 21*(x[1]*x[1])*(x[2]*x[2]) - 12*(x[2]*x[2]*x[2]) - 56*x[0]*(x[2]*x[2]*x[2]) - 16*x[1]*(x[2]*x[2]*x[2]) - 50*(x[2]*x[2]*x[2]*x[2]), -3 - 2*x[0] + 3*(x[0]*x[0]) + 9*(x[0]*x[0]*x[0]) + 5*(x[0]*x[0]*x[0]*x[0]) + 14*x[1] + 8*x[0]*x[1] - 6*(x[0]*x[0])*x[1] + 12*(x[0]*x[0]*x[0])*x[1] + 30*(x[1]*x[1]) - 18*x[0]*(x[1]*x[1]) + 3*(x[0]*x[0])*(x[1]*x[1]) - 32*(x[1]*x[1]*x[1]) + 40*x[0]*(x[1]*x[1]*x[1]) + 50*(x[1]*x[1]*x[1]*x[1]) - 4*x[2] + 14*x[0]*x[2] - 4*(x[0]*x[0])*x[2] + 10*(x[0]*x[0]*x[0])*x[2] - 4*x[1]*x[2] - 28*x[0]*x[1]*x[2] - 36*(x[0]*x[0])*x[1]*x[2] - 6*(x[1]*x[1])*x[2] - 60*x[0]*(x[1]*x[1])*x[2] + 24*(x[1]*x[1]*x[1])*x[2] - 3*(x[2]*x[2]) + 9*x[0]*(x[2]*x[2]) - 24*(x[0]*x[0])*(x[2]*x[2]) - 18*x[1]*(x[2]*x[2]) - 42*x[0]*x[1]*(x[2]*x[2]) + 72*(x[1]*x[1])*(x[2]*x[2]) + 24*(x[2]*x[2]*x[2]) - 16*x[0]*(x[2]*x[2]*x[2]) - 32*x[1]*(x[2]*x[2]*x[2]) - 5*(x[2]*x[2]*x[2]*x[2]), 2 - 18*x[0] - 2*(x[0]*x[0]) - 14*(x[0]*x[0]*x[0]) + 16*(x[0]*x[0]*x[0]*x[0]) - 4*x[1] + 14*x[0]*x[1] - 4*(x[0]*x[0])*x[1] + 10*(x[0]*x[0]*x[0])*x[1] - 2*(x[1]*x[1]) - 14*x[0]*(x[1]*x[1]) - 18*(x[0]*x[0])*(x[1]*x[1]) - 2*(x[1]*x[1]*x[1]) - 20*x[0]*(x[1]*x[1]*x[1]) + 6*(x[1]*x[1]*x[1]*x[1]) - 6*x[2] + 54*x[0]*x[2] - 60*(x[0]*x[0])*x[2] + 48*(x[0]*x[0]*x[0])*x[2] - 6*x[1]*x[2] + 18*x[0]*x[1]*x[2] - 48*(x[0]*x[0])*x[1]*x[2] - 18*(x[1]*x[1])*x[2] - 42*x[0]*(x[1]*x[1])*x[2] + 48*(x[1]*x[1]*x[1])*x[2] + 48*(x[2]*x[2]) - 36*x[0]*(x[2]*x[2]) - 84*(x[0]*x[0])*(x[2]*x[2]) + 72*x[1]*(x[2]*x[2]) - 48*x[0]*x[1]*(x[2]*x[2]) - 48*(x[1]*x[1])*(x[2]*x[2]) + 180*(x[2]*x[2]*x[2]) - 200*x[0]*(x[2]*x[2]*x[2]) - 20*x[1]*(x[2]*x[2]*x[2]) - 150*(x[2]*x[2]*x[2]*x[2])]]
elif funcType == "septic":
    if dimension == 1:
        def func(x):
            return 3 - 7*x[0] + 8*(x[0]*x[0]) + 7*(x[0]*x[0]*x[0]) - 5*(x[0]*x[0]*x[0]*x[0]) - 8*(x[0]*x[0]*x[0]*x[0]*x[0]) - 9*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0]) - 6*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0])
        def dfunc(x):
            return [-7 + 16*x[0] + 21*(x[0]*x[0]) - 20*(x[0]*x[0]*x[0]) - 40*(x[0]*x[0]*x[0]*x[0]) - 54*(x[0]*x[0]*x[0]*x[0]*x[0]) - 42*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0])]
        def ddfunc(x):
            return [[16 + 42*x[0] - 60*(x[0]*x[0]) - 160*(x[0]*x[0]*x[0]) - 270*(x[0]*x[0]*x[0]*x[0]) - 252*(x[0]*x[0]*x[0]*x[0]*x[0])]]
    elif dimension == 2:
        def func(x):
            return -9 - 10*x[0] - 3*(x[0]*x[0]) - 10*(x[0]*x[0]*x[0]) - 6*(x[0]*x[0]*x[0]*x[0]) - 6*(x[0]*x[0]*x[0]*x[0]*x[0]) + 2*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0]) + 8*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]) - 6*x[1] + 8*x[0]*x[1] - 4*(x[0]*x[0])*x[1] - 9*(x[0]*x[0]*x[0])*x[1] - x[0]*x[0]*x[0]*x[0]*x[1] + 4*(x[0]*x[0]*x[0]*x[0]*x[0])*x[1] - 4*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0])*x[1] - 4*(x[1]*x[1]) + 10*x[0]*(x[1]*x[1]) - 8*(x[0]*x[0])*(x[1]*x[1]) + 8*(x[0]*x[0]*x[0])*(x[1]*x[1]) + 9*(x[0]*x[0]*x[0]*x[0])*(x[1]*x[1]) + 8*(x[0]*x[0]*x[0]*x[0]*x[0])*(x[1]*x[1]) + 8*(x[1]*x[1]*x[1]) - 4*x[0]*(x[1]*x[1]*x[1]) + 9*(x[0]*x[0])*(x[1]*x[1]*x[1]) + 6*(x[0]*x[0]*x[0])*(x[1]*x[1]*x[1]) + 6*(x[0]*x[0]*x[0]*x[0])*(x[1]*x[1]*x[1]) + 8*(x[1]*x[1]*x[1]*x[1]) - 4*x[0]*(x[1]*x[1]*x[1]*x[1]) + 6*(x[0]*x[0])*(x[1]*x[1]*x[1]*x[1]) + 3*(x[0]*x[0]*x[0])*(x[1]*x[1]*x[1]*x[1]) + 9*(x[1]*x[1]*x[1]*x[1]*x[1]) + 9*x[0]*(x[1]*x[1]*x[1]*x[1]*x[1]) + 5*(x[0]*x[0])*(x[1]*x[1]*x[1]*x[1]*x[1]) + 2*(x[1]*x[1]*x[1]*x[1]*x[1]*x[1]) - 5*x[0]*(x[1]*x[1]*x[1]*x[1]*x[1]*x[1]) - 10*(x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1])
        def dfunc(x):
            return [-10 - 6*x[0] - 30*(x[0]*x[0]) - 24*(x[0]*x[0]*x[0]) - 30*(x[0]*x[0]*x[0]*x[0]) + 12*(x[0]*x[0]*x[0]*x[0]*x[0]) + 56*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0]) + 8*x[1] - 8*x[0]*x[1] - 27*(x[0]*x[0])*x[1] - 4*(x[0]*x[0]*x[0])*x[1] + 20*(x[0]*x[0]*x[0]*x[0])*x[1] - 24*(x[0]*x[0]*x[0]*x[0]*x[0])*x[1] + 10*(x[1]*x[1]) - 16*x[0]*(x[1]*x[1]) + 24*(x[0]*x[0])*(x[1]*x[1]) + 36*(x[0]*x[0]*x[0])*(x[1]*x[1]) + 40*(x[0]*x[0]*x[0]*x[0])*(x[1]*x[1]) - 4*(x[1]*x[1]*x[1]) + 18*x[0]*(x[1]*x[1]*x[1]) + 18*(x[0]*x[0])*(x[1]*x[1]*x[1]) + 24*(x[0]*x[0]*x[0])*(x[1]*x[1]*x[1]) - 4*(x[1]*x[1]*x[1]*x[1]) + 12*x[0]*(x[1]*x[1]*x[1]*x[1]) + 9*(x[0]*x[0])*(x[1]*x[1]*x[1]*x[1]) + 9*(x[1]*x[1]*x[1]*x[1]*x[1]) + 10*x[0]*(x[1]*x[1]*x[1]*x[1]*x[1]) - 5*(x[1]*x[1]*x[1]*x[1]*x[1]*x[1]), -6 + 8*x[0] - 4*(x[0]*x[0]) - 9*(x[0]*x[0]*x[0]) - x[0]*x[0]*x[0]*x[0] + 4*(x[0]*x[0]*x[0]*x[0]*x[0]) - 4*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0]) - 8*x[1] + 20*x[0]*x[1] - 16*(x[0]*x[0])*x[1] + 16*(x[0]*x[0]*x[0])*x[1] + 18*(x[0]*x[0]*x[0]*x[0])*x[1] + 16*(x[0]*x[0]*x[0]*x[0]*x[0])*x[1] + 24*(x[1]*x[1]) - 12*x[0]*(x[1]*x[1]) + 27*(x[0]*x[0])*(x[1]*x[1]) + 18*(x[0]*x[0]*x[0])*(x[1]*x[1]) + 18*(x[0]*x[0]*x[0]*x[0])*(x[1]*x[1]) + 32*(x[1]*x[1]*x[1]) - 16*x[0]*(x[1]*x[1]*x[1]) + 24*(x[0]*x[0])*(x[1]*x[1]*x[1]) + 12*(x[0]*x[0]*x[0])*(x[1]*x[1]*x[1]) + 45*(x[1]*x[1]*x[1]*x[1]) + 45*x[0]*(x[1]*x[1]*x[1]*x[1]) + 25*(x[0]*x[0])*(x[1]*x[1]*x[1]*x[1]) + 12*(x[1]*x[1]*x[1]*x[1]*x[1]) - 30*x[0]*(x[1]*x[1]*x[1]*x[1]*x[1]) - 70*(x[1]*x[1]*x[1]*x[1]*x[1]*x[1])]
        def ddfunc(x):
            return [[-6 - 60*x[0] - 72*(x[0]*x[0]) - 120*(x[0]*x[0]*x[0]) + 60*(x[0]*x[0]*x[0]*x[0]) + 336*(x[0]*x[0]*x[0]*x[0]*x[0]) - 8*x[1] - 54*x[0]*x[1] - 12*(x[0]*x[0])*x[1] + 80*(x[0]*x[0]*x[0])*x[1] - 120*(x[0]*x[0]*x[0]*x[0])*x[1] - 16*(x[1]*x[1]) + 48*x[0]*(x[1]*x[1]) + 108*(x[0]*x[0])*(x[1]*x[1]) + 160*(x[0]*x[0]*x[0])*(x[1]*x[1]) + 18*(x[1]*x[1]*x[1]) + 36*x[0]*(x[1]*x[1]*x[1]) + 72*(x[0]*x[0])*(x[1]*x[1]*x[1]) + 12*(x[1]*x[1]*x[1]*x[1]) + 18*x[0]*(x[1]*x[1]*x[1]*x[1]) + 10*(x[1]*x[1]*x[1]*x[1]*x[1]), 8 - 8*x[0] - 27*(x[0]*x[0]) - 4*(x[0]*x[0]*x[0]) + 20*(x[0]*x[0]*x[0]*x[0]) - 24*(x[0]*x[0]*x[0]*x[0]*x[0]) + 20*x[1] - 32*x[0]*x[1] + 48*(x[0]*x[0])*x[1] + 72*(x[0]*x[0]*x[0])*x[1] + 80*(x[0]*x[0]*x[0]*x[0])*x[1] - 12*(x[1]*x[1]) + 54*x[0]*(x[1]*x[1]) + 54*(x[0]*x[0])*(x[1]*x[1]) + 72*(x[0]*x[0]*x[0])*(x[1]*x[1]) - 16*(x[1]*x[1]*x[1]) + 48*x[0]*(x[1]*x[1]*x[1]) + 36*(x[0]*x[0])*(x[1]*x[1]*x[1]) + 45*(x[1]*x[1]*x[1]*x[1]) + 50*x[0]*(x[1]*x[1]*x[1]*x[1]) - 30*(x[1]*x[1]*x[1]*x[1]*x[1])], [8 - 8*x[0] - 27*(x[0]*x[0]) - 4*(x[0]*x[0]*x[0]) + 20*(x[0]*x[0]*x[0]*x[0]) - 24*(x[0]*x[0]*x[0]*x[0]*x[0]) + 20*x[1] - 32*x[0]*x[1] + 48*(x[0]*x[0])*x[1] + 72*(x[0]*x[0]*x[0])*x[1] + 80*(x[0]*x[0]*x[0]*x[0])*x[1] - 12*(x[1]*x[1]) + 54*x[0]*(x[1]*x[1]) + 54*(x[0]*x[0])*(x[1]*x[1]) + 72*(x[0]*x[0]*x[0])*(x[1]*x[1]) - 16*(x[1]*x[1]*x[1]) + 48*x[0]*(x[1]*x[1]*x[1]) + 36*(x[0]*x[0])*(x[1]*x[1]*x[1]) + 45*(x[1]*x[1]*x[1]*x[1]) + 50*x[0]*(x[1]*x[1]*x[1]*x[1]) - 30*(x[1]*x[1]*x[1]*x[1]*x[1]), -8 + 20*x[0] - 16*(x[0]*x[0]) + 16*(x[0]*x[0]*x[0]) + 18*(x[0]*x[0]*x[0]*x[0]) + 16*(x[0]*x[0]*x[0]*x[0]*x[0]) + 48*x[1] - 24*x[0]*x[1] + 54*(x[0]*x[0])*x[1] + 36*(x[0]*x[0]*x[0])*x[1] + 36*(x[0]*x[0]*x[0]*x[0])*x[1] + 96*(x[1]*x[1]) - 48*x[0]*(x[1]*x[1]) + 72*(x[0]*x[0])*(x[1]*x[1]) + 36*(x[0]*x[0]*x[0])*(x[1]*x[1]) + 180*(x[1]*x[1]*x[1]) + 180*x[0]*(x[1]*x[1]*x[1]) + 100*(x[0]*x[0])*(x[1]*x[1]*x[1]) + 60*(x[1]*x[1]*x[1]*x[1]) - 150*x[0]*(x[1]*x[1]*x[1]*x[1]) - 420*(x[1]*x[1]*x[1]*x[1]*x[1])]]
    else:
        def func(x):
            return -10 - 4*x[0] + 10*(x[0]*x[0]) + 7*(x[0]*x[0]*x[0]) - 5*(x[0]*x[0]*x[0]*x[0]) - x[0]*x[0]*x[0]*x[0]*x[0] + 8*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0]) + 8*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]) - 9*x[1] + 9*x[0]*x[1] - x[0]*x[0]*x[1] - 7*(x[0]*x[0]*x[0])*x[1] + 5*(x[0]*x[0]*x[0]*x[0])*x[1] + 5*(x[0]*x[0]*x[0]*x[0]*x[0])*x[1] + 9*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0])*x[1] + 8*(x[1]*x[1]) - 10*x[0]*(x[1]*x[1]) - 10*(x[0]*x[0])*(x[1]*x[1]) - 7*(x[0]*x[0]*x[0])*(x[1]*x[1]) + 5*(x[0]*x[0]*x[0]*x[0])*(x[1]*x[1]) - 10*(x[0]*x[0]*x[0]*x[0]*x[0])*(x[1]*x[1]) - 8*(x[1]*x[1]*x[1]) - 9*x[0]*(x[1]*x[1]*x[1]) - 3*(x[0]*x[0])*(x[1]*x[1]*x[1]) - 3*(x[0]*x[0]*x[0])*(x[1]*x[1]*x[1]) + 10*(x[0]*x[0]*x[0]*x[0])*(x[1]*x[1]*x[1]) + x[1]*x[1]*x[1]*x[1] + x[0]*(x[1]*x[1]*x[1]*x[1]) - 8*(x[0]*x[0])*(x[1]*x[1]*x[1]*x[1]) + 9*(x[0]*x[0]*x[0])*(x[1]*x[1]*x[1]*x[1]) + 10*(x[1]*x[1]*x[1]*x[1]*x[1]) + 6*x[0]*(x[1]*x[1]*x[1]*x[1]*x[1]) - 2*(x[0]*x[0])*(x[1]*x[1]*x[1]*x[1]*x[1]) + x[1]*x[1]*x[1]*x[1]*x[1]*x[1] - 6*x[0]*(x[1]*x[1]*x[1]*x[1]*x[1]*x[1]) - 8*(x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]) + 3*x[2] - 3*x[0]*x[2] - 10*(x[0]*x[0])*x[2] - 6*(x[0]*x[0]*x[0])*x[2] + 10*(x[0]*x[0]*x[0]*x[0])*x[2] + 9*(x[0]*x[0]*x[0]*x[0]*x[0])*x[2] + 2*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0])*x[2] - 3*x[1]*x[2] + 7*x[0]*x[1]*x[2] + 7*(x[0]*x[0])*x[1]*x[2] - 4*(x[0]*x[0]*x[0])*x[1]*x[2] + 4*(x[0]*x[0]*x[0]*x[0])*x[1]*x[2] - 8*(x[0]*x[0]*x[0]*x[0]*x[0])*x[1]*x[2] + x[1]*x[1]*x[2] + 8*x[0]*(x[1]*x[1])*x[2] + 2*(x[0]*x[0])*(x[1]*x[1])*x[2] + 10*(x[0]*x[0]*x[0])*(x[1]*x[1])*x[2] - 8*(x[0]*x[0]*x[0]*x[0])*(x[1]*x[1])*x[2] + 2*(x[1]*x[1]*x[1])*x[2] - 8*x[0]*(x[1]*x[1]*x[1])*x[2] - x[0]*x[0]*(x[1]*x[1]*x[1])*x[2] + 7*(x[0]*x[0]*x[0])*(x[1]*x[1]*x[1])*x[2] + 5*(x[1]*x[1]*x[1]*x[1])*x[2] - 4*x[0]*(x[1]*x[1]*x[1]*x[1])*x[2] + 9*(x[0]*x[0])*(x[1]*x[1]*x[1]*x[1])*x[2] - 2*(x[1]*x[1]*x[1]*x[1]*x[1])*x[2] + x[0]*(x[1]*x[1]*x[1]*x[1]*x[1])*x[2] - 9*(x[1]*x[1]*x[1]*x[1]*x[1]*x[1])*x[2] + 10*(x[2]*x[2]) - 2*x[0]*(x[2]*x[2]) - 10*(x[0]*x[0])*(x[2]*x[2]) - 4*(x[0]*x[0]*x[0])*(x[2]*x[2]) + 6*(x[0]*x[0]*x[0]*x[0])*(x[2]*x[2]) + 10*(x[0]*x[0]*x[0]*x[0]*x[0])*(x[2]*x[2]) - 2*x[1]*(x[2]*x[2]) - 8*x[0]*x[1]*(x[2]*x[2]) - x[0]*x[0]*x[1]*(x[2]*x[2]) - 2*(x[0]*x[0]*x[0])*x[1]*(x[2]*x[2]) - 4*(x[0]*x[0]*x[0]*x[0])*x[1]*(x[2]*x[2]) + 8*(x[1]*x[1])*(x[2]*x[2]) + 2*x[0]*(x[1]*x[1])*(x[2]*x[2]) - 10*(x[0]*x[0])*(x[1]*x[1])*(x[2]*x[2]) + 4*(x[0]*x[0]*x[0])*(x[1]*x[1])*(x[2]*x[2]) - 4*(x[1]*x[1]*x[1])*(x[2]*x[2]) - 10*x[0]*(x[1]*x[1]*x[1])*(x[2]*x[2]) + 8*(x[0]*x[0])*(x[1]*x[1]*x[1])*(x[2]*x[2]) - 4*(x[1]*x[1]*x[1]*x[1])*(x[2]*x[2]) + 9*x[0]*(x[1]*x[1]*x[1]*x[1])*(x[2]*x[2]) - 10*(x[1]*x[1]*x[1]*x[1]*x[1])*(x[2]*x[2]) + 3*(x[2]*x[2]*x[2]) - 5*x[0]*(x[2]*x[2]*x[2]) - x[0]*x[0]*(x[2]*x[2]*x[2]) - 9*(x[0]*x[0]*x[0])*(x[2]*x[2]*x[2]) - x[0]*x[0]*x[0]*x[0]*(x[2]*x[2]*x[2]) + 3*x[1]*(x[2]*x[2]*x[2]) + 2*x[0]*x[1]*(x[2]*x[2]*x[2]) - 8*(x[0]*x[0])*x[1]*(x[2]*x[2]*x[2]) + 2*(x[0]*x[0]*x[0])*x[1]*(x[2]*x[2]*x[2]) - 2*(x[1]*x[1])*(x[2]*x[2]*x[2]) + 10*x[0]*(x[1]*x[1])*(x[2]*x[2]*x[2]) - x[0]*x[0]*(x[1]*x[1])*(x[2]*x[2]*x[2]) - 6*(x[1]*x[1]*x[1])*(x[2]*x[2]*x[2]) + 6*x[0]*(x[1]*x[1]*x[1])*(x[2]*x[2]*x[2]) - 6*(x[1]*x[1]*x[1]*x[1])*(x[2]*x[2]*x[2]) + 7*(x[2]*x[2]*x[2]*x[2]) - 2*x[0]*(x[2]*x[2]*x[2]*x[2]) - 8*(x[0]*x[0])*(x[2]*x[2]*x[2]*x[2]) - 8*(x[0]*x[0]*x[0])*(x[2]*x[2]*x[2]*x[2]) + x[1]*(x[2]*x[2]*x[2]*x[2]) + x[0]*x[1]*(x[2]*x[2]*x[2]*x[2]) - x[0]*x[0]*x[1]*(x[2]*x[2]*x[2]*x[2]) + 8*(x[1]*x[1])*(x[2]*x[2]*x[2]*x[2]) - 9*x[0]*(x[1]*x[1])*(x[2]*x[2]*x[2]*x[2]) + 6*(x[1]*x[1]*x[1])*(x[2]*x[2]*x[2]*x[2]) + 6*(x[2]*x[2]*x[2]*x[2]*x[2]) + 10*x[0]*(x[2]*x[2]*x[2]*x[2]*x[2]) + 5*(x[0]*x[0])*(x[2]*x[2]*x[2]*x[2]*x[2]) - x[1]*(x[2]*x[2]*x[2]*x[2]*x[2]) + 5*x[0]*x[1]*(x[2]*x[2]*x[2]*x[2]*x[2]) - x[1]*x[1]*(x[2]*x[2]*x[2]*x[2]*x[2]) + 6*(x[2]*x[2]*x[2]*x[2]*x[2]*x[2]) - 2*x[0]*(x[2]*x[2]*x[2]*x[2]*x[2]*x[2]) + 7*x[1]*(x[2]*x[2]*x[2]*x[2]*x[2]*x[2]) - 2*(x[2]*x[2]*x[2]*x[2]*x[2]*x[2]*x[2])
        def dfunc(x):
            return [-4 + 20*x[0] + 21*(x[0]*x[0]) - 20*(x[0]*x[0]*x[0]) - 5*(x[0]*x[0]*x[0]*x[0]) + 48*(x[0]*x[0]*x[0]*x[0]*x[0]) + 56*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0]) + 9*x[1] - 2*x[0]*x[1] - 21*(x[0]*x[0])*x[1] + 20*(x[0]*x[0]*x[0])*x[1] + 25*(x[0]*x[0]*x[0]*x[0])*x[1] + 54*(x[0]*x[0]*x[0]*x[0]*x[0])*x[1] - 10*(x[1]*x[1]) - 20*x[0]*(x[1]*x[1]) - 21*(x[0]*x[0])*(x[1]*x[1]) + 20*(x[0]*x[0]*x[0])*(x[1]*x[1]) - 50*(x[0]*x[0]*x[0]*x[0])*(x[1]*x[1]) - 9*(x[1]*x[1]*x[1]) - 6*x[0]*(x[1]*x[1]*x[1]) - 9*(x[0]*x[0])*(x[1]*x[1]*x[1]) + 40*(x[0]*x[0]*x[0])*(x[1]*x[1]*x[1]) + x[1]*x[1]*x[1]*x[1] - 16*x[0]*(x[1]*x[1]*x[1]*x[1]) + 27*(x[0]*x[0])*(x[1]*x[1]*x[1]*x[1]) + 6*(x[1]*x[1]*x[1]*x[1]*x[1]) - 4*x[0]*(x[1]*x[1]*x[1]*x[1]*x[1]) - 6*(x[1]*x[1]*x[1]*x[1]*x[1]*x[1]) - 3*x[2] - 20*x[0]*x[2] - 18*(x[0]*x[0])*x[2] + 40*(x[0]*x[0]*x[0])*x[2] + 45*(x[0]*x[0]*x[0]*x[0])*x[2] + 12*(x[0]*x[0]*x[0]*x[0]*x[0])*x[2] + 7*x[1]*x[2] + 14*x[0]*x[1]*x[2] - 12*(x[0]*x[0])*x[1]*x[2] + 16*(x[0]*x[0]*x[0])*x[1]*x[2] - 40*(x[0]*x[0]*x[0]*x[0])*x[1]*x[2] + 8*(x[1]*x[1])*x[2] + 4*x[0]*(x[1]*x[1])*x[2] + 30*(x[0]*x[0])*(x[1]*x[1])*x[2] - 32*(x[0]*x[0]*x[0])*(x[1]*x[1])*x[2] - 8*(x[1]*x[1]*x[1])*x[2] - 2*x[0]*(x[1]*x[1]*x[1])*x[2] + 21*(x[0]*x[0])*(x[1]*x[1]*x[1])*x[2] - 4*(x[1]*x[1]*x[1]*x[1])*x[2] + 18*x[0]*(x[1]*x[1]*x[1]*x[1])*x[2] + x[1]*x[1]*x[1]*x[1]*x[1]*x[2] - 2*(x[2]*x[2]) - 20*x[0]*(x[2]*x[2]) - 12*(x[0]*x[0])*(x[2]*x[2]) + 24*(x[0]*x[0]*x[0])*(x[2]*x[2]) + 50*(x[0]*x[0]*x[0]*x[0])*(x[2]*x[2]) - 8*x[1]*(x[2]*x[2]) - 2*x[0]*x[1]*(x[2]*x[2]) - 6*(x[0]*x[0])*x[1]*(x[2]*x[2]) - 16*(x[0]*x[0]*x[0])*x[1]*(x[2]*x[2]) + 2*(x[1]*x[1])*(x[2]*x[2]) - 20*x[0]*(x[1]*x[1])*(x[2]*x[2]) + 12*(x[0]*x[0])*(x[1]*x[1])*(x[2]*x[2]) - 10*(x[1]*x[1]*x[1])*(x[2]*x[2]) + 16*x[0]*(x[1]*x[1]*x[1])*(x[2]*x[2]) + 9*(x[1]*x[1]*x[1]*x[1])*(x[2]*x[2]) - 5*(x[2]*x[2]*x[2]) - 2*x[0]*(x[2]*x[2]*x[2]) - 27*(x[0]*x[0])*(x[2]*x[2]*x[2]) - 4*(x[0]*x[0]*x[0])*(x[2]*x[2]*x[2]) + 2*x[1]*(x[2]*x[2]*x[2]) - 16*x[0]*x[1]*(x[2]*x[2]*x[2]) + 6*(x[0]*x[0])*x[1]*(x[2]*x[2]*x[2]) + 10*(x[1]*x[1])*(x[2]*x[2]*x[2]) - 2*x[0]*(x[1]*x[1])*(x[2]*x[2]*x[2]) + 6*(x[1]*x[1]*x[1])*(x[2]*x[2]*x[2]) - 2*(x[2]*x[2]*x[2]*x[2]) - 16*x[0]*(x[2]*x[2]*x[2]*x[2]) - 24*(x[0]*x[0])*(x[2]*x[2]*x[2]*x[2]) + x[1]*(x[2]*x[2]*x[2]*x[2]) - 2*x[0]*x[1]*(x[2]*x[2]*x[2]*x[2]) - 9*(x[1]*x[1])*(x[2]*x[2]*x[2]*x[2]) + 10*(x[2]*x[2]*x[2]*x[2]*x[2]) + 10*x[0]*(x[2]*x[2]*x[2]*x[2]*x[2]) + 5*x[1]*(x[2]*x[2]*x[2]*x[2]*x[2]) - 2*(x[2]*x[2]*x[2]*x[2]*x[2]*x[2]), -9 + 9*x[0] - x[0]*x[0] - 7*(x[0]*x[0]*x[0]) + 5*(x[0]*x[0]*x[0]*x[0]) + 5*(x[0]*x[0]*x[0]*x[0]*x[0]) + 9*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0]) + 16*x[1] - 20*x[0]*x[1] - 20*(x[0]*x[0])*x[1] - 14*(x[0]*x[0]*x[0])*x[1] + 10*(x[0]*x[0]*x[0]*x[0])*x[1] - 20*(x[0]*x[0]*x[0]*x[0]*x[0])*x[1] - 24*(x[1]*x[1]) - 27*x[0]*(x[1]*x[1]) - 9*(x[0]*x[0])*(x[1]*x[1]) - 9*(x[0]*x[0]*x[0])*(x[1]*x[1]) + 30*(x[0]*x[0]*x[0]*x[0])*(x[1]*x[1]) + 4*(x[1]*x[1]*x[1]) + 4*x[0]*(x[1]*x[1]*x[1]) - 32*(x[0]*x[0])*(x[1]*x[1]*x[1]) + 36*(x[0]*x[0]*x[0])*(x[1]*x[1]*x[1]) + 50*(x[1]*x[1]*x[1]*x[1]) + 30*x[0]*(x[1]*x[1]*x[1]*x[1]) - 10*(x[0]*x[0])*(x[1]*x[1]*x[1]*x[1]) + 6*(x[1]*x[1]*x[1]*x[1]*x[1]) - 36*x[0]*(x[1]*x[1]*x[1]*x[1]*x[1]) - 56*(x[1]*x[1]*x[1]*x[1]*x[1]*x[1]) - 3*x[2] + 7*x[0]*x[2] + 7*(x[0]*x[0])*x[2] - 4*(x[0]*x[0]*x[0])*x[2] + 4*(x[0]*x[0]*x[0]*x[0])*x[2] - 8*(x[0]*x[0]*x[0]*x[0]*x[0])*x[2] + 2*x[1]*x[2] + 16*x[0]*x[1]*x[2] + 4*(x[0]*x[0])*x[1]*x[2] + 20*(x[0]*x[0]*x[0])*x[1]*x[2] - 16*(x[0]*x[0]*x[0]*x[0])*x[1]*x[2] + 6*(x[1]*x[1])*x[2] - 24*x[0]*(x[1]*x[1])*x[2] - 3*(x[0]*x[0])*(x[1]*x[1])*x[2] + 21*(x[0]*x[0]*x[0])*(x[1]*x[1])*x[2] + 20*(x[1]*x[1]*x[1])*x[2] - 16*x[0]*(x[1]*x[1]*x[1])*x[2] + 36*(x[0]*x[0])*(x[1]*x[1]*x[1])*x[2] - 10*(x[1]*x[1]*x[1]*x[1])*x[2] + 5*x[0]*(x[1]*x[1]*x[1]*x[1])*x[2] - 54*(x[1]*x[1]*x[1]*x[1]*x[1])*x[2] - 2*(x[2]*x[2]) - 8*x[0]*(x[2]*x[2]) - x[0]*x[0]*(x[2]*x[2]) - 2*(x[0]*x[0]*x[0])*(x[2]*x[2]) - 4*(x[0]*x[0]*x[0]*x[0])*(x[2]*x[2]) + 16*x[1]*(x[2]*x[2]) + 4*x[0]*x[1]*(x[2]*x[2]) - 20*(x[0]*x[0])*x[1]*(x[2]*x[2]) + 8*(x[0]*x[0]*x[0])*x[1]*(x[2]*x[2]) - 12*(x[1]*x[1])*(x[2]*x[2]) - 30*x[0]*(x[1]*x[1])*(x[2]*x[2]) + 24*(x[0]*x[0])*(x[1]*x[1])*(x[2]*x[2]) - 16*(x[1]*x[1]*x[1])*(x[2]*x[2]) + 36*x[0]*(x[1]*x[1]*x[1])*(x[2]*x[2]) - 50*(x[1]*x[1]*x[1]*x[1])*(x[2]*x[2]) + 3*(x[2]*x[2]*x[2]) + 2*x[0]*(x[2]*x[2]*x[2]) - 8*(x[0]*x[0])*(x[2]*x[2]*x[2]) + 2*(x[0]*x[0]*x[0])*(x[2]*x[2]*x[2]) - 4*x[1]*(x[2]*x[2]*x[2]) + 20*x[0]*x[1]*(x[2]*x[2]*x[2]) - 2*(x[0]*x[0])*x[1]*(x[2]*x[2]*x[2]) - 18*(x[1]*x[1])*(x[2]*x[2]*x[2]) + 18*x[0]*(x[1]*x[1])*(x[2]*x[2]*x[2]) - 24*(x[1]*x[1]*x[1])*(x[2]*x[2]*x[2]) + x[2]*x[2]*x[2]*x[2] + x[0]*(x[2]*x[2]*x[2]*x[2]) - x[0]*x[0]*(x[2]*x[2]*x[2]*x[2]) + 16*x[1]*(x[2]*x[2]*x[2]*x[2]) - 18*x[0]*x[1]*(x[2]*x[2]*x[2]*x[2]) + 18*(x[1]*x[1])*(x[2]*x[2]*x[2]*x[2]) - x[2]*x[2]*x[2]*x[2]*x[2] + 5*x[0]*(x[2]*x[2]*x[2]*x[2]*x[2]) - 2*x[1]*(x[2]*x[2]*x[2]*x[2]*x[2]) + 7*(x[2]*x[2]*x[2]*x[2]*x[2]*x[2]), 3 - 3*x[0] - 10*(x[0]*x[0]) - 6*(x[0]*x[0]*x[0]) + 10*(x[0]*x[0]*x[0]*x[0]) + 9*(x[0]*x[0]*x[0]*x[0]*x[0]) + 2*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0]) - 3*x[1] + 7*x[0]*x[1] + 7*(x[0]*x[0])*x[1] - 4*(x[0]*x[0]*x[0])*x[1] + 4*(x[0]*x[0]*x[0]*x[0])*x[1] - 8*(x[0]*x[0]*x[0]*x[0]*x[0])*x[1] + x[1]*x[1] + 8*x[0]*(x[1]*x[1]) + 2*(x[0]*x[0])*(x[1]*x[1]) + 10*(x[0]*x[0]*x[0])*(x[1]*x[1]) - 8*(x[0]*x[0]*x[0]*x[0])*(x[1]*x[1]) + 2*(x[1]*x[1]*x[1]) - 8*x[0]*(x[1]*x[1]*x[1]) - x[0]*x[0]*(x[1]*x[1]*x[1]) + 7*(x[0]*x[0]*x[0])*(x[1]*x[1]*x[1]) + 5*(x[1]*x[1]*x[1]*x[1]) - 4*x[0]*(x[1]*x[1]*x[1]*x[1]) + 9*(x[0]*x[0])*(x[1]*x[1]*x[1]*x[1]) - 2*(x[1]*x[1]*x[1]*x[1]*x[1]) + x[0]*(x[1]*x[1]*x[1]*x[1]*x[1]) - 9*(x[1]*x[1]*x[1]*x[1]*x[1]*x[1]) + 20*x[2] - 4*x[0]*x[2] - 20*(x[0]*x[0])*x[2] - 8*(x[0]*x[0]*x[0])*x[2] + 12*(x[0]*x[0]*x[0]*x[0])*x[2] + 20*(x[0]*x[0]*x[0]*x[0]*x[0])*x[2] - 4*x[1]*x[2] - 16*x[0]*x[1]*x[2] - 2*(x[0]*x[0])*x[1]*x[2] - 4*(x[0]*x[0]*x[0])*x[1]*x[2] - 8*(x[0]*x[0]*x[0]*x[0])*x[1]*x[2] + 16*(x[1]*x[1])*x[2] + 4*x[0]*(x[1]*x[1])*x[2] - 20*(x[0]*x[0])*(x[1]*x[1])*x[2] + 8*(x[0]*x[0]*x[0])*(x[1]*x[1])*x[2] - 8*(x[1]*x[1]*x[1])*x[2] - 20*x[0]*(x[1]*x[1]*x[1])*x[2] + 16*(x[0]*x[0])*(x[1]*x[1]*x[1])*x[2] - 8*(x[1]*x[1]*x[1]*x[1])*x[2] + 18*x[0]*(x[1]*x[1]*x[1]*x[1])*x[2] - 20*(x[1]*x[1]*x[1]*x[1]*x[1])*x[2] + 9*(x[2]*x[2]) - 15*x[0]*(x[2]*x[2]) - 3*(x[0]*x[0])*(x[2]*x[2]) - 27*(x[0]*x[0]*x[0])*(x[2]*x[2]) - 3*(x[0]*x[0]*x[0]*x[0])*(x[2]*x[2]) + 9*x[1]*(x[2]*x[2]) + 6*x[0]*x[1]*(x[2]*x[2]) - 24*(x[0]*x[0])*x[1]*(x[2]*x[2]) + 6*(x[0]*x[0]*x[0])*x[1]*(x[2]*x[2]) - 6*(x[1]*x[1])*(x[2]*x[2]) + 30*x[0]*(x[1]*x[1])*(x[2]*x[2]) - 3*(x[0]*x[0])*(x[1]*x[1])*(x[2]*x[2]) - 18*(x[1]*x[1]*x[1])*(x[2]*x[2]) + 18*x[0]*(x[1]*x[1]*x[1])*(x[2]*x[2]) - 18*(x[1]*x[1]*x[1]*x[1])*(x[2]*x[2]) + 28*(x[2]*x[2]*x[2]) - 8*x[0]*(x[2]*x[2]*x[2]) - 32*(x[0]*x[0])*(x[2]*x[2]*x[2]) - 32*(x[0]*x[0]*x[0])*(x[2]*x[2]*x[2]) + 4*x[1]*(x[2]*x[2]*x[2]) + 4*x[0]*x[1]*(x[2]*x[2]*x[2]) - 4*(x[0]*x[0])*x[1]*(x[2]*x[2]*x[2]) + 32*(x[1]*x[1])*(x[2]*x[2]*x[2]) - 36*x[0]*(x[1]*x[1])*(x[2]*x[2]*x[2]) + 24*(x[1]*x[1]*x[1])*(x[2]*x[2]*x[2]) + 30*(x[2]*x[2]*x[2]*x[2]) + 50*x[0]*(x[2]*x[2]*x[2]*x[2]) + 25*(x[0]*x[0])*(x[2]*x[2]*x[2]*x[2]) - 5*x[1]*(x[2]*x[2]*x[2]*x[2]) + 25*x[0]*x[1]*(x[2]*x[2]*x[2]*x[2]) - 5*(x[1]*x[1])*(x[2]*x[2]*x[2]*x[2]) + 36*(x[2]*x[2]*x[2]*x[2]*x[2]) - 12*x[0]*(x[2]*x[2]*x[2]*x[2]*x[2]) + 42*x[1]*(x[2]*x[2]*x[2]*x[2]*x[2]) - 14*(x[2]*x[2]*x[2]*x[2]*x[2]*x[2])]
        def ddfunc(x):
            return [[20 + 42*x[0] - 60*(x[0]*x[0]) - 20*(x[0]*x[0]*x[0]) + 240*(x[0]*x[0]*x[0]*x[0]) + 336*(x[0]*x[0]*x[0]*x[0]*x[0]) - 2*x[1] - 42*x[0]*x[1] + 60*(x[0]*x[0])*x[1] + 100*(x[0]*x[0]*x[0])*x[1] + 270*(x[0]*x[0]*x[0]*x[0])*x[1] - 20*(x[1]*x[1]) - 42*x[0]*(x[1]*x[1]) + 60*(x[0]*x[0])*(x[1]*x[1]) - 200*(x[0]*x[0]*x[0])*(x[1]*x[1]) - 6*(x[1]*x[1]*x[1]) - 18*x[0]*(x[1]*x[1]*x[1]) + 120*(x[0]*x[0])*(x[1]*x[1]*x[1]) - 16*(x[1]*x[1]*x[1]*x[1]) + 54*x[0]*(x[1]*x[1]*x[1]*x[1]) - 4*(x[1]*x[1]*x[1]*x[1]*x[1]) - 20*x[2] - 36*x[0]*x[2] + 120*(x[0]*x[0])*x[2] + 180*(x[0]*x[0]*x[0])*x[2] + 60*(x[0]*x[0]*x[0]*x[0])*x[2] + 14*x[1]*x[2] - 24*x[0]*x[1]*x[2] + 48*(x[0]*x[0])*x[1]*x[2] - 160*(x[0]*x[0]*x[0])*x[1]*x[2] + 4*(x[1]*x[1])*x[2] + 60*x[0]*(x[1]*x[1])*x[2] - 96*(x[0]*x[0])*(x[1]*x[1])*x[2] - 2*(x[1]*x[1]*x[1])*x[2] + 42*x[0]*(x[1]*x[1]*x[1])*x[2] + 18*(x[1]*x[1]*x[1]*x[1])*x[2] - 20*(x[2]*x[2]) - 24*x[0]*(x[2]*x[2]) + 72*(x[0]*x[0])*(x[2]*x[2]) + 200*(x[0]*x[0]*x[0])*(x[2]*x[2]) - 2*x[1]*(x[2]*x[2]) - 12*x[0]*x[1]*(x[2]*x[2]) - 48*(x[0]*x[0])*x[1]*(x[2]*x[2]) - 20*(x[1]*x[1])*(x[2]*x[2]) + 24*x[0]*(x[1]*x[1])*(x[2]*x[2]) + 16*(x[1]*x[1]*x[1])*(x[2]*x[2]) - 2*(x[2]*x[2]*x[2]) - 54*x[0]*(x[2]*x[2]*x[2]) - 12*(x[0]*x[0])*(x[2]*x[2]*x[2]) - 16*x[1]*(x[2]*x[2]*x[2]) + 12*x[0]*x[1]*(x[2]*x[2]*x[2]) - 2*(x[1]*x[1])*(x[2]*x[2]*x[2]) - 16*(x[2]*x[2]*x[2]*x[2]) - 48*x[0]*(x[2]*x[2]*x[2]*x[2]) - 2*x[1]*(x[2]*x[2]*x[2]*x[2]) + 10*(x[2]*x[2]*x[2]*x[2]*x[2]), 9 - 2*x[0] - 21*(x[0]*x[0]) + 20*(x[0]*x[0]*x[0]) + 25*(x[0]*x[0]*x[0]*x[0]) + 54*(x[0]*x[0]*x[0]*x[0]*x[0]) - 20*x[1] - 40*x[0]*x[1] - 42*(x[0]*x[0])*x[1] + 40*(x[0]*x[0]*x[0])*x[1] - 100*(x[0]*x[0]*x[0]*x[0])*x[1] - 27*(x[1]*x[1]) - 18*x[0]*(x[1]*x[1]) - 27*(x[0]*x[0])*(x[1]*x[1]) + 120*(x[0]*x[0]*x[0])*(x[1]*x[1]) + 4*(x[1]*x[1]*x[1]) - 64*x[0]*(x[1]*x[1]*x[1]) + 108*(x[0]*x[0])*(x[1]*x[1]*x[1]) + 30*(x[1]*x[1]*x[1]*x[1]) - 20*x[0]*(x[1]*x[1]*x[1]*x[1]) - 36*(x[1]*x[1]*x[1]*x[1]*x[1]) + 7*x[2] + 14*x[0]*x[2] - 12*(x[0]*x[0])*x[2] + 16*(x[0]*x[0]*x[0])*x[2] - 40*(x[0]*x[0]*x[0]*x[0])*x[2] + 16*x[1]*x[2] + 8*x[0]*x[1]*x[2] + 60*(x[0]*x[0])*x[1]*x[2] - 64*(x[0]*x[0]*x[0])*x[1]*x[2] - 24*(x[1]*x[1])*x[2] - 6*x[0]*(x[1]*x[1])*x[2] + 63*(x[0]*x[0])*(x[1]*x[1])*x[2] - 16*(x[1]*x[1]*x[1])*x[2] + 72*x[0]*(x[1]*x[1]*x[1])*x[2] + 5*(x[1]*x[1]*x[1]*x[1])*x[2] - 8*(x[2]*x[2]) - 2*x[0]*(x[2]*x[2]) - 6*(x[0]*x[0])*(x[2]*x[2]) - 16*(x[0]*x[0]*x[0])*(x[2]*x[2]) + 4*x[1]*(x[2]*x[2]) - 40*x[0]*x[1]*(x[2]*x[2]) + 24*(x[0]*x[0])*x[1]*(x[2]*x[2]) - 30*(x[1]*x[1])*(x[2]*x[2]) + 48*x[0]*(x[1]*x[1])*(x[2]*x[2]) + 36*(x[1]*x[1]*x[1])*(x[2]*x[2]) + 2*(x[2]*x[2]*x[2]) - 16*x[0]*(x[2]*x[2]*x[2]) + 6*(x[0]*x[0])*(x[2]*x[2]*x[2]) + 20*x[1]*(x[2]*x[2]*x[2]) - 4*x[0]*x[1]*(x[2]*x[2]*x[2]) + 18*(x[1]*x[1])*(x[2]*x[2]*x[2]) + x[2]*x[2]*x[2]*x[2] - 2*x[0]*(x[2]*x[2]*x[2]*x[2]) - 18*x[1]*(x[2]*x[2]*x[2]*x[2]) + 5*(x[2]*x[2]*x[2]*x[2]*x[2]), -3 - 20*x[0] - 18*(x[0]*x[0]) + 40*(x[0]*x[0]*x[0]) + 45*(x[0]*x[0]*x[0]*x[0]) + 12*(x[0]*x[0]*x[0]*x[0]*x[0]) + 7*x[1] + 14*x[0]*x[1] - 12*(x[0]*x[0])*x[1] + 16*(x[0]*x[0]*x[0])*x[1] - 40*(x[0]*x[0]*x[0]*x[0])*x[1] + 8*(x[1]*x[1]) + 4*x[0]*(x[1]*x[1]) + 30*(x[0]*x[0])*(x[1]*x[1]) - 32*(x[0]*x[0]*x[0])*(x[1]*x[1]) - 8*(x[1]*x[1]*x[1]) - 2*x[0]*(x[1]*x[1]*x[1]) + 21*(x[0]*x[0])*(x[1]*x[1]*x[1]) - 4*(x[1]*x[1]*x[1]*x[1]) + 18*x[0]*(x[1]*x[1]*x[1]*x[1]) + x[1]*x[1]*x[1]*x[1]*x[1] - 4*x[2] - 40*x[0]*x[2] - 24*(x[0]*x[0])*x[2] + 48*(x[0]*x[0]*x[0])*x[2] + 100*(x[0]*x[0]*x[0]*x[0])*x[2] - 16*x[1]*x[2] - 4*x[0]*x[1]*x[2] - 12*(x[0]*x[0])*x[1]*x[2] - 32*(x[0]*x[0]*x[0])*x[1]*x[2] + 4*(x[1]*x[1])*x[2] - 40*x[0]*(x[1]*x[1])*x[2] + 24*(x[0]*x[0])*(x[1]*x[1])*x[2] - 20*(x[1]*x[1]*x[1])*x[2] + 32*x[0]*(x[1]*x[1]*x[1])*x[2] + 18*(x[1]*x[1]*x[1]*x[1])*x[2] - 15*(x[2]*x[2]) - 6*x[0]*(x[2]*x[2]) - 81*(x[0]*x[0])*(x[2]*x[2]) - 12*(x[0]*x[0]*x[0])*(x[2]*x[2]) + 6*x[1]*(x[2]*x[2]) - 48*x[0]*x[1]*(x[2]*x[2]) + 18*(x[0]*x[0])*x[1]*(x[2]*x[2]) + 30*(x[1]*x[1])*(x[2]*x[2]) - 6*x[0]*(x[1]*x[1])*(x[2]*x[2]) + 18*(x[1]*x[1]*x[1])*(x[2]*x[2]) - 8*(x[2]*x[2]*x[2]) - 64*x[0]*(x[2]*x[2]*x[2]) - 96*(x[0]*x[0])*(x[2]*x[2]*x[2]) + 4*x[1]*(x[2]*x[2]*x[2]) - 8*x[0]*x[1]*(x[2]*x[2]*x[2]) - 36*(x[1]*x[1])*(x[2]*x[2]*x[2]) + 50*(x[2]*x[2]*x[2]*x[2]) + 50*x[0]*(x[2]*x[2]*x[2]*x[2]) + 25*x[1]*(x[2]*x[2]*x[2]*x[2]) - 12*(x[2]*x[2]*x[2]*x[2]*x[2])], [9 - 2*x[0] - 21*(x[0]*x[0]) + 20*(x[0]*x[0]*x[0]) + 25*(x[0]*x[0]*x[0]*x[0]) + 54*(x[0]*x[0]*x[0]*x[0]*x[0]) - 20*x[1] - 40*x[0]*x[1] - 42*(x[0]*x[0])*x[1] + 40*(x[0]*x[0]*x[0])*x[1] - 100*(x[0]*x[0]*x[0]*x[0])*x[1] - 27*(x[1]*x[1]) - 18*x[0]*(x[1]*x[1]) - 27*(x[0]*x[0])*(x[1]*x[1]) + 120*(x[0]*x[0]*x[0])*(x[1]*x[1]) + 4*(x[1]*x[1]*x[1]) - 64*x[0]*(x[1]*x[1]*x[1]) + 108*(x[0]*x[0])*(x[1]*x[1]*x[1]) + 30*(x[1]*x[1]*x[1]*x[1]) - 20*x[0]*(x[1]*x[1]*x[1]*x[1]) - 36*(x[1]*x[1]*x[1]*x[1]*x[1]) + 7*x[2] + 14*x[0]*x[2] - 12*(x[0]*x[0])*x[2] + 16*(x[0]*x[0]*x[0])*x[2] - 40*(x[0]*x[0]*x[0]*x[0])*x[2] + 16*x[1]*x[2] + 8*x[0]*x[1]*x[2] + 60*(x[0]*x[0])*x[1]*x[2] - 64*(x[0]*x[0]*x[0])*x[1]*x[2] - 24*(x[1]*x[1])*x[2] - 6*x[0]*(x[1]*x[1])*x[2] + 63*(x[0]*x[0])*(x[1]*x[1])*x[2] - 16*(x[1]*x[1]*x[1])*x[2] + 72*x[0]*(x[1]*x[1]*x[1])*x[2] + 5*(x[1]*x[1]*x[1]*x[1])*x[2] - 8*(x[2]*x[2]) - 2*x[0]*(x[2]*x[2]) - 6*(x[0]*x[0])*(x[2]*x[2]) - 16*(x[0]*x[0]*x[0])*(x[2]*x[2]) + 4*x[1]*(x[2]*x[2]) - 40*x[0]*x[1]*(x[2]*x[2]) + 24*(x[0]*x[0])*x[1]*(x[2]*x[2]) - 30*(x[1]*x[1])*(x[2]*x[2]) + 48*x[0]*(x[1]*x[1])*(x[2]*x[2]) + 36*(x[1]*x[1]*x[1])*(x[2]*x[2]) + 2*(x[2]*x[2]*x[2]) - 16*x[0]*(x[2]*x[2]*x[2]) + 6*(x[0]*x[0])*(x[2]*x[2]*x[2]) + 20*x[1]*(x[2]*x[2]*x[2]) - 4*x[0]*x[1]*(x[2]*x[2]*x[2]) + 18*(x[1]*x[1])*(x[2]*x[2]*x[2]) + x[2]*x[2]*x[2]*x[2] - 2*x[0]*(x[2]*x[2]*x[2]*x[2]) - 18*x[1]*(x[2]*x[2]*x[2]*x[2]) + 5*(x[2]*x[2]*x[2]*x[2]*x[2]), 16 - 20*x[0] - 20*(x[0]*x[0]) - 14*(x[0]*x[0]*x[0]) + 10*(x[0]*x[0]*x[0]*x[0]) - 20*(x[0]*x[0]*x[0]*x[0]*x[0]) - 48*x[1] - 54*x[0]*x[1] - 18*(x[0]*x[0])*x[1] - 18*(x[0]*x[0]*x[0])*x[1] + 60*(x[0]*x[0]*x[0]*x[0])*x[1] + 12*(x[1]*x[1]) + 12*x[0]*(x[1]*x[1]) - 96*(x[0]*x[0])*(x[1]*x[1]) + 108*(x[0]*x[0]*x[0])*(x[1]*x[1]) + 200*(x[1]*x[1]*x[1]) + 120*x[0]*(x[1]*x[1]*x[1]) - 40*(x[0]*x[0])*(x[1]*x[1]*x[1]) + 30*(x[1]*x[1]*x[1]*x[1]) - 180*x[0]*(x[1]*x[1]*x[1]*x[1]) - 336*(x[1]*x[1]*x[1]*x[1]*x[1]) + 2*x[2] + 16*x[0]*x[2] + 4*(x[0]*x[0])*x[2] + 20*(x[0]*x[0]*x[0])*x[2] - 16*(x[0]*x[0]*x[0]*x[0])*x[2] + 12*x[1]*x[2] - 48*x[0]*x[1]*x[2] - 6*(x[0]*x[0])*x[1]*x[2] + 42*(x[0]*x[0]*x[0])*x[1]*x[2] + 60*(x[1]*x[1])*x[2] - 48*x[0]*(x[1]*x[1])*x[2] + 108*(x[0]*x[0])*(x[1]*x[1])*x[2] - 40*(x[1]*x[1]*x[1])*x[2] + 20*x[0]*(x[1]*x[1]*x[1])*x[2] - 270*(x[1]*x[1]*x[1]*x[1])*x[2] + 16*(x[2]*x[2]) + 4*x[0]*(x[2]*x[2]) - 20*(x[0]*x[0])*(x[2]*x[2]) + 8*(x[0]*x[0]*x[0])*(x[2]*x[2]) - 24*x[1]*(x[2]*x[2]) - 60*x[0]*x[1]*(x[2]*x[2]) + 48*(x[0]*x[0])*x[1]*(x[2]*x[2]) - 48*(x[1]*x[1])*(x[2]*x[2]) + 108*x[0]*(x[1]*x[1])*(x[2]*x[2]) - 200*(x[1]*x[1]*x[1])*(x[2]*x[2]) - 4*(x[2]*x[2]*x[2]) + 20*x[0]*(x[2]*x[2]*x[2]) - 2*(x[0]*x[0])*(x[2]*x[2]*x[2]) - 36*x[1]*(x[2]*x[2]*x[2]) + 36*x[0]*x[1]*(x[2]*x[2]*x[2]) - 72*(x[1]*x[1])*(x[2]*x[2]*x[2]) + 16*(x[2]*x[2]*x[2]*x[2]) - 18*x[0]*(x[2]*x[2]*x[2]*x[2]) + 36*x[1]*(x[2]*x[2]*x[2]*x[2]) - 2*(x[2]*x[2]*x[2]*x[2]*x[2]), -3 + 7*x[0] + 7*(x[0]*x[0]) - 4*(x[0]*x[0]*x[0]) + 4*(x[0]*x[0]*x[0]*x[0]) - 8*(x[0]*x[0]*x[0]*x[0]*x[0]) + 2*x[1] + 16*x[0]*x[1] + 4*(x[0]*x[0])*x[1] + 20*(x[0]*x[0]*x[0])*x[1] - 16*(x[0]*x[0]*x[0]*x[0])*x[1] + 6*(x[1]*x[1]) - 24*x[0]*(x[1]*x[1]) - 3*(x[0]*x[0])*(x[1]*x[1]) + 21*(x[0]*x[0]*x[0])*(x[1]*x[1]) + 20*(x[1]*x[1]*x[1]) - 16*x[0]*(x[1]*x[1]*x[1]) + 36*(x[0]*x[0])*(x[1]*x[1]*x[1]) - 10*(x[1]*x[1]*x[1]*x[1]) + 5*x[0]*(x[1]*x[1]*x[1]*x[1]) - 54*(x[1]*x[1]*x[1]*x[1]*x[1]) - 4*x[2] - 16*x[0]*x[2] - 2*(x[0]*x[0])*x[2] - 4*(x[0]*x[0]*x[0])*x[2] - 8*(x[0]*x[0]*x[0]*x[0])*x[2] + 32*x[1]*x[2] + 8*x[0]*x[1]*x[2] - 40*(x[0]*x[0])*x[1]*x[2] + 16*(x[0]*x[0]*x[0])*x[1]*x[2] - 24*(x[1]*x[1])*x[2] - 60*x[0]*(x[1]*x[1])*x[2] + 48*(x[0]*x[0])*(x[1]*x[1])*x[2] - 32*(x[1]*x[1]*x[1])*x[2] + 72*x[0]*(x[1]*x[1]*x[1])*x[2] - 100*(x[1]*x[1]*x[1]*x[1])*x[2] + 9*(x[2]*x[2]) + 6*x[0]*(x[2]*x[2]) - 24*(x[0]*x[0])*(x[2]*x[2]) + 6*(x[0]*x[0]*x[0])*(x[2]*x[2]) - 12*x[1]*(x[2]*x[2]) + 60*x[0]*x[1]*(x[2]*x[2]) - 6*(x[0]*x[0])*x[1]*(x[2]*x[2]) - 54*(x[1]*x[1])*(x[2]*x[2]) + 54*x[0]*(x[1]*x[1])*(x[2]*x[2]) - 72*(x[1]*x[1]*x[1])*(x[2]*x[2]) + 4*(x[2]*x[2]*x[2]) + 4*x[0]*(x[2]*x[2]*x[2]) - 4*(x[0]*x[0])*(x[2]*x[2]*x[2]) + 64*x[1]*(x[2]*x[2]*x[2]) - 72*x[0]*x[1]*(x[2]*x[2]*x[2]) + 72*(x[1]*x[1])*(x[2]*x[2]*x[2]) - 5*(x[2]*x[2]*x[2]*x[2]) + 25*x[0]*(x[2]*x[2]*x[2]*x[2]) - 10*x[1]*(x[2]*x[2]*x[2]*x[2]) + 42*(x[2]*x[2]*x[2]*x[2]*x[2])], [-3 - 20*x[0] - 18*(x[0]*x[0]) + 40*(x[0]*x[0]*x[0]) + 45*(x[0]*x[0]*x[0]*x[0]) + 12*(x[0]*x[0]*x[0]*x[0]*x[0]) + 7*x[1] + 14*x[0]*x[1] - 12*(x[0]*x[0])*x[1] + 16*(x[0]*x[0]*x[0])*x[1] - 40*(x[0]*x[0]*x[0]*x[0])*x[1] + 8*(x[1]*x[1]) + 4*x[0]*(x[1]*x[1]) + 30*(x[0]*x[0])*(x[1]*x[1]) - 32*(x[0]*x[0]*x[0])*(x[1]*x[1]) - 8*(x[1]*x[1]*x[1]) - 2*x[0]*(x[1]*x[1]*x[1]) + 21*(x[0]*x[0])*(x[1]*x[1]*x[1]) - 4*(x[1]*x[1]*x[1]*x[1]) + 18*x[0]*(x[1]*x[1]*x[1]*x[1]) + x[1]*x[1]*x[1]*x[1]*x[1] - 4*x[2] - 40*x[0]*x[2] - 24*(x[0]*x[0])*x[2] + 48*(x[0]*x[0]*x[0])*x[2] + 100*(x[0]*x[0]*x[0]*x[0])*x[2] - 16*x[1]*x[2] - 4*x[0]*x[1]*x[2] - 12*(x[0]*x[0])*x[1]*x[2] - 32*(x[0]*x[0]*x[0])*x[1]*x[2] + 4*(x[1]*x[1])*x[2] - 40*x[0]*(x[1]*x[1])*x[2] + 24*(x[0]*x[0])*(x[1]*x[1])*x[2] - 20*(x[1]*x[1]*x[1])*x[2] + 32*x[0]*(x[1]*x[1]*x[1])*x[2] + 18*(x[1]*x[1]*x[1]*x[1])*x[2] - 15*(x[2]*x[2]) - 6*x[0]*(x[2]*x[2]) - 81*(x[0]*x[0])*(x[2]*x[2]) - 12*(x[0]*x[0]*x[0])*(x[2]*x[2]) + 6*x[1]*(x[2]*x[2]) - 48*x[0]*x[1]*(x[2]*x[2]) + 18*(x[0]*x[0])*x[1]*(x[2]*x[2]) + 30*(x[1]*x[1])*(x[2]*x[2]) - 6*x[0]*(x[1]*x[1])*(x[2]*x[2]) + 18*(x[1]*x[1]*x[1])*(x[2]*x[2]) - 8*(x[2]*x[2]*x[2]) - 64*x[0]*(x[2]*x[2]*x[2]) - 96*(x[0]*x[0])*(x[2]*x[2]*x[2]) + 4*x[1]*(x[2]*x[2]*x[2]) - 8*x[0]*x[1]*(x[2]*x[2]*x[2]) - 36*(x[1]*x[1])*(x[2]*x[2]*x[2]) + 50*(x[2]*x[2]*x[2]*x[2]) + 50*x[0]*(x[2]*x[2]*x[2]*x[2]) + 25*x[1]*(x[2]*x[2]*x[2]*x[2]) - 12*(x[2]*x[2]*x[2]*x[2]*x[2]), -3 + 7*x[0] + 7*(x[0]*x[0]) - 4*(x[0]*x[0]*x[0]) + 4*(x[0]*x[0]*x[0]*x[0]) - 8*(x[0]*x[0]*x[0]*x[0]*x[0]) + 2*x[1] + 16*x[0]*x[1] + 4*(x[0]*x[0])*x[1] + 20*(x[0]*x[0]*x[0])*x[1] - 16*(x[0]*x[0]*x[0]*x[0])*x[1] + 6*(x[1]*x[1]) - 24*x[0]*(x[1]*x[1]) - 3*(x[0]*x[0])*(x[1]*x[1]) + 21*(x[0]*x[0]*x[0])*(x[1]*x[1]) + 20*(x[1]*x[1]*x[1]) - 16*x[0]*(x[1]*x[1]*x[1]) + 36*(x[0]*x[0])*(x[1]*x[1]*x[1]) - 10*(x[1]*x[1]*x[1]*x[1]) + 5*x[0]*(x[1]*x[1]*x[1]*x[1]) - 54*(x[1]*x[1]*x[1]*x[1]*x[1]) - 4*x[2] - 16*x[0]*x[2] - 2*(x[0]*x[0])*x[2] - 4*(x[0]*x[0]*x[0])*x[2] - 8*(x[0]*x[0]*x[0]*x[0])*x[2] + 32*x[1]*x[2] + 8*x[0]*x[1]*x[2] - 40*(x[0]*x[0])*x[1]*x[2] + 16*(x[0]*x[0]*x[0])*x[1]*x[2] - 24*(x[1]*x[1])*x[2] - 60*x[0]*(x[1]*x[1])*x[2] + 48*(x[0]*x[0])*(x[1]*x[1])*x[2] - 32*(x[1]*x[1]*x[1])*x[2] + 72*x[0]*(x[1]*x[1]*x[1])*x[2] - 100*(x[1]*x[1]*x[1]*x[1])*x[2] + 9*(x[2]*x[2]) + 6*x[0]*(x[2]*x[2]) - 24*(x[0]*x[0])*(x[2]*x[2]) + 6*(x[0]*x[0]*x[0])*(x[2]*x[2]) - 12*x[1]*(x[2]*x[2]) + 60*x[0]*x[1]*(x[2]*x[2]) - 6*(x[0]*x[0])*x[1]*(x[2]*x[2]) - 54*(x[1]*x[1])*(x[2]*x[2]) + 54*x[0]*(x[1]*x[1])*(x[2]*x[2]) - 72*(x[1]*x[1]*x[1])*(x[2]*x[2]) + 4*(x[2]*x[2]*x[2]) + 4*x[0]*(x[2]*x[2]*x[2]) - 4*(x[0]*x[0])*(x[2]*x[2]*x[2]) + 64*x[1]*(x[2]*x[2]*x[2]) - 72*x[0]*x[1]*(x[2]*x[2]*x[2]) + 72*(x[1]*x[1])*(x[2]*x[2]*x[2]) - 5*(x[2]*x[2]*x[2]*x[2]) + 25*x[0]*(x[2]*x[2]*x[2]*x[2]) - 10*x[1]*(x[2]*x[2]*x[2]*x[2]) + 42*(x[2]*x[2]*x[2]*x[2]*x[2]), 20 - 4*x[0] - 20*(x[0]*x[0]) - 8*(x[0]*x[0]*x[0]) + 12*(x[0]*x[0]*x[0]*x[0]) + 20*(x[0]*x[0]*x[0]*x[0]*x[0]) - 4*x[1] - 16*x[0]*x[1] - 2*(x[0]*x[0])*x[1] - 4*(x[0]*x[0]*x[0])*x[1] - 8*(x[0]*x[0]*x[0]*x[0])*x[1] + 16*(x[1]*x[1]) + 4*x[0]*(x[1]*x[1]) - 20*(x[0]*x[0])*(x[1]*x[1]) + 8*(x[0]*x[0]*x[0])*(x[1]*x[1]) - 8*(x[1]*x[1]*x[1]) - 20*x[0]*(x[1]*x[1]*x[1]) + 16*(x[0]*x[0])*(x[1]*x[1]*x[1]) - 8*(x[1]*x[1]*x[1]*x[1]) + 18*x[0]*(x[1]*x[1]*x[1]*x[1]) - 20*(x[1]*x[1]*x[1]*x[1]*x[1]) + 18*x[2] - 30*x[0]*x[2] - 6*(x[0]*x[0])*x[2] - 54*(x[0]*x[0]*x[0])*x[2] - 6*(x[0]*x[0]*x[0]*x[0])*x[2] + 18*x[1]*x[2] + 12*x[0]*x[1]*x[2] - 48*(x[0]*x[0])*x[1]*x[2] + 12*(x[0]*x[0]*x[0])*x[1]*x[2] - 12*(x[1]*x[1])*x[2] + 60*x[0]*(x[1]*x[1])*x[2] - 6*(x[0]*x[0])*(x[1]*x[1])*x[2] - 36*(x[1]*x[1]*x[1])*x[2] + 36*x[0]*(x[1]*x[1]*x[1])*x[2] - 36*(x[1]*x[1]*x[1]*x[1])*x[2] + 84*(x[2]*x[2]) - 24*x[0]*(x[2]*x[2]) - 96*(x[0]*x[0])*(x[2]*x[2]) - 96*(x[0]*x[0]*x[0])*(x[2]*x[2]) + 12*x[1]*(x[2]*x[2]) + 12*x[0]*x[1]*(x[2]*x[2]) - 12*(x[0]*x[0])*x[1]*(x[2]*x[2]) + 96*(x[1]*x[1])*(x[2]*x[2]) - 108*x[0]*(x[1]*x[1])*(x[2]*x[2]) + 72*(x[1]*x[1]*x[1])*(x[2]*x[2]) + 120*(x[2]*x[2]*x[2]) + 200*x[0]*(x[2]*x[2]*x[2]) + 100*(x[0]*x[0])*(x[2]*x[2]*x[2]) - 20*x[1]*(x[2]*x[2]*x[2]) + 100*x[0]*x[1]*(x[2]*x[2]*x[2]) - 20*(x[1]*x[1])*(x[2]*x[2]*x[2]) + 180*(x[2]*x[2]*x[2]*x[2]) - 60*x[0]*(x[2]*x[2]*x[2]*x[2]) + 210*x[1]*(x[2]*x[2]*x[2]*x[2]) - 84*(x[2]*x[2]*x[2]*x[2]*x[2])]]
elif funcType == "sinusoidal":
    if dimension == 1:
        def func(x):
            return 3. + 2.*np.sin(5*x[0])
        def dfunc(x):
            return [10.*np.cos(5*x[0])]
        def ddfunc(x):
            return [[-50.*np.sin(5*x[0])]]
    elif dimension == 2:
        def func(x):
            return 3. + 2.*(np.sin(5*x[0]) + np.sin(5*x[1]))
        def dfunc(x):
            return [10.*np.cos(5*x[0]), 10.*np.cos(5*x[1])]
        def ddfunc(x):
            return [[-50.*np.sin(5*x[0]), 0], [0, -50.*np.sin(5*x[1])]]
    else:
        def func(x):
            return 3. + 2.*(np.sin(5*x[0]) + np.sin(5*x[1]) + np.sin(5*x[2]))
        def dfunc(x):
            return [10.*np.cos(5*x[0]), 10.*np.cos(5*x[1]), 10.*np.cos(5*x[2])]
        def ddfunc(x):
            return [[-50.*np.sin(5*x[0]), 0, 0], [0, -50.*np.sin(5*x[1]), 0], [0, 0, -50.*np.sin(5*x[2])]]
else:
    raise ValueError("function type {} not found".format(funcType))

#-------------------------------------------------------------------------------
# Create RK object
#-------------------------------------------------------------------------------
rk = RKCorrections(orders = set([ZerothOrder, correctionOrder]),
                   dataBase = dataBase,
                   W = WT,
                   volumeType = volumeType,
                   needHessian = testHessian,
                   updateInFinalize = False)
output("rk")
output("rk.correctionOrders")
WR = rk.WR(correctionOrder)
output("WR")
packages = [rk]

#-------------------------------------------------------------------------------
# Create a state directly and initialize physics package
#-------------------------------------------------------------------------------
connectivity = dataBase.connectivityMap()
state = State(dataBase, packages)
derivs = StateDerivatives(dataBase, packages)
rk.initializeProblemStartup(dataBase)
rk.registerState(dataBase, state)
rk.registerDerivatives(dataBase, derivs)
rk.preStepInitialize(dataBase, state, derivs)

#-------------------------------------------------------------------------------
# Get data from state
#-------------------------------------------------------------------------------
position = state.vectorFields(HydroFieldNames.position)
H = state.symTensorFields(HydroFieldNames.H)
volume = state.scalarFields(HydroFieldNames.volume)
zerothCorrections = state.RKCoefficientsFields(RKFieldNames.rkCorrections(ZerothOrder))
corrections = state.RKCoefficientsFields(RKFieldNames.rkCorrections(correctionOrder))
# normal = state.vectorFields(HydroFieldNames.normal)

#-------------------------------------------------------------------------------
# Compute corrections
#-------------------------------------------------------------------------------
rk_time = time.time()
WR.computeCorrections(connectivity, volume, position, H, testHessian,
                      zerothCorrections, corrections)
rk_time = time.time() - rk_time
output("rk_time")

if quitAfterTiming:
    quit()

#-------------------------------------------------------------------------------
# Try interpolation
#-------------------------------------------------------------------------------
# Get the nodes to check
if numToCheck == -1:
    nodesToCheck = list(range(nodes.numNodes))
elif numToCheck > 0:
    nodesToCheck = random.sample(list(range(nodes.numNodes)), numToCheck)
else:
    raise ValueError("numToCheck must be -1 or positive")

# Fill the FieldList we're interpolating from
fill_time = time.time()
answer_vals = dataBase.newGlobalScalarFieldList(0.0, name="initial values")
for i in range(nodes.numNodes):
    answer_vals[0][i] = func(position(0,i))
fill_time = time.time() - fill_time
output("fill_time")

interp_time = time.time()
interp_vals = interpolateRK(answer_vals, position, volume, H, connectivity, WR, corrections)
grad_vals = gradientRK(answer_vals, position, volume, H, connectivity, WR, corrections)
if testHessian:
    hess_vals = hessianRK(answer_vals, position, volume, H, connectivity, WR, corrections)
interp_time = time.time() - interp_time
output("interp_time")

# Check the results
check_time = time.time()
vals = np.zeros((nodes.numNodes, 2))
dvals = np.zeros((nodes.numNodes, dimension, 2))
ddvals = np.zeros((nodes.numNodes, dimension, dimension, 2))
for i in nodesToCheck:
    vals[i,0] = interp_vals(0,i)
    vals[i,1] = answer_vals(0,i)
    dvals[i,:,0] = grad_vals(0,i)
    dvals[i,:,1] =  dfunc(position(0,i))
    if testHessian:
        for irow in range(dimension):
            ddvals[i,irow,:,0] = hess_vals(0,i).getRow(irow)
        ddvals[i,:,:,1] = ddfunc(position(0,i))
check_time = time.time() - check_time
output("check_time")

#     xi = position(0, i)
#     fi = func(xi)
#     def addToValues(nj, j):
#         xj = position(nj, j)
#         if type(xj) is not type(xi):
#             raise TypeError, "error in xj, i = {}, j = {}".format(i, j)
#         fj = func(xj)
#         xij = xi - xj
#         vj = volume(nj, j)
#         w, dw, ddw = getKernel(ni, i, nj, j)
#         vals[i,0] += vj * w * fj
#         dvals[i,:,0] += vj * dw * fj
#         if testHessian:
#             ddvals[i,:,:,0] += vj * ddw * fj
#     connectivityi = connectivity.connectivityForNode(ni, i)
#     for nj, neighbors in enumerate(connectivityi):
#         for j in neighbors:
#             addToValues(nj, j)
#     addToValues(ni, i)
#     vals[i,1] = fi
#     dvals[i,:,1] = dfunc(xi)
#     if testHessian:
#         ddvals[i,:,:,1] = ddfunc(xi)

# interp_time = time.time() - interp_time
# output("interp_time")

#-------------------------------------------------------------------------------
# Optionally plot results
#-------------------------------------------------------------------------------
def plotThings(numLocal, anaLocal, title):
    num = mpi.gather(numLocal)
    ana = mpi.gather(anaLocal)
    if mpi.rank == 0:
        num = np.concatenate(num)
        ana = np.concatenate(ana)
        colors = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00']
        fig, ax1 = plt.subplots()
        ax1.plot(num, label="num", color=colors[0], linestyle="-")
        ax1.plot(ana, label="ana", color=colors[1], linestyle="--")
        plt.legend()
        ax2 = ax1.twinx()
        err = np.divide(num - ana, np.mean(np.abs(ana)) + 1.e-15)
        ax2.semilogy(err, label="err", color=colors[2], linestyle="-.")
        plt.title(title)
        plt.legend()
if plot:
    plotThings(vals[:,0], vals[:,1], "vals")
    for d in range(dimension):
        plotThings(dvals[:,d,0], dvals[:,d,1], "dvals{}".format(d))
    if testHessian:
        for d1 in range(dimension):
            for d2 in range(dimension):
                plotThings(ddvals[:,d1,d2,0], ddvals[:,d1,d2,1], "ddvals{}-{}".format(d1, d2))
    plt.show()
    
#-------------------------------------------------------------------------------
# Check error
#-------------------------------------------------------------------------------
def getError(num, ana):
    err = np.sum(np.abs(num - ana))
    tot = np.sum(np.abs(ana))
    err = mpi.allreduce(err, mpi.SUM)
    tot = mpi.allreduce(tot, mpi.SUM)
    return err / (tot + 1.e-11)
error = getError(vals[:,0], vals[:,1])
output("error")
derror = [getError(dvals[:,d,0], dvals[:,d,1]) for d in range(dimension)]
output("derror")
if testHessian:
    dderror = [getError(ddvals[:,d1,d2,0], ddvals[:,d1,d2,1]) for d1 in range(dimension) for d2 in range(dimension)]
    output("dderror")

if correctionOrder == RKOrder.ZerothOrder:
    ni = 0
    zerothErr = np.zeros((nodes.numNodes))
    for i in range(nodes.numNodes):
        zerothErr[i] = getError(np.array(corrections(ni, i)), np.array(zerothCorrections(ni, i)))
    output("np.amax(zerothErr)")
    if checkConditions and any([e0 > tolerance for e0 in zerothErr]):
        raise ValueError("zeroth corrections do not agree")

if checkConditions:
    if error > tolerance:
        raise ValueError("error is greater than tolerance")
    if funcType != "constant":
        if any([de > 10*tolerance for de in derror]):
            raise ValueError("gradient error is greater than tolerance")
    if testHessian and funcType != "constant" and funcType != "linear":
        if any([dde > 100*tolerance for dde in dderror]):
            raise ValueError("hessian error is greater than tolerance")
        
