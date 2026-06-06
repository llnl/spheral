#-------------------------------------------------------------------------------
# Test doing integrals using RK
#-------------------------------------------------------------------------------
from Spheral import *
from SpheralTestUtilities import *
import os, shutil
import numpy as np
import time
from matplotlib import pyplot as plt
title("Reproducing kernel integral test")

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

    # CRK options
    correctionOrder = LinearOrder,
    volumeType = RKMassOverDensity,
    needHessian = False,

    # Testing options
    randomizeNodes = False,
    ranfrac = 0.2,
    useBaseKernel = False, # Test using standard SPH kernel
    computeCorrectionsDirectly = False,
    
    # Manufactured parameters
    funcType = "linear",
    
    # Interpolation kernel choice
    nPerh = 4.01,
    hminmult = 1.e-3,
    hmaxmult = 1.e3,

    # Material parameters
    rho0 = 2.5e-7,
    gamma = 5.0/3.0,
    mu = 1.0,

    # Data dir
    dataDirBase = "dumps-RKInterpolation",

    # Plotting
    plot = False,
)

if nPerh < int(correctionOrder):
    print("nPerh is not large enough for correction order: {} < {}".format(nPerh, int(correctionOrder)))
    
if mpi.procs > 1:
    raise ValueError("parallel node generation not working")
    
#-------------------------------------------------------------------------------
# Choose correct corrections
#-------------------------------------------------------------------------------
if dimension == 1:
    from Spheral1d import *
    ZerothRKUtilities = RKUtilities1d0
    if correctionOrder == ZerothOrder:
        RKCorrections = RKCorrections1d0
        RKUtilities = RKUtilities1d0
    elif correctionOrder == LinearOrder:
        RKCorrections = RKCorrections1d1
        RKUtilities = RKUtilities1d1
    elif correctionOrder == QuadraticOrder:
        RKCorrections = RKCorrections1d2
        RKUtilities = RKUtilities1d2
    elif correctionOrder == CubicOrder:
        RKCorrections = RKCorrections1d3
        RKUtilities = RKUtilities1d3
    elif correctionOrder == QuarticOrder:
        RKCorrections = RKCorrections1d4
        RKUtilities = RKUtilities1d4
    elif correctionOrder == QuinticOrder:
        RKCorrections = RKCorrections1d5
        RKUtilities = RKUtilities1d5
    elif correctionOrder == SexticOrder:
        RKCorrections = RKCorrections1d6
        RKUtilities = RKUtilities1d6
    elif correctionOrder == SepticOrder:
        RKCorrections = RKCorrections1d7
        RKUtilities = RKUtilities1d7
    else:
        raise ValueError("correction order \"{}\" not found".format(correctionOrder))
elif dimension == 2:
    from Spheral2d import *
    ZerothRKUtilities = RKUtilities2d0
    if correctionOrder == ZerothOrder:
        RKCorrections = RKCorrections2d0
        RKUtilities = RKUtilities2d0
    elif correctionOrder == LinearOrder:
        RKCorrections = RKCorrections2d1
        RKUtilities = RKUtilities2d1
    elif correctionOrder == QuadraticOrder:
        RKCorrections = RKCorrections2d2
        RKUtilities = RKUtilities2d2
    elif correctionOrder == CubicOrder:
        RKCorrections = RKCorrections2d3
        RKUtilities = RKUtilities2d3
    elif correctionOrder == QuarticOrder:
        RKCorrections = RKCorrections2d4
        RKUtilities = RKUtilities2d4
    elif correctionOrder == QuinticOrder:
        RKCorrections = RKCorrections2d5
        RKUtilities = RKUtilities2d5
    elif correctionOrder == SexticOrder:
        RKCorrections = RKCorrections2d6
        RKUtilities = RKUtilities2d6
    elif correctionOrder == SepticOrder:
        RKCorrections = RKCorrections2d7
        RKUtilities = RKUtilities2d7
    else:
        raise ValueError("correction order \"{}\" not found".format(correctionOrder))
else:
    from Spheral3d import *
    ZerothRKUtilities = RKUtilities3d0
    if correctionOrder == ZerothOrder:
        RKCorrections = RKCorrections3d0
        RKUtilities = RKUtilities3d0
    elif correctionOrder == LinearOrder:
        RKCorrections = RKCorrections3d1
        RKUtilities = RKUtilities3d1
    elif correctionOrder == QuadraticOrder:
        RKCorrections = RKCorrections3d2
        RKUtilities = RKUtilities3d2
    elif correctionOrder == CubicOrder:
        RKCorrections = RKCorrections3d3
        RKUtilities = RKUtilities3d3
    elif correctionOrder == QuarticOrder:
        RKCorrections = RKCorrections3d4
        RKUtilities = RKUtilities3d4
    elif correctionOrder == QuinticOrder:
        RKCorrections = RKCorrections3d5
        RKUtilities = RKUtilities3d5
    elif correctionOrder == SexticOrder:
        RKCorrections = RKCorrections3d6
        RKUtilities = RKUtilities3d6
    elif correctionOrder == SepticOrder:
        RKCorrections = RKCorrections3d7
        RKUtilities = RKUtilities3d7
    else:
        raise ValueError("correction order \"{}\" not found".format(correctionOrder))

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
seed = 4898201204
random.seed(seed)

if randomizeNodes:
    print("randomizing nodes")
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
bounds = vector_of_Boundary()
method = SPHSmoothingScale()
iterateIdealH(dataBase,
              bounds,
              WT,
              method,
              100, # max h iterations
              1.e-4) # h tolerance
dataBase.updateConnectivityMap(True)

#-------------------------------------------------------------------------------
# Create RK object
#-------------------------------------------------------------------------------
rk = RKCorrections(dataBase = dataBase,
                   W = WT,
                   volumeType = volumeType,
                   needHessian = needHessian)
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
surfaceArea = state.scalarFields(HydroFieldNames.surfaceArea)
normal = state.vectorFields(HydroFieldNames.normal)
zerothCorrections = state.vector_of_doubleFields(HydroFieldNames.rkZerothCorrections)
corrections = state.vector_of_doubleFields(HydroFieldNames.rkCorrections)

#-------------------------------------------------------------------------------
# Compute corrections
#-------------------------------------------------------------------------------
if computeCorrectionsDirectly:
    rk_time = time.time()
    RKUtilities.computeCorrections(connectivity, WT, volume, position, H, needHessian,
                                   zerothCorrections, corrections)
    rk_time = time.time() - rk_time
    # RKUtilities.computeNormal(connectivity, WT, volume, position, H, corrections,
    #                           surfaceArea, normal)
    ZerothRKUtilities.computeNormal(connectivity, WT, volume, position, H, zerothCorrections,
                                    surfaceArea, normal)
else:
    rk_time = time.time()
    rk.initialize(0.0, 0.0, dataBase, state, derivs)
    rk_time = time.time() - rk_time
output("rk_time")

#-------------------------------------------------------------------------------
# Set up a simple method to calculate the kernel
#-------------------------------------------------------------------------------
# Base kernel
def getBaseKernel(ni, i, nj, j):
    xi = position(ni, i)
    xj = position(nj, j)
    xij = xi - xj
    Hj = H(nj, j)
    w = RKUtilities.evaluateBaseKernel(WT, xij, Hj)
    dw = RKUtilities.evaluateBaseGradient(WT, xij, Hj)
    if needHessian:
        ddw = RKUtilities.evaluateBaseHessian(WT, xij, Hj)
        ddw = np.reshape(ddw, (dimension, dimension))
    else:
        ddw = np.zeros((dimension, dimension))
    return w, dw, ddw

# RK kernel
def getNewKernel(ni, i, nj, j):
    xi = position(ni, i)
    xj = position(nj, j)
    xij = xi - xj
    Hj = H(nj, j)
    c = corrections(ni, i)
    w = RKUtilities.evaluateKernel(WT, xij, Hj, c)
    dw = RKUtilities.evaluateGradient(WT, xij, Hj, c)
    if needHessian:
        ddw = RKUtilities.evaluateHessian(WT, xij, Hj, c)
        ddw = np.reshape(ddw, (dimension, dimension))
    else:
        ddw = np.zeros((dimension, dimension))
    return w, dw, ddw

# Kernel to test with
getKernel = getBaseKernel if useBaseKernel else getNewKernel

#-------------------------------------------------------------------------------
# Get a simple function for testing surface integrals
#-------------------------------------------------------------------------------
def func(x):
    return 4 + 0.5 * x[0] * x[0]

#-------------------------------------------------------------------------------
# Check whether normals I'm calculating make any sense
#-------------------------------------------------------------------------------
# Get expected surface area
expectedNS = np.zeros((dimension))
if dimension == 1:
    expectedS = 2.0
elif dimension == 2:
    expectedS = 2*(x1 - x0 + y1 - y0)
else:
    expectedS = 2*(x1-x0)*(y1-y0) + 2*(x1-x0)*(z1-z0) + 2*(y1-y0)*(z1-z0)

# Get total surface area
totalNS = np.zeros((dimension))
totalS = 0.0
ni = 0
for i in range(nodes.numNodes):
    xi = position(ni, i)
    normali = normal(ni, i)
    si = surfaceArea(ni, i)
    totalNS += normali * si
    totalS += si
    print(xi, si, normali)

output("expectedNS")
output("totalNS")
output("expectedS")
output("totalS")

if dimension == 1:
    ni = 0
    xmid = 0.5 * (x0 + x1)
    totalBNS = np.zeros((2))
    for i in range(nodes.numNodes):
        xi = position(ni, i)
        normali = normal(ni, i)
        si = surfaceArea(ni, i)
        if xi[0] < xmid:
            totalBNS[0] += normali[0] * si
        else:
            totalBNS[1] += normali[0] * si
    output("totalBNS")



# #-------------------------------------------------------------------------------
# # Test out Dilts method for surface integrals
# #-------------------------------------------------------------------------------
# # Check whether point is on boundary
# def checkOnBoundary(x, h):
#     radius = 1. / h
#     for d in range(dimension):
#         for k in range(2):
#             if np.abs(x[d] - lims[d][k]) <= radius:
#                 return True
#     return False
# # Get the +- direction for each component of the vector
# def quadDir(x):
#     return x / np.abs(x)

# # Check which points are on the boundary
# onBoundary = np.zeros((nodes.numNodes), dtype=bool)
# ni = 0
# for i in range(nodes.numNodes):
#     xi = position(ni, i)
#     Hi = H(ni, i)
#     onBoundary[i] = checkOnBoundary(xi, Hi.Determinant())
# ni = nb = 0
# for i in range(nodes.numNodes):
#     if onBoundary[i]:
#         normali = np.zeros((dimension))
#         connectivityi = connectivity.connectivityForNode(ni, i)
#         xi = position(ni, i)
#         vi = volume(ni, i)
#         pii, dpii, ddpii = getKernel(ni, i, ni, i)
#         Aii = vi * vi * dpii
#         Bij = 2 * Aii
#         normali += Bij
#         for nj, neighbors in enumerate(connectivityi):
#             for j in neighbors:
#                 if onBoundary[j]:
#                     vj = volume(nj, j)
#                     pji, dpji, ddpji = getKernel(ni, i, nj, j)
#                     pij, dpij, ddpij = getKernel(nj, j, ni, i)
#                     Aij = vi * vj * dpji
#                     Aji = vi * vj * dpij
#                     Bij = Aij + Aji
#                     normali += Bij
#         # print xi, normali / np.linalg.norm(normali), quadDir(xi) - quadDir(normali)

# #-------------------------------------------------------------------------------
# # Test out surface integral for transport
# #-------------------------------------------------------------------------------
# ni = nb = 0
# for i in range(nodes.numNodes):
#     xi = position(ni, i)
#     b = 0 if xi[0] < 0 else nodes.numNodes - 1
#     normb = -1 if xi[0] < 0 else 1
#     xb = position(nb, b)
#     pib, dpib, ddpib = getKernel(nb, b, ni, i)
#     fb = func(xb)
#     num1 = pib * fb
#     num2 = 0.
#     connectivityi = connectivity.connectivityForNode(ni, i)
#     for nj, neighbors in enumerate(connectivityi):
#         for j in neighbors:
#             xj = position(nj, j)
#             vj = volume(nj, j)
#             pij, dpij, ddpij = getKernel(nj, j, ni, i)
#             pji, dpji, ddpji = getKernel(ni, i, nj, j)
#             fj = func(xj)
#             num2 += vj * (dpij[0] + dpji[0]) * fj
#     vi = volume(ni, i)
#     pii, dpii, ddpii = getKernel(ni, i, ni, i)
#     fi = func(xi)
#     num2 += vi * (dpii[0] + dpii[0]) * fi
#     # if num1 != 0:
#     print normb * num1, num2#, num2 / num1

#-------------------------------------------------------------------------------
# Test out surface integrals
#-------------------------------------------------------------------------------
# def getSurfaceIntegral(i, j):
#     ni = nj = 0
#     xi = position(ni, i)
#     xj = position(nj, j)
#     xij = xi - xj
#     vi = volume(ni, i)
#     vj = volume(nj, j)
#     wj, dwj, ddwj = getKernel(ni, i, nj, j)
#     wi, dwi, ddwi = getKernel(nj, j, ni, i)
#     # ns = 2 * (dwj + dwi) / (wj / vi + wi / vj)
#     ns = (dwj + dwi) / (wi * wj)
#     nsnorm = np.linalg.norm(ns)
#     n = ns / nsnorm if nsnorm > 0 else np.zeros_like(ns)
#     print xi, xj, n
#     return n

# norm = np.zeros(nodes.numNodes)
# test = np.zeros(nodes.numNodes)
# def testSurfaceIntegral1d(i):
#     ni = nj = 0
#     xi = position(ni, i)
#     vi = volume(ni, i)
#     normb = 1 if xi[0] > 0 else -1
#     nb = nodes.numNodes - 1 if xi[0] > 0 else 0
#     wib, dwib, ddwib = getKernel(ni, i, nj, nb)
#     connectivityi = connectivity.connectivityForNode(ni, i)
#     def addThing(ni, i, nj, j):
#         vj = volume(nj, j)
#         wij, dwij, ddwij = getKernel(ni, i, nj, j)
#         wji, dwji, ddwji = getKernel(nj, j, ni, i)
#         norm[i] += vj * (dwij[0] + dwji[0])
#         test[i] += vj * wij
#     for nj, neighbors in enumerate(connectivityi):
#         for j in neighbors:
#             addThing(ni, i, nj, j)
#     addThing(ni, i, ni, i)
#     # print test[i]
#     if wib > 0:
#         norm[i] = norm[i] / wib
#         print xi[0], norm[i]

# norm = np.zeros(2)
# denom = np.zeros(2)
# normi = np.zeros((nodes.numNodes))
# denomi = np.zeros((nodes.numNodes))
# normi2 = np.zeros((nodes.numNodes))
# denomi2 = np.zeros((nodes.numNodes))
# def testOnePointSurface1d(i, j):
#     ni = nj = nb = 0
#     vi = volume(ni, i)
#     vj = volume(nj, j)
#     xi = position(ni, i)
#     xj = position(nj, j)
#     xmid = 0.5 * (xi[0] + xj[0])
#     b = 0 if xmid < 0. else nodes.numNodes - 1
#     wji, dwji, ddwji = getKernel(nj, j, ni, i)
#     wij, dwij, ddwij = getKernel(ni, i, nj, j)
#     wbi, dwbi, ddwbi = getKernel(nb, b, ni, i)
#     wbj, dwbj, ddwbj = getKernel(nb, b, nj, j)
#     if wbi != 0 and wbj != 0:
#         normi[i] += vj * (dwji[0] + dwij[0])
#         denomi[i] += vj * wbi * wbj
#         normi2[i] += vj * dwji[0]
#         denomi2[i] = wbi
#         if xmid < 0.:
#             normb = -1.
#             norm[0] += vi * vj * (dwji[0] + dwij[0])
#             denom[0] += vi * vj * wbi * wbj
#             # print "left", dwij[0], dwji[0]
#         else:
#             normb = 1.
#             norm[1] += vi * vj * (dwji[0] + dwij[0])
#             denom[1] += vi * vj * wbi * wbj
#             # print "right", dwij[0], dwji[0]
#         # print normb * wib * wjb, dwij[0] + dwji[0]
        
# ni = 0
# for i in range(nodes.numNodes):
#     connectivityi = connectivity.connectivityForNode(ni, i)
#     for nj, neighbors in enumerate(connectivityi):
#         for j in neighbors:
#             testOnePointSurface1d(i, j)
#     testOnePointSurface1d(i, i)
#     # testSurfaceIntegral1d(i)
#     # connectivityi = connectivity.connectivityForNode(ni, i)
#     # for nj, neighbors in enumerate(connectivityi):
#     #     for j in neighbors:
#     #         getSurfaceIntegral(i, j)
# # normi = normi / (denomi + 1.e-15)
# # normi2 = normi2 / (denomi2 + 1.e-15)
# print norm
# print normi, denomi
# print normi2, denomi2
#-------------------------------------------------------------------------------
# Quit, as I wanted to keep below code for reference but it isn't used
#-------------------------------------------------------------------------------
quit()

#-------------------------------------------------------------------------------
# Try interpolation
#-------------------------------------------------------------------------------
interp_time = time.time()
vals = np.zeros((nodes.numNodes, 2))
dvals = np.zeros((nodes.numNodes, dimension, 2))
ddvals = np.zeros((nodes.numNodes, dimension, dimension, 2))
ni = 0
for i in range(nodes.numNodes):
    xi = position(ni, i)
    fi = func(xi)
    def addToValues(nj, j):
        xj = position(nj, j)
        if type(xj) is not type(xi):
            raise TypeError("error in xj, i = {}, j = {}".format(i, j))
        fj = func(xj)
        vj = volume(nj, j)
        w, dw, ddw = getKernel(ni, i, nj, j)
        vals[i,0] += vj * w * fj
        dvals[i,:,0] += vj * dw * fj
        if needHessian:
            ddvals[i,:,:,0] += vj * ddw * fj
    connectivityi = connectivity.connectivityForNode(ni, i)
    for nj, neighbors in enumerate(connectivityi):
        for j in neighbors:
            addToValues(nj, j)
    addToValues(ni, i)
    vals[i,1] = fi
    dvals[i,:,1] = dfunc(xi)
    if needHessian:
        ddvals[i,:,:,1] = ddfunc(xi)

interp_time = time.time() - interp_time
output("interp_time")

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
    if needHessian:
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
if needHessian:
    dderror = [getError(ddvals[:,d1,d2,0], ddvals[:,d1,d2,1]) for d1 in range(dimension) for d2 in range(dimension)]
    output("dderror")
    
if error > tolerance:
    raise ValueError("error is greater than tolerance")
if funcType != "constant":
    if any([de > tolerance for de in derror]):
        raise ValueError("gradient error is greater than tolerance")
if needHessian and funcType != "constant" and funcType != "linear":
    if any([dde > tolerance for dde in dderror]):
        raise ValueError("hessian error is greater than tolerance")
        
