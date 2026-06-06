#ATS:t1 = test(SELF, "--dimension 1 --order 100 --tolerance 3.0e-5", label="integration, 1d", np=1)
#ATS:t2 = test(SELF, "--dimension 2 --nx 10 --ny 10 --order 10 --tolerance 2.0e-5", label="integration, 2d", np=1)
#ATS:t3 = test(SELF, "--dimension 3 --nx 5 --ny 5 --nz 5 --order 6 --tolerance 1.0e-5", label="integration, 3d", np=1)
#ATS:r1 = test(SELF, "--dimension 1 --nx 20 --order 100 --correctionOrderIntegration 1", label="integration, 1d, rk1", np=1)
#ATS:r1 = test(SELF, "--dimension 1 --nx 20 --nPerh 10.01 --order 100 --correctionOrderIntegration 4", label="integration, 1d, rk4", np=1)
#ATS:r2 = test(SELF, "--dimension 2 --nx 20 --ny 20 --order 10 --correctionOrderIntegration 1", label="integration, 2d, rk1", np=1)
#ATS:t3 = test(SELF, "--dimension 3 --nx 5 --ny 5 --nz 5 --order 6 --correctionOrderIntegration 1", label="integration, 3d, rk1", np=1)

#-------------------------------------------------------------------------------
# Test the kernel integration methods
#-------------------------------------------------------------------------------
from Spheral import *
from SpheralTestUtilities import *
import os, shutil
import numpy as np
import time
from matplotlib import pyplot as plt
title("Kernel integration test")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(
    # Spatial parameters
    dimension = 1,
    nx = 20,
    ny = 20,
    nz = 20,
    x0 = -2.0,
    x1 = 2.0,
    y0 = -2.0,
    y1 = 2.0,
    z0 = -2.0,
    z1 = 2.0,

    # RK options
    useRK = False,
    correctionOrderIntegration = -1,
    volumeType = RKVoronoiVolume,

    # Integration parameters
    clipBoundaries = True,
    
    # Testing options
    randomizeNodes = False,
    ranfrac = 0.2,
    testHessian = False,
    printErrors = False,
    outputAllSurfaceVolume = False,
    checkAllSurfaceVolume = False,
    
    # Manufactured parameters
    funcType = "linear",
    
    # Integrator parameters
    order = 4,
    hMult = 1.0,
    
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
    dataDirBase = "dumps-KernelIntegration",

    # Error
    tolerance = 1.e-4,
    
    # Plotting
    plot = False,
)

if nPerh < correctionOrderIntegration:
    print("nPerh is not large enough for correction order: {} < {}".format(nPerh, correctionOrderIntegration))
    
if mpi.procs > 1:
    raise ValueError("test is not written for parallel")
    
#-------------------------------------------------------------------------------
# Choose correct dimension aliases
#-------------------------------------------------------------------------------
exec("from Spheral%id import *" % dimension, globals())

#-------------------------------------------------------------------------------
# Set up data
#-------------------------------------------------------------------------------
# Set limits
lims = [[x0, x1], [y0, y1], [z0, z1]]
delta = [(x1 - x0)/nx, (y1 - y0) / ny, (z1 - z0) / nz]
length = [x1 - x0, y1 - y0, z1 - z0]
midpoints = [0.5*(x0+x1), 0.5*(y0+y1), 0.5*(z0+z1)]
analytic_volume = 1.
for d in range(dimension):
    analytic_volume *= lims[d][1] - lims[d][0]
if dimension == 1:
    analytic_surface = 2.
elif dimension == 2:
    analytic_surface = 2 * (length[0] + length[1])
elif dimension == 3:
    analytic_surface = 2 * (length[0] * length[1] +
                            length[0] * length[2] +
                            length[1] * length[2])
# Get point spacing
units = MKS()
deltaMin = min(delta[:dimension])
deltaMax = max(delta[:dimension])
hmin = deltaMin * nPerh * hminmult
hmax = deltaMax * nPerh * hmaxmult

output("delta")

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
dataDir = os.path.join(dataDirBase,
                       "correctionOrder={}".format(correctionOrderIntegration),
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
WT = TableKernel(WendlandC2Kernel(), 1000)
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
seed = 4587592729430
random.seed(seed)

if randomizeNodes:
    dx = delta[0]
    dy = delta[1]
    dz = delta[2]
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
              [method],
              [],
              100, # max h iterations
              1.e-4) # h tolerance
dataBase.updateConnectivityMap(True, False) # need ghost and overlap connectivity

#-------------------------------------------------------------------------------
# Create RK object
#-------------------------------------------------------------------------------
correctionOrder = LinearOrder # We don't actually use this
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
# Add faceted boundaries
#-------------------------------------------------------------------------------
if clipBoundaries:
    if dimension == 1:
        points = vector_of_Vector([Vector(x0), Vector(x1)])
        facetedBoundary = Box1d(points)
    elif dimension == 2:
        points = vector_of_Vector([Vector(x0, y0), Vector(x0, y1),
                                   Vector(x1, y0), Vector(x1, y1)])
        facetedBoundary = Polygon(points)
    else:
        points = vector_of_Vector([Vector(x0, y0, z0), Vector(x0, y0, z1),
                                   Vector(x0, y1, z0), Vector(x0, y1, z1),
                                   Vector(x1, y0, z0), Vector(x1, y0, z1),
                                   Vector(x1, y1, z0), Vector(x1, y1, z1)])
        facetedBoundary = Polyhedron(points)
    output("facetedBoundary")
    output("facetedBoundary.volume")
    rk.addFacetedBoundary(facetedBoundary)

#-------------------------------------------------------------------------------
# Create a state directly and initialize physics package
#-------------------------------------------------------------------------------
connectivity = dataBase.connectivityMap(True, False)
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
cells = state.facetedVolumeFields(HydroFieldNames.cells)
flags = state.vector_of_CellFaceFlagFields(HydroFieldNames.cellFaceFlags)
# normal = state.vectorFields(HydroFieldNames.normal)

#-------------------------------------------------------------------------------
# Compute corrections
#-------------------------------------------------------------------------------
WR.computeCorrections(connectivity, volume, position, H, testHessian,
                      zerothCorrections, corrections)

#-------------------------------------------------------------------------------
# Create flat connectivity
#-------------------------------------------------------------------------------
connectivity_time = time.time()
flatConnectivity = FlatConnectivity()
flatConnectivity.computeIndices(dataBase)
flatConnectivity.computeSurfaceIndices(dataBase, state)
connectivity_time = time.time() - connectivity_time
output("connectivity_time")

#-------------------------------------------------------------------------------
# Create integrator, add integrals, and perform integration
#-------------------------------------------------------------------------------
if correctionOrderIntegration < 0:
    integrationKernel = SPHIntegrationKernel(WT)
else:
    integrationKernel = eval("RKIntegrationKernel{}d{}(WT)".format(dimension, correctionOrderIntegration))
output("integrationKernel")
integrator = KernelIntegrator(order, integrationKernel, dataBase, flatConnectivity)
vlK_f = LinearKernel()
vlG_f = LinearGrad() # This and slKn_f should be equal
vbKK_f = BilinearKernelKernel()
vbGK_f = BilinearGradKernel()
vbKG_f = BilinearKernelGrad()
vbGdG_f = BilinearGradDotGrad()
vbGpG_f = BilinearGradProdGrad()
slKn_f = LinearSurfaceNormalKernel()
sbKKn_f = BilinearSurfaceNormalKernelKernel() # This and sbKKn2_f should be equal
sbKGdn_f = BilinearSurfaceNormalKernelDotGrad()
sbKKn2_f = BilinearSurfaceNormalKernelKernelFromGrad()
vcc_f = CellCoefficient() # Calculate the volume directly
scn_f = SurfaceNormalCoefficient() # Calculate the surface area directly
volumeIntegrals = [vlK_f, vlG_f, vbKK_f, vbGK_f, vbKG_f, vbGdG_f, vbGpG_f, sbKKn2_f, vcc_f]
surfaceIntegrals = [slKn_f, sbKKn_f, sbKGdn_f, scn_f]
integrals = volumeIntegrals + surfaceIntegrals

for integral in integrals:
    integrator.addIntegral(integral)

integrator.setState(0.0, # time
                    state)

integration_time = time.time()
integrator.performIntegration()
integration_time = time.time() - integration_time
output("integration_time")
output("integrator.totalNumSubcells()")
output("integrator.totalNumSubfacets()")

# Get the integrals
vlK = list(vlK_f.values())
vlG = list(vlG_f.values())
vbKK = list(vbKK_f.values())
vbGK = list(vbGK_f.values())
vbKG = list(vbKG_f.values())
vbGdG = list(vbGdG_f.values())
vbGpG = list(vbGpG_f.values())
slKn = list(slKn_f.values())
sbKKn = list(sbKKn_f.values())
sbKGdn = list(sbKGdn_f.values())
sbKKn2 = list(sbKKn2_f.values())
slKn2 = list(vlG_f.values()) # The surface integral is the linear integral for constant coefficients
vcc = list(vcc_f.values())
scn = list(scn_f.values())

#-------------------------------------------------------------------------------
# Verify that volumes calculated in different ways are equal:
# 1. Analytic volume
# 2. From the Voronoi cell volumes
# 3. From the calculated Spheral volumes
# 4. From the integral package
#-------------------------------------------------------------------------------
checksum = 0

totcell = 0.
totvol = 0.
totint = 0.
for i in range(nodes.numNodes):
    totcell += cells[0][i].volume
    totvol += volume[0][i]
    totint += vcc[i]
    # print cells[0][i].xmin.x, cells[0][i].xmax.x
volumes = [analytic_volume, totcell, totvol, totint]
output("volumes")
for v in volumes:
    if np.abs(v - analytic_volume) > tolerance:
        print("volumes not correct")
        checksum += 1

#-------------------------------------------------------------------------------
# Verify the areas in the same way:
# 1. Analytic surface area
# 2. Integrated surface area
#-------------------------------------------------------------------------------
totarea = 0.
for i in range(nodes.numNodes):
    totarea += np.sum(np.abs(scn[i]))
print("areas: ", totarea, analytic_surface)
if np.abs(totarea - analytic_surface) > tolerance:
    print("areas not correct")
    checksum += 1

#-------------------------------------------------------------------------------
# Check numerical integration. These integrals are calculated in the Mathematica
# notebook TestIntegrator.nb and only apply to specific cases due to dependence
# on xi/j and Hi/j. If these are changed, the tests need to be updated.
#-------------------------------------------------------------------------------
nPerhTest = 4.01
if (nx == 20) and (dimension == 1) and (not useRK) and (nPerh == nPerhTest) and (not randomizeNodes) and (correctionOrderIntegration < 0):
    indi = 10
    indj = 11
    indij = flatConnectivity.localToFlat(indi, indj)
    vals = [["xi", position(0, indi).x, 0.1],
            ["xj", position(0, indj).x, 0.3],
            ["Hi", H(0, indi).xx, 1.2468380523035534],
            ["Hj", H(0, indj).xx, 1.2468380523035534],
            ["vlK",  vlK[indi], 1.0],
            ["vlG",  vlG[indi].x, 0.0],
            ["vbKK",  vbKK[indi][indij], 1.21521457112],
            ["vbGK",  vbGK[indi][indij].x, -7.4859553716],
            ["vbKG",  vbKG[indi][indij].x, 7.4859553716],
            ["vbGdG",  vbGdG[indi][indij], -5.82735340663],
            ["vbGpG",  vbGpG[indi][indij].xx, -5.82735340663]]
    print("i = {}, j = {}".format(indi, indj))
    print("\tdelta: ", delta[0])
    for val in vals:
        err = val[1] - val[2]
        print("\t{}\t{}\t{}\t{}".format(val[0], val[1], val[2], err))
        if np.abs(err) > tolerance:
            print("\ttolerance fail")
            checksum += 1

    indi = 0
    indj = 1
    indij = flatConnectivity.localToFlat(indi, indj)
    numSurfaces = flatConnectivity.numSurfaces(indi)
    print("i = {}, j = {}".format(indi, indj))
    print("\tdelta: ", 2*delta[0])
    vals = [["xi", position(0, indi).x, -1.9],
            ["xj", position(0, indj).x, -1.7],
            ["Hi", H(0, indi).xx, 0.6538380071103822],
            ["Hj", H(0, indj).xx, 0.7890618854483368],
            ["slKn1",  slKn[indi][0].x, -1.49469156773],
            ["slKn2",  slKn[indj][0].x, -0.697064575841],
            ["slKKn",  sbKKn[indi][0 + numSurfaces * indij].x, -1.04189654368],
            ["vlK1",  vlK[indi], 0.658571492971],
            ["vlK2",  vlK[indj], 0.934263499614],
            ["vlG1",  vlG[indi].x, -1.49469156773],
            ["vlG2",  vlG[indj].x, -0.697064575841],
            ["vbKK",  vbKK[indi][indij], 0.962348197392],
            ["vbGK",  vbGK[indi][indij].x, -2.26201532957],
            ["vbKG",  vbKG[indi][indij].x, 1.22011878589],
            ["vbGdG",  vbGdG[indi][indij], 4.06549020154],
            ["vbGpG",  vbGpG[indi][indij].xx, 4.06549020154]]
    for val in vals:
        err = val[1] - val[2]
        print("\t{}\t{}\t{}\t{}".format(val[0], val[1], val[2], err))
        if np.abs(err) > tolerance:
            print("tolerance fail")
            checksum += 1

if (nx == 10) and (ny == 10) and (dimension == 2) and (not useRK) and (nPerh == nPerhTest) and (not randomizeNodes) and (correctionOrderIntegration < 0):
    indi = 5
    indj = 14
    print("i = {}, j = {}".format(indi, indj))
    indij = flatConnectivity.localToFlat(indi, indj)
    normali = Vector(0.0, -1.0)
    inds = flatConnectivity.surfaceIndex(indi, normali)
    output("inds")
    numSurfaces = flatConnectivity.numSurfaces(indi)
    vals = [["xix", position(0, indi).x, 0.2],
            ["xiy", position(0, indi).y, -1.8],
            ["xjx", position(0, indj).x, -0.2],
            ["xjy", position(0, indj).y, -1.4],
            ["Hixx", H(0, indi).xx, 0.474869572654424],
            ["Hixy", H(0, indi).xy, 0.0],
            ["Hiyy", H(0, indi).yy, 0.474869572654424],
            ["Hjxx", H(0, indj).xx, 0.5625194593752989],
            ["Hjxy", H(0, indj).xy, 0.0],
            ["Hjyy", H(0, indj).yy, 0.5625194593752989],
            ["slKn1x",  slKn[indi][0].x, 0.0],
            ["slKn1y",  slKn[indi][0].y, -1.0962491605],
            ["slKn2x",  slKn[indj][0].x, 0.0],
            ["slKn2y",  slKn[indj][0].y, -0.0593986825441],
            ["slKKnx",  sbKKn[indi][0 + numSurfaces * indij].x, 0.0],
            ["slKKny",  sbKKn[indi][0 + numSurfaces * indij].y, -0.0390433300028],
            ["vlK1",  vlK[indi], 0.759716083761],
            ["vlK2",  vlK[indj], 0.996634444907],
            ["vlG1x",  vlG[indi].x, 0.0],
            ["vlG1y",  vlG[indi].y, -1.09624915461],
            ["vlG2x",  vlG[indj].x, 0.0],
            ["vlG2y",  vlG[indj].y, -0.059398683218],
            ["vbKK",  vbKK[indi][indij], 0.366764699432],
            ["vbGKx",  vbGK[indi][indij].x, 1.06895893859],
            ["vbGKy",  vbGK[indi][indij].y, -1.08118265739],
            ["vbKGx",  vbKG[indi][indij].x, -1.06895893419],
            ["vbKGy",  vbKG[indi][indij].y, 1.04213935049],
            ["vbGdG",  vbGdG[indi][indij], -0.758141358187],
            ["vbGpGxx",  vbGpG[indi][indij].xx, -0.326603680722],
            ["vbGpGxy",  vbGpG[indi][indij].xy, 2.91099707839],
            ["vbGpGyx",  vbGpG[indi][indij].yx, 3.03962223718],
            ["vbGpGyy",  vbGpG[indi][indij].yy, -0.431537657039]]
    for val in vals:
        err = val[1] - val[2]
        print("\t{}\t{}\t{}\t{}".format(val[0], val[1], val[2], err))
        if np.abs(err) > tolerance:
            print("tolerance fail")
            checksum += 1

if (nx == 5) and (ny == 5) and (nz == 5) and (dimension == 3) and (not useRK) and (nPerh == nPerhTest) and (not randomizeNodes) and (correctionOrderIntegration < 0):
    indi = 30
    indj = 31
    print("i = {}, j = {}".format(indi, indj))
    indij = flatConnectivity.localToFlat(indi, indj)
    normali1 = Vector(-1.0, 0.0, 0.0)
    normali2 = Vector(0.0, -1.0, 0.0)
    normali3 = Vector(0.0, 0.0, -1.0)
    inds1 = flatConnectivity.surfaceIndex(indi, normali1)
    inds2 = flatConnectivity.surfaceIndex(indi, normali2)
    inds3 = flatConnectivity.surfaceIndex(indi, normali3)
    numSurfaces = flatConnectivity.numSurfaces(indi)
    vals =  [["xix", position(0, indi).x, -1.6],
             ["xiy", position(0, indi).y, -0.8],
             ["xiz", position(0, indi).z, -0.8],
             ["xjx", position(0, indj).x, -0.8],
             ["xjy", position(0, indj).y, -0.8],
             ["xjz", position(0, indj).z, -0.8],             
             ["Hixx", H(0, indi).xx, 0.2193721769093367],
             ["Hixy", H(0, indi).xy, 0.0],
             ["Hixz", H(0, indi).xz, 0.0],
             ["Hiyy", H(0, indi).yy, 0.2193721769093367],
             ["Hiyz", H(0, indi).yz, 0.0],
             ["Hizz", H(0, indi).zz, 0.2193721769093367],
             ["Hjxx", H(0, indj).xx, 0.25496138893896947],
             ["Hjxy", H(0, indj).xy, 0.0],
             ["Hjxz", H(0, indj).xz, 0.0],
             ["Hjyy", H(0, indj).yy, 0.25496138893896947],
             ["Hjyz", H(0, indj).yz, 0.0],
             ["Hjzz", H(0, indj).zz, 0.25496138893896947],
             ["slKn1x",  slKn[indi][inds1].x, -0.510394848431],
             ["slKn2y",  slKn[indi][inds2].y, -0.0691607900558],
             ["slKn3z",  slKn[indi][inds3].z, -0.0691607898569],
             ["slKKn1x",  sbKKn[indi][inds1 + numSurfaces * indij].x, -0.00700584509725],   # These seem a bit flipped from Mathematica script
             ["slKKn2y",  sbKKn[indi][inds2 + numSurfaces * indij].y, -0.000747543627853],
             ["slKKn3z",  sbKKn[indi][inds3 + numSurfaces * indij].z, -0.000747543627853],
             ["vlK1",  vlK[indi], 0.715947297279],
             ["vlK2",  vlK[indj], 0.979715710562],
             ["vlG1x",  vlG[indi].x, -0.510394848431],
             ["vlG1y",  vlG[indi].y, -0.0691607900558],
             ["vlG1z",  vlG[indi].z, -0.0691607898569],
             ["vbKK",  vbKK[indi][indij], 0.0760195704411],
             ["vbGKx",  vbGK[indi][indij].x, -0.0942685322422],
             ["vbGKy",  vbGK[indi][indij].y, -0.000287654687878],
             ["vbGKz",  vbGK[indi][indij].y, -0.000287654735128],
             ["vbKGx",  vbKG[indi][indij].x, 0.0872626882115],
             ["vbKGy",  vbKG[indi][indij].y, -0.000459888928958],
             ["vbKGy",  vbKG[indi][indij].z, -0.00045988895322],
             ["vbGdG",  vbGdG[indi][indij], 0.226244752765],
             ["vbGpGxx",  vbGpG[indi][indij].xx, -0.0000209266431491],
             ["vbGpGxy",  vbGpG[indi][indij].xy, 0.00071708214423],
             ["vbGpGxz",  vbGpG[indi][indij].xz, 0.000717093764359],
             ["vbGpGyx",  vbGpG[indi][indij].yx, -0.000423206021051],
             ["vbGpGyy",  vbGpG[indi][indij].yy, 0.113132837973],
             ["vbGpGyz",  vbGpG[indi][indij].yz, 1.70817198392e-7],
             ["vbGpGzx",  vbGpG[indi][indij].zx, -0.000423206014427],
             ["vbGpGzy",  vbGpG[indi][indij].zy, 1.70817198392e-7],
             ["vbGpGzy",  vbGpG[indi][indij].zz, 0.113132839503]]
    for val in vals:
        err = val[1] - val[2]
        print("\t{}\t{}\t{}\t{}".format(val[0], val[1], val[2], err))
        if np.abs(err) > tolerance:
            print("tolerance fail")
            checksum += 1
            
#-------------------------------------------------------------------------------
# Check whether surface and volume integrals agree for the integrals that can
# be written either way.
# - Bilinear surface integral:
#   \int_{S}n^{\alpha}u_{i}u_{j}=\int_{V}\partial^{\alpha}u_{i}u_{j}+\int_{V}u_{i}\partial^{\alpha}u_{j}
# - Linear surface integral:
#   \int_{S}n^{\alpha}u_{i}=\int_{V}\partial^{\alpha}u_{i}
#-------------------------------------------------------------------------------
print("surface-volume equivalence")
av_neighbors = 0.
av_surfaces = 0.
num_overlaps = 0
lin_err = [0., 0.]
bil_err = [0., 0.]
lin_sum = [0., 0.]
bil_sum = [0., 0.]
sv_time = time.time()
if checkAllSurfaceVolume:
    num_tested = nodes.numNodes
    ivals = list(range(nodes.numNodes))
else:
    num_tested = 20
    ivals = np.random.choice(list(range(nodes.numNodes)), size=num_tested)
    # num_tested = 2
    # ivals = [0, 1]

for i in ivals:
    numElements = flatConnectivity.numNeighbors(i)
    numSurfaces = flatConnectivity.numSurfaces(i)
    av_neighbors += numElements
    av_surfaces += numSurfaces
    sbKKni = sbKKn[i]
    sbKKn2i = sbKKn2[i]
    # print i, numElements, numSurfaces, position(0, i), flags(0, i), flatConnectivity.normal([0, i])
    # print sbKKni
    for j in range(numElements):
        t1 = sbKKn2i[j]
        t2 = Vector()
        for s in range(numSurfaces):
            t2 += sbKKni[s + numSurfaces * j]
        bil_err[0] += (t2 - t1).magnitude()
        bil_err[1] += (t2 - t1).magnitude2()
        bil_sum[0] += t2.magnitude()
        bil_sum[1] += t2.magnitude2()
        if outputAllSurfaceVolume:
            print("\tsbKKn\t{}\t{}\t{}\t{}\t{}".format(i, j, t1, t2, t2 - t1))
        num_overlaps += 1
    t1 = slKn2[i]
    t2 = Vector()
    slKni = slKn[i]
    for s in range(numSurfaces):
        t2 += slKni[s]
    lin_err[0] += (t2 - t1).magnitude()
    lin_err[1] += (t2 - t1).magnitude2()
    lin_sum[0] += t2.magnitude()
    lin_sum[1] += t2.magnitude2()
    if outputAllSurfaceVolume:
        print("\tslKn\ti\t{}\t{}\t{}".format(i, t1, t2, t2 - t1))
sv_time = time.time() - sv_time
av_neighbors /= num_tested
av_surfaces /= num_tested
for i in range(2):
    bil_err[i] /= (bil_sum[i] + 1.e-16)
    lin_err[i] /= (lin_sum[i] + 1.e-16)
bil_err[1] = np.sqrt(bil_err[1])
lin_err[1] = np.sqrt(lin_err[1])
output("av_neighbors")
output("av_surfaces")
for err in bil_err:
    if err > tolerance * 10:
        checksum += 1
        print("bilinear error too high")
output("bil_err")
for err in lin_err:
    if err > tolerance * 10:
        checksum += 1
        print("linear error too high")
output("lin_err")
output("sv_time")
    
#-------------------------------------------------------------------------------
# Output how many things have failed
#-------------------------------------------------------------------------------
output("checksum")
if checksum > 0:
    raise ValueError("too many errors")
