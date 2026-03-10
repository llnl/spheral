#-------------------------------------------------------------------------------
# A set of tests to compare how the CRK formalism compares to the RK formalisml; The former being inconsistent and conservative, the latter being consistent and non-conservative.
#-------------------------------------------------------------------------------
from Spheral import *
from SpheralTestUtilities import *

import os, sys
import shutil
from SolidSpheral1d import *
from SpheralTestUtilities import *



title("Interpolation tests")

#-------------------------------------------------------------------------------
# Generic problem parameters
#-------------------------------------------------------------------------------
commandLine(
    # Parameters for seeding nodes.
    nx1 = 50,     
    nx2 = 50,
    rho1 = 1.0,
    rho2 = 1.0,
    eps1 = 0.0,
    eps2 = 0.0,
    x0 = 0.0,
    x1 = 0.5,
    x2 = 1.0,
    nPerh = 2.01,
    hmin = 0.0001, 
    hmax = 10.0,

    # What order of reproducing kernel should we use (0,1,2)?
    correctionOrder = LinearOrder,

    #Kernel Type
    KernelConstructor = NBSplineKernel,
    order = 7,
    
    # Should we randomly perturb the positions?
    ranfrac = 0.2,
    seed = 14892042,

    # What test problem are we doing?
    testDim = "1d",
    testCase = "linear",

    # The fields we're going to interpolate.
    # Linear coefficients: y = y0 + m0*x
    y0 = 1.0,
    m0 = 1.0,

    # Quadratic coefficients: y = y2 + m2*x^2
    y2 = 1.0,
    m2 = 0.5,

    # Cubic coefficients: y = y3 + m3*x^3
    y3 = 1.0,
    m3 = 3.0,

    gamma = 5.0/3.0,
    mu = 1.0,

    # Parameters for iterating H.
    iterateH = True,
    maxHIterations = 200,
    Htolerance = 1.0e-4,

    # Parameters for passing the test
    interpolationTolerance = 5.0e-7,
    derivativeTolerance = 5.0e-5,

    graphics = True,
    graphBij = False,
    plotKernels = False,
    outputFile = None,
    plotSPH = True,
    dataCut = False,
    dataCutMin = 0.0,
    dataCutMax = 1.0,
)

assert testCase in ("linear", "quadratic", "cubic", "step")
assert testDim in ("1d", "2d", "3d")
assert correctionOrder == LinearOrder
if not dataCut:
 dataCutMax = x2
 dataCutMin = x0
  
nDim = 1
if testDim == "2d":
  nDim = 2
elif testDim == "3d":
  nDim = 3

FacetedVolume = {"1d" : Box1d,
                 "2d" : Polygon,
                 "3d" : Polyhedron}[testDim]

#-------------------------------------------------------------------------------
# Appropriately set generic object names based on the test dimensionality.
#-------------------------------------------------------------------------------
exec("from Spheral%s import *" % testDim, globals())

## import Spheral
## for name in [x for x in Spheral.__dict__ if testDim in x]:
##     exec("%s = Spheral.__dict__['%s']" % (name.replace(testDim, ""), name))

#-------------------------------------------------------------------------------
# Create a random number generator.
#-------------------------------------------------------------------------------
import random
random.seed(seed)

#-------------------------------------------------------------------------------
# Material properties.
#-------------------------------------------------------------------------------
eos = GammaLawGasMKS(gamma, mu)

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
#WT = TableKernel(BSplineKernel(), 1000)
#output("WT")
#kernelExtent = WT.kernelExtent

if KernelConstructor==NBSplineKernel:
  WT = TableKernel(NBSplineKernel(order), 1000)
else:
  WT = TableKernel(KernelConstructor(), 1000)
output("WT")
kernelExtent = WT.kernelExtent

#-------------------------------------------------------------------------------
# Make the NodeList.
#-------------------------------------------------------------------------------
nodes1 = makeFluidNodeList("nodes1", eos,
                           hmin = hmin,
                           hmax = hmax,
                           kernelExtent = kernelExtent,
                           nPerh = nPerh)
output("nodes1")
output("nodes1.hmin")
output("nodes1.hmax")
output("nodes1.nodesPerSmoothingScale")

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
if testDim == "1d":
    from DistributeNodes import distributeNodesInRange1d
    distributeNodesInRange1d([(nodes1, [(nx1, rho1, (x0, x1)),
                                        (nx2, rho2, (x1, x2))])], nPerh = nPerh)
elif testDim == "2d":
    from DistributeNodes import distributeNodes2d
    from GenerateNodeDistribution2d import GenerateNodeDistribution2d
    from CompositeNodeDistribution import CompositeNodeDistribution
    gen1 = GenerateNodeDistribution2d(nx1, nx1 + nx2, rho1,
                                      distributionType = "lattice",
                                      xmin = (x0, x0),
                                      xmax = (x1, x2),
                                      nNodePerh = nPerh,
                                      SPH = True)
    gen2 = GenerateNodeDistribution2d(nx2, nx1 + nx2, rho2,
                                      distributionType = "lattice",
                                      xmin = (x1, x0),
                                      xmax = (x2, x2),
                                      nNodePerh = nPerh,
                                      SPH = True)
    gen = CompositeNodeDistribution(gen1, gen2)
    distributeNodes2d((nodes1, gen))

elif testDim == "3d":
    from DistributeNodes import distributeNodes3d
    from GenerateNodeDistribution3d import GenerateNodeDistribution3d
    from CompositeNodeDistribution import CompositeNodeDistribution
    gen1 = GenerateNodeDistribution3d(nx1, nx1 + nx2, nx1 + nx2, rho1,
                                      distributionType = "lattice",
                                      xmin = (x0, x0, x0),
                                      xmax = (x1, x2, x2),
                                      nNodePerh = nPerh,
                                      SPH = True)
    gen2 = GenerateNodeDistribution3d(nx2, nx1 + nx2, nx1 + nx2, rho2,
                                      distributionType = "lattice",
                                      xmin = (x1, x0, x0),
                                      xmax = (x2, x2, x2),
                                      nNodePerh = nPerh,
                                      SPH = True)
    gen = CompositeNodeDistribution(gen1, gen2)
    distributeNodes3d((nodes1, gen))

else:
    raise ValueError("Only tests cases for 1d,2d and 3d.") 

output("nodes1.numNodes")

# Set node properties.
eps = nodes1.specificThermalEnergy()
for i in range(nx1):
    eps[i] = eps1
for i in range(nx2):
    eps[i + nx1] = eps2

#-------------------------------------------------------------------------------
# Optionally randomly jitter the node positions.
#-------------------------------------------------------------------------------
dx1 = (x1 - x0)/nx1
dx2 = (x2 - x1)/nx2
dy = (x2 - x0)/(nx1 + nx2)
dz = (x2 - x0)/(nx1 + nx2)
pos = nodes1.positions()
import numpy as np
np.random.seed(18)
rposx = np.random.uniform(x0,x2,nodes1.numInternalNodes)
rposy = np.random.uniform(x0,x2,nodes1.numInternalNodes)
rposz = np.random.uniform(x0,x2,nodes1.numInternalNodes)
rposx = np.sort(rposx)
rposy = np.sort(rposy)
rposz = np.sort(rposz)
for i in range(nodes1.numInternalNodes):
    if pos[i] < x1:
        dx = dx1
    else:
        dx = dx2
    if testDim == "1d":
        pos[i].x += ranfrac * dx * random.uniform(-1.0, 1.0)
        #pos[i].x = rposx[i]
    elif testDim == "2d":
        pos[i].x += ranfrac * dx * random.uniform(-1.0, 1.0)
        pos[i].y += ranfrac * dy * random.uniform(-1.0, 1.0)
        #pos[i].x = rposx[i]
        #pos[i].y = rposy[i]
    elif testDim == "3d":
        pos[i].x += ranfrac * dx * random.uniform(-1.0, 1.0)
        pos[i].y += ranfrac * dy * random.uniform(-1.0, 1.0)
        pos[i].z += ranfrac * dz * random.uniform(-1.0, 1.0)
        #pos[i].x = rposx[i]
        #pos[i].y = rposy[i]
        #pos[i].z = rposz[i]

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node list
#-------------------------------------------------------------------------------
db = DataBase()
db.appendNodeList(nodes1)
output("db")
output("db.appendNodeList(nodes1)")
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Iterate the h to convergence if requested.
#-------------------------------------------------------------------------------
if iterateH:
    bounds = vector_of_Boundary()
    method = SPHSmoothingScale()
    iterateIdealH(db,
                  bounds,
                  WT,
                  method,
                  maxHIterations,
                  Htolerance)

#-------------------------------------------------------------------------------
# Initialize our field.
#-------------------------------------------------------------------------------
f = ScalarField("test field", nodes1)
gf = ScalarField("gradient test field", nodes1)
for i in range(nodes1.numInternalNodes):
    x = nodes1.positions()[i].x
    if testCase == "linear":
        f[i] = y0 + m0*x
        gf[i] = m0
    elif testCase == "quadratic":
        f[i] = y2 + m2*x*x
        gf[i] = 2.0*m2*x
    elif testCase == "cubic":
        f[i] = y3 + m3*x*x*x
        gf[i] = 3.0*m3*x*x
    elif testCase == "step":
        if x < x1:
            f[i] = y0
            gf[i] = 0.0
        else:
            f[i] = 2*y0
            gf[i] = 0.0

#-------------------------------------------------------------------------------
# Prepare variables to accumulate the test values.
#-------------------------------------------------------------------------------
accCRKSPH  = VectorField("CRKSPH type III interpolated acceleration values", nodes1)       #Locally Conservative and Inconsistant 
accRKSPHI  = VectorField("RKSPH type I interpolated acceleration values", nodes1)          #Non-conservative and Consistent 
accRKSPHII = VectorField("RKSPH type II interpolated acceleration values", nodes1)         #Globally Conservative and Inconsistant 
accRKSPHIV = VectorField("RKSPH type IV interpolated acceleration values", nodes1)         #Locally Conservative and Inconsistant 
accRKSPHV  = VectorField("RKSPH type V interpolated acceleration values", nodes1)          #Non-conservative and Consistent
accSPH     = VectorField("SPH interpolated acceleration values", nodes1)                   #It is what it is

accBCRKSPH  = VectorField("CRKSPH type III interpolated acceleration values with Bij terms", nodes1)     
accBRKSPHII = VectorField("RKSPH type II interpolated acceleration values with Bij terms", nodes1)     
accBRKSPHIV = VectorField("RKSPH type IV interpolated acceleration values with Bij terms", nodes1)     

fRK = ScalarField("RK interpolated values", nodes1)
gfRK = VectorField("RK interpolated derivative values", nodes1)
fSPH = ScalarField("SPH interpolated values", nodes1)
gfSPH = VectorField("SPH interpolated derivative values", nodes1)

A_fl = db.newFluidScalarFieldList(0.0, "A")
B_fl = db.newFluidVectorFieldList(Vector.zero, "B")
C_fl = db.newFluidTensorFieldList(Tensor.zero, "C")
gradA_fl = db.newFluidVectorFieldList(Vector.zero, "gradA")
gradB_fl = db.newFluidTensorFieldList(Tensor.zero, "gradB")
gradC_fl = db.newFluidThirdRankTensorFieldList(ThirdRankTensor.zero, "gradB")

m0_fl = db.newFluidScalarFieldList(0.0, "m0")
m1_fl = db.newFluidVectorFieldList(Vector.zero, "m1")
m2_fl = db.newFluidSymTensorFieldList(SymTensor.zero, "m2")
m3_fl = db.newFluidThirdRankTensorFieldList(ThirdRankTensor.zero, "m3")
m4_fl = db.newFluidFourthRankTensorFieldList(FourthRankTensor.zero, "m4")
gradm0_fl = db.newFluidVectorFieldList(Vector.zero, "grad m0")
gradm1_fl = db.newFluidTensorFieldList(Tensor.zero, "grad m1")
gradm2_fl = db.newFluidThirdRankTensorFieldList(ThirdRankTensor.zero, "grad m2")
gradm3_fl = db.newFluidFourthRankTensorFieldList(FourthRankTensor.zero, "grad m3")
gradm4_fl = db.newFluidFifthRankTensorFieldList(FifthRankTensor.zero, "grad m4")

db.updateConnectivityMap(True)
cm = db.connectivityMap()
position_fl = db.fluidPosition
mass_fl = db.fluidMass
H_fl = db.fluidHfield

# Compute the volumes to use as weighting.
weight_fl = db.fluidMass
#weight_fl = db.newFluidScalarFieldList(1.0, "volume")
#computeCRKSPHSumVolume(cm, WT, position_fl, mass_fl, H_fl, weight_fl)

# Now the moments and corrections
computeCRKSPHMoments(cm, WT, weight_fl, position_fl, H_fl, correctionOrder, NodeCoupling(),
                     m0_fl, m1_fl, m2_fl, m3_fl, m4_fl, gradm0_fl, gradm1_fl, gradm2_fl, gradm3_fl, gradm4_fl)
computeCRKSPHCorrections(m0_fl, m1_fl, m2_fl, m3_fl, m4_fl, gradm0_fl, gradm1_fl, gradm2_fl, gradm3_fl, gradm4_fl,
                         correctionOrder,
                         A_fl, B_fl, C_fl, gradA_fl, gradB_fl, gradC_fl)

# # Since we use the density for the SPH case, it's worth finding it's sum value.
# rho_fl = db.fluidMassDensity
# computeSPHSumMassDensity(cm, WT, True, position_fl, mass_fl, H_fl, rho_fl)

# Extract the field state for the following calculations.
positions = position_fl[0]
weight = weight_fl[0]
mass = mass_fl[0]
rho = nodes1.massDensity()
H = H_fl[0]
A = A_fl[0]
B = B_fl[0]
C = C_fl[0]
gradA = gradA_fl[0]
gradB = gradB_fl[0]
gradC = gradC_fl[0]

#-------------------------------------------------------------------------------
# Find the boundary Particles (Assumes Lattice in the generators)
#-------------------------------------------------------------------------------

isBound = [False]*nodes1.numInternalNodes
boundIndx = []
if testDim == "1d":
  boundIndx = [0,nodes1.numInternalNodes - 1]
elif testDim == "2d":
  for iglobal in range(nodes1.numInternalNodes):
     nx = nx1 + nx2
     ny = nx1 + nx2
     i = iglobal % nx
     j = (iglobal // nx) 
     if i == 0 or i == nx - 1 or j == 0 or j == ny - 1:
       boundIndx.append(iglobal)
elif testDim == "3d":
  for iglobal in range(nodes1.numInternalNodes):
       nx = nx1 + nx2
       ny = nx1 + nx2
       nz = nx1 + nx2
       i = iglobal % nx
       j = (iglobal // nx) % ny
       k = iglobal // (nx*ny)
       if i == 0 or i == nx - 1 or j == 0 or j == ny - 1 or k == 0 or k == nz - 1:
           boundIndx.append(iglobal)

for i in range(nodes1.numInternalNodes):
  if i in boundIndx:
    isBound[i] = True
  else:
    neighbors = cm.connectivityForNode(nodes1, i)
    assert len(neighbors) == 1
    for j in neighbors[0]:
      if j in boundIndx:
        isBound[i] = True
        break
  
#print isBound  
#-------------------------------------------------------------------------------
# Calculate the Acceration for all the RK schemes
#-------------------------------------------------------------------------------
for i in range(nodes1.numInternalNodes):
    ri = positions[i]
    Hi = H[i]
    Hdeti = H[i].Determinant()
    wi = weight[i]
    mi = mass[i]
    rhoi = rho[i]
    fi = f[i]
    Ai = A[i]
    Bi = B[i]
    Ci = C[i]
    gradAi = gradA[i]
    gradBi = gradB[i]
    gradCi = gradC[i]

    # Self contribution for interpolation
    W0 = WT.kernelValue(0.0, Hdeti)
    fRK[i] = wi*fi*W0*Ai;
    gfRK[i] = wi*fi*W0*(Ai*Bi+gradAi);

    fSPH[i] = wi*fi*W0

    # Self contribution for acceleration. Gradient at zero is not zero for some RK types.
    accRKSPHI[i]   = -wi*fi*W0*(Ai*Bi+gradAi)/rhoi;
    accRKSPHII[i]  =  wi*fi*W0*(Ai*Bi+gradAi)/rhoi;
    accBRKSPHII[i] =  wi*fi*W0*(Ai*Bi+gradAi)/rhoi;
    if isBound[i]:
          accBCRKSPH[i]  = -wi*(fi+fi)*W0*(Ai*Bi+gradAi)/rhoi
          accBRKSPHIV[i] = -wi*(fi+fi)*W0*(Ai*Bi+gradAi)/rhoi
          accBRKSPHII[i] = -wi*fi*W0*(Ai*Bi+gradAi)/rhoi;
      

    # Go over them neighbors.
    neighbors = cm.connectivityForNode(nodes1, i)
    assert len(neighbors) == 1
    for j in neighbors[0]:
        rj = positions[j]
        Hj = H[j]
        Hdetj = H[j].Determinant()
        wj = weight[j]
        fj = f[j]
        mj = mass[j]
        rhoj = rho[j]
        Aj = A[j]
        Bj = B[j]
        Cj = C[j]
        gradAj = gradA[j]
        gradBj = gradB[j]
        gradCj = gradC[j]

        # The standard SPH kernel and it's gradient.
        rij = ri - rj
        etai = Hi*rij
        etaj = Hj*rij
       
        #SPH Kernels and Gradients
        Wj = WT.kernelValue(etaj.magnitude(), Hdetj) 
        Wi = WT.kernelValue(etai.magnitude(), Hdeti) 
        Wij = 0.5*(Wi+Wj)
        gradWj = Hj*etaj.unitVector() * WT.gradValue(etaj.magnitude(), Hdetj)
        gradWi = Hi*etai.unitVector() * WT.gradValue(etai.magnitude(), Hdeti)
        gradWij = 0.5*(gradWj+gradWi)


        # Increment our interpolated values.
        fSPH[i] += fj * wj*Wj

        # Increment the derivatives.
        gfSPH[i] += fj * wj*gradWj



        #RK Kernels and Gradients
        rkWj   = CRKSPHKernel(WT, correctionOrder,  rij,  etai, Hdeti,  etaj, Hdetj, Ai, Bi, Ci)
        rkWi   = CRKSPHKernel(WT, correctionOrder, -rij, -etaj, Hdetj, -etai, Hdeti, Aj, Bj, Cj);
        gradrkWi  = Vector.zero
        gradrkWj  = Vector.zero
        CRKSPHKernelAndGradient(WT, correctionOrder,  rij,  etai, Hi, Hdeti,  etaj, Hj, Hdetj, Ai, Bi, Ci, gradAi, gradBi, gradCi, gradrkWj)
        CRKSPHKernelAndGradient(WT, correctionOrder, -rij, -etaj, Hj, Hdetj, -etai, Hi, Hdeti, Aj, Bj, Cj, gradAj, gradBj, gradCj, gradrkWi)
        deltagrad = gradrkWj - gradrkWi

        #accCRKSPH[i]  -= wj*(0.5*(fi+fj)*deltagrad)/rhoi
        accCRKSPH[i]  -= (wj*0.5*(fi+fj)*gradrkWj/rhoi - wi*0.5*(fi+fj)*gradrkWi*mj/(mi*rhoj)) #For some reason this is exactly the same as the commented out form to machine precision, keeping things consistent using this form
        accRKSPHI[i]  -= wj*fj*gradrkWj/rhoi;
        accRKSPHII[i] += wj*fj*gradrkWj/rhoi;
        #accRKSPHIV[i] -= wj*(fj*gradrkWj - fi*gradrkWi)/rhoi;
        accRKSPHIV[i] -= (wj*fj*gradrkWj/rhoi - wi*fi*gradrkWi*mj/(rhoj*mi)) #For some reason this is exactly the same as the commented out form to machine precision, keeping things consistent using this form
        accRKSPHV[i]  -= wj*(fj-fi)*gradrkWj/rhoi;
 
        # And of course SPH.
        accSPH[i] -= mj*(fi/(rhoi*rhoi) + fj/(rhoj*rhoj))*gradWij

        #Bij Terms 
        if isBound[i] and isBound[j]:
          accBCRKSPH[i]  -= wj*(fi+fj)*gradrkWj/rhoi
          accBRKSPHII[i] -= wi*fj*gradrkWi*mj/(mi*rhoj);
          accBRKSPHIV[i] -= wj*(fi+fj)*gradrkWj/rhoi
        else:
          #accBCRKSPH[i] -= wj*(0.5*(fi+fj)*deltagrad)/rhoi
          accBCRKSPH[i]  -= (wj*0.5*(fi+fj)*gradrkWj/rhoi - wi*0.5*(fi+fj)*gradrkWi*mj/(mi*rhoj)) 
          accBRKSPHII[i] += wj*fj*gradrkWj/rhoi;
          #accBRKSPHIV[i] -= wj*(fj*gradrkWj - fi*gradrkWi)/rhoi;
          accBRKSPHIV[i] -= (wj*fj*gradrkWj/rhoi - wi*fi*gradrkWi*mj/(rhoj*mi)) 

        #Check RK interpolation
        fRK[i] += wj*fj*rkWj
        gfRK[i] += wj*fj*gradrkWj

#-------------------------------------------------------------------------------
# Prepare the answer to check against.
#-------------------------------------------------------------------------------
xans = [positions[i].x for i in range(nodes1.numInternalNodes)]
axans = ScalarField("accelertion answer x component", nodes1)
ayans = ScalarField("accelertion answer y component", nodes1)
azans = ScalarField("accelertion answer z component", nodes1)
for i in range(nodes1.numInternalNodes):
    if testCase == "linear":
        axans[i] = -m0/rho1
        ayans[i] = 0.0
        azans[i] = 0.0
    elif testCase == "quadratic":
        axans[i] = -2.0*m2*xans[i]/rho1
        ayans[i] = 0.0
        azans[i] = 0.0
    elif testCase == "cubic":
        axans[i] = -3.0*m3*xans[i]*xans[i]/rho1
        ayans[i] = 0.0
        azans[i] = 0.0
    elif testCase == "step":
        if i < nx1:
            axans[i] = 0.0
            ayans[i] = 0.0
            azans[i] = 0.0
        else:
            axans[i] = 0.0
            ayans[i] = 0.0
            azans[i] = 0.0

#-------------------------------------------------------------------------------
# Check our answers accuracy.
#-------------------------------------------------------------------------------
errfRK      =  ScalarField("RK interpolation error", nodes1)
errgfRK     =  ScalarField("RK interpolation error", nodes1)
errxCRKSPH  =  ScalarField("CRKSPH consistency error", nodes1)
errxRKSPHI  =  ScalarField("RKSPH Type I consistency error", nodes1)
errxRKSPHII =  ScalarField("RKSPH Type II consistency error", nodes1)
errxRKSPHIV =  ScalarField("RKSPH Type IV consistency error", nodes1)
errxRKSPHV  =  ScalarField("RKSPH Type V consistency error", nodes1)
errxSPH     =  ScalarField("SPH consistency error", nodes1)

errxBCRKSPH  =  ScalarField("BCRKSPH consistency error", nodes1)
errxBRKSPHII =  ScalarField("BRKSPHII consistency error", nodes1)
errxBRKSPHIV =  ScalarField("BRKSPHIV consistency error", nodes1)
for i in range(nodes1.numInternalNodes):
    errfRK[i]      =  fRK[i] - f[i]
    errgfRK[i]     =  gfRK[i][0] - gf[i]

    errxBCRKSPH[i]   =  accBCRKSPH[i][0] - axans[i]
    errxBRKSPHII[i]  =  accBRKSPHII[i][0] - axans[i]
    errxBRKSPHIV[i]  =  accBRKSPHIV[i][0] - axans[i]

    errxCRKSPH[i]  =  accCRKSPH[i][0] - axans[i]
    errxRKSPHI[i]  =  accRKSPHI[i][0] - axans[i]
    errxRKSPHII[i] =  accRKSPHII[i][0] - axans[i]
    errxRKSPHIV[i] =  accRKSPHIV[i][0] - axans[i]
    errxRKSPHV[i]  =  accRKSPHV[i][0] - axans[i]
    errxSPH[i]     =  accSPH[i][0] - axans[i]

maxfRKerror = max([abs(x) for x in errfRK])
maxgfRKerror = max([abs(x) for x in errgfRK])
maxaxCRKSPHerror = max([abs(x) for x in errxCRKSPH])
maxaxRKSPHIerror = max([abs(x) for x in errxRKSPHI])
maxaxRKSPHIIerror = max([abs(x) for x in errxRKSPHII])
maxaxRKSPHIVerror = max([abs(x) for x in errxRKSPHIV])
maxaxRKSPHVerror = max([abs(x) for x in errxRKSPHV])
maxaxSPHerror = max([abs(x) for x in errxSPH])

maxaxBCRKSPHerror  = max([abs(x) for x in errxBCRKSPH])
maxaxBRKSPHIIerror = max([abs(x) for x in errxBRKSPHII])
maxaxBRKSPHIVerror = max([abs(x) for x in errxBRKSPHIV])
xans = [positions[i].x for i in range(nodes1.numInternalNodes)]
xidx = [i for i in range(len(xans)) if xans[i] <= dataCutMax and xans[i] >= dataCutMin]
xans = [xans[i] for i in xidx]
gfRK = [gfRK.internalValues()[i].x for i in xidx]
fRK = [fRK.internalValues()[i] for i in xidx]
gfSPH = [gfSPH.internalValues()[i].x for i in xidx]
fSPH = [fSPH.internalValues()[i] for i in xidx]
f = [f.internalValues()[i] for i in xidx]
gf = [gf.internalValues()[i] for i in xidx]
axans = [axans.internalValues()[i] for i in xidx]
accCRKSPH    = [accCRKSPH.internalValues()[i].x for i in xidx]
accRKSPHI    = [accRKSPHI.internalValues()[i].x for i in xidx]
accRKSPHII   = [accRKSPHII.internalValues()[i].x for i in xidx]
accRKSPHIV   = [accRKSPHIV.internalValues()[i].x for i in xidx]
accRKSPHV    = [accRKSPHV.internalValues()[i].x for i in xidx]
accSPH       = [accSPH.internalValues()[i].x for i in xidx]
accBCRKSPH   = [accBCRKSPH.internalValues()[i].x for i in xidx]
accBRKSPHII  = [accBRKSPHII.internalValues()[i].x for i in xidx]
accBRKSPHIV  = [accBRKSPHIV.internalValues()[i].x for i in xidx]
errfRK       = [errfRK.internalValues()[i] for i in xidx]
errgfRK      = [errgfRK.internalValues()[i] for i in xidx]
errxCRKSPH   = [errxCRKSPH.internalValues()[i] for i in xidx]
errxRKSPHI   = [errxRKSPHI.internalValues()[i] for i in xidx]
errxRKSPHII  = [errxRKSPHII.internalValues()[i] for i in xidx]
errxRKSPHIV  = [errxRKSPHIV.internalValues()[i] for i in xidx]
errxRKSPHV   = [errxRKSPHV.internalValues()[i] for i in xidx]
errxSPH      = [errxSPH.internalValues()[i] for i in xidx]
errxBCRKSPH  = [errxBCRKSPH.internalValues()[i] for i in xidx]
errxBRKSPHII = [errxBRKSPHII.internalValues()[i] for i in xidx]
errxBRKSPHIV = [errxBRKSPHIV.internalValues()[i] for i in xidx]

#print fRK.internalValues()[5], " FVAL\n"
#print gfRK.internalValues()[5].x, " FVAL\n"
#-------------------------------------------------------------------------------
# Plot the things.
#-------------------------------------------------------------------------------
if graphics:
    from SpheralGnuPlotUtilities import *
    import Gnuplot
#Data Cut

    #Initial Pressure Filed
    initdata = Gnuplot.Data(xans, f,
                           with_ = "lines",
                           title = "Answer",
                           inline = True)
    interpdata = Gnuplot.Data(xans, fRK,
                           with_ = "points",
                           title = "RK Interp",
                           inline = True)
    initDervdata = Gnuplot.Data(xans, gf,
                           with_ = "lines",
                           title = "Answer",
                           inline = True)
    interpDervdata = Gnuplot.Data(xans, gfRK,
                           with_ = "points",
                           title = "RK Interp Derv",
                           inline = True)
    ansdata = Gnuplot.Data(xans, axans,
                           with_ = "lines",
                           title = "Answer",
                           inline = True)
    #Bij fixed RK data
    BCRKSPHdata  = Gnuplot.Data(xans, accBCRKSPH,
                            with_ = "points",
                            title = "Bij-CRKSPH",
                            inline = True)
    BRKSPHIIdata = Gnuplot.Data(xans, accBRKSPHII,
                            with_ = "points",
                            title = "Bij-RKSPHII",
                            inline = True)
    BRKSPHIVdata = Gnuplot.Data(xans, accBRKSPHIV,
                            with_ = "points",
                            title = "Bij-RKSPHIV",
                            inline = True)
    #RK Schemes Data
    CRKSPHdata  = Gnuplot.Data(xans, accCRKSPH,
                            with_ = "points",
                            title = "CRKSPH",
                            inline = True)
    RKSPHIdata  = Gnuplot.Data(xans, accRKSPHI,
                            with_ = "points",
                            title = "RKSPH I",
                            inline = True)
    RKSPHIIdata = Gnuplot.Data(xans, accRKSPHII,
                            with_ = "points",
                            title = "RKSPH II",
                            inline = True)

    RKSPHIVdata = Gnuplot.Data(xans, accRKSPHIV,
                            with_ = "points",
                            title = "RKSPH IV",
                            inline = True)

    RKSPHVdata  = Gnuplot.Data(xans, accRKSPHV,
                            with_ = "points",
                            title = "RKSPH V",
                            inline = True)

    SPHdata  = Gnuplot.Data(xans, accSPH,
                            with_ = "points",
                            title = "SPH",
                            inline = True)
    #Bij fixed Error Data
    errBCRKSPHdata  = Gnuplot.Data(xans, errxBCRKSPH,
                            with_ = "points",
                            title = "Bij-CRKSPH",
                            inline = True)
    errBRKSPHIIdata = Gnuplot.Data(xans, errxBRKSPHII,
                            with_ = "points",
                            title = "Bij-RKSPHII",
                            inline = True)
    errBRKSPHIVdata = Gnuplot.Data(xans, errxBRKSPHIV,
                            with_ = "points",
                            title = "Bij-RKSPHIV",
                            inline = True)
    #Error Data
    errCRKSPHdata  = Gnuplot.Data(xans, errxCRKSPH,
                            with_ = "points",
                            title = "CRKSPH",
                            inline = True)
    errRKSPHIdata  = Gnuplot.Data(xans, errxRKSPHI,
                            with_ = "points",
                            title = "RKSPH I",
                            inline = True)
    errRKSPHIIdata = Gnuplot.Data(xans, errxRKSPHII,
                            with_ = "points",
                            title = "RKSPH II",
                            inline = True)

    errRKSPHIVdata = Gnuplot.Data(xans, errxRKSPHIV,
                            with_ = "points",
                            title = "RKSPH IV",
                            inline = True)

    errRKSPHVdata  = Gnuplot.Data(xans, errxRKSPHV,
                            with_ = "points",
                            title = "RKSPH V",
                            inline = True)

    errSPHdata  = Gnuplot.Data(xans, errxSPH,
                               with_ = "points",
                               title = "SPH",
                               inline = True)

    p0 = generateNewGnuPlot()
    p0.plot(initdata)
    p0.replot(interpdata)
    p0.title("Initial P Field and Interpolated RK")
    p0.refresh()

    p01 = generateNewGnuPlot()
    p01.plot(initDervdata)
    p01.replot(interpDervdata)
    p01.title("Initial grad P Field and Interpolated RK")
    p01.refresh()

    p1 = generateNewGnuPlot()
    p1.plot(ansdata)
    p1.replot(CRKSPHdata)
    p1.replot(RKSPHIdata)
    #p1.replot(RKSPHIIdata)
    p1.replot(RKSPHIVdata)
    p1.replot(RKSPHVdata)
    p1.replot(SPHdata)
    p1("set key top left")
    p1.title("x acceleration values")
    p1.refresh()

    p11 = generateNewGnuPlot()
    p11.plot(ansdata)
    p11.replot(CRKSPHdata)
    p11.title("x CRK acceleration values")
    p11.refresh()

    p2 = generateNewGnuPlot()
    p2.plot(errCRKSPHdata)
    p2.replot(errRKSPHIdata)
    #p2.replot(errRKSPHIIdata)
    p2.replot(errRKSPHIVdata)
    p2.replot(errRKSPHVdata)
    p2.replot(errSPHdata)
    p2("set key top left")
    p2.title("Error in acceleration")
    p2.refresh()

    p21 = generateNewGnuPlot()
    p21.plot(errCRKSPHdata)
    p21("set key top left")
    p21.title("CRK Error in acceleration")
    p21.refresh()

    #p22 = generateNewGnuPlot()
    #p22.replot(errSPHdata)
    #p22("set key top left")
    #p22.title("SPH Error in acceleration")
    #p22.refresh()

    #p23 = generateNewGnuPlot()
    #p23.replot(errRKSPHIdata)
    #p23("set key top left")
    #p23.title("RK1 Error in acceleration")
    #p23.refresh()


    if graphBij:
      p3 = generateNewGnuPlot()
      p3.plot(ansdata)
      p3.replot(CRKSPHdata)
      p3.replot(BCRKSPHdata)
      p3("set key top left")
      p3.title("x acceleration values For Testing CRK Bij")
      p3.refresh()

      p4 = generateNewGnuPlot()
      p4.plot(ansdata)
      p4.replot(RKSPHIIdata)
      p4.replot(BRKSPHIIdata)
      p4("set key top left")
      p4.title("x acceleration values For Testing RKSPH II Bij")
      p4.refresh()
     
      p5 = generateNewGnuPlot()
      p5.plot(ansdata)
      p5.replot(RKSPHIVdata)
      p5.replot(BRKSPHIVdata)
      p5("set key top left")
      p5.title("x acceleration values For Testing RKSPH IV Bij")
      p5.refresh()


      p6 = generateNewGnuPlot()
      p6.plot(errBCRKSPHdata)
      p6("set key top left")
      p6.title("Error in acceleration for Testing CRK Bij")
      p6.refresh()

      p7 = generateNewGnuPlot()
      p7.plot(errBRKSPHIIdata)
      p7("set key top left")
      p7.title("Error in acceleration for Testing RKSPH II Bij")
      p7.refresh()

      p8 = generateNewGnuPlot()
      p8.plot(errBRKSPHIVdata)
      p8("set key top left")
      p8.title("Error in acceleration for Testing RKSPH IV Bij")
      p8.refresh()


from Pnorm import Pnorm
print("L1 errors: CRKSPH = %g, RKSPH I = %g, RKSPH II = %g, RKSPH IV = %g, RKSPH V = %g, SPH = %g, BCRKSPH = %g, BRKSPHII = %g, BRKSPHIV = %g," % (Pnorm(errxCRKSPH, xans).pnorm(1),
                                                                                                  Pnorm(errxRKSPHI, xans).pnorm(1),
                                                                                                  Pnorm(errxRKSPHII, xans).pnorm(1),
                                                                                                  Pnorm(errxRKSPHIV, xans).pnorm(1),
                                                                                                  Pnorm(errxRKSPHV, xans).pnorm(1),
                                                                                                  Pnorm(errxSPH, xans).pnorm(1),
            											  Pnorm(errxBCRKSPH, xans).pnorm(1),
            											  Pnorm(errxBRKSPHII, xans).pnorm(1),
            											  Pnorm(errxBRKSPHIV, xans).pnorm(1)))
print("Maximum errors: CRKSPH = %g, RKSPH I = %g, RKSPH II = %g, RKSPH IV = %g, RKSPH V = %g, SPH = %g, BCRKSPH = %g, BRKSPHII = %g, BRKSPHIV = %g" % (maxaxCRKSPHerror, maxaxRKSPHIerror, maxaxRKSPHIIerror, maxaxRKSPHIVerror, maxaxRKSPHVerror, maxaxSPHerror, maxaxBCRKSPHerror, maxaxBRKSPHIIerror, maxaxBRKSPHIVerror))
print("L1 Interpolation Error RK = %g, Max err = %g, L1 Derivative Error Rk = %g, Max err = %g" % (Pnorm(errfRK, xans).pnorm(1),maxfRKerror, Pnorm(errgfRK, xans).pnorm(1),maxgfRKerror))
if outputFile: 
  fl = open(outputFile+".ascii", "w") 
  fl.write(("# " + 8*"%15s \t " + "\n") % ("x", "Dv/Dt", "CRK", "RK Type 1", "SPH", "CRK Error (difference)", "RK Error", "SPH Error"))
  for i in range(len(xidx)):
    fl.write((8*"%16.12e " + "\n") % (xans[i], axans[i], accCRKSPH[i], accRKSPHI[i], accSPH[i], errxCRKSPH[i], errxRKSPHI[i], errxSPH[i]))
  fl.close()
  fl = open(outputFile+"_interpolate.ascii", "w") 
  fl.write(("# " + 7*"%15s \t" + "\n") % ("x", "P", "grad P", "RK (P estimate)", "RK (grad P estimate)", "SPH (P estimate)", "SPH (grad P estimate)"))
  for i in range(len(xidx)):
    fl.write((7*"%16.12e " + "\n") % (xans[i], f[i], gf[i], fRK[i], gfRK[i], fSPH[i], gfSPH[i]))
  fl.close()
