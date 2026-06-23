#-------------------------------------------------------------------------------
# Test interpolation of a constant field using the TableKernel
#-------------------------------------------------------------------------------
from Spheral import *
from GenerateNodeProfile import *
from SortAndDivideRedistributeNodes import *
import matplotlib.pyplot as plt

# Command line arguments
commandLine(n = 100,
            nPerh = 4.0,
            F0 = 2.0,
            a = 0.0,
            b = 0.0,
            r0 = 0.0,
            r1 = 1.0)

def F(r):
    return F0 + a*r + b*r*r

def gradF(r):
    return a + 2.0*b*r
    
# Build the SphericalKernel
W1 = TableKernel1d(BSplineKernel1d(), 500)
W2 = TableKernel1d(WendlandC4Kernel1d(), 500)

for W in (W1,):

    # Generate some points
    eos = GammaLawGasMKS1d(2.0, 1.0)
    nodes = makeFluidNodeList1d("nodes", eos, 
                                hmin = 1e-10,
                                hmax = 1e10,
                                nPerh = nPerh,
                                kernelExtent = W.kernelExtent)
    gen = GenerateNodeProfile1d(nx = n,
                                rho = 1.0,
                                xmin = r0,
                                xmax = r1,
                                nNodePerh = nPerh)
    distributeNodes1d((nodes, gen))
    pos = nodes.positions()
    mass = nodes.mass()
    rho = nodes.massDensity()
    H = nodes.Hfield()

    # Build the DataBase and connectivity
    db = DataBase1d()
    db.appendNodeList(nodes)
    db.updateConnectivityMap(False, False, False)
    cm = db.connectivityMap()
    pairs = cm.nodePairList

    # Create our constant Field
    field0 = ScalarField1d("Initial field", nodes, F0)
    for i in range(nodes.numInternalNodes):
        field0[i] = F(pos[i].x)

    # Interpolate to the new field
    field1 = ScalarField1d("Interpolated field", nodes)
    grad_field1 = VectorField1d("Interpolated gradient", nodes)
    rho1 = ScalarField1d("Interpolated rho (gather)", nodes)

    def pair_sum(i, j):
        ri, mi, rhoi, Hi = pos[i], mass[i], rho[i], H[i]
        rj, mj, rhoj, Hj = pos[j], mass[j], rho[j], H[j]
        Wijj, gradWijj, deltaWsum = W.kernelAndGrad(Hj*rj, Hj*ri, Hj)
        field1[i] += mj/rhoj * field0[j] * Wijj
        grad_field1[i] += mj/rhoj * field0[j] * gradWijj
        rho1[i] += mj * Wijj

    for pair in pairs:
        pair_sum(pair.i_node, pair.j_node)
        pair_sum(pair.j_node, pair.i_node)
    for i in range(nodes.numInternalNodes):
        pair_sum(i, i)

    # Plot the results
    fig1 = plt.figure()
    rvals = [posi.x for posi in pos.internalValues()]
    plt.plot(rvals, field0.internalValues(), label="Analytic field")
    plt.plot(rvals, field1.internalValues(), label="Scatter interpolated")
    plt.legend()

    fig2 = plt.figure()
    plt.plot(rvals, [gradF(pos[i].x) for i in range(nodes.numInternalNodes)], label="Analytic gradient")
    plt.plot(rvals, [grad_field1[i].x for i in range(nodes.numInternalNodes)], label="Numerical gradient")
    plt.legend()

plt.show()

    
