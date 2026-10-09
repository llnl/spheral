#ATS:t0 = test(SELF, label="Field reductions of geometric types (serial)")
#ATS:t1 = test(SELF, np=3, label="Field reductions of geometric types (parallel)")
#-------------------------------------------------------------------------------
# Check the global (allReduce based) Field sumElements/min/max for the
# geometric types (Vector, Tensor, SymTensor) as well as scalars.
#-------------------------------------------------------------------------------
import unittest
import mpi
from Spheral import *

class TestFieldReductions(unittest.TestCase):

    def check(self, nDim):
        Vector = eval("Vector%id" % nDim)
        Tensor = eval("Tensor%id" % nDim)
        SymTensor = eval("SymTensor%id" % nDim)
        eos = eval("GammaLawGasMKS%id(5.0/3.0, 1.0)" % nDim)
        n = 5
        nodes = eval("makeFluidNodeList%id('nodes%id', eos, numInternal=%i)" % (nDim, nDim, n))

        # Every rank holds n points; values depend on the rank so the global
        # min and max differ from the local ones.
        r = mpi.rank
        nprocs = mpi.procs
        S = eval("ScalarField%id" % nDim)("S", nodes)
        V = eval("VectorField%id" % nDim)("V", nodes)
        T = eval("TensorField%id" % nDim)("T", nodes)
        H = eval("SymTensorField%id" % nDim)("H", nodes)
        for i in range(n):
            S[i] = float(r + i)
            V[i] = Vector.one * float(r + i)
            T[i] = Tensor.one * float(r + i)
            H[i] = SymTensor.one * float(r + i)

        # Sum over i of (r + i), then over ranks.
        expected = sum(n*rr + n*(n - 1)/2 for rr in range(nprocs))
        self.assertAlmostEqual(S.sumElements(), expected)
        self.assertEqual(V.sumElements(), Vector.one*expected)
        self.assertEqual(T.sumElements(), Tensor.one*expected)
        self.assertEqual(H.sumElements(), SymTensor.one*expected)

        # Global min and max.
        self.assertAlmostEqual(S.min(), 0.0)
        self.assertAlmostEqual(S.max(), float(nprocs - 1 + n - 1))
        self.assertEqual(V.min(), Vector.zero)
        self.assertEqual(V.max(), Vector.one*float(nprocs - 1 + n - 1))

    def test1d(self):
        self.check(1)

    def test2d(self):
        self.check(2)

    def test3d(self):
        self.check(3)

if __name__ == "__main__":
    unittest.main()
