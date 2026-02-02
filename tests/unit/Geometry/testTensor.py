#ATS:test(SELF, label="GeomTensor unit tests")
# Unit tests for the Spheral++ GeomTensor classes.
#
# Note these test exercise *almost* all of the functionality of the GeomTensor
# classes, with the exception of rotationalTransform and the eigen value/vector
# routines.  Those are more extensively tested in the separate unit tests in
# testEigen.py.

from SpheralTestUtilities import fuzzyEqual
from math import *
import unittest

from Spheral import *

# Create a global random number generator.
import random
random.seed(690)
nrandom = 10000

#-------------------------------------------------------------------------------
# Helper method to compute determinants.
#-------------------------------------------------------------------------------
def computeDeterminant(tensor, nDimensions):
    assert nDimensions in (1, 2, 3)
    if nDimensions == 1:
        return tensor(0,0)
    elif nDimensions == 2:
        return tensor(0,0)*tensor(1,1) - tensor(0,1)*tensor(1,0)
    else:
        return (tensor(0,0)*tensor(1,1)*tensor(2,2) +
                tensor(0,1)*tensor(1,2)*tensor(2,0) +
                tensor(0,2)*tensor(1,0)*tensor(2,1) -
                tensor(0,0)*tensor(1,2)*tensor(2,1) -
                tensor(0,1)*tensor(1,0)*tensor(2,2) -
                tensor(0,2)*tensor(1,1)*tensor(2,0))

#-------------------------------------------------------------------------------
# Helper method to compute a random rotation matrix.
#-------------------------------------------------------------------------------
def randomRotationMatrix(n):
    theta = random.uniform(0.0, 2.0*pi)
    if n == 1:
        return Tensor1d(1.0)
    elif n == 2:
        return Tensor2d(cos(theta), sin(theta),
                        -sin(theta), cos(theta))
    elif n == 3:
        phi = random.uniform(0.0, pi)
        psi = random.uniform(0.0, pi)
        return Tensor3d(cos(psi)*cos(phi) - cos(theta)*sin(phi)*sin(psi),
                        -sin(psi)*cos(phi) - cos(theta)*sin(phi)*cos(psi),
                        sin(theta)*sin(phi),
                        cos(psi)*sin(phi) + cos(theta)*cos(phi)*sin(psi),
                        -sin(psi)*sin(phi) + cos(theta)*cos(phi)*cos(psi),
                        -sin(theta)*cos(phi),
                        sin(theta)*sin(psi),
                        sin(theta)*cos(psi),
                        cos(theta))
    else:
        raise "randomRotionMatrix: %s is an invalid number of dimensions" % str(n)

#-------------------------------------------------------------------------------
# Random symmetric matrix.
#-------------------------------------------------------------------------------
def randomSymmetricMatrix(n):
    assert n == 1 or n == 2 or n == 3
    minValue, maxValue = -1.0, 1.0
    R = randomRotationMatrix(n)
    if n == 1:
        t = SymTensor1d(random.uniform(minValue, maxValue))
    elif n == 2:
        t = SymTensor2d(random.uniform(minValue, maxValue), 0.0,
                        0.0, random.uniform(minValue, maxValue))
    elif n == 3:
        t = SymTensor3d(random.uniform(minValue, maxValue), 0.0, 0.0,
                        0.0, random.uniform(minValue, maxValue), 0.0,
                        0.0, 0.0, random.uniform(minValue, maxValue))
    t.rotationalTransform(R)
    return t

#-------------------------------------------------------------------------------
# Random positive symmetric matrix (all positive eigenvalues).
#-------------------------------------------------------------------------------
def randomPositiveSymmetricMatrix(n):
    assert n == 1 or n == 2 or n == 3
    minValue, maxValue = 0.0, 1.0
    R = randomRotationMatrix(n)
    if n == 1:
        t = SymTensor1d(random.uniform(minValue, maxValue))
    elif n == 2:
        t = SymTensor2d(random.uniform(minValue, maxValue), 0.0,
                        0.0, random.uniform(minValue, maxValue))
    elif n == 3:
        t = SymTensor3d(random.uniform(minValue, maxValue), 0.0, 0.0,
                        0.0, random.uniform(minValue, maxValue), 0.0,
                        0.0, 0.0, random.uniform(minValue, maxValue))
    t.rotationalTransform(R)
    return t

#-------------------------------------------------------------------------------
# Generic tensor tests.
#-------------------------------------------------------------------------------
class TensorTestBase:

    def testCopy(self):
        t = self.TensorType(self.lhs)
        for row in range(self.TensorType.nDimensions):
            for col in range(self.TensorType.nDimensions):
                assert t(row, col) == self.lhs(row, col)
        return

    def testCopyFullTensor(self):
        n = self.TensorType.nDimensions
        t = self.FullTensorType(*[0.5*((i//n)*n + (i%n) + (i%n)*n + (i//n))
                                  for i in range(n*n)])
        st = self.SymmetricTensorType(t)
        for row in range(n):
            for col in range(n):
                assert st(row, col) == 0.5*(row*n + col +
                                            col*n + row)
        return

    def testCopySymTensor(self):
        n = self.TensorType.nDimensions
        st = self.SymmetricTensorType(*[0.5*((i//n)*n + (i%n) + (i%n)*n + (i//n))
                                        for i in range(n*n)])
        t = self.FullTensorType(st)
        for row in range(n):
            for col in range(n):
                assert t(row, col) == 0.5*(row*n + col +
                                           col*n + row)
        return

    def testGetRow(self):
        for row in range(self.TensorType.nDimensions):
            vec = self.lhs.getRow(row)
            for col in range(self.TensorType.nDimensions):
                assert vec(col) == self.lhs(row, col)
        return

    def testGetColumn(self):
        for col in range(self.TensorType.nDimensions):
            vec = self.lhs.getColumn(col)
            for row in range(self.TensorType.nDimensions):
                assert vec(row) == self.lhs(row, col)
        return

    def testZero(self):
        self.lhs.Zero()
        for row in range(self.TensorType.nDimensions):
            for col in range(self.TensorType.nDimensions):
                assert self.lhs(row, col) == 0.0
        return


    def testIdentity(self):
        self.lhs.Identity()
        for row in range(self.TensorType.nDimensions):
            for col in range(self.TensorType.nDimensions):
                if row == col:
                    assert self.lhs(row, col) == 1.0
                else:
                    assert self.lhs(row, col) == 0.0
        return

    def testNegative(self):
        ten = -(self.lhs)
        for row in range(self.TensorType.nDimensions):
            for col in range(self.TensorType.nDimensions):
                assert ten(row, col) == -(self.lhs(row, col))
        return

    def testTensorAddition(self):
        result = self.lhs + self.rhs
        assert isinstance(result, self.TensorType)
        for row in range(self.TensorType.nDimensions):
            for col in range(self.TensorType.nDimensions):
                assert result(row, col) == self.lhs(row, col) + self.rhs(row, col)
        return

    def testTensorSubtraction(self):
        result = self.lhs - self.rhs
        assert isinstance(result, self.TensorType)
        for row in range(self.TensorType.nDimensions):
            for col in range(self.TensorType.nDimensions):
                assert result(row, col) == self.lhs(row, col) - self.rhs(row, col)
        return

    def testTensorMultiplication(self):
        result = self.lhs * self.rhs
        assert isinstance(result, self.FullTensorType)
        for row in range(self.TensorType.nDimensions):
            for col in range(self.TensorType.nDimensions):
                check = 0.0
                for i in range(self.TensorType.nDimensions):
                    check += self.lhs(row, i) * self.rhs(i, col)
                assert fuzzyEqual(result(row, col), check)
        return

    def testOtherTensorAddition(self):
        result = self.lhs + self.other
        assert isinstance(result, self.FullTensorType)
        for row in range(self.TensorType.nDimensions):
            for col in range(self.TensorType.nDimensions):
                assert result(row, col) == self.lhs(row, col) + self.other(row, col)
        return

    def testOtherTensorSubtraction(self):
        result = self.lhs - self.other
        assert isinstance(result, self.FullTensorType)
        for row in range(self.TensorType.nDimensions):
            for col in range(self.TensorType.nDimensions):
                assert result(row, col) == self.lhs(row, col) - self.other(row, col)
        return

    def testOtherTensorMultiplication(self):
        result = self.lhs * self.other
        assert isinstance(result, self.FullTensorType)
        for row in range(self.TensorType.nDimensions):
            for col in range(self.TensorType.nDimensions):
                check = 0.0
                for i in range(self.TensorType.nDimensions):
                    check += self.lhs(row, i) * self.other(i, col)
                assert fuzzyEqual(result(row, col), check)
        return
     
    def testVectorMultiplication(self):
        result = self.lhs * self.vec
        assert isinstance(result, self.VectorType)
        for row in range(self.TensorType.nDimensions):
            check = 0.0
            for col in range(self.TensorType.nDimensions):
                check += self.lhs(row, col) * self.vec(col)
            assert result(row) == check
        return

    # def testScalarAddition(self):
    #     val = 44.0
    #     result = self.lhs + val
    #     self.assertTrue(isinstance(result, self.TensorType),
    #                     "%s is not instance of %s" % (str(type(result)),
    #                                                   str(type(self.TensorType))))
    #     for row in xrange(self.TensorType.nDimensions):
    #         for col in xrange(self.TensorType.nDimensions):
    #             assert result(row, col) == self.lhs(row, col) + val
    #     return

    # def testScalarSubtraction(self):
    #     val = 44.0
    #     result = self.lhs - val
    #     assert isinstance(result, self.TensorType)
    #     for row in xrange(self.TensorType.nDimensions):
    #         for col in xrange(self.TensorType.nDimensions):
    #             assert result(row, col) == self.lhs(row, col) - val
    #     return

    def testScalarMultiplication(self):
        val = 44.0
        result = self.lhs * val
        assert isinstance(result, self.TensorType)
        for row in range(self.TensorType.nDimensions):
            for col in range(self.TensorType.nDimensions):
                assert result(row, col) == self.lhs(row, col) * val
        return

    def testScalarDivision(self):
        val = 44.0
        result = self.lhs / val
        assert isinstance(result, self.TensorType)
        for row in range(self.TensorType.nDimensions):
            for col in range(self.TensorType.nDimensions):
                assert fuzzyEqual(result(row, col), self.lhs(row, col) / val)
        return

    def testInPlaceTensorAddition(self):
        result = self.TensorType(self.lhs)
        result += self.rhs
        assert isinstance(result, self.TensorType)
        for row in range(self.TensorType.nDimensions):
            for col in range(self.TensorType.nDimensions):
                assert result(row, col) == self.lhs(row, col) + self.rhs(row, col)
        return

    def testInPlaceTensorSubtraction(self):
        result = self.TensorType(self.lhs)
        result -= self.rhs
        assert isinstance(result, self.TensorType)
        for row in range(self.TensorType.nDimensions):
            for col in range(self.TensorType.nDimensions):
                assert result(row, col) == self.lhs(row, col) - self.rhs(row, col)
        return

    # def testInPlaceScalarAddition(self):
    #     val = 44.0
    #     result = self.TensorType(self.lhs)
    #     result += val
    #     assert isinstance(result, self.TensorType)
    #     for row in xrange(self.TensorType.nDimensions):
    #         for col in xrange(self.TensorType.nDimensions):
    #             assert result(row, col) == self.lhs(row, col) + val
    #     return

    # def testInPlaceScalarSubtraction(self):
    #     val = 44.0
    #     result = self.TensorType(self.lhs)
    #     result -= val
    #     assert isinstance(result, self.TensorType)
    #     for row in xrange(self.TensorType.nDimensions):
    #         for col in xrange(self.TensorType.nDimensions):
    #             assert result(row, col) == self.lhs(row, col) - val
    #     return

    def testInPlaceScalarMultiplication(self):
        val = 44.0
        result = self.TensorType(self.lhs)
        result *= val
        assert isinstance(result, self.TensorType)
        for row in range(self.TensorType.nDimensions):
            for col in range(self.TensorType.nDimensions):
                assert result(row, col) == self.lhs(row, col) * val
        return

    def testInPlaceScalarDivision(self):
        val = 44.0
        result = self.TensorType(self.lhs)
        result /= val
        assert isinstance(result, self.TensorType)
        for row in range(self.TensorType.nDimensions):
            for col in range(self.TensorType.nDimensions):
                self.assertTrue(fuzzyEqual(result(row, col), self.lhs(row, col) / val),
                                "Tensor inplace division : %g != %g" % (result(row, col),
                                                                        self.lhs(row, col) / val))
        return

    def testEqual(self):
        x = self.TensorType(self.lhs)
        assert x == self.lhs
        assert not (self.lhs == self.rhs)
        return

    def testNotEqual(self):
        x = self.TensorType(self.lhs)
        assert self.lhs != self.rhs
        assert not (self.lhs != x)
        return

    def testLessThan(self):
        assert self.rhs < self.lhs
        assert not self.lhs < self.rhs
        return

    def testGreaterThan(self):
        assert self.lhs > self.rhs
        assert not self.rhs > self.lhs
        return

    def testLessThanOrEqual(self):
        assert self.rhs <= self.lhs
        assert self.lhs <= self.lhs
        assert not self.lhs <= self.rhs
        return

    def testGreaterThanOrEqual(self):
        assert self.lhs >= self.rhs
        assert self.lhs >= self.lhs
        assert not self.rhs >= self.lhs
        return

    def testSymmetric(self):
        result = self.lhs.Symmetric()
        assert isinstance(result, self.SymmetricTensorType)
        for row in range(self.TensorType.nDimensions):
            for col in range(self.TensorType.nDimensions):
                assert result(row, col) == 0.5*(self.lhs(row, col) +
                                                self.lhs(col, row))
        return

    def testSkewSymmetric(self):
        result = self.lhs.SkewSymmetric()
        assert isinstance(result, self.FullTensorType)
        for row in range(self.TensorType.nDimensions):
            for col in range(self.TensorType.nDimensions):
                self.assertTrue(result(row, col) == 0.5*(self.lhs(row, col) - self.lhs(col, row)),
                                "SkewSymmetric failed: %s %s" % (str(self.lhs), result))
        return

    def testTranspose(self):
        result = self.lhs.Transpose()
        for row in range(self.TensorType.nDimensions):
            for col in range(self.TensorType.nDimensions):
                assert result(row, col) == self.lhs(col, row)
        return

    def testInverse(self):
        result = self.lhs.Inverse()
        check = result*self.lhs
        for row in range(self.TensorType.nDimensions):
            for col in range(self.TensorType.nDimensions):
                if row == col:
                    assert fuzzyEqual(check(row, col), 1.0)
                else:
                    self.assertTrue(fuzzyEqual(check(row, col), 0.0, 1.0e-5),
                                    "Off diagonal not zero: %s" % str(check))
        return

    def testDiagonal(self):
        result = self.lhs.diagonalElements()
        assert isinstance(result, self.VectorType)
        for i in range(self.TensorType.nDimensions):
            assert result(i) == self.lhs(i,i)
        return

    def testTrace(self):
        result = self.lhs.Trace()
        check = 0.0
        for i in range(self.TensorType.nDimensions):
            check += self.lhs(i,i)
        assert result == check
        return

    def testDeterminant(self):
        result = self.lhs.Determinant()
        check = computeDeterminant(self.lhs, self.TensorType.nDimensions)
        assert fuzzyEqual(result, check)
        return

    def testDotUnitSym(self):
        one = self.SymmetricTensorType.one
        resultr = one * self.rhs
        resultl = self.lhs * one
        for row in range(self.TensorType.nDimensions):
            for col in range(self.TensorType.nDimensions):
                assert fuzzyEqual(resultl(row, col), self.lhs(row, col))
                assert fuzzyEqual(resultr(row, col), self.rhs(row, col))
        return

    def testDotUnitFull(self):
        one = self.FullTensorType.one
        resultr = one * self.rhs
        resultl = self.lhs * one
        for row in range(self.TensorType.nDimensions):
            for col in range(self.TensorType.nDimensions):
                assert fuzzyEqual(resultl(row, col), self.lhs(row, col))
                assert fuzzyEqual(resultr(row, col), self.rhs(row, col))
        return
                
    def testDoubleDot(self):
        result = self.lhs.doubledot(self.rhs)
        check = 0
        for row in range(self.TensorType.nDimensions):
            for col in range(self.TensorType.nDimensions):
                check += self.lhs(row, col)*self.rhs(col, row)
        self.assertTrue(fuzzyEqual(result, check),
                        "Doubledot check failure: %f != %f" % (result, check))
        return

    def testSelfDoubleDot(self):
        result = self.lhs.selfDoubledot()
        check = self.lhs.doubledot(self.lhs)
        self.assertTrue(fuzzyEqual(result, check),
                        "selfDoubledot check failure: %f != %f" % (result, check))
        return

    def testSquare(self):
        result = self.lhs.square()
        check = self.lhs * self.lhs
        for row in range(self.TensorType.nDimensions):
            for col in range(self.TensorType.nDimensions):
                self.assertTrue(fuzzyEqual(result(row, col), check(row, col)),
                                "Bad value: (%i,%i), %g != %g" % (row, col, result(row, col), check(row, col)))
 
    def testSquareElements(self):
        result = self.lhs.squareElements()
        for row in range(self.TensorType.nDimensions):
            for col in range(self.TensorType.nDimensions):
                self.assertTrue(fuzzyEqual(result(row, col), self.lhs(row, col)**2),
                                "Square matrix elements failure: %g != %g" % (result(row, col), self.lhs(row, col)**2))
        return

    def testRotation(self):
        R = randomRotationMatrix(self.TensorType.nDimensions)
        assert fuzzyEqual(R.Determinant(), 1.0)
        check = R*self.lhs*R.Transpose()
        self.lhs.rotationalTransform(R)
        for row in range(self.TensorType.nDimensions):
            for col in range(self.TensorType.nDimensions):
                assert fuzzyEqual(self.lhs(row, col), check(row, col))
        return

#-------------------------------------------------------------------------------
# Tensor tests that only the FullTensor type supports.
#-------------------------------------------------------------------------------
class FullTensorTestBase:

    def testInPlaceOtherTensorAddition(self):
        result = self.TensorType(self.lhs)
        result += self.other
        assert isinstance(result, self.TensorType)
        for row in range(self.TensorType.nDimensions):
            for col in range(self.TensorType.nDimensions):
                assert result(row, col) == self.lhs(row, col) + self.other(row, col)
        return

    def testInPlaceOtherTensorSubtraction(self):
        result = self.TensorType(self.lhs)
        result -= self.other
        assert isinstance(result, self.TensorType)
        for row in range(self.TensorType.nDimensions):
            for col in range(self.TensorType.nDimensions):
                assert result(row, col) == self.lhs(row, col) - self.other(row, col)
        return

##    def testInPlaceOtherTensorMultiplication(self):
##        result = self.TensorType(self.lhs)
##        result *= self.other
##        assert isinstance(result, self.FullTensorType)
##        for row in xrange(self.TensorType.nDimensions):
##            for col in xrange(self.TensorType.nDimensions):
##                check = 0.0
##                for i in xrange(self.TensorType.nDimensions):
##                    check += self.lhs(row, i) * self.other(i, col)
##                assert fuzzyEqual(result(row, col), check)
##        return

#-------------------------------------------------------------------------------
# Tensor tests that only the SymmetricTensor type supports.
#-------------------------------------------------------------------------------
class SymmetricTensorTestBase:

    def testSqrt(self):
        tol = {1: 1e-5,
               2: 1e-5,
               3: 1e-3}[self.TensorType.nDimensions]
        for t in range(nrandom):
            st = randomPositiveSymmetricMatrix(self.TensorType.nDimensions)
            st12 = st.sqrt()
            st12squared = st12*st12
            diff = st12squared - st
            check = diff.doubledot(diff)
            self.assertTrue(fuzzyEqual(check, 0.0, tol),
                            "SQRT failure %s != %s, eigen values=%s, %g" % (str(st12squared), str(st),
                                                                            str(st.eigenValues()), check))
        return

    def testCube(self):
        for t in range(nrandom):
            st = randomSymmetricMatrix(self.TensorType.nDimensions)
            st3 = st.cube()
            check = st*st*st
            for row in range(self.TensorType.nDimensions):
                for col in range(self.TensorType.nDimensions):
                    self.assertTrue(fuzzyEqual(st3(row, col), check(row, col)),
                                    "CUBE failure %s != %s" % (str(st3), str(check)))
        return

    def testCuberoot(self):
        tol = {1: 1e-5,
               2: 1e-5,
               3: 5e-4}[self.TensorType.nDimensions]
        for t in range(nrandom):
            st = randomSymmetricMatrix(self.TensorType.nDimensions)
            st13 = st.cuberoot()
            st13cube = st13*st13*st13
            diff = st13cube - st
            check = diff.doubledot(diff)
            self.assertTrue(fuzzyEqual(check, 0.0, tol),
                            "CUBEROOT failure %s != %s, eigen values=%s, check=%g" % (str(st13cube), str(st),
                                                                                      str(st.eigenValues()), check))
        return

##     def testPow(self):
##         tol = {1: 1e-5,
##                2: 1e-5,
##                3: 1e-3}[self.TensorType.nDimensions]
##         for t in xrange(nrandom):
##             st = randomSymmetricMatrix(self.TensorType.nDimensions)
##             p = abs(random.choice([random.uniform(-5.0, -0.5), random.uniform(0.5, 5.0)]))
##             stp = st.pow(p)
##             stpi = stp.pow(1.0/p)
##             diff = stpi - st
##             check = diff.doubledot(diff)
##             self.assertTrue(fuzzyEqual(check, 0.0, tol),
##                             "POW failure %s !=\n                            %s,\n power = %f,\n eigen values=%s\n              %s,\n check=%g" % (str(stpi), str(st), p,
##                                                                                                                                                   str(stpi.eigenValues()),
##                                                                                                                                                   str(st.eigenValues()),
##                                                                                                                                                   check))
        return

#-------------------------------------------------------------------------------
# Full tensor (1-D)
#-------------------------------------------------------------------------------
class Tensor1dTest(TensorTestBase, FullTensorTestBase, unittest.TestCase):

    def setUp(self):
        self.TensorType = Tensor1d
        self.FullTensorType = Tensor1d
        self.SymmetricTensorType = SymTensor1d
        self.VectorType = Vector1d
        self.rhs = Tensor1d(1.0)
        self.lhs = Tensor1d(10.0)
        self.other = SymTensor1d(14.0)
        self.vec = Vector1d(-2.5)
        return

    def tearDown(self):
        return

#-------------------------------------------------------------------------------
# Full tensor (2-D)
#-------------------------------------------------------------------------------
class Tensor2dTest(TensorTestBase, FullTensorTestBase, unittest.TestCase):

    def setUp(self):
        self.TensorType = Tensor2d
        self.FullTensorType = Tensor2d
        self.SymmetricTensorType = SymTensor2d
        self.VectorType = Vector2d
        self.rhs = Tensor2d(1.0, 2.0,
                            3.0, 4.0)
        self.lhs = Tensor2d(10.0, 5.0,
                            1.0,  3.0)
        self.other = SymTensor2d(14.0, 7.0,
                                 7.0, 10.0)
        self.vec = Vector2d(-2.5, -1.0)
        return

    def tearDown(self):
        return

#-------------------------------------------------------------------------------
# Full tensor (3-D)
#-------------------------------------------------------------------------------
class Tensor3dTest(TensorTestBase, FullTensorTestBase, unittest.TestCase):

    def setUp(self):
        self.TensorType = Tensor3d
        self.FullTensorType = Tensor3d
        self.SymmetricTensorType = SymTensor3d
        self.VectorType = Vector3d
        self.rhs = Tensor3d(1.0, 2.0, 3.0,
                            4.0, 5.0, 6.0,
                            7.0, 8.0, 9.0)
        self.lhs = Tensor3d(-4.0, -7.0, -3.0,
                            -9.0, -8.0, -7.3,
                            -9.9, -8.3, -9.9)
        self.other = SymTensor3d(14.0, 7.0, 3.0,
                                 7.0, 10.0, -4.0,
                                 3.0, -4.0, 6.0)
        self.vec = Vector3d(-2.5, -1.0, 89.0)
        return

    def tearDown(self):
        return

#-------------------------------------------------------------------------------
# Symmetric tensor (1-D)
#-------------------------------------------------------------------------------
class SymTensor1dTest(TensorTestBase, SymmetricTensorTestBase, unittest.TestCase):

    def setUp(self):
        self.TensorType = SymTensor1d
        self.FullTensorType = Tensor1d
        self.SymmetricTensorType = SymTensor1d
        self.VectorType = Vector1d
        self.rhs = SymTensor1d(1.0)
        self.lhs = SymTensor1d(10.0)
        self.other = Tensor1d(14.0)
        self.vec = Vector1d(-2.5)
        return

    def tearDown(self):
        return

#-------------------------------------------------------------------------------
# Symmetric tensor (2-D)
#-------------------------------------------------------------------------------
class SymTensor2dTest(TensorTestBase, SymmetricTensorTestBase, unittest.TestCase):

    def setUp(self):
        self.TensorType = SymTensor2d
        self.FullTensorType = Tensor2d
        self.SymmetricTensorType = SymTensor2d
        self.VectorType = Vector2d
        self.rhs = SymTensor2d(1.0, 2.0,
                               2.0, 4.0)
        self.lhs = SymTensor2d(10.0, 5.0,
                               5.0, 3.0)
        self.other = Tensor2d(14.0, 3.0,
                              3.0, 10.0)
        self.vec = Vector2d(-2.5, -1.0)
        return

    def tearDown(self):
        return

#-------------------------------------------------------------------------------
# Symmetric tensor (3-D)
#-------------------------------------------------------------------------------
class SymTensor3dTest(TensorTestBase, SymmetricTensorTestBase, unittest.TestCase):

    def setUp(self):
        self.TensorType = SymTensor3d
        self.FullTensorType = Tensor3d
        self.SymmetricTensorType = SymTensor3d
        self.VectorType = Vector3d
        self.rhs = SymTensor3d(1.0, 2.0, 3.0,
                               2.0, 5.0, 6.0,
                               3.0, 6.0, 9.0)
        self.lhs = SymTensor3d(-4.0, -7.0, -3.0,
                               -7.0, -8.0, -7.3,
                               -3.0, -7.3, -9.9)
        self.other = Tensor3d(14.0, 7.0, 3.0,
                              7.0, 10.0, -4.0,
                              3.0, -4.0, 6.0)
        self.vec = Vector3d(-2.5, -1.0, 89.0)
        return

    def tearDown(self):
        return

if __name__ == "__main__":
    unittest.main()
