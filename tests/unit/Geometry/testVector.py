#ATS:test(SELF, label="GeomVector unit tests")
# Unit tests for the Spheral++ GeomVector

from SpheralTestUtilities import fuzzyEqual
from math import *
import unittest

from Spheral import *

# Create a global random number generator.
import random
random.seed(889)

#-------------------------------------------------------------------------------
# Generic vector tests.
#-------------------------------------------------------------------------------
class VectorTestBase:

    def testCopy(self):
        v = self.VectorType(self.lhs)
        for i in range(self.VectorType.nDimensions):
            assert v(i) == self.lhs(i)
        return

    def testGetX(self):
        assert self.lhs.x == 10.0

    def testSetX(self):
        check = random.uniform(-1e10, 1e10)
        self.lhs.x = check
        assert self.lhs.x == check
        assert self.lhs(0) == check

    def testZero(self):
        self.lhs.Zero()
        for i in range(self.VectorType.nDimensions):
            assert self.lhs(i) == 0.0
        return

    def testNegative(self):
        v = -(self.lhs)
        for i in range(self.VectorType.nDimensions):
            assert v(i) == -(self.lhs(i))
        return

    def testVectorAddition(self):
        result = self.lhs + self.rhs
        assert isinstance(result, self.VectorType)
        for i in range(self.VectorType.nDimensions):
            assert result(i) == self.lhs(i) + self.rhs(i)
        return

    def testVectorSubtraction(self):
        result = self.lhs - self.rhs
        assert isinstance(result, self.VectorType)
        for i in range(self.VectorType.nDimensions):
            assert result(i) == self.lhs(i) - self.rhs(i)
        return

    def testScalarMultiplication(self):
        val = 44.0
        result = self.lhs * val
        assert isinstance(result, self.VectorType)
        for i in range(self.VectorType.nDimensions):
            assert result(i) == self.lhs(i) * val
        return

    def testScalarDivision(self):
        val = 44.0
        result = self.lhs / val
        assert isinstance(result, self.VectorType)
        for i in range(self.VectorType.nDimensions):
            assert fuzzyEqual(result(i), self.lhs(i) / val)
        return

    def testInPlaceVectorAddition(self):
        result = self.VectorType(self.lhs)
        result += self.rhs
        assert isinstance(result, self.VectorType)
        for i in range(self.VectorType.nDimensions):
            assert result(i) == self.lhs(i) + self.rhs(i)
        return

    def testInPlaceVectorSubtraction(self):
        result = self.VectorType(self.lhs)
        result -= self.rhs
        assert isinstance(result, self.VectorType)
        for i in range(self.VectorType.nDimensions):
            assert result(i) == self.lhs(i) - self.rhs(i)
        return

    def testInPlaceScalarMultiplication(self):
        val = 44.0
        result = self.VectorType(self.lhs)
        result *= val
        assert isinstance(result, self.VectorType)
        for i in range(self.VectorType.nDimensions):
            assert result(i) == self.lhs(i) * val
        return

    def testInPlaceScalarDivision(self):
        val = 44.0
        result = self.VectorType(self.lhs)
        result /= val
        assert isinstance(result, self.VectorType)
        for i in range(self.VectorType.nDimensions):
            assert fuzzyEqual(result(i), self.lhs(i) / val)
        return

    def testEqual(self):
        x = self.VectorType(self.lhs)
        assert x == self.lhs
        assert not (self.lhs == self.rhs)
        return

    def testNotEqual(self):
        x = self.VectorType(self.lhs)
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

    def testDot(self):
        result = self.lhs.dot(self.rhs)
        check = 0.0
        for i in range(self.VectorType.nDimensions):
            check += self.lhs(i) * self.rhs(i)
        assert result == check

    def testSelfDyad(self):
        result = self.lhs.selfdyad()
        check = self.lhs.dyad(self.lhs)
        assert isinstance(result, self.SymTensorType)
        for i in range(self.VectorType.nDimensions):
            for j in range(self.VectorType.nDimensions):
                assert result(i,j) == check(i,j)

    def testVectorMultiplication(self):
        result = self.lhs * self.rhs
        assert isinstance(result, self.TensorType)
        check = self.lhs.dyad(self.rhs)
        assert result == check

    def testUnitVector(self):
        result = self.lhs.unitVector()
        assert isinstance(result, self.VectorType)
        assert fuzzyEqual(result.magnitude(), 1.0)
        assert fuzzyEqual(self.lhs.dot(result), self.lhs.magnitude())

    def testMagnitude(self):
        result = self.lhs.magnitude()
        check = 0.0
        for i in range(self.VectorType.nDimensions):
            check += (self.lhs(i))**2
        check = sqrt(check)
        assert fuzzyEqual(result, check)

    def testMagnitude2(self):
        result = self.lhs.magnitude2()
        check = 0.0
        for i in range(self.VectorType.nDimensions):
            check += (self.lhs(i))**2
        assert fuzzyEqual(result, check)

    def testMinElement(self):
        result = self.lhs.minElement()
        check = min([self.lhs(i) for i in range(self.VectorType.nDimensions)])
        assert result == check

    def testMaxElement(self):
        result = self.lhs.maxElement()
        check = max([self.lhs(i) for i in range(self.VectorType.nDimensions)])
        assert result == check

    def testSumElements(self):
        result = self.lhs.sumElements()
        check = sum([self.lhs(i) for i in range(self.VectorType.nDimensions)])
        assert result == check

#-------------------------------------------------------------------------------
# 1-D
#-------------------------------------------------------------------------------
class Vector1dTest(VectorTestBase, unittest.TestCase):

    def setUp(self):
        self.VectorType = Vector1d
        self.TensorType = Tensor1d
        self.SymTensorType = SymTensor1d
        self.lhs = Vector1d(10.0)
        self.rhs = Vector1d(-1.0)
        return

    def tearDown(self):
        return

    def testCross(self):
        result = self.lhs.cross(self.rhs)
        assert isinstance(result, Vector3d)
        assert ((result.x == 0.0) and
                (result.y == 0.0) and
                (result.z == 0.0))

    def testDyad(self):
        result = self.lhs.dyad(self.rhs)
        assert isinstance(result, self.TensorType)
        assert result.xx == self.lhs.x * self.rhs.x

#-------------------------------------------------------------------------------
# 2-D
#-------------------------------------------------------------------------------
class Vector2dTest(VectorTestBase, unittest.TestCase):

    def setUp(self):
        self.VectorType = Vector2d
        self.TensorType = Tensor2d
        self.SymTensorType = SymTensor2d
        self.lhs = Vector2d(10.0, 431.0)
        self.rhs = Vector2d(-1.0, -10.0)
        return

    def tearDown(self):
        return

    def testGetY(self):
        assert self.lhs.y == 431.0

    def testSetY(self):
        check = random.uniform(-1e10, 1e10)
        self.lhs.y = check
        assert self.lhs.y == check
        assert self.lhs(1) == check

    def testCross(self):
        result = self.lhs.cross(self.rhs)
        assert isinstance(result, Vector3d)
        check = Vector3d(0.0,
                         0.0,
                         self.lhs.x * self.rhs.y - self.lhs.y * self.rhs.x)
        self.assertTrue(result == check,
                        "cross product failure: %s != %s" % (str(result), str(check)))

    def testDyad(self):
        result = self.lhs.dyad(self.rhs)
        assert isinstance(result, self.TensorType)
        assert result == self.TensorType(self.lhs.x * self.rhs.x,
                                         self.lhs.x * self.rhs.y,
                                         self.lhs.y * self.rhs.x,
                                         self.lhs.y * self.rhs.y)

#-------------------------------------------------------------------------------
# 3-D
#-------------------------------------------------------------------------------
class Vector3dTest(VectorTestBase, unittest.TestCase):

    def setUp(self):
        self.VectorType = Vector3d
        self.TensorType = Tensor3d
        self.SymTensorType = SymTensor3d
        self.lhs = Vector3d(10.0, 431.0, 945.5)
        self.rhs = Vector3d(-1.0, -10.0, -208.0)
        return

    def tearDown(self):
        return

    def testGetY(self):
        assert self.lhs.y == 431.0

    def testSetY(self):
        check = random.uniform(-1e10, 1e10)
        self.lhs.y = check
        assert self.lhs.y == check
        assert self.lhs(1) == check

    def testGetZ(self):
        assert self.lhs.z == 945.5

    def testSetY(self):
        check = random.uniform(-1e10, 1e10)
        self.lhs.z = check
        assert self.lhs.z == check
        assert self.lhs(2) == check

    def testCross(self):
        result = self.lhs.cross(self.rhs)
        assert isinstance(result, Vector3d)
        check = Vector3d(self.lhs.y * self.rhs.z - self.lhs.z * self.rhs.y,
                         self.lhs.z * self.rhs.x - self.lhs.x * self.rhs.z,
                         self.lhs.x * self.rhs.y - self.lhs.y * self.rhs.x)
        self.assertTrue(result == check,
                        "cross product failure: %s != %s" % (str(result), str(check)))

    def testDyad(self):
        result = self.lhs.dyad(self.rhs)
        assert isinstance(result, self.TensorType)
        assert result == self.TensorType(self.lhs.x * self.rhs.x,
                                         self.lhs.x * self.rhs.y,
                                         self.lhs.x * self.rhs.z,
                                         self.lhs.y * self.rhs.x,
                                         self.lhs.y * self.rhs.y,
                                         self.lhs.y * self.rhs.z,
                                         self.lhs.z * self.rhs.x,
                                         self.lhs.z * self.rhs.y,
                                         self.lhs.z * self.rhs.z)

#-------------------------------------------------------------------------------
# Run those tests.
#-------------------------------------------------------------------------------
if __name__ == "__main__":
    unittest.main()
