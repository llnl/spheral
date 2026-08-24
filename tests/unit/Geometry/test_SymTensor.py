"""
Unit tests for GeomSymmetricTensor classes (SymTensor1d, SymTensor2d, SymTensor3d)
Tests the Python bindings for the Spheral geometric symmetric tensor types.
"""

import unittest
import math
from Spheral1d import SymTensor1d, Vector1d, Tensor1d
from Spheral2d import SymTensor2d, Vector2d, Tensor2d, EigenStruct2d
from Spheral3d import SymTensor3d, Vector3d, Tensor3d, EigenStruct3d


class TestSymTensor1d(unittest.TestCase):
    """Test suite for SymTensor1d"""

    def setUp(self):
        """Set up test fixtures"""
        self.zero = SymTensor1d()
        self.identity = SymTensor1d(1.0, 1.0, 1.0)
        self.st1 = SymTensor1d(2.0, 3.0, 4.0)
        self.st2 = SymTensor1d(5.0, 6.0, 7.0)

    def test_default_constructor(self):
        """Test default constructor creates zero tensor"""
        st = SymTensor1d()
        self.assertEqual(st.xx, 0.0)
        self.assertEqual(st.yy, 0.0)
        self.assertEqual(st.zz, 0.0)

    def test_scalar_constructor(self):
        """Test construction with explicit values"""
        st = SymTensor1d(1.0, 2.0, 3.0)
        self.assertEqual(st.xx, 1.0)
        self.assertEqual(st.yy, 2.0)
        self.assertEqual(st.zz, 3.0)

    def test_single_value_constructor(self):
        """Test construction with single value (diagonal)"""
        st = SymTensor1d(5.0)
        self.assertEqual(st.xx, 5.0)

    def test_copy_constructor(self):
        """Test copy constructor"""
        st = SymTensor1d(self.st1)
        self.assertEqual(st.xx, self.st1.xx)
        self.assertEqual(st.yy, self.st1.yy)
        self.assertEqual(st.zz, self.st1.zz)

    def test_properties(self):
        """Test property access"""
        self.assertEqual(self.st1.xx, 2.0)
        self.assertEqual(self.st1.yy, 3.0)
        self.assertEqual(self.st1.zz, 4.0)

        # Test property setting
        st = SymTensor1d()
        st.xx = 5.0
        st.yy = 6.0
        st.zz = 7.0
        self.assertEqual(st.xx, 5.0)
        self.assertEqual(st.yy, 6.0)
        self.assertEqual(st.zz, 7.0)

    def test_indexing(self):
        """Test element access via indexing"""
        self.assertEqual(self.st1[0], 2.0)
        self.assertEqual(self.st1[1], 3.0)
        self.assertEqual(self.st1[2], 4.0)

        # Test index assignment
        st = SymTensor1d()
        st[0] = 1.0
        st[1] = 2.0
        st[2] = 3.0
        self.assertEqual(st[0], 1.0)
        self.assertEqual(st[1], 2.0)
        self.assertEqual(st[2], 3.0)

    def test_call_operator(self):
        """Test () operator for element access"""
        st = SymTensor1d(2.0, 3.0, 4.0)
        self.assertEqual(st(0, 0), 2.0)

        # Test setting via call operator
        st(0, 0, 5.0)
        self.assertEqual(st(0, 0), 5.0)

    def test_length(self):
        """Test __len__ returns correct number of elements"""
        self.assertEqual(len(self.st1), 3)

    def test_iteration(self):
        """Test iteration over tensor elements"""
        elements = list(self.st1)
        self.assertEqual(len(elements), 3)
        self.assertEqual(elements[0], 2.0)
        self.assertEqual(elements[1], 3.0)
        self.assertEqual(elements[2], 4.0)

    def test_negation(self):
        """Test unary negation"""
        st = -self.st1
        self.assertEqual(st.xx, -2.0)
        self.assertEqual(st.yy, -3.0)
        self.assertEqual(st.zz, -4.0)

    def test_addition(self):
        """Test symmetric tensor addition"""
        st = self.st1 + self.st2
        self.assertEqual(st.xx, 7.0)
        self.assertEqual(st.yy, 9.0)
        self.assertEqual(st.zz, 11.0)

    def test_subtraction(self):
        """Test symmetric tensor subtraction"""
        st = self.st1 - self.st2
        self.assertEqual(st.xx, -3.0)
        self.assertEqual(st.yy, -3.0)
        self.assertEqual(st.zz, -3.0)

    def test_scalar_multiplication(self):
        """Test scalar multiplication"""
        st = self.st1 * 2.0
        self.assertEqual(st.xx, 4.0)
        self.assertEqual(st.yy, 6.0)
        self.assertEqual(st.zz, 8.0)

        # Test reverse multiplication
        st = 2.0 * self.st1
        self.assertEqual(st.xx, 4.0)
        self.assertEqual(st.yy, 6.0)
        self.assertEqual(st.zz, 8.0)

    def test_scalar_division(self):
        """Test scalar division"""
        st = self.st1 / 2.0
        self.assertEqual(st.xx, 1.0)
        self.assertEqual(st.yy, 1.5)
        self.assertEqual(st.zz, 2.0)

    def test_in_place_addition(self):
        """Test in-place addition"""
        st = SymTensor1d(self.st1)
        st += self.st2
        self.assertEqual(st.xx, 7.0)
        self.assertEqual(st.yy, 9.0)
        self.assertEqual(st.zz, 11.0)

    def test_in_place_subtraction(self):
        """Test in-place subtraction"""
        st = SymTensor1d(self.st1)
        st -= self.st2
        self.assertEqual(st.xx, -3.0)
        self.assertEqual(st.yy, -3.0)
        self.assertEqual(st.zz, -3.0)

    def test_in_place_scalar_multiplication(self):
        """Test in-place scalar multiplication"""
        st = SymTensor1d(self.st1)
        st *= 2.0
        self.assertEqual(st.xx, 4.0)
        self.assertEqual(st.yy, 6.0)
        self.assertEqual(st.zz, 8.0)

    def test_in_place_scalar_division(self):
        """Test in-place scalar division"""
        st = SymTensor1d(self.st1)
        st /= 2.0
        self.assertEqual(st.xx, 1.0)
        self.assertEqual(st.yy, 1.5)
        self.assertEqual(st.zz, 2.0)

    def test_equality(self):
        """Test equality operator"""
        st1 = SymTensor1d(1.0, 2.0, 3.0)
        st2 = SymTensor1d(1.0, 2.0, 3.0)
        st3 = SymTensor1d(1.0, 2.0, 4.0)
        self.assertTrue(st1 == st2)
        self.assertFalse(st1 == st3)

    def test_inequality(self):
        """Test inequality operator"""
        st1 = SymTensor1d(1.0, 2.0, 3.0)
        st2 = SymTensor1d(1.0, 2.0, 3.0)
        st3 = SymTensor1d(1.0, 2.0, 4.0)
        self.assertFalse(st1 != st2)
        self.assertTrue(st1 != st3)

    def test_zero_method(self):
        """Test Zero() method"""
        st = SymTensor1d(self.st1)
        st.Zero()
        self.assertEqual(st.xx, 0.0)
        self.assertEqual(st.yy, 0.0)
        self.assertEqual(st.zz, 0.0)

    def test_identity_method(self):
        """Test Identity() method"""
        st = SymTensor1d()
        st.Identity()
        self.assertEqual(st.xx, 1.0)
        self.assertEqual(st.yy, 1.0)
        self.assertEqual(st.zz, 1.0)

    def test_trace(self):
        """Test Trace() method"""
        trace = self.st1.Trace()
        self.assertEqual(trace, 2.0)  # Only xx counts in 1D

    def test_determinant(self):
        """Test Determinant() method"""
        det = self.st1.Determinant()
        self.assertEqual(det, 2.0)  # Only xx counts in 1D

    def test_transpose(self):
        """Test Transpose() method (same as self for symmetric)"""
        st = self.st1.Transpose()
        self.assertEqual(st.xx, self.st1.xx)
        self.assertEqual(st.yy, self.st1.yy)
        self.assertEqual(st.zz, self.st1.zz)

    def test_symmetric(self):
        """Test Symmetric() method (returns self for symmetric tensor)"""
        st = self.st1.Symmetric()
        self.assertEqual(st.xx, self.st1.xx)
        self.assertEqual(st.yy, self.st1.yy)
        self.assertEqual(st.zz, self.st1.zz)

    def test_diagonal_elements(self):
        """Test diagonalElements() method"""
        diag = self.st1.diagonalElements()
        self.assertEqual(diag.x, self.st1.xx)

    def test_max_abs_element(self):
        """Test maxAbsElement() method"""
        st = SymTensor1d(-5.0, 3.0, -1.0)
        self.assertEqual(st.maxAbsElement(), 5.0)

    def test_square(self):
        """Test square() method (matrix multiplication)"""
        st = SymTensor1d(2.0, 1.0, 1.0)
        st2 = st.square()
        self.assertEqual(st2.xx, 4.0)  # 2*2

    def test_square_elements(self):
        """Test squareElements() method"""
        st = self.st1.squareElements()
        self.assertEqual(st.xx, 4.0)   # 2^2
        self.assertEqual(st.yy, 9.0)   # 3^2
        self.assertEqual(st.zz, 16.0)  # 4^2

    def test_sqrt(self):
        """Test sqrt() method"""
        st = SymTensor1d(4.0, 1.0, 1.0)
        st_sqrt = st.sqrt()
        # Verify sqrt^2 = original
        st_squared = st_sqrt.square()
        self.assertAlmostEqual(st_squared.xx, 4.0, places=10)

    def test_pow(self):
        """Test pow() method"""
        st = SymTensor1d(2.0, 1.0, 1.0)
        st_pow = st.pow(3.0)
        # For 1D, this should be approximately 2^3 = 8
        self.assertAlmostEqual(st_pow.xx, 8.0, places=10)

    def test_static_properties(self):
        """Test static constexpr properties"""
        self.assertEqual(SymTensor1d.nDimensions, 1)
        self.assertEqual(SymTensor1d.numElements, 3)

        # Test zero and one
        z = SymTensor1d.zero
        self.assertEqual(z.xx, 0.0)
        self.assertEqual(z.yy, 0.0)
        self.assertEqual(z.zz, 0.0)

        o = SymTensor1d.one
        self.assertEqual(o.xx, 1.0)
        self.assertEqual(o.yy, 1.0)
        self.assertEqual(o.zz, 1.0)


class TestSymTensor2d(unittest.TestCase):
    """Test suite for SymTensor2d"""

    def setUp(self):
        """Set up test fixtures"""
        self.zero = SymTensor2d()
        self.identity = SymTensor2d(1, 0, 0, 1, 1)
        self.st1 = SymTensor2d(1, 2, 2, 4, 5)
        self.st2 = SymTensor2d(2, 1, 1, 3, 2)

    def test_constructor(self):
        """Test 2D constructor"""
        st = SymTensor2d(1, 2, 2, 4, 5)
        self.assertEqual(st.xx, 1)
        self.assertEqual(st.xy, 2)
        self.assertEqual(st.yx, 2)  # Must be symmetric
        self.assertEqual(st.yy, 4)
        self.assertEqual(st.zz, 5)

    def test_symmetry(self):
        """Test that xy == yx for symmetric tensor"""
        st = SymTensor2d(1, 2, 2, 4, 5)
        self.assertEqual(st.xy, st.yx)

    def test_properties_2d(self):
        """Test 2D-specific properties"""
        self.assertEqual(self.st1.xx, 1)
        self.assertEqual(self.st1.xy, 2)
        self.assertEqual(self.st1.yy, 4)
        self.assertEqual(self.st1.zz, 5)

        # Test setting
        st = SymTensor2d()
        st.xy = 5.0
        self.assertEqual(st.xy, 5.0)
        self.assertEqual(st.yx, 5.0)  # Should be symmetric

    def test_length_2d(self):
        """Test __len__ for 2D"""
        self.assertEqual(len(self.st1), 4)

    def test_trace_2d(self):
        """Test Trace() for 2D"""
        trace = self.st1.Trace()
        self.assertEqual(trace, 5.0)  # xx + yy

    def test_determinant_2d(self):
        """Test Determinant() for 2D"""
        det = self.st1.Determinant()
        expected = 1*4 - 2*2  # xx*yy - xy*yx
        self.assertEqual(det, expected)

    def test_transpose_2d(self):
        """Test Transpose() for 2D (same as original)"""
        st = self.st1.Transpose()
        self.assertEqual(st.xx, self.st1.xx)
        self.assertEqual(st.xy, self.st1.xy)
        self.assertEqual(st.yy, self.st1.yy)

    def test_vector_multiplication_2d(self):
        """Test symmetric tensor-vector multiplication"""
        v = Vector2d(1.0, 2.0)
        result = self.st1 * v
        # [1 2] [1]   [1*1 + 2*2]   [5]
        # [2 4] [2] = [2*1 + 4*2] = [10]
        self.assertAlmostEqual(result.x, 5.0)
        self.assertAlmostEqual(result.y, 10.0)

    def test_symtensor_multiplication_2d(self):
        """Test symmetric tensor-symmetric tensor multiplication"""
        result = self.st1 * self.st2
        # Returns a Tensor (may not be symmetric)
        # [1 2] [2 1]   [1*2+2*1  1*1+2*3]   [4  7]
        # [2 4] [1 3] = [2*2+4*1  2*1+4*3] = [8 14]
        self.assertAlmostEqual(result.xx, 4.0)
        self.assertAlmostEqual(result.xy, 7.0)
        self.assertAlmostEqual(result.yx, 8.0)
        self.assertAlmostEqual(result.yy, 14.0)

    def test_dot_vector_2d(self):
        """Test dot product with vector"""
        v = Vector2d(1.0, 2.0)
        result = self.st1.dot(v)
        self.assertAlmostEqual(result.x, 5.0)
        self.assertAlmostEqual(result.y, 10.0)

    def test_dot_tensor_2d(self):
        """Test dot product with tensor"""
        t = Tensor2d(1, 0, 0, 1, 1)
        result = self.st1.dot(t)
        # Result is a Tensor
        self.assertAlmostEqual(result.xx, 1.0)
        self.assertAlmostEqual(result.yy, 4.0)

    def test_doubledot_2d(self):
        """Test double dot product"""
        result = self.st1.doubledot(self.st2)
        # a_ij * b_ji = sum of element-wise products
        expected = 1*2 + 2*1 + 2*1 + 4*3
        self.assertAlmostEqual(result, expected)

    def test_self_doubledot_2d(self):
        """Test self double dot"""
        result = self.st1.selfDoubledot()
        # a_ij * a_ji: xx*xx + xy*yx + yx*xy + yy*yy
        # st1 = [1 2] (symmetric, so yx=xy=2)
        #       [2 4]
        # selfDoubledot = 1*1 + 2*2 + 2*2 + 4*4 = 1 + 4 + 4 + 16 = 25
        expected = 1*1 + 2*2 + 2*2 + 4*4
        self.assertAlmostEqual(result, expected)

    def test_inverse_2d(self):
        """Test Inverse() method"""
        st = SymTensor2d(4, 2, 2, 5, 1)
        inv = st.Inverse()
        # Verify A * A^-1 = I
        identity = st * inv
        self.assertAlmostEqual(identity.xx, 1.0, places=10)
        self.assertAlmostEqual(identity.xy, 0.0, places=10)
        self.assertAlmostEqual(identity.yy, 1.0, places=10)

    def test_get_row_2d(self):
        """Test getRow() method"""
        row0 = self.st1.getRow(0)
        self.assertEqual(row0.x, 1)
        self.assertEqual(row0.y, 2)

        row1 = self.st1.getRow(1)
        self.assertEqual(row1.x, 2)
        self.assertEqual(row1.y, 4)

    def test_get_column_2d(self):
        """Test getColumn() method"""
        col0 = self.st1.getColumn(0)
        self.assertEqual(col0.x, 1)
        self.assertEqual(col0.y, 2)

        col1 = self.st1.getColumn(1)
        self.assertEqual(col1.x, 2)
        self.assertEqual(col1.y, 4)

    def test_diagonal_elements_2d(self):
        """Test diagonalElements() for 2D"""
        diag = self.st1.diagonalElements()
        self.assertEqual(diag.x, 1)
        self.assertEqual(diag.y, 4)

    def test_square_2d(self):
        """Test square() returns SymTensor"""
        result = self.st1.square()
        # Result should be symmetric
        self.assertAlmostEqual(result.xy, result.yx)

    def test_sqrt_2d(self):
        """Test sqrt() for 2D"""
        st = SymTensor2d(4, 0, 0, 9, 1)
        st_sqrt = st.sqrt()
        st_squared = st_sqrt.square()
        self.assertAlmostEqual(st_squared.xx, 4.0, places=10)
        self.assertAlmostEqual(st_squared.yy, 9.0, places=10)

    def test_pow_2d(self):
        """Test pow() for 2D"""
        st = SymTensor2d(2, 0, 0, 3, 1)
        st_pow = st.pow(2.0)
        # Should be same as square for power 2
        st_squared = st.square()
        self.assertAlmostEqual(st_pow.xx, st_squared.xx, places=10)
        self.assertAlmostEqual(st_pow.yy, st_squared.yy, places=10)

    def test_eigenvalues_2d(self):
        """Test eigenValues() method"""
        # Use a diagonal matrix for predictable eigenvalues
        st = SymTensor2d(3, 0, 0, 2, 1)
        eigs = st.eigenValues()
        vals = sorted([eigs.x, eigs.y])
        expected = sorted([3.0, 2.0])
        self.assertAlmostEqual(vals[0], expected[0], places=10)
        self.assertAlmostEqual(vals[1], expected[1], places=10)

    def test_eigenvectors_2d(self):
        """Test eigenVectors() method"""
        st = SymTensor2d(3, 0, 0, 2, 1)
        eigen_struct = st.eigenVectors()
        # Should return an EigenStruct with eigenvalues and eigenvectors
        eigs = eigen_struct.eigenValues
        self.assertIsNotNone(eigs)

    def test_static_properties_2d(self):
        """Test static properties for 2D"""
        self.assertEqual(SymTensor2d.nDimensions, 2)
        self.assertEqual(SymTensor2d.numElements, 4)


class TestSymTensor3d(unittest.TestCase):
    """Test suite for SymTensor3d"""

    def setUp(self):
        """Set up test fixtures"""
        self.zero = SymTensor3d()
        self.identity = SymTensor3d(1, 0, 0, 0, 1, 0, 0, 0, 1)
        self.st1 = SymTensor3d(1, 2, 3, 2, 5, 6, 3, 6, 9)
        self.st2 = SymTensor3d(2, 0, 0, 0, 3, 0, 0, 0, 4)

    def test_constructor_3d(self):
        """Test 3D constructor"""
        st = SymTensor3d(1, 2, 3, 2, 5, 6, 3, 6, 9)
        self.assertEqual(st.xx, 1)
        self.assertEqual(st.xy, 2)
        self.assertEqual(st.xz, 3)
        self.assertEqual(st.yx, 2)  # Must be symmetric
        self.assertEqual(st.yy, 5)
        self.assertEqual(st.yz, 6)
        self.assertEqual(st.zx, 3)  # Must be symmetric
        self.assertEqual(st.zy, 6)  # Must be symmetric
        self.assertEqual(st.zz, 9)

    def test_symmetry_3d(self):
        """Test symmetry in 3D"""
        st = SymTensor3d(1, 2, 3, 2, 5, 6, 3, 6, 9)
        self.assertEqual(st.xy, st.yx)
        self.assertEqual(st.xz, st.zx)
        self.assertEqual(st.yz, st.zy)

    def test_properties_3d(self):
        """Test 3D-specific properties"""
        self.assertEqual(self.st1.xz, 3)
        self.assertEqual(self.st1.yz, 6)

        # Test setting
        st = SymTensor3d()
        st.xz = 1.0
        st.yz = 2.0
        self.assertEqual(st.xz, 1.0)
        self.assertEqual(st.zx, 1.0)
        self.assertEqual(st.yz, 2.0)
        self.assertEqual(st.zy, 2.0)

    def test_length_3d(self):
        """Test __len__ for 3D"""
        self.assertEqual(len(self.st1), 6)

    def test_trace_3d(self):
        """Test Trace() for 3D"""
        trace = self.st1.Trace()
        self.assertEqual(trace, 15.0)  # xx + yy + zz = 1 + 5 + 9

    def test_determinant_3d(self):
        """Test Determinant() for 3D"""
        st = SymTensor3d(2, 0, 0, 0, 3, 0, 0, 0, 4)
        det = st.Determinant()
        self.assertEqual(det, 24.0)  # Diagonal matrix: 2*3*4

    def test_transpose_3d(self):
        """Test Transpose() for 3D (same as original)"""
        st = self.st1.Transpose()
        self.assertEqual(st.xx, self.st1.xx)
        self.assertEqual(st.xy, self.st1.xy)
        self.assertEqual(st.xz, self.st1.xz)
        self.assertEqual(st.yy, self.st1.yy)
        self.assertEqual(st.yz, self.st1.yz)
        self.assertEqual(st.zz, self.st1.zz)

    def test_vector_multiplication_3d(self):
        """Test symmetric tensor-vector multiplication in 3D"""
        v = Vector3d(1.0, 2.0, 3.0)
        result = self.st1 * v
        # [1 2 3] [1]   [1+4+9]    [14]
        # [2 5 6] [2] = [2+10+18] = [30]
        # [3 6 9] [3]   [3+12+27]  [42]
        self.assertAlmostEqual(result.x, 14.0)
        self.assertAlmostEqual(result.y, 30.0)
        self.assertAlmostEqual(result.z, 42.0)

    def test_inverse_3d(self):
        """Test Inverse() for 3D"""
        st = SymTensor3d(2, 0, 0, 0, 3, 0, 0, 0, 4)
        inv = st.Inverse()
        # Verify A * A^-1 = I
        identity = st * inv
        self.assertAlmostEqual(identity.xx, 1.0, places=10)
        self.assertAlmostEqual(identity.yy, 1.0, places=10)
        self.assertAlmostEqual(identity.zz, 1.0, places=10)
        self.assertAlmostEqual(identity.xy, 0.0, places=10)
        self.assertAlmostEqual(identity.xz, 0.0, places=10)
        self.assertAlmostEqual(identity.yz, 0.0, places=10)

    def test_get_row_3d(self):
        """Test getRow() for 3D"""
        row0 = self.st1.getRow(0)
        self.assertEqual(row0.x, 1)
        self.assertEqual(row0.y, 2)
        self.assertEqual(row0.z, 3)

        row1 = self.st1.getRow(1)
        self.assertEqual(row1.x, 2)
        self.assertEqual(row1.y, 5)
        self.assertEqual(row1.z, 6)

        row2 = self.st1.getRow(2)
        self.assertEqual(row2.x, 3)
        self.assertEqual(row2.y, 6)
        self.assertEqual(row2.z, 9)

    def test_get_column_3d(self):
        """Test getColumn() for 3D"""
        col0 = self.st1.getColumn(0)
        self.assertEqual(col0.x, 1)
        self.assertEqual(col0.y, 2)
        self.assertEqual(col0.z, 3)

        col1 = self.st1.getColumn(1)
        self.assertEqual(col1.x, 2)
        self.assertEqual(col1.y, 5)
        self.assertEqual(col1.z, 6)

        col2 = self.st1.getColumn(2)
        self.assertEqual(col2.x, 3)
        self.assertEqual(col2.y, 6)
        self.assertEqual(col2.z, 9)

    def test_diagonal_elements_3d(self):
        """Test diagonalElements() for 3D"""
        diag = self.st1.diagonalElements()
        self.assertEqual(diag.x, 1)
        self.assertEqual(diag.y, 5)
        self.assertEqual(diag.z, 9)

    def test_eigenvalues_3d(self):
        """Test eigenValues() for 3D"""
        # Use a diagonal matrix for predictable eigenvalues
        st = SymTensor3d(2, 0, 0, 0, 3, 0, 0, 0, 5)
        eigs = st.eigenValues()
        vals = sorted([eigs.x, eigs.y, eigs.z])
        expected = sorted([2.0, 3.0, 5.0])
        for v1, v2 in zip(vals, expected):
            self.assertAlmostEqual(v1, v2, places=10)

    def test_eigenvectors_3d(self):
        """Test eigenVectors() for 3D"""
        st = SymTensor3d(2, 0, 0, 0, 3, 0, 0, 0, 5)
        eigen_struct = st.eigenVectors()
        eigs = eigen_struct.eigenValues
        self.assertIsNotNone(eigs)
        # For a diagonal matrix, eigenvectors should be basis vectors

    def test_sqrt_3d(self):
        """Test sqrt() for 3D"""
        st = SymTensor3d(4, 0, 0, 0, 9, 0, 0, 0, 16)
        st_sqrt = st.sqrt()
        st_squared = st_sqrt.square()
        self.assertAlmostEqual(st_squared.xx, 4.0, places=10)
        self.assertAlmostEqual(st_squared.yy, 9.0, places=10)
        self.assertAlmostEqual(st_squared.zz, 16.0, places=10)

    def test_pow_3d(self):
        """Test pow() for 3D"""
        st = SymTensor3d(2, 0, 0, 0, 3, 0, 0, 0, 4)
        st_cubed = st.pow(3.0)
        # For diagonal matrix: eigenvalues raised to power
        self.assertAlmostEqual(st_cubed.xx, 8.0, places=10)
        self.assertAlmostEqual(st_cubed.yy, 27.0, places=10)
        self.assertAlmostEqual(st_cubed.zz, 64.0, places=10)

    def test_static_properties_3d(self):
        """Test static properties for 3D"""
        self.assertEqual(SymTensor3d.nDimensions, 3)
        self.assertEqual(SymTensor3d.numElements, 6)

    def test_3d_methods(self):
        """Test 3D-specific methods"""
        diag3d = self.st1.diagonalElements3D()
        self.assertEqual(diag3d.x, 1)
        self.assertEqual(diag3d.y, 5)
        self.assertEqual(diag3d.z, 9)

        trace3d = self.st1.Trace3D()
        self.assertEqual(trace3d, 15.0)

        eigen3d = self.st2.eigenVectors3D()
        self.assertIsNotNone(eigen3d)


class TestSymTensorTensorInteraction(unittest.TestCase):
    """Test interaction between SymTensor and Tensor"""

    def test_symtensor_from_tensor_2d(self):
        """Test constructing SymTensor from Tensor (extracts symmetric part)"""
        t = Tensor2d(1, 2, 3, 4, 5)
        st = SymTensor2d(t)
        # Should extract symmetric part
        self.assertEqual(st.xx, 1)
        self.assertEqual(st.yy, 4)

    def test_symtensor_tensor_addition_2d(self):
        """Test adding SymTensor and Tensor"""
        st = SymTensor2d(1, 1, 1, 2, 3)
        t = Tensor2d(1, 2, 3, 4, 5)
        result = st + t
        # Result is a Tensor
        self.assertEqual(result.xx, 2)
        self.assertEqual(result.xy, 3)
        self.assertEqual(result.yx, 4)
        self.assertEqual(result.yy, 6)

    def test_symtensor_tensor_multiplication_2d(self):
        """Test multiplying SymTensor and Tensor"""
        st = SymTensor2d(1, 0, 0, 1, 1)
        t = Tensor2d(2, 0, 0, 2, 1)
        result = st * t
        # Result is a Tensor
        self.assertEqual(result.xx, 2)
        self.assertEqual(result.yy, 2)


class TestSymTensorRepr(unittest.TestCase):
    """Test string representation"""

    def test_repr_1d(self):
        """Test __repr__ for 1D"""
        st = SymTensor1d(1.0, 2.0, 3.0)
        repr_str = repr(st)
        self.assertIn("SymTensor1d", repr_str)

    def test_repr_2d(self):
        """Test __repr__ for 2D"""
        st = SymTensor2d(1, 2, 2, 4, 5)
        repr_str = repr(st)
        self.assertIn("SymTensor2d", repr_str)

    def test_repr_3d(self):
        """Test __repr__ for 3D"""
        st = SymTensor3d(1, 2, 3, 2, 5, 6, 3, 6, 9)
        repr_str = repr(st)
        self.assertIn("SymTensor3d", repr_str)


class TestSymTensorAdvanced(unittest.TestCase):
    """Test advanced SymTensor functionality"""

    def test_rotational_transform_3d(self):
        """Test rotational transform"""
        st = SymTensor3d(2, 0, 0, 0, 3, 0, 0, 0, 4)
        # Identity rotation - should not change
        rotation = Tensor3d(1, 0, 0, 0, 1, 0, 0, 0, 1)
        st.rotationalTransform(rotation)
        self.assertAlmostEqual(st.xx, 2.0, places=10)
        self.assertAlmostEqual(st.yy, 3.0, places=10)
        self.assertAlmostEqual(st.zz, 4.0, places=10)

    def test_skew_symmetric_returns_zero_tensor(self):
        """Test that SkewSymmetric of a symmetric tensor returns zero"""
        st = SymTensor3d(1, 2, 3, 2, 5, 6, 3, 6, 9)
        skew = st.SkewSymmetric()
        # All elements should be near zero
        self.assertAlmostEqual(skew.xy, 0.0, places=10)
        self.assertAlmostEqual(skew.xz, 0.0, places=10)
        self.assertAlmostEqual(skew.yz, 0.0, places=10)

    def test_comparison_operators(self):
        """Test comparison operators"""
        st1 = SymTensor3d(1, 0, 0, 0, 2, 0, 0, 0, 3)
        st2 = SymTensor3d(2, 0, 0, 0, 3, 0, 0, 0, 4)
        st3 = SymTensor3d(1, 0, 0, 0, 2, 0, 0, 0, 3)

        self.assertTrue(st1 == st3)
        self.assertFalse(st1 == st2)
        self.assertTrue(st1 != st2)
        self.assertTrue(st1 < st2)
        self.assertTrue(st2 > st1)
        self.assertTrue(st1 <= st3)
        self.assertTrue(st1 >= st3)


if __name__ == '__main__':
    unittest.main()
