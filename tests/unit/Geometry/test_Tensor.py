"""
Unit tests for GeomTensor classes (Tensor1d, Tensor2d, Tensor3d)
Tests the Python bindings for the Spheral geometric tensor types.
"""

import unittest
import math
from Spheral1d import Tensor1d, Vector1d, SymTensor1d
from Spheral2d import Tensor2d, Vector2d, SymTensor2d
from Spheral3d import Tensor3d, Vector3d, SymTensor3d


class TestTensor1d(unittest.TestCase):
    """Test suite for Tensor1d"""

    def setUp(self):
        """Set up test fixtures"""
        self.zero = Tensor1d()
        self.identity = Tensor1d(1.0, 1.0, 1.0)
        self.t1 = Tensor1d(2.0, 3.0, 4.0)
        self.t2 = Tensor1d(5.0, 6.0, 7.0)

    def test_default_constructor(self):
        """Test default constructor creates zero tensor"""
        t = Tensor1d()
        self.assertEqual(t.xx, 0.0)
        self.assertEqual(t.yy, 0.0)
        self.assertEqual(t.zz, 0.0)

    def test_scalar_constructor(self):
        """Test construction with explicit values"""
        t = Tensor1d(1.0, 2.0, 3.0)
        self.assertEqual(t.xx, 1.0)
        self.assertEqual(t.yy, 2.0)
        self.assertEqual(t.zz, 3.0)

    def test_copy_constructor(self):
        """Test copy constructor"""
        t = Tensor1d(self.t1)
        self.assertEqual(t.xx, self.t1.xx)
        self.assertEqual(t.yy, self.t1.yy)
        self.assertEqual(t.zz, self.t1.zz)

    def test_properties(self):
        """Test property access"""
        self.assertEqual(self.t1.xx, 2.0)
        self.assertEqual(self.t1.yy, 3.0)
        self.assertEqual(self.t1.zz, 4.0)

        # Test property setting
        t = Tensor1d()
        t.xx = 5.0
        t.yy = 6.0
        t.zz = 7.0
        self.assertEqual(t.xx, 5.0)
        self.assertEqual(t.yy, 6.0)
        self.assertEqual(t.zz, 7.0)

    def test_indexing(self):
        """Test element access via indexing"""
        self.assertEqual(self.t1[0], 2.0)
        self.assertEqual(self.t1[1], 3.0)
        self.assertEqual(self.t1[2], 4.0)

        # Test index assignment
        t = Tensor1d()
        t[0] = 1.0
        t[1] = 2.0
        t[2] = 3.0
        self.assertEqual(t[0], 1.0)
        self.assertEqual(t[1], 2.0)
        self.assertEqual(t[2], 3.0)

    def test_call_operator(self):
        """Test () operator for element access"""
        t = Tensor1d(2.0, 3.0, 4.0)
        self.assertEqual(t(0, 0), 2.0)

        # Test setting via call operator
        t(0, 0, 5.0)
        self.assertEqual(t(0, 0), 5.0)

    def test_length(self):
        """Test __len__ returns correct number of elements"""
        self.assertEqual(len(self.t1), 3)

    def test_iteration(self):
        """Test iteration over tensor elements"""
        elements = list(self.t1)
        self.assertEqual(len(elements), 3)
        self.assertEqual(elements[0], 2.0)
        self.assertEqual(elements[1], 3.0)
        self.assertEqual(elements[2], 4.0)

    def test_negation(self):
        """Test unary negation"""
        t = -self.t1
        self.assertEqual(t.xx, -2.0)
        self.assertEqual(t.yy, -3.0)
        self.assertEqual(t.zz, -4.0)

    def test_addition(self):
        """Test tensor addition"""
        t = self.t1 + self.t2
        self.assertEqual(t.xx, 7.0)
        self.assertEqual(t.yy, 9.0)
        self.assertEqual(t.zz, 11.0)

    def test_subtraction(self):
        """Test tensor subtraction"""
        t = self.t1 - self.t2
        self.assertEqual(t.xx, -3.0)
        self.assertEqual(t.yy, -3.0)
        self.assertEqual(t.zz, -3.0)

    def test_scalar_multiplication(self):
        """Test scalar multiplication"""
        t = self.t1 * 2.0
        self.assertEqual(t.xx, 4.0)
        self.assertEqual(t.yy, 6.0)
        self.assertEqual(t.zz, 8.0)

        # Test reverse multiplication
        t = 2.0 * self.t1
        self.assertEqual(t.xx, 4.0)
        self.assertEqual(t.yy, 6.0)
        self.assertEqual(t.zz, 8.0)

    def test_scalar_division(self):
        """Test scalar division"""
        t = self.t1 / 2.0
        self.assertEqual(t.xx, 1.0)
        self.assertEqual(t.yy, 1.5)
        self.assertEqual(t.zz, 2.0)

    def test_in_place_addition(self):
        """Test in-place addition"""
        t = Tensor1d(self.t1)
        t += self.t2
        self.assertEqual(t.xx, 7.0)
        self.assertEqual(t.yy, 9.0)
        self.assertEqual(t.zz, 11.0)

    def test_in_place_subtraction(self):
        """Test in-place subtraction"""
        t = Tensor1d(self.t1)
        t -= self.t2
        self.assertEqual(t.xx, -3.0)
        self.assertEqual(t.yy, -3.0)
        self.assertEqual(t.zz, -3.0)

    def test_in_place_scalar_multiplication(self):
        """Test in-place scalar multiplication"""
        t = Tensor1d(self.t1)
        t *= 2.0
        self.assertEqual(t.xx, 4.0)
        self.assertEqual(t.yy, 6.0)
        self.assertEqual(t.zz, 8.0)

    def test_in_place_scalar_division(self):
        """Test in-place scalar division"""
        t = Tensor1d(self.t1)
        t /= 2.0
        self.assertEqual(t.xx, 1.0)
        self.assertEqual(t.yy, 1.5)
        self.assertEqual(t.zz, 2.0)

    def test_equality(self):
        """Test equality operator"""
        t1 = Tensor1d(1.0, 2.0, 3.0)
        t2 = Tensor1d(1.0, 2.0, 3.0)
        t3 = Tensor1d(1.0, 2.0, 4.0)
        self.assertTrue(t1 == t2)
        self.assertFalse(t1 == t3)

    def test_inequality(self):
        """Test inequality operator"""
        t1 = Tensor1d(1.0, 2.0, 3.0)
        t2 = Tensor1d(1.0, 2.0, 3.0)
        t3 = Tensor1d(1.0, 2.0, 4.0)
        self.assertFalse(t1 != t2)
        self.assertTrue(t1 != t3)

    def test_zero_method(self):
        """Test Zero() method"""
        t = Tensor1d(self.t1)
        t.Zero()
        self.assertEqual(t.xx, 0.0)
        self.assertEqual(t.yy, 0.0)
        self.assertEqual(t.zz, 0.0)

    def test_identity_method(self):
        """Test Identity() method"""
        t = Tensor1d()
        t.Identity()
        self.assertEqual(t.xx, 1.0)
        self.assertEqual(t.yy, 1.0)
        self.assertEqual(t.zz, 1.0)

    def test_trace(self):
        """Test Trace() method"""
        trace = self.t1.Trace()
        self.assertEqual(trace, 2.0)  # Only xx counts in 1D

    def test_determinant(self):
        """Test Determinant() method"""
        det = self.t1.Determinant()
        self.assertEqual(det, 2.0)  # Only xx counts in 1D

    def test_transpose(self):
        """Test Transpose() method"""
        t = self.t1.Transpose()
        # In 1D, transpose is same as original for diagonal
        self.assertEqual(t.xx, self.t1.xx)
        self.assertEqual(t.yy, self.t1.yy)
        self.assertEqual(t.zz, self.t1.zz)

    def test_diagonal_elements(self):
        """Test diagonalElements() method"""
        diag = self.t1.diagonalElements()
        self.assertEqual(diag.x, self.t1.xx)

    def test_max_abs_element(self):
        """Test maxAbsElement() method"""
        t = Tensor1d(-5.0, 3.0, -1.0)
        self.assertEqual(t.maxAbsElement(), 5.0)

    def test_square(self):
        """Test square() method (matrix multiplication)"""
        t = Tensor1d(2.0, 1.0, 1.0)
        t2 = t.square()
        self.assertEqual(t2.xx, 4.0)  # 2*2

    def test_square_elements(self):
        """Test squareElements() method"""
        t = self.t1.squareElements()
        self.assertEqual(t.xx, 4.0)   # 2^2
        self.assertEqual(t.yy, 9.0)   # 3^2
        self.assertEqual(t.zz, 16.0)  # 4^2

    def test_static_properties(self):
        """Test static constexpr properties"""
        self.assertEqual(Tensor1d.nDimensions, 1)
        self.assertEqual(Tensor1d.numElements, 3)

        # Test zero and one
        z = Tensor1d.zero
        self.assertEqual(z.xx, 0.0)
        self.assertEqual(z.yy, 0.0)
        self.assertEqual(z.zz, 0.0)

        o = Tensor1d.one
        self.assertEqual(o.xx, 1.0)
        self.assertEqual(o.yy, 1.0)
        self.assertEqual(o.zz, 1.0)


class TestTensor2d(unittest.TestCase):
    """Test suite for Tensor2d"""

    def setUp(self):
        """Set up test fixtures"""
        self.zero = Tensor2d()
        self.identity = Tensor2d(1, 0, 0, 1, 1)
        self.t1 = Tensor2d(1, 2, 3, 4, 5)
        self.t2 = Tensor2d(2, 1, 1, 2, 3)

    def test_constructor(self):
        """Test 2D constructor"""
        t = Tensor2d(1, 2, 3, 4, 5)
        self.assertEqual(t.xx, 1)
        self.assertEqual(t.xy, 2)
        self.assertEqual(t.yx, 3)
        self.assertEqual(t.yy, 4)
        self.assertEqual(t.zz, 5)

    def test_properties_2d(self):
        """Test 2D-specific properties"""
        self.assertEqual(self.t1.xx, 1)
        self.assertEqual(self.t1.xy, 2)
        self.assertEqual(self.t1.yx, 3)
        self.assertEqual(self.t1.yy, 4)
        self.assertEqual(self.t1.zz, 5)

        # Test setting
        t = Tensor2d()
        t.xy = 5.0
        t.yx = 6.0
        self.assertEqual(t.xy, 5.0)
        self.assertEqual(t.yx, 6.0)

    def test_length_2d(self):
        """Test __len__ for 2D"""
        self.assertEqual(len(self.t1), 5)

    def test_trace_2d(self):
        """Test Trace() for 2D"""
        trace = self.t1.Trace()
        self.assertEqual(trace, 5.0)  # xx + yy

    def test_determinant_2d(self):
        """Test Determinant() for 2D"""
        det = self.t1.Determinant()
        expected = 1*4 - 2*3  # xx*yy - xy*yx
        self.assertEqual(det, expected)

    def test_transpose_2d(self):
        """Test Transpose() for 2D"""
        t = self.t1.Transpose()
        self.assertEqual(t.xx, self.t1.xx)
        self.assertEqual(t.xy, self.t1.yx)
        self.assertEqual(t.yx, self.t1.xy)
        self.assertEqual(t.yy, self.t1.yy)

    def test_vector_multiplication_2d(self):
        """Test tensor-vector multiplication"""
        v = Vector2d(1.0, 2.0)
        result = self.t1 * v
        # [1 2] [1]   [1*1 + 2*2]   [5]
        # [3 4] [2] = [3*1 + 4*2] = [11]
        self.assertAlmostEqual(result.x, 5.0)
        self.assertAlmostEqual(result.y, 11.0)

    def test_tensor_multiplication_2d(self):
        """Test tensor-tensor multiplication"""
        result = self.t1 * self.t2
        # Matrix multiplication
        # [1 2] [2 1]   [1*2+2*1  1*1+2*2]   [4  5]
        # [3 4] [1 2] = [3*2+4*1  3*1+4*2] = [10 11]
        self.assertAlmostEqual(result.xx, 4.0)
        self.assertAlmostEqual(result.xy, 5.0)
        self.assertAlmostEqual(result.yx, 10.0)
        self.assertAlmostEqual(result.yy, 11.0)

    def test_dot_vector_2d(self):
        """Test dot product with vector"""
        v = Vector2d(1.0, 2.0)
        result = self.t1.dot(v)
        self.assertAlmostEqual(result.x, 5.0)
        self.assertAlmostEqual(result.y, 11.0)

    def test_dot_tensor_2d(self):
        """Test dot product with tensor"""
        result = self.t1.dot(self.t2)
        self.assertAlmostEqual(result.xx, 4.0)
        self.assertAlmostEqual(result.xy, 5.0)
        self.assertAlmostEqual(result.yx, 10.0)
        self.assertAlmostEqual(result.yy, 11.0)

    def test_doubledot_2d(self):
        """Test double dot product"""
        # a_ij * b_ji
        result = self.t1.doubledot(self.t2)
        # Sum of element-wise products: xx*xx + xy*yx + yx*xy + yy*yy
        expected = 1*2 + 2*1 + 3*1 + 4*2
        self.assertAlmostEqual(result, expected)

    def test_self_doubledot_2d(self):
        """Test self double dot"""
        result = self.t1.selfDoubledot()
        # a_ij * a_ji: xx*xx + xy*yx + yx*xy + yy*yy
        # t1 = [1 2]
        #      [3 4]
        # selfDoubledot = 1*1 + 2*3 + 3*2 + 4*4 = 1 + 6 + 6 + 16 = 29
        expected = 1*1 + 2*3 + 3*2 + 4*4
        self.assertAlmostEqual(result, expected)

    def test_symmetric_2d(self):
        """Test Symmetric() extracts symmetric part"""
        sym = self.t1.Symmetric()
        # Symmetric part: (A + A^T)/2
        self.assertAlmostEqual(sym.xx, 1.0)
        self.assertAlmostEqual(sym.xy, 2.5)  # (2+3)/2
        self.assertAlmostEqual(sym.yy, 4.0)

    def test_skew_symmetric_2d(self):
        """Test SkewSymmetric() extracts skew-symmetric part"""
        skew = self.t1.SkewSymmetric()
        # Skew-symmetric part: (A - A^T)/2
        self.assertAlmostEqual(skew.xx, 0.0)
        self.assertAlmostEqual(skew.xy, -0.5)  # (2-3)/2
        self.assertAlmostEqual(skew.yx, 0.5)   # (3-2)/2
        self.assertAlmostEqual(skew.yy, 0.0)

    def test_inverse_2d(self):
        """Test Inverse() method"""
        t = Tensor2d(4, 7, 2, 6, 1)  # Use a matrix with non-zero det
        inv = t.Inverse()
        # Verify A * A^-1 = I
        identity = t * inv
        self.assertAlmostEqual(identity.xx, 1.0, places=10)
        self.assertAlmostEqual(identity.xy, 0.0, places=10)
        self.assertAlmostEqual(identity.yx, 0.0, places=10)
        self.assertAlmostEqual(identity.yy, 1.0, places=10)

    def test_get_row_2d(self):
        """Test getRow() method"""
        row0 = self.t1.getRow(0)
        self.assertEqual(row0.x, 1)
        self.assertEqual(row0.y, 2)

        row1 = self.t1.getRow(1)
        self.assertEqual(row1.x, 3)
        self.assertEqual(row1.y, 4)

    def test_get_column_2d(self):
        """Test getColumn() method"""
        col0 = self.t1.getColumn(0)
        self.assertEqual(col0.x, 1)
        self.assertEqual(col0.y, 3)

        col1 = self.t1.getColumn(1)
        self.assertEqual(col1.x, 2)
        self.assertEqual(col1.y, 4)

    def test_diagonal_elements_2d(self):
        """Test diagonalElements() for 2D"""
        diag = self.t1.diagonalElements()
        self.assertEqual(diag.x, 1)
        self.assertEqual(diag.y, 4)

    def test_static_properties_2d(self):
        """Test static properties for 2D"""
        self.assertEqual(Tensor2d.nDimensions, 2)
        self.assertEqual(Tensor2d.numElements, 5)


class TestTensor3d(unittest.TestCase):
    """Test suite for Tensor3d"""

    def setUp(self):
        """Set up test fixtures"""
        self.zero = Tensor3d()
        self.identity = Tensor3d(1, 0, 0, 0, 1, 0, 0, 0, 1)
        self.t1 = Tensor3d(1, 2, 3, 4, 5, 6, 7, 8, 9)
        self.t2 = Tensor3d(1, 0, 0, 0, 1, 0, 0, 0, 1)

    def test_constructor_3d(self):
        """Test 3D constructor"""
        t = Tensor3d(1, 2, 3, 4, 5, 6, 7, 8, 9)
        self.assertEqual(t.xx, 1)
        self.assertEqual(t.xy, 2)
        self.assertEqual(t.xz, 3)
        self.assertEqual(t.yx, 4)
        self.assertEqual(t.yy, 5)
        self.assertEqual(t.yz, 6)
        self.assertEqual(t.zx, 7)
        self.assertEqual(t.zy, 8)
        self.assertEqual(t.zz, 9)

    def test_properties_3d(self):
        """Test 3D-specific properties"""
        self.assertEqual(self.t1.xz, 3)
        self.assertEqual(self.t1.yz, 6)
        self.assertEqual(self.t1.zx, 7)
        self.assertEqual(self.t1.zy, 8)

        # Test setting
        t = Tensor3d()
        t.xz = 1.0
        t.yz = 2.0
        t.zx = 3.0
        t.zy = 4.0
        self.assertEqual(t.xz, 1.0)
        self.assertEqual(t.yz, 2.0)
        self.assertEqual(t.zx, 3.0)
        self.assertEqual(t.zy, 4.0)

    def test_length_3d(self):
        """Test __len__ for 3D"""
        self.assertEqual(len(self.t1), 9)

    def test_trace_3d(self):
        """Test Trace() for 3D"""
        trace = self.t1.Trace()
        self.assertEqual(trace, 15.0)  # xx + yy + zz = 1 + 5 + 9

    def test_determinant_3d(self):
        """Test Determinant() for 3D"""
        t = Tensor3d(1, 0, 0, 0, 1, 0, 0, 0, 2)
        det = t.Determinant()
        self.assertEqual(det, 2.0)  # Diagonal matrix

    def test_transpose_3d(self):
        """Test Transpose() for 3D"""
        t = self.t1.Transpose()
        self.assertEqual(t.xx, self.t1.xx)
        self.assertEqual(t.xy, self.t1.yx)
        self.assertEqual(t.xz, self.t1.zx)
        self.assertEqual(t.yx, self.t1.xy)
        self.assertEqual(t.yy, self.t1.yy)
        self.assertEqual(t.yz, self.t1.zy)
        self.assertEqual(t.zx, self.t1.xz)
        self.assertEqual(t.zy, self.t1.yz)
        self.assertEqual(t.zz, self.t1.zz)

    def test_vector_multiplication_3d(self):
        """Test tensor-vector multiplication in 3D"""
        v = Vector3d(1.0, 2.0, 3.0)
        result = self.t1 * v
        # [1 2 3] [1]   [1+4+9]    [14]
        # [4 5 6] [2] = [4+10+18] = [32]
        # [7 8 9] [3]   [7+16+27]  [50]
        self.assertAlmostEqual(result.x, 14.0)
        self.assertAlmostEqual(result.y, 32.0)
        self.assertAlmostEqual(result.z, 50.0)

    def test_inverse_3d(self):
        """Test Inverse() for 3D"""
        t = Tensor3d(1, 0, 0, 0, 2, 0, 0, 0, 3)  # Diagonal matrix
        inv = t.Inverse()
        # Verify A * A^-1 = I
        identity = t * inv
        self.assertAlmostEqual(identity.xx, 1.0, places=10)
        self.assertAlmostEqual(identity.yy, 1.0, places=10)
        self.assertAlmostEqual(identity.zz, 1.0, places=10)
        self.assertAlmostEqual(identity.xy, 0.0, places=10)
        self.assertAlmostEqual(identity.xz, 0.0, places=10)
        self.assertAlmostEqual(identity.yz, 0.0, places=10)

    def test_get_row_3d(self):
        """Test getRow() for 3D"""
        row0 = self.t1.getRow(0)
        self.assertEqual(row0.x, 1)
        self.assertEqual(row0.y, 2)
        self.assertEqual(row0.z, 3)

        row1 = self.t1.getRow(1)
        self.assertEqual(row1.x, 4)
        self.assertEqual(row1.y, 5)
        self.assertEqual(row1.z, 6)

        row2 = self.t1.getRow(2)
        self.assertEqual(row2.x, 7)
        self.assertEqual(row2.y, 8)
        self.assertEqual(row2.z, 9)

    def test_get_column_3d(self):
        """Test getColumn() for 3D"""
        col0 = self.t1.getColumn(0)
        self.assertEqual(col0.x, 1)
        self.assertEqual(col0.y, 4)
        self.assertEqual(col0.z, 7)

        col1 = self.t1.getColumn(1)
        self.assertEqual(col1.x, 2)
        self.assertEqual(col1.y, 5)
        self.assertEqual(col1.z, 8)

        col2 = self.t1.getColumn(2)
        self.assertEqual(col2.x, 3)
        self.assertEqual(col2.y, 6)
        self.assertEqual(col2.z, 9)

    def test_diagonal_elements_3d(self):
        """Test diagonalElements() for 3D"""
        diag = self.t1.diagonalElements()
        self.assertEqual(diag.x, 1)
        self.assertEqual(diag.y, 5)
        self.assertEqual(diag.z, 9)

    def test_eigenvalues_3d(self):
        """Test eigenValues() for symmetric matrix"""
        # Use a symmetric matrix for predictable eigenvalues
        t = Tensor3d(2, 0, 0, 0, 3, 0, 0, 0, 1)
        eigs = t.eigenValues()
        # For diagonal matrix, eigenvalues are the diagonal elements
        vals = sorted([eigs.x, eigs.y, eigs.z])
        expected = sorted([2.0, 3.0, 1.0])
        for v1, v2 in zip(vals, expected):
            self.assertAlmostEqual(v1, v2, places=10)

    def test_static_properties_3d(self):
        """Test static properties for 3D"""
        self.assertEqual(Tensor3d.nDimensions, 3)
        self.assertEqual(Tensor3d.numElements, 9)

    def test_3d_methods(self):
        """Test 3D-specific methods"""
        diag3d = self.t1.diagonalElements3D()
        self.assertEqual(diag3d.x, 1)
        self.assertEqual(diag3d.y, 5)
        self.assertEqual(diag3d.z, 9)

        trace3d = self.t1.Trace3D()
        self.assertEqual(trace3d, 15.0)


class TestTensorSymTensorInteraction(unittest.TestCase):
    """Test interaction between Tensor and SymTensor"""

    def test_tensor_from_symtensor_2d(self):
        """Test constructing Tensor from SymTensor"""
        st = SymTensor2d(1, 2, 2, 3, 4)
        t = Tensor2d(st)
        self.assertEqual(t.xx, 1)
        self.assertEqual(t.xy, 2)
        self.assertEqual(t.yx, 2)  # Should be symmetric
        self.assertEqual(t.yy, 3)

    def test_tensor_symtensor_addition_2d(self):
        """Test adding Tensor and SymTensor"""
        t = Tensor2d(1, 2, 3, 4, 5)
        st = SymTensor2d(1, 1, 1, 2, 3)
        result = t + st
        self.assertEqual(result.xx, 2)
        self.assertEqual(result.xy, 3)
        self.assertEqual(result.yx, 4)
        self.assertEqual(result.yy, 6)

    def test_tensor_symtensor_multiplication_2d(self):
        """Test multiplying Tensor and SymTensor"""
        t = Tensor2d(1, 0, 0, 1, 1)
        st = SymTensor2d(2, 0, 0, 2, 1)
        result = t * st
        self.assertEqual(result.xx, 2)
        self.assertEqual(result.yy, 2)


class TestTensorRepr(unittest.TestCase):
    """Test string representation"""

    def test_repr_1d(self):
        """Test __repr__ for 1D"""
        t = Tensor1d(1.0, 2.0, 3.0)
        repr_str = repr(t)
        self.assertIn("Tensor1d", repr_str)

    def test_repr_2d(self):
        """Test __repr__ for 2D"""
        t = Tensor2d(1, 2, 3, 4, 5)
        repr_str = repr(t)
        self.assertIn("Tensor2d", repr_str)

    def test_repr_3d(self):
        """Test __repr__ for 3D"""
        t = Tensor3d(1, 2, 3, 4, 5, 6, 7, 8, 9)
        repr_str = repr(t)
        self.assertIn("Tensor3d", repr_str)


if __name__ == '__main__':
    unittest.main()
