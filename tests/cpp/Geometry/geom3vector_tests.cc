#include "test-basic-exec-policies.hh"
#include "test-utilities.hh"

#include "Geometry/Geom3Vector.hh"

using Vector = Spheral::Geom3Vector;

class Geom3VectorTest : public ::testing::Test {};

// Setting up Typed Test Suite for Geom3Vector
TYPED_TEST_SUITE_P(Geom3VectorTypedTest);
template <typename T> class Geom3VectorTypedTest : public Geom3VectorTest {};

/*
 * Test case for the Geom3Vector constructor.
 * This test verifies that the constructor correctly initializes the vector's
 * components. A vector `v` is created with components (1.0, 2.0, 3.0), and the
 * test asserts that v.x(), v.y(), and v.z() return these values.
 */
GPU_TYPED_TEST_P(Geom3VectorTypedTest, Ctor) {
  using WORK_EXEC_POLICY = TypeParam;
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    Vector v(1.0, 2.0, 3.0);
    SPHERAL_ASSERT_EQ(v.x(), 1.0);
    SPHERAL_ASSERT_EQ(v.y(), 2.0);
    SPHERAL_ASSERT_EQ(v.z(), 3.0);
  EXEC_IN_SPACE_END()
}

/*
 * Test case for the Geom3Vector assignment operators.
 * This test verifies that vector assignment works correctly.
 * It first assigns a scalar value to all components of a vector `v` and asserts
 * the result. It then assigns another vector `u` to `v` and asserts that `v`'s
 * components match `u`'s.
 */
GPU_TYPED_TEST_P(Geom3VectorTypedTest, Assignment) {
  using WORK_EXEC_POLICY = TypeParam;
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    double s = 5.0;
    Vector v(s, s, s);
    SPHERAL_ASSERT_EQ(v.x(), 5.0);
    SPHERAL_ASSERT_EQ(v.y(), 5.0);
    SPHERAL_ASSERT_EQ(v.z(), 5.0);

    Vector u(1.0, 2.0, 3.0);
    v = u;
    SPHERAL_ASSERT_EQ(v.x(), 1.0);
    SPHERAL_ASSERT_EQ(v.y(), 2.0);
    SPHERAL_ASSERT_EQ(v.z(), 3.0);
  EXEC_IN_SPACE_END()
}

/*
 * Test case for the Geom3Vector element accessors.
 * This test verifies that the parenthesis operators
 * correctly access and modify vector components. It asserts that `v(i)`
 * return the correct values, and then modifies components using these
 * accessors and asserts the new values.
 */
GPU_TYPED_TEST_P(Geom3VectorTypedTest, Accessors) {
  using WORK_EXEC_POLICY = TypeParam;
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    Vector v(1.0, 2.0, 3.0);
    SPHERAL_ASSERT_EQ(v(0), 1.0);
    SPHERAL_ASSERT_EQ(v(1), 2.0);
    SPHERAL_ASSERT_EQ(v(2), 3.0);

    v(0) = 4.0;
    SPHERAL_ASSERT_EQ(v.x(), 4.0);
  EXEC_IN_SPACE_END()
}

/*
 * Test case for the Geom3Vector component setters.
 * This test verifies that the `x(val)`, `y(val)`, and `z(val)` methods
 * correctly set the vector's components. A vector `v` is modified using these
 * setters, and the test asserts that the components have been updated.
 */
GPU_TYPED_TEST_P(Geom3VectorTypedTest, Setters) {
  using WORK_EXEC_POLICY = TypeParam;
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    Vector v;
    v.x(1.0);
    v.y(2.0);
    v.z(3.0);
    SPHERAL_ASSERT_EQ(v.x(), 1.0);
    SPHERAL_ASSERT_EQ(v.y(), 2.0);
    SPHERAL_ASSERT_EQ(v.z(), 3.0);
  EXEC_IN_SPACE_END()
}

/*
 * Test case for the Geom3Vector's Zero() method.
 * This test verifies that the Zero() method correctly sets all components of
 * the vector to 0.0. A vector `v` is initialized with non-zero values, Zero()
 * is called, and the test asserts that all components are now 0.0.
 */
GPU_TYPED_TEST_P(Geom3VectorTypedTest, Zero) {
  using WORK_EXEC_POLICY = TypeParam;
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    Vector v(1.0, 2.0, 3.0);
    v.Zero();
    SPHERAL_ASSERT_EQ(v.x(), 0.0);
    SPHERAL_ASSERT_EQ(v.y(), 0.0);
    SPHERAL_ASSERT_EQ(v.z(), 0.0);
  EXEC_IN_SPACE_END()
}
/*
 * Test case for the Geom3Vector unary minus operator.
 * This test verifies that the unary minus operator correctly negates all
 * components of the vector. A vector `v2` is created by negating `v1`, and the
 * test asserts that each component of `v2` is the negated value of the
 * corresponding component in `v1`.
 */
GPU_TYPED_TEST_P(Geom3VectorTypedTest, UnaryMinus) {
  using WORK_EXEC_POLICY = TypeParam;
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    Vector v1(1.0, -2.0, 3.0);
    Vector v2 = -v1;
    SPHERAL_ASSERT_EQ(v2.x(), -1.0);
    SPHERAL_ASSERT_EQ(v2.y(), 2.0);
    SPHERAL_ASSERT_EQ(v2.z(), -3.0);
  EXEC_IN_SPACE_END()
}

/*
 * Test case for Geom3Vector addition and subtraction.
 * This test verifies that the addition and subtraction operators produce the
 * correct results. It asserts that the components of `v1 + v2` and `v2 - v1`
 * are the sums and differences of the components of `v1` and `v2`.
 */
GPU_TYPED_TEST_P(Geom3VectorTypedTest, AddSub) {
  using WORK_EXEC_POLICY = TypeParam;
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    Vector v1(1.0, 2.0, 3.0);
    Vector v2(4.0, 5.0, 6.0);

    Vector v_add = v1 + v2;
    SPHERAL_ASSERT_EQ(v_add.x(), 5.0);
    SPHERAL_ASSERT_EQ(v_add.y(), 7.0);
    SPHERAL_ASSERT_EQ(v_add.z(), 9.0);

    Vector v_sub = v2 - v1;
    SPHERAL_ASSERT_EQ(v_sub.x(), 3.0);
    SPHERAL_ASSERT_EQ(v_sub.y(), 3.0);
    SPHERAL_ASSERT_EQ(v_sub.z(), 3.0);
  EXEC_IN_SPACE_END()
}

/*
 * Test case for Geom3Vector scalar multiplication and division.
 * This test verifies that scalar multiplication and division operators produce
 * the correct results. It asserts that the components of `v * s` and `v / s`
 * are the product and quotient of the components of `v` and a scalar `s`.
 */
GPU_TYPED_TEST_P(Geom3VectorTypedTest, ScalarMulDiv) {
  using WORK_EXEC_POLICY = TypeParam;
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    Vector v(1.0, 2.0, 3.0);
    double s = 2.0;

    Vector v_mul = v * s;
    SPHERAL_ASSERT_EQ(v_mul.x(), 2.0);
    SPHERAL_ASSERT_EQ(v_mul.y(), 4.0);
    SPHERAL_ASSERT_EQ(v_mul.z(), 6.0);

    Vector v_div = v / s;
    SPHERAL_ASSERT_EQ(v_div.x(), 0.5);
    SPHERAL_ASSERT_EQ(v_div.y(), 1.0);
    SPHERAL_ASSERT_EQ(v_div.z(), 1.5);
  EXEC_IN_SPACE_END()
}

/*
 * Test case for Geom3Vector in-place addition and subtraction.
 * This test verifies that the `+=` and `-=` operators produce the correct
 * results. It asserts that `v1 += v2` and `v1 -= v2` correctly modify the
 * components of `v1`.
 */
GPU_TYPED_TEST_P(Geom3VectorTypedTest, InPlaceAddSub) {
  using WORK_EXEC_POLICY = TypeParam;
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    Vector v1(1.0, 2.0, 3.0);
    Vector v2(4.0, 5.0, 6.0);
    
    v1 += v2;
    SPHERAL_ASSERT_EQ(v1.x(), 5.0);
    SPHERAL_ASSERT_EQ(v1.y(), 7.0);
    SPHERAL_ASSERT_EQ(v1.z(), 9.0);

    v1 -= v2;
    SPHERAL_ASSERT_EQ(v1.x(), 1.0);
    SPHERAL_ASSERT_EQ(v1.y(), 2.0);
    SPHERAL_ASSERT_EQ(v1.z(), 3.0);
  EXEC_IN_SPACE_END()
}

/*
 * Test case for Geom3Vector in-place scalar multiplication and division.
 * This test verifies that the `*=` and `/=` operators produce the correct
 * results. It asserts that `v *= s` and `v /= s` correctly modify the
 * components of `v`.
 */
GPU_TYPED_TEST_P(Geom3VectorTypedTest, InPlaceScalarMulDiv) {
  using WORK_EXEC_POLICY = TypeParam;
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    Vector v(1.0, 2.0, 3.0);
    double s = 2.0;

    v *= s;
    SPHERAL_ASSERT_EQ(v.x(), 2.0);
    SPHERAL_ASSERT_EQ(v.y(), 4.0);
    SPHERAL_ASSERT_EQ(v.z(), 6.0);

    v /= s;
    SPHERAL_ASSERT_EQ(v.x(), 1.0);
    SPHERAL_ASSERT_EQ(v.y(), 2.0);
    SPHERAL_ASSERT_EQ(v.z(), 3.0);
  EXEC_IN_SPACE_END()
}

/*
 * Test case for Geom3Vector comparison operators.
 * This test verifies that the `==`, `!=`, `<`, `>`, `<=`, and `>=` operators
 * produce the correct boolean results. It asserts the outcomes of comparing
 * several vectors.
 */
GPU_TYPED_TEST_P(Geom3VectorTypedTest, Comparison) {
  using WORK_EXEC_POLICY = TypeParam;
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    Vector v1(1.0, 2.0, 3.0);
    Vector v2(1.0, 2.0, 3.0);
    Vector v3(4.0, 5.0, 6.0);
    Vector v4(0.0, 0.0, 0.0);

    SPHERAL_ASSERT_TRUE(v1 == v2);
    SPHERAL_ASSERT_TRUE(v1 != v3);
    SPHERAL_ASSERT_TRUE(v1 < v3);
    SPHERAL_ASSERT_TRUE(v3 > v1);
    SPHERAL_ASSERT_TRUE(v1 <= v2);
    SPHERAL_ASSERT_TRUE(v1 >= v2);
    SPHERAL_ASSERT_TRUE(v4 < v1);
  EXEC_IN_SPACE_END()
}

/*
 * Test case for the Geom3Vector dot product.
 * This test verifies that the `dot()` method calculates the correct dot product
 * of two vectors. It asserts that `v1.dot(v2)` equals the sum of the products
 * of their corresponding components.
 */
GPU_TYPED_TEST_P(Geom3VectorTypedTest, Dot) {
  using WORK_EXEC_POLICY = TypeParam;
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    Vector v1(1.0, 2.0, 3.0);
    Vector v2(4.0, 5.0, 6.0);
    double dot_product = v1.dot(v2);
    SPHERAL_ASSERT_EQ(dot_product, 1.0*4.0 + 2.0*5.0 + 3.0*6.0); // 4 + 10 + 18 = 32
  EXEC_IN_SPACE_END()
}

/*
 * Test case for the Geom3Vector cross product.
 * This test verifies that the `cross()` method calculates the correct cross
 * product of two vectors. It asserts that the cross product of the x-axis and
 * y-axis vectors results in the z-axis vector.
 */
GPU_TYPED_TEST_P(Geom3VectorTypedTest, Cross) {
  using WORK_EXEC_POLICY = TypeParam;
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    Vector v1(1.0, 0.0, 0.0);
    Vector v2(0.0, 1.0, 0.0);
    Vector v_cross = v1.cross(v2);
    SPHERAL_ASSERT_EQ(v_cross.x(), 0.0);
    SPHERAL_ASSERT_EQ(v_cross.y(), 0.0);
    SPHERAL_ASSERT_EQ(v_cross.z(), 1.0);
  EXEC_IN_SPACE_END()
}

/*
 * Test case for the Geom3Vector magnitude.
 * This test verifies that the `magnitude()` and `magnitude2()` methods
 * calculate the correct values. It asserts that `v.magnitude2()` (squared
 * magnitude) and `v.magnitude()` are correct for a 3-4-5 right triangle vector.
 */
GPU_TYPED_TEST_P(Geom3VectorTypedTest, Magnitude) {
  using WORK_EXEC_POLICY = TypeParam;
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    Vector v(3.0, 4.0, 0.0);
    SPHERAL_ASSERT_EQ(v.magnitude2(), 25.0);
    SPHERAL_ASSERT_EQ(v.magnitude(), 5.0);
  EXEC_IN_SPACE_END()
}

/*
 * Test case for the Geom3Vector unit vector.
 * This test verifies that the `unitVector()` method calculates the correct unit
 * vector. It asserts that the components of the unit vector are correct and
 * that its magnitude is 1.0.
 */
GPU_TYPED_TEST_P(Geom3VectorTypedTest, UnitVector) {
  using WORK_EXEC_POLICY = TypeParam;
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    Vector v(3.0, 4.0, 0.0);
    Vector u = v.unitVector();
    SPHERAL_ASSERT_FLOAT_EQ(u.x(), 0.6);
    SPHERAL_ASSERT_FLOAT_EQ(u.y(), 0.8);
    SPHERAL_ASSERT_FLOAT_EQ(u.z(), 0.0);
    SPHERAL_ASSERT_FLOAT_EQ(u.magnitude(), 1.0);
  EXEC_IN_SPACE_END()
}

/*
 * Test case for the Geom3Vector element-wise statistics.
 * This test verifies that the `minElement`, `maxElement`, and
 * `sumElements` methods calculate the correct values. It asserts the results of
 * these methods for a vector with positive and negative components.
 */
GPU_TYPED_TEST_P(Geom3VectorTypedTest, ElementStats) {
  using WORK_EXEC_POLICY = TypeParam;
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    Vector v(-1.0, 2.0, -3.0);
    SPHERAL_ASSERT_EQ(v.minElement(), -3.0);
    SPHERAL_ASSERT_EQ(v.maxElement(), 2.0);
    SPHERAL_ASSERT_EQ(v.sumElements(), -2.0);
  EXEC_IN_SPACE_END()
}

REGISTER_TYPED_TEST_SUITE_P(Geom3VectorTypedTest, Ctor, Assignment, Accessors,
                            Setters, Zero, UnaryMinus, AddSub, ScalarMulDiv,
                            InPlaceAddSub, InPlaceScalarMulDiv, Comparison, Dot,
                            Cross, Magnitude, UnitVector, ElementStats);


INSTANTIATE_TYPED_TEST_SUITE_P(Geom3Vector, Geom3VectorTypedTest,
                               typename Spheral::Test<EXEC_TYPES>::Types, );
