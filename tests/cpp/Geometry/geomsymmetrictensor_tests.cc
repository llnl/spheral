#include "test-basic-exec-policies.hh"
#include "test-utilities.hh"

#include "Geometry/GeomSymmetricTensor.hh"
#include "Geometry/GeomTensor.hh"
#include "Geometry/GeomVector.hh"

using SymTensor = Spheral::GeomSymmetricTensor<3>;
using Tensor = Spheral::GeomTensor<3>;
using Vector = Spheral::GeomVector<3>;

class GeomSymmetricTensorTest : public ::testing::Test {};

// Setting up Typed Test Suite for GeomSymmetricTensor
TYPED_TEST_SUITE_P(GeomSymmetricTensorTypedTest);
template <typename T> class GeomSymmetricTensorTypedTest : public GeomSymmetricTensorTest {};

GPU_TYPED_TEST_P(GeomSymmetricTensorTypedTest, Ctor) {
  using WORK_EXEC_POLICY = TypeParam;
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    SymTensor T1;
    SPHERAL_ASSERT_EQ(T1.xx(), 0.0);
    SPHERAL_ASSERT_EQ(T1.xy(), 0.0);
    SPHERAL_ASSERT_EQ(T1.xz(), 0.0);
    SPHERAL_ASSERT_EQ(T1.yy(), 0.0);
    SPHERAL_ASSERT_EQ(T1.yz(), 0.0);
    SPHERAL_ASSERT_EQ(T1.zz(), 0.0);

    SymTensor T2(1.0, 2.0, 3.0, 2.0, 5.0, 6.0, 3.0, 6.0, 9.0);
    SPHERAL_ASSERT_EQ(T2.xx(), 1.0);
    SPHERAL_ASSERT_EQ(T2.xy(), 2.0);
    SPHERAL_ASSERT_EQ(T2.xz(), 3.0);
    SPHERAL_ASSERT_EQ(T2.yy(), 5.0);
    SPHERAL_ASSERT_EQ(T2.yz(), 6.0);
    SPHERAL_ASSERT_EQ(T2.zz(), 9.0);

    SymTensor T3(T2);
    SPHERAL_ASSERT_EQ(T3.xx(), 1.0);
    SPHERAL_ASSERT_EQ(T3.xy(), 2.0);
    SPHERAL_ASSERT_EQ(T3.xz(), 3.0);
    SPHERAL_ASSERT_EQ(T3.yy(), 5.0);
    SPHERAL_ASSERT_EQ(T3.yz(), 6.0);
    SPHERAL_ASSERT_EQ(T3.zz(), 9.0);

    Tensor GT(1.0, 2.0, 3.0, 2.0, 5.0, 6.0, 3.0, 6.0, 9.0);
    SymTensor T4(GT);
    SPHERAL_ASSERT_EQ(T4.xx(), 1.0);
    SPHERAL_ASSERT_EQ(T4.xy(), 2.0);
    SPHERAL_ASSERT_EQ(T4.xz(), 3.0);
    SPHERAL_ASSERT_EQ(T4.yy(), 5.0);
    SPHERAL_ASSERT_EQ(T4.yz(), 6.0);
    SPHERAL_ASSERT_EQ(T4.zz(), 9.0);
  EXEC_IN_SPACE_END()
}

GPU_TYPED_TEST_P(GeomSymmetricTensorTypedTest, Assignment) {
  using WORK_EXEC_POLICY = TypeParam;
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    SymTensor T1(1.0, 2.0, 3.0, 2.0, 5.0, 6.0, 3.0, 6.0, 9.0);
    SymTensor T2;
    T2 = T1;
    SPHERAL_ASSERT_EQ(T2.xx(), 1.0);
    SPHERAL_ASSERT_EQ(T2.xy(), 2.0);
    SPHERAL_ASSERT_EQ(T2.xz(), 3.0);
    SPHERAL_ASSERT_EQ(T2.yy(), 5.0);
    SPHERAL_ASSERT_EQ(T2.yz(), 6.0);
    SPHERAL_ASSERT_EQ(T2.zz(), 9.0);

    Tensor GT(1.0, 2.0, 3.0, 2.0, 5.0, 6.0, 3.0, 6.0, 9.0);
    T2 = GT;
    SPHERAL_ASSERT_EQ(T2.xx(), 1.0);
    SPHERAL_ASSERT_EQ(T2.xy(), 2.0);
    SPHERAL_ASSERT_EQ(T2.xz(), 3.0);
    SPHERAL_ASSERT_EQ(T2.yy(), 5.0);
    SPHERAL_ASSERT_EQ(T2.yz(), 6.0);
    SPHERAL_ASSERT_EQ(T2.zz(), 9.0);
  EXEC_IN_SPACE_END()
}

GPU_TYPED_TEST_P(GeomSymmetricTensorTypedTest, Accessors) {
  using WORK_EXEC_POLICY = TypeParam;
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    SymTensor T1(1.0, 2.0, 3.0, 2.0, 5.0, 6.0, 3.0, 6.0, 9.0);
    SPHERAL_ASSERT_EQ(T1(0, 0), 1.0);
    SPHERAL_ASSERT_EQ(T1(0, 1), 2.0);
    SPHERAL_ASSERT_EQ(T1(0, 2), 3.0);
    SPHERAL_ASSERT_EQ(T1(1, 0), 2.0);
    SPHERAL_ASSERT_EQ(T1(1, 1), 5.0);
    SPHERAL_ASSERT_EQ(T1(1, 2), 6.0);
    SPHERAL_ASSERT_EQ(T1(2, 0), 3.0);
    SPHERAL_ASSERT_EQ(T1(2, 1), 6.0);
    SPHERAL_ASSERT_EQ(T1(2, 2), 9.0);

    SPHERAL_ASSERT_EQ(T1[0], 1.0);
    SPHERAL_ASSERT_EQ(T1[1], 2.0);
    SPHERAL_ASSERT_EQ(T1[5], 9.0);

    T1(0,0) = 10.0;
    SPHERAL_ASSERT_EQ(T1.xx(), 10.0);
    T1[1] = 20.0;
    SPHERAL_ASSERT_EQ(T1.xy(), 20.0);
  EXEC_IN_SPACE_END()
}

GPU_TYPED_TEST_P(GeomSymmetricTensorTypedTest, Setters) {
  using WORK_EXEC_POLICY = TypeParam;
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    SymTensor T1;
    T1.xx(1.0);
    T1.xy(2.0);
    T1.xz(3.0);
    T1.yy(5.0);
    T1.yz(6.0);
    T1.zz(9.0);
    SPHERAL_ASSERT_EQ(T1.xx(), 1.0);
    SPHERAL_ASSERT_EQ(T1.xy(), 2.0);
    SPHERAL_ASSERT_EQ(T1.xz(), 3.0);
    SPHERAL_ASSERT_EQ(T1.yy(), 5.0);
    SPHERAL_ASSERT_EQ(T1.yz(), 6.0);
    SPHERAL_ASSERT_EQ(T1.zz(), 9.0);
  EXEC_IN_SPACE_END()
}

GPU_TYPED_TEST_P(GeomSymmetricTensorTypedTest, GetSetRowsColumns) {
  using WORK_EXEC_POLICY = TypeParam;
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    SymTensor T1(1.0, 2.0, 3.0, 2.0, 5.0, 6.0, 3.0, 6.0, 9.0);
    Vector row0 = T1.getRow(0);
    SPHERAL_ASSERT_EQ(row0.x(), 1.0);
    SPHERAL_ASSERT_EQ(row0.y(), 2.0);
    SPHERAL_ASSERT_EQ(row0.z(), 3.0);

    Vector col1 = T1.getColumn(1);
    SPHERAL_ASSERT_EQ(col1.x(), 2.0);
    SPHERAL_ASSERT_EQ(col1.y(), 5.0);
    SPHERAL_ASSERT_EQ(col1.z(), 6.0);

    Vector newRow(10.0, 11.0, 12.0);
    T1.setRow(0, newRow);
    SPHERAL_ASSERT_EQ(T1.xx(), 10.0);
    SPHERAL_ASSERT_EQ(T1.xy(), 11.0);
    SPHERAL_ASSERT_EQ(T1.xz(), 12.0);

    Vector newCol(13.0, 14.0, 15.0);
    T1.setColumn(1, newCol);
    SPHERAL_ASSERT_EQ(T1.xy(), 13.0);
    SPHERAL_ASSERT_EQ(T1.yy(), 14.0);
    SPHERAL_ASSERT_EQ(T1.yz(), 15.0);
  EXEC_IN_SPACE_END()
}

GPU_TYPED_TEST_P(GeomSymmetricTensorTypedTest, ZeroIdentity) {
  using WORK_EXEC_POLICY = TypeParam;
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    SymTensor T1(1.0, 2.0, 3.0, 2.0, 5.0, 6.0, 3.0, 6.0, 9.0);
    T1.Zero();
    SPHERAL_ASSERT_EQ(T1.xx(), 0.0);
    SPHERAL_ASSERT_EQ(T1.yy(), 0.0);
    SPHERAL_ASSERT_EQ(T1.zz(), 0.0);

    T1.Identity();
    SPHERAL_ASSERT_EQ(T1.xx(), 1.0);
    SPHERAL_ASSERT_EQ(T1.xy(), 0.0);
    SPHERAL_ASSERT_EQ(T1.yy(), 1.0);
    SPHERAL_ASSERT_EQ(T1.zz(), 1.0);
  EXEC_IN_SPACE_END()
}

GPU_TYPED_TEST_P(GeomSymmetricTensorTypedTest, UnaryMinus) {
  using WORK_EXEC_POLICY = TypeParam;
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    SymTensor T1(1.0, 2.0, 3.0, 2.0, 5.0, 6.0, 3.0, 6.0, 9.0);
    SymTensor T2 = -T1;
    SPHERAL_ASSERT_EQ(T2.xx(), -1.0);
    SPHERAL_ASSERT_EQ(T2.xy(), -2.0);
    SPHERAL_ASSERT_EQ(T2.xz(), -3.0);
    SPHERAL_ASSERT_EQ(T2.yy(), -5.0);
    SPHERAL_ASSERT_EQ(T2.yz(), -6.0);
    SPHERAL_ASSERT_EQ(T2.zz(), -9.0);
  EXEC_IN_SPACE_END()
}

GPU_TYPED_TEST_P(GeomSymmetricTensorTypedTest, AddSub) {
  using WORK_EXEC_POLICY = TypeParam;
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    SymTensor T1(1.0, 2.0, 3.0, 2.0, 5.0, 6.0, 3.0, 6.0, 9.0);
    SymTensor T2(9.0, 8.0, 7.0, 8.0, 5.0, 4.0, 7.0, 4.0, 1.0);

    SymTensor T_add = T1 + T2;
    SPHERAL_ASSERT_EQ(T_add.xx(), 10.0);
    SPHERAL_ASSERT_EQ(T_add.xy(), 10.0);
    SPHERAL_ASSERT_EQ(T_add.xz(), 10.0);
    SPHERAL_ASSERT_EQ(T_add.yy(), 10.0);
    SPHERAL_ASSERT_EQ(T_add.yz(), 10.0);
    SPHERAL_ASSERT_EQ(T_add.zz(), 10.0);

    SymTensor T_sub = T1 - T2;
    SPHERAL_ASSERT_EQ(T_sub.xx(), -8.0);
    SPHERAL_ASSERT_EQ(T_sub.xy(), -6.0);
    SPHERAL_ASSERT_EQ(T_sub.xz(), -4.0);
    SPHERAL_ASSERT_EQ(T_sub.yy(), 0.0);
    SPHERAL_ASSERT_EQ(T_sub.yz(), 2.0);
    SPHERAL_ASSERT_EQ(T_sub.zz(), 8.0);

    Tensor GT(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
    Tensor T_add_gen = T1 + GT;
    SPHERAL_ASSERT_EQ(T_add_gen.xx(), 2.0);
    SPHERAL_ASSERT_EQ(T_add_gen.xy(), 4.0);
    SPHERAL_ASSERT_EQ(T_add_gen.xz(), 6.0);
    SPHERAL_ASSERT_EQ(T_add_gen.yx(), 6.0);
    SPHERAL_ASSERT_EQ(T_add_gen.yy(), 10.0);
    SPHERAL_ASSERT_EQ(T_add_gen.yz(), 12.0);
    SPHERAL_ASSERT_EQ(T_add_gen.zx(), 10.0);
    SPHERAL_ASSERT_EQ(T_add_gen.zy(), 14.0);
    SPHERAL_ASSERT_EQ(T_add_gen.zz(), 18.0);

    Tensor T_sub_gen = T1 - GT;
    SPHERAL_ASSERT_EQ(T_sub_gen.xx(), 0.0);
    SPHERAL_ASSERT_EQ(T_sub_gen.xy(), 0.0);
    SPHERAL_ASSERT_EQ(T_sub_gen.xz(), 0.0);
    SPHERAL_ASSERT_EQ(T_sub_gen.yx(), -2.0);
    SPHERAL_ASSERT_EQ(T_sub_gen.yy(), 0.0);
    SPHERAL_ASSERT_EQ(T_sub_gen.yz(), 0.0);
    SPHERAL_ASSERT_EQ(T_sub_gen.zx(), -4.0);
    SPHERAL_ASSERT_EQ(T_sub_gen.zy(), -2.0);
    SPHERAL_ASSERT_EQ(T_sub_gen.zz(), 0.0);
  EXEC_IN_SPACE_END()
}

GPU_TYPED_TEST_P(GeomSymmetricTensorTypedTest, ScalarMulDiv) {
  using WORK_EXEC_POLICY = TypeParam;
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    SymTensor T1(1.0, 2.0, 3.0, 2.0, 5.0, 6.0, 3.0, 6.0, 9.0);
    double s = 2.0;

    SymTensor T_mul = T1 * s;
    SPHERAL_ASSERT_EQ(T_mul.xx(), 2.0);
    SPHERAL_ASSERT_EQ(T_mul.xy(), 4.0);
    SPHERAL_ASSERT_EQ(T_mul.xz(), 6.0);
    SPHERAL_ASSERT_EQ(T_mul.yy(), 10.0);
    SPHERAL_ASSERT_EQ(T_mul.yz(), 12.0);
    SPHERAL_ASSERT_EQ(T_mul.zz(), 18.0);

    SymTensor T_div = T1 / s;
    SPHERAL_ASSERT_EQ(T_div.xx(), 0.5);
    SPHERAL_ASSERT_EQ(T_div.xy(), 1.0);
    SPHERAL_ASSERT_EQ(T_div.xz(), 1.5);
    SPHERAL_ASSERT_EQ(T_div.yy(), 2.5);
    SPHERAL_ASSERT_EQ(T_div.yz(), 3.0);
    SPHERAL_ASSERT_EQ(T_div.zz(), 4.5);
  EXEC_IN_SPACE_END()
}

GPU_TYPED_TEST_P(GeomSymmetricTensorTypedTest, InPlaceAddSub) {
  using WORK_EXEC_POLICY = TypeParam;
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    SymTensor T1(1.0, 2.0, 3.0, 2.0, 5.0, 6.0, 3.0, 6.0, 9.0);
    SymTensor T2(9.0, 8.0, 7.0, 8.0, 5.0, 4.0, 7.0, 4.0, 1.0);
    
    T1 += T2;
    SPHERAL_ASSERT_EQ(T1.xx(), 10.0);
    SPHERAL_ASSERT_EQ(T1.xy(), 10.0);
    SPHERAL_ASSERT_EQ(T1.xz(), 10.0);
    SPHERAL_ASSERT_EQ(T1.yy(), 10.0);
    SPHERAL_ASSERT_EQ(T1.yz(), 10.0);
    SPHERAL_ASSERT_EQ(T1.zz(), 10.0);

    T1 -= T2;
    SPHERAL_ASSERT_EQ(T1.xx(), 1.0);
    SPHERAL_ASSERT_EQ(T1.xy(), 2.0);
    SPHERAL_ASSERT_EQ(T1.xz(), 3.0);
    SPHERAL_ASSERT_EQ(T1.yy(), 5.0);
    SPHERAL_ASSERT_EQ(T1.yz(), 6.0);
    SPHERAL_ASSERT_EQ(T1.zz(), 9.0);
  EXEC_IN_SPACE_END()
}

GPU_TYPED_TEST_P(GeomSymmetricTensorTypedTest, InPlaceScalarMulDiv) {
  using WORK_EXEC_POLICY = TypeParam;
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    SymTensor T1(1.0, 2.0, 3.0, 2.0, 5.0, 6.0, 3.0, 6.0, 9.0);
    double s = 2.0;

    T1 *= s;
    SPHERAL_ASSERT_EQ(T1.xx(), 2.0);
    SPHERAL_ASSERT_EQ(T1.xy(), 4.0);
    SPHERAL_ASSERT_EQ(T1.xz(), 6.0);
    SPHERAL_ASSERT_EQ(T1.yy(), 10.0);
    SPHERAL_ASSERT_EQ(T1.yz(), 12.0);
    SPHERAL_ASSERT_EQ(T1.zz(), 18.0);

    T1 /= s;
    SPHERAL_ASSERT_EQ(T1.xx(), 1.0);
    SPHERAL_ASSERT_EQ(T1.xy(), 2.0);
    SPHERAL_ASSERT_EQ(T1.xz(), 3.0);
    SPHERAL_ASSERT_EQ(T1.yy(), 5.0);
    SPHERAL_ASSERT_EQ(T1.yz(), 6.0);
    SPHERAL_ASSERT_EQ(T1.zz(), 9.0);
  EXEC_IN_SPACE_END()
}

GPU_TYPED_TEST_P(GeomSymmetricTensorTypedTest, Comparison) {
  using WORK_EXEC_POLICY = TypeParam;
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    SymTensor T1(1.0, 2.0, 3.0, 2.0, 5.0, 6.0, 3.0, 6.0, 9.0);
    SymTensor T2(1.0, 2.0, 3.0, 2.0, 5.0, 6.0, 3.0, 6.0, 9.0);
    SymTensor T3(9.0, 8.0, 7.0, 8.0, 5.0, 4.0, 7.0, 4.0, 1.0);

    SPHERAL_ASSERT_TRUE(T1 == T2);
    SPHERAL_ASSERT_TRUE(T1 != T3);
    SPHERAL_ASSERT_TRUE(T1 < T3);
    SPHERAL_ASSERT_TRUE(T3 > T1);
    SPHERAL_ASSERT_TRUE(T1 <= T2);
    SPHERAL_ASSERT_TRUE(T1 >= T2);
  EXEC_IN_SPACE_END()
}

GPU_TYPED_TEST_P(GeomSymmetricTensorTypedTest, Transpose) {
  using WORK_EXEC_POLICY = TypeParam;
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    SymTensor T1(1.0, 2.0, 3.0, 2.0, 5.0, 6.0, 3.0, 6.0, 9.0);
    SymTensor T_trans = T1.Transpose();
    SPHERAL_ASSERT_EQ(T_trans.xx(), 1.0);
    SPHERAL_ASSERT_EQ(T_trans.xy(), 2.0);
    SPHERAL_ASSERT_EQ(T_trans.xz(), 3.0);
    SPHERAL_ASSERT_EQ(T_trans.yy(), 5.0);
    SPHERAL_ASSERT_EQ(T_trans.yz(), 6.0);
    SPHERAL_ASSERT_EQ(T_trans.zz(), 9.0);
  EXEC_IN_SPACE_END()
}

GPU_TYPED_TEST_P(GeomSymmetricTensorTypedTest, TraceDeterminant) {
  using WORK_EXEC_POLICY = TypeParam;
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    SymTensor T1(1.0, 2.0, 3.0, 2.0, 5.0, 6.0, 3.0, 6.0, 9.0);
    SPHERAL_ASSERT_EQ(T1.Trace(), 15.0);
    SPHERAL_ASSERT_EQ(T1.Determinant(), 0.0);
  EXEC_IN_SPACE_END()
}

GPU_TYPED_TEST_P(GeomSymmetricTensorTypedTest, DotProduct) {
  using WORK_EXEC_POLICY = TypeParam;
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    SymTensor T1(1.0, 2.0, 3.0, 2.0, 5.0, 6.0, 3.0, 6.0, 9.0);
    SymTensor T2(9.0, 8.0, 7.0, 8.0, 5.0, 4.0, 7.0, 4.0, 1.0);
    Tensor T_dot = T1.dot(T2);
    // Expected values calculated manually
    SPHERAL_ASSERT_EQ(T_dot.xx(), 46.0);
    SPHERAL_ASSERT_EQ(T_dot.xy(), 30.0);
    SPHERAL_ASSERT_EQ(T_dot.xz(), 18.0);
    SPHERAL_ASSERT_EQ(T_dot.yx(), 100.0);
    SPHERAL_ASSERT_EQ(T_dot.yy(), 65.0);
    SPHERAL_ASSERT_EQ(T_dot.yz(), 40.0);
    SPHERAL_ASSERT_EQ(T_dot.zx(), 138.0);
    SPHERAL_ASSERT_EQ(T_dot.zy(), 90.0);
    SPHERAL_ASSERT_EQ(T_dot.zz(), 54.0);

    Vector V1(1.0, 2.0, 3.0);
    Vector V_dot = T1.dot(V1);
    SPHERAL_ASSERT_EQ(V_dot.x(), 14.0);
    SPHERAL_ASSERT_EQ(V_dot.y(), 30.0);
    SPHERAL_ASSERT_EQ(V_dot.z(), 42.0);
  EXEC_IN_SPACE_END()
}

GPU_TYPED_TEST_P(GeomSymmetricTensorTypedTest, DoubleDotProduct) {
  using WORK_EXEC_POLICY = TypeParam;
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    SymTensor T1(1.0, 2.0, 3.0, 2.0, 5.0, 6.0, 3.0, 6.0, 9.0);
    SymTensor T2(9.0, 8.0, 7.0, 8.0, 5.0, 4.0, 7.0, 4.0, 1.0);
    double ddot = T1.doubledot(T2);
    // Expected value calculated manually
    SPHERAL_ASSERT_EQ(ddot, 1.0*9.0 + 2.0*8.0 + 3.0*7.0 +
                            2.0*8.0 + 5.0*5.0 + 6.0*4.0 +
                            3.0*7.0 + 6.0*4.0 + 9.0*1.0);
  EXEC_IN_SPACE_END()
}

GPU_TYPED_TEST_P(GeomSymmetricTensorTypedTest, Square) {
  using WORK_EXEC_POLICY = TypeParam;
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    SymTensor T1(1.0, 2.0, 3.0, 2.0, 5.0, 6.0, 3.0, 6.0, 9.0);
    SymTensor T_sq = T1.square();
    // Expected values calculated manually
    SPHERAL_ASSERT_EQ(T_sq.xx(), 14.0);
    SPHERAL_ASSERT_EQ(T_sq.xy(), 30.0);
    SPHERAL_ASSERT_EQ(T_sq.xz(), 42.0);
    SPHERAL_ASSERT_EQ(T_sq.yy(), 65.0);
    SPHERAL_ASSERT_EQ(T_sq.yz(), 90.0);
    SPHERAL_ASSERT_EQ(T_sq.zz(), 126.0);
  EXEC_IN_SPACE_END()
}

GPU_TYPED_TEST_P(GeomSymmetricTensorTypedTest, SquareElements) {
  using WORK_EXEC_POLICY = TypeParam;
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    SymTensor T1(1.0, 2.0, 3.0, 2.0, 5.0, 6.0, 3.0, 6.0, 9.0);
    SymTensor T_sq_elem = T1.squareElements();
    SPHERAL_ASSERT_EQ(T_sq_elem.xx(), 1.0);
    SPHERAL_ASSERT_EQ(T_sq_elem.xy(), 4.0);
    SPHERAL_ASSERT_EQ(T_sq_elem.xz(), 9.0);
    SPHERAL_ASSERT_EQ(T_sq_elem.yy(), 25.0);
    SPHERAL_ASSERT_EQ(T_sq_elem.yz(), 36.0);
    SPHERAL_ASSERT_EQ(T_sq_elem.zz(), 81.0);
  EXEC_IN_SPACE_END()
}

GPU_TYPED_TEST_P(GeomSymmetricTensorTypedTest, EigenVector) {
  using WORK_EXEC_POLICY = TypeParam;
  const double inv_sqt = 1./std::sqrt(2.);
  const Vector ref_vals(-3., -1., 1.);
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    SymTensor T1(0.0, 0.0, 1.0, 0.0, -3.0, 0.0, 1.0, 0.0, 0.0);
    Spheral::EigenStruct<3> result = T1.eigenVectors();
    Tensor vecs = result.eigenVectors;
    Vector vals = result.eigenValues;
    int found = 0;
    for (auto& ref : ref_vals) {
      for (auto& v : vals) {
        if (std::abs(ref - v) < 1.E-12) {
          found += 1;
        }
      }
    }
    SPHERAL_ASSERT_EQ(found, 3);
  EXEC_IN_SPACE_END()
}

REGISTER_TYPED_TEST_SUITE_P(GeomSymmetricTensorTypedTest, Ctor, Assignment, Accessors,
                            Setters, GetSetRowsColumns, ZeroIdentity, UnaryMinus, AddSub,
                            ScalarMulDiv, InPlaceAddSub, InPlaceScalarMulDiv, Comparison,
                            Transpose, TraceDeterminant, DotProduct, DoubleDotProduct,
                            Square, SquareElements, EigenVector);

INSTANTIATE_TYPED_TEST_SUITE_P(GeomSymmetricTensor, GeomSymmetricTensorTypedTest,
                               typename Spheral::Test<EXEC_TYPES>::Types, );
