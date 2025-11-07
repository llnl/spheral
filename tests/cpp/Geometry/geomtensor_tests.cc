#include "test-basic-exec-policies.hh"
#include "test-utilities.hh"

#include "Geometry/GeomTensor.hh"
#include "Geometry/GeomSymmetricTensor.hh"
#include "Geometry/GeomVector.hh"

using Tensor = Spheral::GeomTensor<3>;
using SymTensor = Spheral::GeomSymmetricTensor<3>;
using Vector = Spheral::GeomVector<3>;

class GeomTensorTest : public ::testing::Test {};

// Setting up Typed Test Suite for GeomTensor
TYPED_TEST_SUITE_P(GeomTensorTypedTest);
template <typename T> class GeomTensorTypedTest : public GeomTensorTest {};

GPU_TYPED_TEST_P(GeomTensorTypedTest, Ctor) {
  using WORK_EXEC_POLICY = TypeParam;
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    Tensor T1;
    SPHERAL_ASSERT_EQ(T1.xx(), 0.0);
    SPHERAL_ASSERT_EQ(T1.xy(), 0.0);
    SPHERAL_ASSERT_EQ(T1.xz(), 0.0);
    SPHERAL_ASSERT_EQ(T1.yx(), 0.0);
    SPHERAL_ASSERT_EQ(T1.yy(), 0.0);
    SPHERAL_ASSERT_EQ(T1.yz(), 0.0);
    SPHERAL_ASSERT_EQ(T1.zx(), 0.0);
    SPHERAL_ASSERT_EQ(T1.zy(), 0.0);
    SPHERAL_ASSERT_EQ(T1.zz(), 0.0);

    Tensor T2(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
    SPHERAL_ASSERT_EQ(T2.xx(), 1.0);
    SPHERAL_ASSERT_EQ(T2.xy(), 2.0);
    SPHERAL_ASSERT_EQ(T2.xz(), 3.0);
    SPHERAL_ASSERT_EQ(T2.yx(), 4.0);
    SPHERAL_ASSERT_EQ(T2.yy(), 5.0);
    SPHERAL_ASSERT_EQ(T2.yz(), 6.0);
    SPHERAL_ASSERT_EQ(T2.zx(), 7.0);
    SPHERAL_ASSERT_EQ(T2.zy(), 8.0);
    SPHERAL_ASSERT_EQ(T2.zz(), 9.0);

    Tensor T3(T2);
    SPHERAL_ASSERT_EQ(T3.xx(), 1.0);
    SPHERAL_ASSERT_EQ(T3.xy(), 2.0);
    SPHERAL_ASSERT_EQ(T3.xz(), 3.0);
    SPHERAL_ASSERT_EQ(T3.yx(), 4.0);
    SPHERAL_ASSERT_EQ(T3.yy(), 5.0);
    SPHERAL_ASSERT_EQ(T3.yz(), 6.0);
    SPHERAL_ASSERT_EQ(T3.zx(), 7.0);
    SPHERAL_ASSERT_EQ(T3.zy(), 8.0);
    SPHERAL_ASSERT_EQ(T3.zz(), 9.0);

    SymTensor ST(1.0, 2.0, 3.0,
                 2.0, 5.0, 6.0,
                 3.0, 6.0, 9.0);
    Tensor T4(ST);
    SPHERAL_ASSERT_EQ(T4.xx(), 1.0);
    SPHERAL_ASSERT_EQ(T4.xy(), 2.0);
    SPHERAL_ASSERT_EQ(T4.xz(), 3.0);
    SPHERAL_ASSERT_EQ(T4.yx(), 2.0);
    SPHERAL_ASSERT_EQ(T4.yy(), 5.0);
    SPHERAL_ASSERT_EQ(T4.yz(), 6.0);
    SPHERAL_ASSERT_EQ(T4.zx(), 3.0);
    SPHERAL_ASSERT_EQ(T4.zy(), 6.0);
    SPHERAL_ASSERT_EQ(T4.zz(), 9.0);
  EXEC_IN_SPACE_END()
}

GPU_TYPED_TEST_P(GeomTensorTypedTest, Assignment) {
  using WORK_EXEC_POLICY = TypeParam;
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    Tensor T1(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
    Tensor T2;
    T2 = T1;
    SPHERAL_ASSERT_EQ(T2.xx(), 1.0);
    SPHERAL_ASSERT_EQ(T2.xy(), 2.0);
    SPHERAL_ASSERT_EQ(T2.xz(), 3.0);
    SPHERAL_ASSERT_EQ(T2.yx(), 4.0);
    SPHERAL_ASSERT_EQ(T2.yy(), 5.0);
    SPHERAL_ASSERT_EQ(T2.yz(), 6.0);
    SPHERAL_ASSERT_EQ(T2.zx(), 7.0);
    SPHERAL_ASSERT_EQ(T2.zy(), 8.0);
    SPHERAL_ASSERT_EQ(T2.zz(), 9.0);

    SymTensor ST(1.0, 2.0, 3.0,
                 2.0, 5.0, 6.0,
                 3.0, 6.0, 9.0);
    T2 = ST;
    SPHERAL_ASSERT_EQ(T2.xx(), 1.0);
    SPHERAL_ASSERT_EQ(T2.xy(), 2.0);
    SPHERAL_ASSERT_EQ(T2.xz(), 3.0);
    SPHERAL_ASSERT_EQ(T2.yx(), 2.0);
    SPHERAL_ASSERT_EQ(T2.yy(), 5.0);
    SPHERAL_ASSERT_EQ(T2.yz(), 6.0);
    SPHERAL_ASSERT_EQ(T2.zx(), 3.0);
    SPHERAL_ASSERT_EQ(T2.zy(), 6.0);
    SPHERAL_ASSERT_EQ(T2.zz(), 9.0);
  EXEC_IN_SPACE_END()
}

GPU_TYPED_TEST_P(GeomTensorTypedTest, Accessors) {
  using WORK_EXEC_POLICY = TypeParam;
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    Tensor T1(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
    SPHERAL_ASSERT_EQ(T1(0, 0), 1.0);
    SPHERAL_ASSERT_EQ(T1(0, 1), 2.0);
    SPHERAL_ASSERT_EQ(T1(0, 2), 3.0);
    SPHERAL_ASSERT_EQ(T1(1, 0), 4.0);
    SPHERAL_ASSERT_EQ(T1(1, 1), 5.0);
    SPHERAL_ASSERT_EQ(T1(1, 2), 6.0);
    SPHERAL_ASSERT_EQ(T1(2, 0), 7.0);
    SPHERAL_ASSERT_EQ(T1(2, 1), 8.0);
    SPHERAL_ASSERT_EQ(T1(2, 2), 9.0);

    SPHERAL_ASSERT_EQ(T1[0], 1.0);
    SPHERAL_ASSERT_EQ(T1[1], 2.0);
    SPHERAL_ASSERT_EQ(T1[8], 9.0);

    T1(0,0) = 10.0;
    SPHERAL_ASSERT_EQ(T1.xx(), 10.0);
    T1[1] = 20.0;
    SPHERAL_ASSERT_EQ(T1.xy(), 20.0);
  EXEC_IN_SPACE_END()
}

GPU_TYPED_TEST_P(GeomTensorTypedTest, Setters) {
  using WORK_EXEC_POLICY = TypeParam;
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    Tensor T1;
    T1.xx(1.0);
    T1.xy(2.0);
    T1.xz(3.0);
    T1.yx(4.0);
    T1.yy(5.0);
    T1.yz(6.0);
    T1.zx(7.0);
    T1.zy(8.0);
    T1.zz(9.0);
    SPHERAL_ASSERT_EQ(T1.xx(), 1.0);
    SPHERAL_ASSERT_EQ(T1.xy(), 2.0);
    SPHERAL_ASSERT_EQ(T1.xz(), 3.0);
    SPHERAL_ASSERT_EQ(T1.yx(), 4.0);
    SPHERAL_ASSERT_EQ(T1.yy(), 5.0);
    SPHERAL_ASSERT_EQ(T1.yz(), 6.0);
    SPHERAL_ASSERT_EQ(T1.zx(), 7.0);
    SPHERAL_ASSERT_EQ(T1.zy(), 8.0);
    SPHERAL_ASSERT_EQ(T1.zz(), 9.0);
  EXEC_IN_SPACE_END()
}

GPU_TYPED_TEST_P(GeomTensorTypedTest, GetSetRowsColumns) {
  using WORK_EXEC_POLICY = TypeParam;
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    Tensor T1(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
    Vector row0 = T1.getRow(0);
    SPHERAL_ASSERT_EQ(row0.x(), 1.0);
    SPHERAL_ASSERT_EQ(row0.y(), 2.0);
    SPHERAL_ASSERT_EQ(row0.z(), 3.0);

    Vector col1 = T1.getColumn(1);
    SPHERAL_ASSERT_EQ(col1.x(), 2.0);
    SPHERAL_ASSERT_EQ(col1.y(), 5.0);
    SPHERAL_ASSERT_EQ(col1.z(), 8.0);

    Vector newRow(10.0, 11.0, 12.0);
    T1.setRow(0, newRow);
    SPHERAL_ASSERT_EQ(T1.xx(), 10.0);
    SPHERAL_ASSERT_EQ(T1.xy(), 11.0);
    SPHERAL_ASSERT_EQ(T1.xz(), 12.0);

    Vector newCol(13.0, 14.0, 15.0);
    T1.setColumn(1, newCol);
    SPHERAL_ASSERT_EQ(T1.xy(), 13.0);
    SPHERAL_ASSERT_EQ(T1.yy(), 14.0);
    SPHERAL_ASSERT_EQ(T1.zy(), 15.0);
  EXEC_IN_SPACE_END()
}

GPU_TYPED_TEST_P(GeomTensorTypedTest, ZeroIdentity) {
  using WORK_EXEC_POLICY = TypeParam;
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    Tensor T1(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
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

GPU_TYPED_TEST_P(GeomTensorTypedTest, UnaryMinus) {
  using WORK_EXEC_POLICY = TypeParam;
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    Tensor T1(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
    Tensor T2 = -T1;
    SPHERAL_ASSERT_EQ(T2.xx(), -1.0);
    SPHERAL_ASSERT_EQ(T2.xy(), -2.0);
    SPHERAL_ASSERT_EQ(T2.xz(), -3.0);
    SPHERAL_ASSERT_EQ(T2.yx(), -4.0);
    SPHERAL_ASSERT_EQ(T2.yy(), -5.0);
    SPHERAL_ASSERT_EQ(T2.yz(), -6.0);
    SPHERAL_ASSERT_EQ(T2.zx(), -7.0);
    SPHERAL_ASSERT_EQ(T2.zy(), -8.0);
    SPHERAL_ASSERT_EQ(T2.zz(), -9.0);
  EXEC_IN_SPACE_END()
}

GPU_TYPED_TEST_P(GeomTensorTypedTest, AddSub) {
  using WORK_EXEC_POLICY = TypeParam;
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    Tensor T1(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
    Tensor T2(9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0);

    Tensor T_add = T1 + T2;
    SPHERAL_ASSERT_EQ(T_add.xx(), 10.0);
    SPHERAL_ASSERT_EQ(T_add.xy(), 10.0);
    SPHERAL_ASSERT_EQ(T_add.xz(), 10.0);
    SPHERAL_ASSERT_EQ(T_add.yx(), 10.0);
    SPHERAL_ASSERT_EQ(T_add.yy(), 10.0);
    SPHERAL_ASSERT_EQ(T_add.yz(), 10.0);
    SPHERAL_ASSERT_EQ(T_add.zx(), 10.0);
    SPHERAL_ASSERT_EQ(T_add.zy(), 10.0);
    SPHERAL_ASSERT_EQ(T_add.zz(), 10.0);

    Tensor T_sub = T1 - T2;
    SPHERAL_ASSERT_EQ(T_sub.xx(), -8.0);
    SPHERAL_ASSERT_EQ(T_sub.xy(), -6.0);
    SPHERAL_ASSERT_EQ(T_sub.xz(), -4.0);
    SPHERAL_ASSERT_EQ(T_sub.yx(), -2.0);
    SPHERAL_ASSERT_EQ(T_sub.yy(), 0.0);
    SPHERAL_ASSERT_EQ(T_sub.yz(), 2.0);
    SPHERAL_ASSERT_EQ(T_sub.zx(), 4.0);
    SPHERAL_ASSERT_EQ(T_sub.zy(), 6.0);
    SPHERAL_ASSERT_EQ(T_sub.zz(), 8.0);

    SymTensor ST(1.0, 2.0, 3.0,
                 2.0, 5.0, 6.0,
                 3.0, 6.0, 9.0);
    T_add = T1 + ST;
    SPHERAL_ASSERT_EQ(T_add.xx(), 2.0);
    SPHERAL_ASSERT_EQ(T_add.xy(), 4.0);
    SPHERAL_ASSERT_EQ(T_add.xz(), 6.0);
    SPHERAL_ASSERT_EQ(T_add.yx(), 6.0);
    SPHERAL_ASSERT_EQ(T_add.yy(), 10.0);
    SPHERAL_ASSERT_EQ(T_add.yz(), 12.0);
    SPHERAL_ASSERT_EQ(T_add.zx(), 10.0);
    SPHERAL_ASSERT_EQ(T_add.zy(), 14.0);
    SPHERAL_ASSERT_EQ(T_add.zz(), 18.0);

    T_sub = T1 - ST;
    SPHERAL_ASSERT_EQ(T_sub.xx(), 0.0);
    SPHERAL_ASSERT_EQ(T_sub.xy(), 0.0);
    SPHERAL_ASSERT_EQ(T_sub.xz(), 0.0);
    SPHERAL_ASSERT_EQ(T_sub.yx(), 2.0);
    SPHERAL_ASSERT_EQ(T_sub.yy(), 0.0);
    SPHERAL_ASSERT_EQ(T_sub.yz(), 0.0);
    SPHERAL_ASSERT_EQ(T_sub.zx(), 4.0);
    SPHERAL_ASSERT_EQ(T_sub.zy(), 2.0);
    SPHERAL_ASSERT_EQ(T_sub.zz(), 0.0);
  EXEC_IN_SPACE_END()
}

GPU_TYPED_TEST_P(GeomTensorTypedTest, ScalarMulDiv) {
  using WORK_EXEC_POLICY = TypeParam;
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    Tensor T1(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
    double s = 2.0;

    Tensor T_mul = T1 * s;
    SPHERAL_ASSERT_EQ(T_mul.xx(), 2.0);
    SPHERAL_ASSERT_EQ(T_mul.xy(), 4.0);
    SPHERAL_ASSERT_EQ(T_mul.xz(), 6.0);
    SPHERAL_ASSERT_EQ(T_mul.yx(), 8.0);
    SPHERAL_ASSERT_EQ(T_mul.yy(), 10.0);
    SPHERAL_ASSERT_EQ(T_mul.yz(), 12.0);
    SPHERAL_ASSERT_EQ(T_mul.zx(), 14.0);
    SPHERAL_ASSERT_EQ(T_mul.zy(), 16.0);
    SPHERAL_ASSERT_EQ(T_mul.zz(), 18.0);

    Tensor T_div = T1 / s;
    SPHERAL_ASSERT_EQ(T_div.xx(), 0.5);
    SPHERAL_ASSERT_EQ(T_div.xy(), 1.0);
    SPHERAL_ASSERT_EQ(T_div.xz(), 1.5);
    SPHERAL_ASSERT_EQ(T_div.yx(), 2.0);
    SPHERAL_ASSERT_EQ(T_div.yy(), 2.5);
    SPHERAL_ASSERT_EQ(T_div.yz(), 3.0);
    SPHERAL_ASSERT_EQ(T_div.zx(), 3.5);
    SPHERAL_ASSERT_EQ(T_div.zy(), 4.0);
    SPHERAL_ASSERT_EQ(T_div.zz(), 4.5);
  EXEC_IN_SPACE_END()
}

GPU_TYPED_TEST_P(GeomTensorTypedTest, InPlaceAddSub) {
  using WORK_EXEC_POLICY = TypeParam;
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    Tensor T1(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
    Tensor T2(9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0);

    T1 += T2;
    SPHERAL_ASSERT_EQ(T1.xx(), 10.0);
    SPHERAL_ASSERT_EQ(T1.xy(), 10.0);
    SPHERAL_ASSERT_EQ(T1.xz(), 10.0);
    SPHERAL_ASSERT_EQ(T1.yx(), 10.0);
    SPHERAL_ASSERT_EQ(T1.yy(), 10.0);
    SPHERAL_ASSERT_EQ(T1.yz(), 10.0);
    SPHERAL_ASSERT_EQ(T1.zx(), 10.0);
    SPHERAL_ASSERT_EQ(T1.zy(), 10.0);
    SPHERAL_ASSERT_EQ(T1.zz(), 10.0);

    T1 -= T2;
    SPHERAL_ASSERT_EQ(T1.xx(), 1.0);
    SPHERAL_ASSERT_EQ(T1.xy(), 2.0);
    SPHERAL_ASSERT_EQ(T1.xz(), 3.0);
    SPHERAL_ASSERT_EQ(T1.yx(), 4.0);
    SPHERAL_ASSERT_EQ(T1.yy(), 5.0);
    SPHERAL_ASSERT_EQ(T1.yz(), 6.0);
    SPHERAL_ASSERT_EQ(T1.zx(), 7.0);
    SPHERAL_ASSERT_EQ(T1.zy(), 8.0);
    SPHERAL_ASSERT_EQ(T1.zz(), 9.0);
  EXEC_IN_SPACE_END()
}

GPU_TYPED_TEST_P(GeomTensorTypedTest, InPlaceScalarMulDiv) {
  using WORK_EXEC_POLICY = TypeParam;
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    Tensor T1(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
    double s = 2.0;

    T1 *= s;
    SPHERAL_ASSERT_EQ(T1.xx(), 2.0);
    SPHERAL_ASSERT_EQ(T1.xy(), 4.0);
    SPHERAL_ASSERT_EQ(T1.xz(), 6.0);
    SPHERAL_ASSERT_EQ(T1.yx(), 8.0);
    SPHERAL_ASSERT_EQ(T1.yy(), 10.0);
    SPHERAL_ASSERT_EQ(T1.yz(), 12.0);
    SPHERAL_ASSERT_EQ(T1.zx(), 14.0);
    SPHERAL_ASSERT_EQ(T1.zy(), 16.0);
    SPHERAL_ASSERT_EQ(T1.zz(), 18.0);

    T1 /= s;
    SPHERAL_ASSERT_EQ(T1.xx(), 1.0);
    SPHERAL_ASSERT_EQ(T1.xy(), 2.0);
    SPHERAL_ASSERT_EQ(T1.xz(), 3.0);
    SPHERAL_ASSERT_EQ(T1.yx(), 4.0);
    SPHERAL_ASSERT_EQ(T1.yy(), 5.0);
    SPHERAL_ASSERT_EQ(T1.yz(), 6.0);
    SPHERAL_ASSERT_EQ(T1.zx(), 7.0);
    SPHERAL_ASSERT_EQ(T1.zy(), 8.0);
    SPHERAL_ASSERT_EQ(T1.zz(), 9.0);
  EXEC_IN_SPACE_END()
}

GPU_TYPED_TEST_P(GeomTensorTypedTest, Comparison) {
  using WORK_EXEC_POLICY = TypeParam;
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    Tensor T1(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
    Tensor T2(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
    Tensor T3(8.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0);

    SPHERAL_ASSERT_TRUE(T1 == T2);
    SPHERAL_ASSERT_TRUE(T1 != T3);
    SPHERAL_ASSERT_TRUE(T1 < T3);
    SPHERAL_ASSERT_TRUE(T3 > T1);
    SPHERAL_ASSERT_TRUE(T1 <= T2);
    SPHERAL_ASSERT_TRUE(T1 >= T2);
  EXEC_IN_SPACE_END()
}

GPU_TYPED_TEST_P(GeomTensorTypedTest, Transpose) {
  using WORK_EXEC_POLICY = TypeParam;
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    Tensor T1(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
    Tensor T_trans = T1.Transpose();
    SPHERAL_ASSERT_EQ(T_trans.xx(), 1.0);
    SPHERAL_ASSERT_EQ(T_trans.xy(), 4.0);
    SPHERAL_ASSERT_EQ(T_trans.xz(), 7.0);
    SPHERAL_ASSERT_EQ(T_trans.yx(), 2.0);
    SPHERAL_ASSERT_EQ(T_trans.yy(), 5.0);
    SPHERAL_ASSERT_EQ(T_trans.yz(), 8.0);
    SPHERAL_ASSERT_EQ(T_trans.zx(), 3.0);
    SPHERAL_ASSERT_EQ(T_trans.zy(), 6.0);
    SPHERAL_ASSERT_EQ(T_trans.zz(), 9.0);
  EXEC_IN_SPACE_END()
}

GPU_TYPED_TEST_P(GeomTensorTypedTest, TraceDeterminant) {
  using WORK_EXEC_POLICY = TypeParam;
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    Tensor T1(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
    SPHERAL_ASSERT_EQ(T1.Trace(), 15.0);
    SPHERAL_ASSERT_EQ(T1.Determinant(), 0.0);
  EXEC_IN_SPACE_END()
}

GPU_TYPED_TEST_P(GeomTensorTypedTest, DotProduct) {
  using WORK_EXEC_POLICY = TypeParam;
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    Tensor T1(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
    Tensor T2(9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0);
    Tensor T_dot = T1.dot(T2);
    // Expected values calculated manually
    SPHERAL_ASSERT_EQ(T_dot.xx(), 30.0);
    SPHERAL_ASSERT_EQ(T_dot.xy(), 24.0);
    SPHERAL_ASSERT_EQ(T_dot.xz(), 18.0);
    SPHERAL_ASSERT_EQ(T_dot.yx(), 84.0);
    SPHERAL_ASSERT_EQ(T_dot.yy(), 69.0);
    SPHERAL_ASSERT_EQ(T_dot.yz(), 54.0);
    SPHERAL_ASSERT_EQ(T_dot.zx(), 138.0);
    SPHERAL_ASSERT_EQ(T_dot.zy(), 114.0);
    SPHERAL_ASSERT_EQ(T_dot.zz(), 90.0);

    Vector V1(1.0, 2.0, 3.0);
    Vector V_dot = T1.dot(V1);
    SPHERAL_ASSERT_EQ(V_dot.x(), 14.0);
    SPHERAL_ASSERT_EQ(V_dot.y(), 32.0);
    SPHERAL_ASSERT_EQ(V_dot.z(), 50.0);
  EXEC_IN_SPACE_END()
}

// TODO: Check this w/ Mike Owen
GPU_TYPED_TEST_P(GeomTensorTypedTest, DoubleDotProduct) {
  //using WORK_EXEC_POLICY = TypeParam;
  //EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
  //  Tensor T1(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
  //  Tensor T2(9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0);
  //  double ddot = T1.doubledot(T2);
  //  // Expected value calculated manually
  //  SPHERAL_ASSERT_EQ(ddot, 1.0*9.0 + 2.0*8.0 + 3.0*7.0 +
  //                          4.0*6.0 + 5.0*5.0 + 6.0*4.0 +
  //                          7.0*3.0 + 8.0*2.0 + 9.0*1.0);
  //EXEC_IN_SPACE_END()
}

GPU_TYPED_TEST_P(GeomTensorTypedTest, Square) {
  using WORK_EXEC_POLICY = TypeParam;
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    Tensor T1(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
    Tensor T_sq = T1.square();
    // Expected values calculated manually
    SPHERAL_ASSERT_EQ(T_sq.xx(), 30.0);
    SPHERAL_ASSERT_EQ(T_sq.xy(), 36.0);
    SPHERAL_ASSERT_EQ(T_sq.xz(), 42.0);
    SPHERAL_ASSERT_EQ(T_sq.yx(), 66.0);
    SPHERAL_ASSERT_EQ(T_sq.yy(), 81.0);
    SPHERAL_ASSERT_EQ(T_sq.yz(), 96.0);
    SPHERAL_ASSERT_EQ(T_sq.zx(), 102.0);
    SPHERAL_ASSERT_EQ(T_sq.zy(), 126.0);
    SPHERAL_ASSERT_EQ(T_sq.zz(), 150.0);
  EXEC_IN_SPACE_END()
}

GPU_TYPED_TEST_P(GeomTensorTypedTest, SquareElements) {
  using WORK_EXEC_POLICY = TypeParam;
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    Tensor T1(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
    Tensor T_sq_elem = T1.squareElements();
    SPHERAL_ASSERT_EQ(T_sq_elem.xx(), 1.0);
    SPHERAL_ASSERT_EQ(T_sq_elem.xy(), 4.0);
    SPHERAL_ASSERT_EQ(T_sq_elem.xz(), 9.0);
    SPHERAL_ASSERT_EQ(T_sq_elem.yx(), 16.0);
    SPHERAL_ASSERT_EQ(T_sq_elem.yy(), 25.0);
    SPHERAL_ASSERT_EQ(T_sq_elem.yz(), 36.0);
    SPHERAL_ASSERT_EQ(T_sq_elem.zx(), 49.0);
    SPHERAL_ASSERT_EQ(T_sq_elem.zy(), 64.0);
    SPHERAL_ASSERT_EQ(T_sq_elem.zz(), 81.0);
  EXEC_IN_SPACE_END()
}

GPU_TYPED_TEST_P(GeomTensorTypedTest, EigenValues) {
  using WORK_EXEC_POLICY = TypeParam;
  const Vector ref_vals(-3., -1., 1.);
  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
    Tensor T1(0.0, 0.0, 1.0, 0.0, -3.0, 0.0, 1.0, 0.0, 0.0);
    Vector vals = T1.eigenValues();
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

REGISTER_TYPED_TEST_SUITE_P(GeomTensorTypedTest, Ctor, Assignment, Accessors,
                            Setters, GetSetRowsColumns, ZeroIdentity, UnaryMinus, AddSub,
                            ScalarMulDiv, InPlaceAddSub, InPlaceScalarMulDiv, Comparison,
                            Transpose, TraceDeterminant, DotProduct, DoubleDotProduct,
                            Square, SquareElements, EigenValues);

INSTANTIATE_TYPED_TEST_SUITE_P(GeomTensor, GeomTensorTypedTest,
                               typename Spheral::Test<EXEC_TYPES>::Types, );
