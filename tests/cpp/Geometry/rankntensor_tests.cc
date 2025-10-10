#include "rankntensor_test_types.hh"
#include "test-utilities.hh"

#include "Geometry/RankNTensor.hh"

/*
 * Test fixture for Rank-N Tensor tests.
 */
class RankNTensorTest : public ::testing::Test {};

/*
 * Typed test suite for Rank-N Tensors.
 */
template <typename T>
class RankNTensorTypedTest : public RankNTensorTest {};
TYPED_TEST_SUITE_P(RankNTensorTypedTest);

/*
 * Test case for the constructor of a Rank-N tensor.
 * Verifies that the constructor correctly initializes the tensor elements.
 * A tensor `t` is created, and the test asserts that all elements are
 * initialized to 0.0.
 */
GPU_TYPED_TEST_P(RankNTensorTypedTest, Ctor) {
  using WORK_EXEC_POLICY = typename camp::at<TypeParam, camp::num<0>>::type;
  using TensorType = typename camp::at<TypeParam, camp::num<1>>::type;

  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
  TensorType t;
  for (typename TensorType::size_type i = 0; i < TensorType::numElements; ++i) {
    SPHERAL_ASSERT_EQ(t[i], 0.0);
  }
  EXEC_IN_SPACE_END()
}

/*
 * Test case for the assignment operators of a Rank-N tensor.
 * Verifies that both scalar and tensor assignment work correctly.
 * It first assigns a scalar value to all elements of a tensor `t` and asserts
 * the result. It then assigns another tensor `u` to `t` and asserts that `t`'s
 * elements match `u`'s.
 */
GPU_TYPED_TEST_P(RankNTensorTypedTest, Assignment) {
  using WORK_EXEC_POLICY = typename camp::at<TypeParam, camp::num<0>>::type;
  using TensorType = typename camp::at<TypeParam, camp::num<1>>::type;

  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
  TensorType t;
  t = 5.0;
  for (typename TensorType::size_type i = 0; i < TensorType::numElements; ++i) {
    SPHERAL_ASSERT_EQ(t[i], 5.0);
  }

  TensorType u(1.0);
  t = u;
  for (typename TensorType::size_type i = 0; i < TensorType::numElements; ++i) {
    SPHERAL_ASSERT_EQ(t[i], 1.0);
  }
  EXEC_IN_SPACE_END()
}

/*
 * Test case for the element accessors of a Rank-N tensor.
 * Verifies that the square bracket operator correctly accesses and modifies
 * tensor elements. It asserts that `t[i]` returns the correct values, and then
 * modifies elements using this accessor and asserts the new values.
 */
GPU_TYPED_TEST_P(RankNTensorTypedTest, Accessors) {
  using WORK_EXEC_POLICY = typename camp::at<TypeParam, camp::num<0>>::type;
  using TensorType = typename camp::at<TypeParam, camp::num<1>>::type;

  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
  TensorType t;
  for (typename TensorType::size_type i = 0; i < TensorType::numElements; ++i) {
    t[i] = i;
  }

  for (typename TensorType::size_type i = 0; i < TensorType::numElements; ++i) {
    SPHERAL_ASSERT_EQ(t[i], i);
  }
  EXEC_IN_SPACE_END()
}

/*
 * Test case for the Zero() method of a Rank-N tensor.
 * Verifies that the Zero() method correctly sets all elements of the tensor to
 * 0.0. A tensor `t` is initialized with non-zero values, Zero() is called,
 * and the test asserts that all elements are now 0.0.
 */
GPU_TYPED_TEST_P(RankNTensorTypedTest, Zero) {
  using WORK_EXEC_POLICY = typename camp::at<TypeParam, camp::num<0>>::type;
  using TensorType = typename camp::at<TypeParam, camp::num<1>>::type;

  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
  TensorType t(1.0);
  t.Zero();
  for (typename TensorType::size_type i = 0; i < TensorType::numElements; ++i) {
    SPHERAL_ASSERT_EQ(t[i], 0.0);
  }
  EXEC_IN_SPACE_END()
}

/*
 * Test case for the unary minus operator of a Rank-N tensor.
 * Verifies that the unary minus operator correctly negates all elements of the
 * tensor. A tensor `t2` is created by negating `t1`, and the test asserts that
 * each element of `t2` is the negated value of the corresponding element in
 * `t1`.
 */
GPU_TYPED_TEST_P(RankNTensorTypedTest, UnaryMinus) {
  using WORK_EXEC_POLICY = typename camp::at<TypeParam, camp::num<0>>::type;
  using TensorType = typename camp::at<TypeParam, camp::num<1>>::type;

  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
  TensorType t1(1.0);
  TensorType t2 = -t1;
  for (typename TensorType::size_type i = 0; i < TensorType::numElements; ++i) {
    SPHERAL_ASSERT_EQ(t2[i], -1.0);
  }
  EXEC_IN_SPACE_END()
}

/*
 * Test case for Rank-N tensor addition and subtraction.
 * Verifies that the addition and subtraction operators produce the correct
 * results. It asserts that the elements of `t1 + t2` and `t2 - t1` are the
 * sums and differences of the elements of `t1` and `t2`.
 */
GPU_TYPED_TEST_P(RankNTensorTypedTest, AddSub) {
  using WORK_EXEC_POLICY = typename camp::at<TypeParam, camp::num<0>>::type;
  using TensorType = typename camp::at<TypeParam, camp::num<1>>::type;

  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
  TensorType t1(1.0);
  TensorType t2(2.0);

  TensorType t_add = t1 + t2;
  for (typename TensorType::size_type i = 0; i < TensorType::numElements; ++i) {
    SPHERAL_ASSERT_EQ(t_add[i], 3.0);
  }

  TensorType t_sub = t2 - t1;
  for (typename TensorType::size_type i = 0; i < TensorType::numElements; ++i) {
    SPHERAL_ASSERT_EQ(t_sub[i], 1.0);
  }
  EXEC_IN_SPACE_END()
}

/*
 * Test case for Rank-N tensor scalar multiplication and division.
 * Verifies that scalar multiplication and division operators produce the
 * correct results. It asserts that the elements of `t * s` and `t / s` are the
 * product and quotient of the elements of `t` and a scalar `s`.
 */
GPU_TYPED_TEST_P(RankNTensorTypedTest, ScalarMulDiv) {
  using WORK_EXEC_POLICY = typename camp::at<TypeParam, camp::num<0>>::type;
  using TensorType = typename camp::at<TypeParam, camp::num<1>>::type;

  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
  TensorType t(1.0);
  double s = 2.0;

  TensorType t_mul = t * s;
  for (typename TensorType::size_type i = 0; i < TensorType::numElements; ++i) {
    SPHERAL_ASSERT_EQ(t_mul[i], 2.0);
  }

  TensorType t_div = t / s;
  for (typename TensorType::size_type i = 0; i < TensorType::numElements; ++i) {
    SPHERAL_ASSERT_EQ(t_div[i], 0.5);
  }
  EXEC_IN_SPACE_END()
}

/*
 * Test case for Rank-N tensor in-place addition and subtraction.
 * Verifies that the `+=` and `-=` operators produce the correct results. It
 * asserts that `t1 += t2` and `t1 -= t2` correctly modify the elements of `t1`.
 */
GPU_TYPED_TEST_P(RankNTensorTypedTest, InPlaceAddSub) {
  using WORK_EXEC_POLICY = typename camp::at<TypeParam, camp::num<0>>::type;
  using TensorType = typename camp::at<TypeParam, camp::num<1>>::type;

  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
  TensorType t1(1.0);
  TensorType t2(2.0);

  t1 += t2;
  for (typename TensorType::size_type i = 0; i < TensorType::numElements; ++i) {
    SPHERAL_ASSERT_EQ(t1[i], 3.0);
  }

  t1 -= t2;
  for (typename TensorType::size_type i = 0; i < TensorType::numElements; ++i) {
    SPHERAL_ASSERT_EQ(t1[i], 1.0);
  }
  EXEC_IN_SPACE_END()
}

/*
 * Test case for Rank-N tensor in-place scalar multiplication and division.
 * Verifies that the `*=` and `/=` operators produce the correct results. It
 * asserts that `t *= s` and `t /= s` correctly modify the elements of `t`.
 */
GPU_TYPED_TEST_P(RankNTensorTypedTest, InPlaceScalarMulDiv) {
  using WORK_EXEC_POLICY = typename camp::at<TypeParam, camp::num<0>>::type;
  using TensorType = typename camp::at<TypeParam, camp::num<1>>::type;

  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
  TensorType t(1.0);
  double s = 2.0;

  t *= s;
  for (typename TensorType::size_type i = 0; i < TensorType::numElements; ++i) {
    SPHERAL_ASSERT_EQ(t[i], 2.0);
  }

  t /= s;
  for (typename TensorType::size_type i = 0; i < TensorType::numElements; ++i) {
    SPHERAL_ASSERT_EQ(t[i], 1.0);
  }
  EXEC_IN_SPACE_END()
}

/*
 * Test case for Rank-N tensor comparison operators.
 * Verifies that the `==`, `!=`, `<`, `>`, `<=`, and `>=` operators produce the
 * correct boolean results. It asserts the outcomes of comparing several
 * tensors.
 */
GPU_TYPED_TEST_P(RankNTensorTypedTest, Comparison) {
  using WORK_EXEC_POLICY = typename camp::at<TypeParam, camp::num<0>>::type;
  using TensorType = typename camp::at<TypeParam, camp::num<1>>::type;

  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
  TensorType t1(1.0);
  TensorType t2(1.0);
  TensorType t3(2.0);
  TensorType t4(0.0);

  SPHERAL_ASSERT_TRUE(t1 == t2);
  SPHERAL_ASSERT_TRUE(t1 != t3);
  SPHERAL_ASSERT_TRUE(t1 < t3);
  SPHERAL_ASSERT_TRUE(t3 > t1);
  SPHERAL_ASSERT_TRUE(t1 <= t2);
  SPHERAL_ASSERT_TRUE(t1 >= t2);
  SPHERAL_ASSERT_TRUE(t4 < t1);
  EXEC_IN_SPACE_END()
}

/*
 * Test case for the Rank-N tensor double dot product.
 * Verifies that the `doubledot()` method calculates the correct double dot
 * product of two tensors. It asserts that `t1.doubledot(t2)` equals the sum of
 * the products of their corresponding elements.
 */
GPU_TYPED_TEST_P(RankNTensorTypedTest, DoubleDot) {
  using WORK_EXEC_POLICY = typename camp::at<TypeParam, camp::num<0>>::type;
  using TensorType = typename camp::at<TypeParam, camp::num<1>>::type;

  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
  TensorType t1(2.0);
  TensorType t2(3.0);
  double dot_product = t1.doubledot(t2);
  SPHERAL_ASSERT_EQ(dot_product, 2.0 * 3.0 * TensorType::numElements);
  EXEC_IN_SPACE_END()
}

/*
 * Test case for the Rank-N tensor element-wise statistics.
 * Verifies that `maxAbsElement` calculates the correct value. It asserts the
 * result for a tensor with positive and negative elements.
 */
GPU_TYPED_TEST_P(RankNTensorTypedTest, ElementStats) {
  using WORK_EXEC_POLICY = typename camp::at<TypeParam, camp::num<0>>::type;
  using TensorType = typename camp::at<TypeParam, camp::num<1>>::type;

  EXEC_IN_SPACE_BEGIN(WORK_EXEC_POLICY)
  TensorType t;
  for (size_t i = 0; i < TensorType::numElements; ++i) {
    t[i] = (i % 2 == 0) ? i : -(double)i;
  }

  SPHERAL_ASSERT_EQ(t.maxAbsElement(), TensorType::numElements - 1);
  EXEC_IN_SPACE_END()
}

/*
 * Register all typed tests for the Rank-N tensor test suite.
 */
REGISTER_TYPED_TEST_SUITE_P(RankNTensorTypedTest, Ctor, Assignment, Accessors,
                            Zero, UnaryMinus, AddSub, ScalarMulDiv,
                            InPlaceAddSub, InPlaceScalarMulDiv, Comparison,
                            DoubleDot, ElementStats);


INSTANTIATE_TYPED_TEST_SUITE_P(RankNTensor, RankNTensorTypedTest,
                               RANK_N_TENSOR_TEST_TYPES, );
