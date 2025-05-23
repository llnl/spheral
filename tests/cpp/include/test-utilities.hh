#ifndef SPHERAL_TEST_UTIILITIES_HH
#define SPHERAL_TEST_UTIILITIES_HH

#include "RAJA/RAJA.hpp"
#include "assert.h"
#include "config.hh"
#include "gtest/gtest.h"

using TRS_UINT = RAJA::TypedRangeSegment<unsigned>;
using LOOP_EXEC_POLICY = RAJA::seq_exec;

#define EXEC_IN_SPACE_BEGIN(POL)                                               \
  RAJA::forall<POL>(TRS_UINT(0,1), [=] SPHERAL_HOST_DEVICE (int) {

#define EXEC_IN_SPACE_END()                                                    \
  });

#if defined(SPHERAL_ENABLE_HIP) && defined(__HIPCC__)
#define SPHERAL_GPU_ACTIVE
#endif // SPHERAL_ENABLE_CUDA && __CUDACC__

#if !defined(SPHERAL_GPU_ACTIVE)

#define SPHERAL_ASSERT_EQ(LHS, RHS) ASSERT_EQ(LHS, RHS);
#define SPHERAL_ASSERT_NE(LHS, RHS) ASSERT_NE(LHS, RHS);
#define SPHERAL_ASSERT_TRUE(VAL) ASSERT_TRUE(VAL);
#define SPHERAL_ASSERT_FALSE(VAL) ASSERT_FALSE(VAL);

#else

#define SPHERAL_ASSERT_EQ(LHS, RHS)                                            \
  if (LHS != RHS) {                                                            \
    printf("ERROR @ cuda_assert\n");                                           \
    assert(0);                                                                 \
  }

#define SPHERAL_ASSERT_NE(LHS, RHS)                                            \
  if (LHS == RHS) {                                                            \
    printf("ERROR @ cuda_assert\n");                                           \
    assert(0);                                                                 \
  }

#define SPHERALSPHERAL_ASSERT_TRUE(VAL)                                        \
  if (!result) {                                                               \
    printf("ERROR @ cuda_assert\n");                                           \
    assert(0);                                                                 \
  }

#define SPHERALSPHERAL_ASSERT_FALSE(VAL)                                       \
  if (result) {                                                                \
    printf("ERROR @ cuda_assert\n");                                           \
    assert(0);                                                                 \
  }

#endif

// Macro used throughout LLNLProjects to get around calling
// HOST_DEVICE lambdas from within the "Testing" class directly
#define GPU_TYPED_TEST(X, Y)                                                   \
  template <typename TypeParam> static void gpu_test_##X##Y();                 \
  TYPED_TEST(X, Y) { gpu_test_##X##Y<TypeParam>(); }                           \
  template <typename TypeParam> static void gpu_test_##X##Y()

#define GPU_TYPED_TEST_P(X, Y)                                                 \
  template <typename TypeParam, typename TestFixture>                          \
  static void gpu_test_##X##Y(TestFixture *gpu_this);                          \
  TYPED_TEST_P(X, Y) { gpu_test_##X##Y<TypeParam>(this); }                     \
  template <typename TypeParam, typename TestFixture>                          \
  static void gpu_test_##X##Y(TestFixture *gpu_this)

#endif // SPHERAL_TEST_UTIILITIES_HH
