// Debug log printing can be quickly enabled for this unit test by uncommenting the
// definition below even if Spheral was not configured w/ SPHERAL_ENABLE_LOGGER=On.
// #define SPHERAL_ENABLE_LOGGER

#include "chai/ExecutionSpaces.hpp"
#include "chai/Types.hpp"
#include "test-basic-exec-policies.hh"
#include "test-utilities.hh"

#include "Utilities/QuadraticInterpolator.hh"
#include <Utilities/Logger.hh>
#include <functional>

using QI = Spheral::QuadraticInterpolator;
using QIView = Spheral::QuadraticInterpolatorView;

class QuadraticInterpolatorTest : public ::testing::Test {
public:
  const size_t NV = 41;
  const double xmin = 10.;
  const double xmax = 100.;
  SPHERAL_HOST_DEVICE static double func(const double x) {
    return 1. + x*(2. + 3.*x);
  }
  std::vector<double> makeVec(size_t N) {
    double xstep = (xmax - xmin)/((double)N - 1.);
    std::vector<double> yvals(N);
    for (size_t i = 0; i < N; ++i) {
      double x = xmin + xstep*(double)i;
      yvals[i] = func(x);
    }
    return yvals;
  }
};

// Setting up G Test for QI
TYPED_TEST_SUITE_P(QuadraticInterpolatorTypedTest);
template <typename T> class QuadraticInterpolatorTypedTest : public QuadraticInterpolatorTest {};

// Test copy and assignment constructors
GPU_TYPED_TEST_P(QuadraticInterpolatorTypedTest, CopyAssign) {
  const size_t NV = gpu_this->NV;
  QI qiref(gpu_this->xmin, gpu_this->xmax, NV, gpu_this->func);
  {
    QI qi1(qiref);
    SPHERAL_ASSERT_EQ(qi1.size(), qiref.size());
    QI qi2 = qiref;
    SPHERAL_ASSERT_EQ(qi2.size(), qiref.size());
    QIView qi1_view = qi1.view();
    QIView qi2_view = qi2.view();
    QIView qiref_view = qiref;
    // Ensure the underlying data pointer is different from the initial QI
    EXEC_IN_SPACE_BEGIN(TypeParam)
      SPHERAL_ASSERT_NE(qi1_view.data(), qiref_view.data());
      SPHERAL_ASSERT_NE(qi2_view.data(), qiref_view.data());
    EXEC_IN_SPACE_END()    
    RAJA::forall<TypeParam>(TRS_UINT(0, NV),
      [=] (size_t i) {
        SPHERAL_ASSERT_EQ(qi1_view[i], qiref_view[i]);
        SPHERAL_ASSERT_EQ(qi2_view[i], qiref_view[i]);
      });
  }
}

// Assignment and reinitialization must release a previously migrated View
// before replacing the vector-backed coefficient storage.
GPU_TYPED_TEST_P(QuadraticInterpolatorTypedTest, AssignAndReinitializeAfterMigration) {
  const size_t NV = gpu_this->NV;
  QI source(gpu_this->xmin, gpu_this->xmax, NV, gpu_this->func);
  QI target(0.0, 1.0, 7u, gpu_this->func);

  {
    const auto oldView = target.view();
    EXEC_IN_SPACE_BEGIN(TypeParam)
      SPHERAL_ASSERT_EQ(oldView.size(), 9u);
    EXEC_IN_SPACE_END()
  }

  target = source;
  const auto assignedSize = source.size();
  {
    const auto sourceView = source.view();
    const auto targetView = target.view();
    EXEC_IN_SPACE_BEGIN(TypeParam)
      SPHERAL_ASSERT_NE(targetView.data(), sourceView.data());
    EXEC_IN_SPACE_END()
    RAJA::forall<TypeParam>(TRS_UINT(0, target.size()),
      [=] SPHERAL_HOST_DEVICE (size_t i) {
        SPHERAL_ASSERT_EQ(targetView[i], sourceView[i]);
      });
  }

  // Migrate the assigned coefficients, then replace them with a differently
  // sized table and verify the newly acquired snapshot.
  {
    const auto migratedView = target.view();
    EXEC_IN_SPACE_BEGIN(TypeParam)
      SPHERAL_ASSERT_EQ(migratedView.size(), assignedSize);
    EXEC_IN_SPACE_END()
  }
  const size_t NV2 = NV + 4u;       // Quadratic tables require odd sample counts.
  target.initialize(-2.0, 2.0, NV2, gpu_this->func);
  const auto targetView = target.view();
  const double dx = 4.0/(NV2 - 1u);
  RAJA::forall<TypeParam>(TRS_UINT(0, NV2),
    [=] SPHERAL_HOST_DEVICE (size_t i) {
      const double x = -2.0 + dx*i;
      SPHERAL_ASSERT_FLOAT_EQ(targetView(x), gpu_this->func(x));
    });
}

// Test initialize using a func
GPU_TYPED_TEST_P(QuadraticInterpolatorTypedTest, FuncCtorTest) {
  const size_t NV = gpu_this->NV;
  const double xmin = gpu_this->xmin;
  const double xmax = gpu_this->xmax;
  QI qih(xmin, xmax, NV, gpu_this->func);
  {
    size_t N = qih.size();
    QIView qi = qih.view();
    EXEC_IN_SPACE_BEGIN(TypeParam)
      SPHERAL_ASSERT_EQ(qi.size(), N);
    EXEC_IN_SPACE_END()
    const double xstep = (xmax - xmin)/((double)NV - 1.);
    RAJA::forall<TypeParam>(TRS_UINT(0, NV),
      [=] (size_t i) {
        double x = xmin + xstep*(double)i;
        double rval = gpu_this->func(x);
        double ival = qi(x);
        SPHERAL_ASSERT_FLOAT_EQ(rval, ival);
      });
  }
}

// Test initialize using a vector
GPU_TYPED_TEST_P(QuadraticInterpolatorTypedTest, VecCtorTest) {
  const size_t NV = gpu_this->NV;
  std::vector<double> yvals = gpu_this->makeVec(NV);
  const double xmin = gpu_this->xmin;
  const double xmax = gpu_this->xmax;
  QI qih(xmin, xmax, yvals);
  size_t N = qih.size();
  QIView qi = qih.view();
  EXEC_IN_SPACE_BEGIN(TypeParam)
    SPHERAL_ASSERT_EQ(qi.size(), N);
  EXEC_IN_SPACE_END()
  const double xstep = (xmax - xmin)/((double)NV - 1.);
  RAJA::forall<TypeParam>(TRS_UINT(0, NV),
    [=] (size_t i) {
      double x = xmin + xstep*(double)i;
      double rval = gpu_this->func(x);
      double ival = qi(x);
      SPHERAL_ASSERT_FLOAT_EQ(rval, ival);
    });
}

REGISTER_TYPED_TEST_SUITE_P(QuadraticInterpolatorTypedTest, CopyAssign,
                            AssignAndReinitializeAfterMigration, FuncCtorTest,
                            VecCtorTest);

INSTANTIATE_TYPED_TEST_SUITE_P(QuadraticInterpolator, QuadraticInterpolatorTypedTest,
                               typename Spheral::Test<EXEC_TYPES>::Types, );
