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
    Spheral::QIView qi1_view = qi1.view();
    Spheral::QIView qi2_view = qi2.view();
    Spheral::QIView qiref_view = qiref;
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

// Test initialize using a func
GPU_TYPED_TEST_P(QuadraticInterpolatorTypedTest, FuncCtorTest) {
  const size_t NV = gpu_this->NV;
  const double xmin = gpu_this->xmin;
  const double xmax = gpu_this->xmax;
  QI qih(xmin, xmax, NV, gpu_this->func);
  {
    size_t N = qih.size();
    Spheral::QIView qi = qih.view();
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
  Spheral::QIView qi = qih.view();
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

REGISTER_TYPED_TEST_SUITE_P(QuadraticInterpolatorTypedTest, CopyAssign, FuncCtorTest,
                            VecCtorTest);

INSTANTIATE_TYPED_TEST_SUITE_P(QuadraticInterpolator, QuadraticInterpolatorTypedTest,
                               typename Spheral::Test<EXEC_TYPES>::Types, );
