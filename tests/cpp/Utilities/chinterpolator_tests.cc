// Debug log printing can be quickly enabled for this unit test by uncommenting the
// definition below even if Spheral was not configured w/ SPHERAL_ENABLE_LOGGER=On.
// #define SPHERAL_ENABLE_LOGGER

#include "chai/ExecutionSpaces.hpp"
#include "chai/Types.hpp"
#include "test-basic-exec-policies.hh"
#include "test-utilities.hh"

#include "Utilities/CubicHermiteInterpolator.hh"
#include <Utilities/Logger.hh>
#include <functional>

using CHI = Spheral::CubicHermiteInterpolator;

class CubicHermiteInterpolatorTest : public ::testing::Test {
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

// Setting up G Test for CHI
TYPED_TEST_SUITE_P(CubicHermiteInterpolatorTypedTest);
template <typename T> class CubicHermiteInterpolatorTypedTest : public CubicHermiteInterpolatorTest {};

// Test copy and assignment constructors
GPU_TYPED_TEST_P(CubicHermiteInterpolatorTypedTest, CopyAssign) {
  const size_t NV = gpu_this->NV;
  CHI chiref(gpu_this->xmin, gpu_this->xmax, NV, gpu_this->func);
  {
    CHI chi1(chiref);
    SPHERAL_ASSERT_EQ(chi1.size(), chiref.size());
    CHI chi2 = chiref;
    SPHERAL_ASSERT_EQ(chi2.size(), chiref.size());
    Spheral::CHIView chi1_view = chi1.view();
    Spheral::CHIView chi2_view = chi2.view();
    Spheral::CHIView chiref_view = chiref;
    // Ensure the underlying data pointer is different from the initial CHI
    EXEC_IN_SPACE_BEGIN(TypeParam)
      SPHERAL_ASSERT_NE(chi1_view.data(), chiref_view.data());
      SPHERAL_ASSERT_NE(chi2_view.data(), chiref_view.data());
    EXEC_IN_SPACE_END()    
    RAJA::forall<TypeParam>(TRS_UINT(0, NV),
      [=] (size_t i) {
        SPHERAL_ASSERT_EQ(chi1_view[i], chiref_view[i]);
        SPHERAL_ASSERT_EQ(chi2_view[i], chiref_view[i]);
      });
  }
}

// Test initialize using a func
GPU_TYPED_TEST_P(CubicHermiteInterpolatorTypedTest, FuncCtorTest) {
  const size_t NV = gpu_this->NV;
  const double xmin = gpu_this->xmin;
  const double xmax = gpu_this->xmax;
  CHI chih(xmin, xmax, NV, gpu_this->func);
  {
    size_t N = chih.size();
    Spheral::CHIView chi = chih.view();
    EXEC_IN_SPACE_BEGIN(TypeParam)
      SPHERAL_ASSERT_EQ(chi.size(), N);
    EXEC_IN_SPACE_END()
    const double xstep = (xmax - xmin)/((double)NV - 1.);
    RAJA::forall<TypeParam>(TRS_UINT(0, NV),
      [=] (size_t i) {
        double x = xmin + xstep*(double)i;
        double rval = gpu_this->func(x);
        double ival = chi(x);
        SPHERAL_ASSERT_FLOAT_EQ(rval, ival);
      });
  }
}

// Test initialize using a vector
GPU_TYPED_TEST_P(CubicHermiteInterpolatorTypedTest, VecCtorTest) {
  const size_t NV = gpu_this->NV;
  std::vector<double> yvals = gpu_this->makeVec(NV);
  const double xmin = gpu_this->xmin;
  const double xmax = gpu_this->xmax;
  CHI chih(xmin, xmax, yvals);
  size_t N = chih.size();
  Spheral::CHIView chi = chih.view();
  EXEC_IN_SPACE_BEGIN(TypeParam)
    SPHERAL_ASSERT_EQ(chi.size(), N);
  EXEC_IN_SPACE_END()
  const double xstep = (xmax - xmin)/((double)NV - 1.);
  RAJA::forall<TypeParam>(TRS_UINT(0, NV),
    [=] (size_t i) {
      double x = xmin + xstep*(double)i;
      double rval = gpu_this->func(x);
      double ival = chi(x);
      SPHERAL_ASSERT_FLOAT_EQ(rval, ival);
    });
}

REGISTER_TYPED_TEST_SUITE_P(CubicHermiteInterpolatorTypedTest, CopyAssign, FuncCtorTest,
                            VecCtorTest);

INSTANTIATE_TYPED_TEST_SUITE_P(CubicHermiteInterpolator, CubicHermiteInterpolatorTypedTest,
                               typename Spheral::Test<EXEC_TYPES>::Types, );
