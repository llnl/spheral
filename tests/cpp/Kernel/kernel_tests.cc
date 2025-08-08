// Debug log printing can be quickly enabled for this unit test by uncommenting the
// definition below even if Spheral was not configured w/ SPHERAL_ENABLE_LOGGER=On.
// #define SPHERAL_ENABLE_LOGGER

#include "chai/ExecutionSpaces.hpp"
#include "chai/Types.hpp"
#include "test-basic-exec-policies.hh"
#include "test-utilities.hh"

#include "Kernel/ExpInvKernel.hh"
#include <Utilities/Logger.hh>

using ExpInvKernel1D = Spheral::ExpInvKernel<Spheral::Dim<1>>;
using ExpInvKernel2D = Spheral::ExpInvKernel<Spheral::Dim<2>>;
using ExpInvKernel3D = Spheral::ExpInvKernel<Spheral::Dim<3>>;

class ExpInvKernelTest : public ::testing::Test {
};

// Setting up G Test for Kernel
TYPED_TEST_SUITE_P(ExpInvKernelTypedTest);
template <typename T> class ExpInvKernelTypedTest : public ExpInvKernelTest {};

// Test copy and assignment constructors
GPU_TYPED_TEST_P(ExpInvKernelTypedTest, CopyAssign) {
  const size_t N = 100;
  const double Hdet = 1.;
  ExpInvKernel1D ref_kernel;
  const double dx = (ref_kernel.kernelExtent() - 0.) / (double)(N - 1);
  const double x_sample = 0.5*dx*(double)N;
  const double ref_val = ref_kernel.kernelValue(x_sample, Hdet);
  EXEC_IN_SPACE_BEGIN(TypeParam)
    ExpInvKernel1D kernel1d;
    const double val = kernel1d.kernelValue(x_sample, Hdet);
    SPHERAL_ASSERT_EQ(val, ref_val);
  EXEC_IN_SPACE_END()
}

REGISTER_TYPED_TEST_SUITE_P(ExpInvKernelTypedTest, CopyAssign);

INSTANTIATE_TYPED_TEST_SUITE_P(ExpInvKernel, ExpInvKernelTypedTest,
                               typename Spheral::Test<EXEC_TYPES>::Types, );
