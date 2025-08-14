// Debug log printing can be quickly enabled for this unit test by uncommenting the
// definition below even if Spheral was not configured w/ SPHERAL_ENABLE_LOGGER=On.
// #define SPHERAL_ENABLE_LOGGER

#include "chai/ExecutionSpaces.hpp"
#include "chai/Types.hpp"
#include "test-basic-exec-policies.hh"
#include "test-utilities.hh"

#include "Kernel/TableKernel.hh"
#include "Kernel/ExpInvKernel.hh"
#include <Utilities/Logger.hh>

using TableKernel1D = Spheral::TableKernel<Spheral::Dim<1>>;
using TableKernel2D = Spheral::TableKernel<Spheral::Dim<2>>;
using TableKernel3D = Spheral::TableKernel<Spheral::Dim<3>>;
using TKV1D = Spheral::TableKernelView<Spheral::Dim<1>>;
using TKV2D = Spheral::TableKernelView<Spheral::Dim<2>>;
using TKV3D = Spheral::TableKernelView<Spheral::Dim<3>>;
using ExpInv1D = Spheral::ExpInvKernel<Spheral::Dim<1>>;
using ExpInv2D = Spheral::ExpInvKernel<Spheral::Dim<2>>;
using ExpInv3D = Spheral::ExpInvKernel<Spheral::Dim<3>>;

class TableKernelTest : public ::testing::Test {
};

// Setting up G Test for Kernel
TYPED_TEST_SUITE_P(TableKernelTypedTest);
template <typename T> class TableKernelTypedTest : public TableKernelTest {};

// Test copy and assignment constructors
GPU_TYPED_TEST_P(TableKernelTypedTest, CopyAssign) {
  const size_t N = 100;
  const double Hdet = 1.;
  ExpInv2D ei_kernel;
  TableKernel2D ref_kernel(ei_kernel);
  TableKernel2D copy_kernel(ref_kernel);
  const double dx = (ref_kernel.kernelExtent() - 0.) / (double)(N - 1);
  const double x_sample = 0.5*dx*(double)N;
  const double ref_val = ref_kernel.kernelValue(x_sample, Hdet);
  TKV2D tkv1 = ref_kernel.view();
  TKV2D tkv2 = ref_kernel;
  TKV2D tkv3(tkv2);
  TKV2D tkv4 = copy_kernel.view();
  EXEC_IN_SPACE_BEGIN(TypeParam)
    const double val1 = tkv1.kernelValue(x_sample, Hdet);
    const double val2 = tkv2.kernelValue(x_sample, Hdet);
    const double val3 = tkv3.kernelValue(x_sample, Hdet);
    const double val4 = tkv4.kernelValue(x_sample, Hdet);
    SPHERAL_ASSERT_FLOAT_EQ(val1, ref_val);
    SPHERAL_ASSERT_FLOAT_EQ(val2, ref_val);
    SPHERAL_ASSERT_FLOAT_EQ(val3, ref_val);
    SPHERAL_ASSERT_FLOAT_EQ(val4, ref_val);
  EXEC_IN_SPACE_END()
}

REGISTER_TYPED_TEST_SUITE_P(TableKernelTypedTest, CopyAssign);

INSTANTIATE_TYPED_TEST_SUITE_P(TableKernel, TableKernelTypedTest,
                               typename Spheral::Test<EXEC_TYPES>::Types, );
