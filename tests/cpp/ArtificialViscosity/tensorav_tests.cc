// Debug log printing can be quickly enabled for this unit test by uncommenting the
// definition below even if Spheral was not configured w/ SPHERAL_ENABLE_LOGGER=On.
// #define SPHERAL_ENABLE_LOGGER

#include "chai/ExecutionSpaces.hpp"
#include "chai/Types.hpp"
#include "test-basic-exec-policies.hh"
#include "test-utilities.hh"

#include "NodeList/NodeList.hh"
#include "ArtificialViscosity/TensorMonaghanGingoldViscosity.hh"
#include "Kernel/TableKernel.hh"
#include "Kernel/WendlandC4Kernel.hh"
#include <Utilities/Logger.hh>

using Scalar= Spheral::Dim<3>::Scalar;
using Tensor = Spheral::Dim<3>::Tensor;
using ArtVisc3D = Spheral::ArtificialViscosity<Spheral::Dim<3>>;
using ArtViscTensorView3D = Spheral::ArtificialViscosityView<Spheral::Dim<3>, Tensor>;
using TableKernel3D = Spheral::TableKernel<Spheral::Dim<3>>;
using WC4Kernel3D = Spheral::WendlandC4Kernel<Spheral::Dim<3>>;
using NodeList_t = Spheral::NodeList<Spheral::Dim<3>>;
using TensorArtVisc = Spheral::TensorMonaghanGingoldViscosity<Spheral::Dim<3>>;
using TensorArtView = Spheral::TensorMonaghanGingoldViscosityView<Spheral::Dim<3>>;

class TensorAVTest : public ::testing::Test {
};

// Setting up G Test for ArtificialViscosity
TYPED_TEST_SUITE_P(TensorAVTypedTest);
template <typename T> class TensorAVTypedTest : public TensorAVTest {};

// Test copy and assignment constructors
GPU_TYPED_TEST_P(TensorAVTypedTest, InitTests) {
  WC4Kernel3D w_kernel;
  TableKernel3D t_kernel(w_kernel, 200);
  Scalar kernelExtent = t_kernel.kernelExtent();
  Scalar Cl = kernelExtent;
  Scalar Cq = 2.0*std::pow(0.5*kernelExtent, 2);
  // Initialize with certain variables
  TensorArtVisc Q(Cl, Cq, t_kernel);
  // This grabs the dynamically casted view, used in the code
  chai::managed_ptr<ArtViscTensorView3D> Qview = Q.getTensorView();

  // Test if initialized variables are set properly
  SPHERAL_ASSERT_FLOAT_EQ(Q.Cl(), Cl);
  SPHERAL_ASSERT_FLOAT_EQ(Qview->Cl(), Cl);
  SPHERAL_ASSERT_FLOAT_EQ(Qview->Cq(), Cq);
  EXEC_IN_SPACE_BEGIN(TypeParam)
    SPHERAL_ASSERT_FLOAT_EQ(Qview->Cl(), Cl);
    SPHERAL_ASSERT_FLOAT_EQ(Qview->Cq(), Cq);
  EXEC_IN_SPACE_END()

  // Modify variables in the value class
  Cl = 11.5; Cq = 12.5;
  Q.Cl(Cl);
  Q.Cq(Cq);

  // Test that the value and view classes are properly updated
  SPHERAL_ASSERT_FLOAT_EQ(Qview->Cl(), Cl);
  EXEC_IN_SPACE_BEGIN(TypeParam)
    SPHERAL_ASSERT_FLOAT_EQ(Qview->Cl(), Cl);
  EXEC_IN_SPACE_END()
}

REGISTER_TYPED_TEST_SUITE_P(TensorAVTypedTest, InitTests);

INSTANTIATE_TYPED_TEST_SUITE_P(TensorArtificialViscosity, TensorAVTypedTest,
                               typename Spheral::Test<EXEC_TYPES>::Types, );
