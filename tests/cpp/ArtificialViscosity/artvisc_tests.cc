// Debug log printing can be quickly enabled for this unit test by uncommenting the
// definition below even if Spheral was not configured w/ SPHERAL_ENABLE_LOGGER=On.
// #define SPHERAL_ENABLE_LOGGER

#include "chai/ExecutionSpaces.hpp"
#include "chai/Types.hpp"
#include "test-basic-exec-policies.hh"
#include "test-utilities.hh"

#include "NodeList/NodeList.hh"
#include "ArtificialViscosity/LimitedMonaghanGingoldViscosity.hh"
#include "Kernel/TableKernel.hh"
#include "Kernel/WendlandC4Kernel.hh"
#include <Utilities/Logger.hh>

using Scalar = Spheral::Dim<3>::Scalar;
using ArtVisc3D = Spheral::ArtificialViscosity<Spheral::Dim<3>>;
using ArtViscScalarView3D = Spheral::ArtificialViscosityView<Spheral::Dim<3>, Scalar>;
using TableKernel3D = Spheral::TableKernel<Spheral::Dim<3>>;
using WC4Kernel3D = Spheral::WendlandC4Kernel<Spheral::Dim<3>>;
using NodeList_t = Spheral::NodeList<Spheral::Dim<3>>;
using LimMonGArtVisc = Spheral::LimitedMonaghanGingoldViscosity<Spheral::Dim<3>>;
using LimMonGArtView = Spheral::LimitedMonaghanGingoldViscosityView<Spheral::Dim<3>>;

class ArtViscTest : public ::testing::Test {
};

// Setting up G Test for ArtificialViscosity
TYPED_TEST_SUITE_P(ArtViscTypedTest);
template <typename T> class ArtViscTypedTest : public ArtViscTest {};

// Test copy and assignment constructors
GPU_TYPED_TEST_P(ArtViscTypedTest, InitTests) {
  WC4Kernel3D w_kernel;
  TableKernel3D t_kernel(w_kernel, 200);
  Scalar kernelExtent = t_kernel.kernelExtent();
  Scalar Cl = kernelExtent;
  Scalar Cq = 2.0*std::pow(0.5*kernelExtent, 2);
  Scalar etaCF = 0.5;
  Scalar etaFF = 0.3;
  bool linear = true;
  bool quad = false;
  // Initialize with certain variables
  LimMonGArtVisc Q(Cl, Cq, t_kernel, linear, quad, etaCF, etaFF);
  // This grabs the dynamically casted view, used in the code
  chai::managed_ptr<ArtViscScalarView3D> Qview = Q.getScalarView();
  // This grabs the view directly, only used for testing
  chai::managed_ptr<LimMonGArtView> QAview = Q.getView();

  // Test if initialized variables are set properly
  SPHERAL_ASSERT_FLOAT_EQ(Q.etaCritFrac(), etaCF);
  SPHERAL_ASSERT_FLOAT_EQ(QAview->etaCritFrac(), etaCF);
  SPHERAL_ASSERT_FLOAT_EQ(Qview->Cl(), Cl);
  SPHERAL_ASSERT_EQ(QAview->linearInExpansion(), linear);
  EXEC_IN_SPACE_BEGIN(TypeParam)
    SPHERAL_ASSERT_FLOAT_EQ(Qview->Cl(), Cl);
    SPHERAL_ASSERT_FLOAT_EQ(QAview->etaCritFrac(), etaCF);
    SPHERAL_ASSERT_EQ(QAview->linearInExpansion(), linear);
  EXEC_IN_SPACE_END()

  // Modify variables in the value class
  Cl = 11.5; Cq = 12.5; linear = false; etaCF = 1.3;
  Q.Cl(Cl);
  Q.Cq(Cq);
  Q.linearInExpansion(linear);
  Q.etaCritFrac(etaCF);

  // Test that the value and view classes are properly updated
  SPHERAL_ASSERT_FLOAT_EQ(Q.etaCritFrac(), etaCF);
  SPHERAL_ASSERT_FLOAT_EQ(QAview->etaCritFrac(), etaCF);
  SPHERAL_ASSERT_FLOAT_EQ(Qview->Cl(), Cl);
  SPHERAL_ASSERT_EQ(QAview->linearInExpansion(), linear);
  EXEC_IN_SPACE_BEGIN(TypeParam)
    SPHERAL_ASSERT_FLOAT_EQ(Qview->Cl(), Cl);
    SPHERAL_ASSERT_FLOAT_EQ(QAview->etaCritFrac(), etaCF);
    SPHERAL_ASSERT_EQ(QAview->linearInExpansion(), linear);
  EXEC_IN_SPACE_END()
}

REGISTER_TYPED_TEST_SUITE_P(ArtViscTypedTest, InitTests);

INSTANTIATE_TYPED_TEST_SUITE_P(ArtificialViscosity, ArtViscTypedTest,
                               typename Spheral::Test<EXEC_TYPES>::Types, );
