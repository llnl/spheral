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
using MonGArtVisc = Spheral::LimitedMonaghanGingoldViscosity<Spheral::Dim<3>>;
using MonGArtView = Spheral::LimitedMonaghanGingoldViscosityView<Spheral::Dim<3>>;

//static constexpr int NL_count = 100;

// class FakeSPH {
// public:
//   FakeSPH(ArtVisc3D& A) :
//     m_ArtVisc(A) {
//     m_nl = NodeList_t("NL", NL_count, 0);
//   }

// private:
//   ArtVisc3D& m_ArtVisc;
//   NodeList_T m_nl;
// };

class ArtViscTest : public ::testing::Test {
};

// Setting up G Test for ArtificialViscosity
TYPED_TEST_SUITE_P(ArtViscTypedTest);
template <typename T> class ArtViscTypedTest : public ArtViscTest {};

// Test copy and assignment constructors
GPU_TYPED_TEST_P(ArtViscTypedTest, DiffInit) {
  WC4Kernel3D w_kernel;
  TableKernel3D t_kernel(w_kernel, 200);
  Scalar kernelExtent = t_kernel.kernelExtent();
  Scalar Cl = kernelExtent;
  Scalar Cq = 2.0*std::pow(0.5*kernelExtent, 2);
  Scalar etaCF = 0.5;
  Scalar etaFF = 0.3;
  bool linear = true;
  bool quad = false;
  MonGArtVisc Q(Cl, Cq, t_kernel, linear, quad, etaCF, etaFF);
  Cl = 11.5; Cq = 12.5;
  Q.Cl(Cl);
  Q.Cq(Cq);
  Q.linearInExpansion(false);
  chai::managed_ptr<ArtViscScalarView3D> Qview = Q.getScalarView();
  SPHERAL_ASSERT_FLOAT_EQ(Qview->Cl(), Cl);
  EXEC_IN_SPACE_BEGIN(TypeParam)
    SPHERAL_ASSERT_FLOAT_EQ(Qview->Cl(), Cl);
  EXEC_IN_SPACE_END()
}

// Test copy and assignment constructors
// GPU_TYPED_TEST_P(TableKernelTypedTest, FillTest) {
// }

REGISTER_TYPED_TEST_SUITE_P(ArtViscTypedTest, DiffInit);//, FillTest);

INSTANTIATE_TYPED_TEST_SUITE_P(ArtificialViscosity, ArtViscTypedTest,
                               typename Spheral::Test<EXEC_TYPES>::Types, );
