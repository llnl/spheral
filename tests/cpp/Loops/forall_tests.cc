// Debug log printing can be quickly enabled for this unit test by uncommenting the
// definition below even if Spheral was not configured w/ SPHERAL_ENABLE_LOGGER=On.
// #define SPHERAL_ENABLE_LOGGER

#include "chai/ExecutionSpaces.hpp"
#include "chai/Types.hpp"
#include "test-basic-exec-policies.hh"
#include "test-utilities.hh"

#include "Field/Field.hh"

#include "NodeList/FluidNodeList.hh"
#include "Neighbor/TreeNeighbor.hh"
#include "Neighbor/NodePairList.hh"
#include "Material/GammaLawGas.hh"
#include "Kernel/TableKernel.hh"
#include "Kernel/BSplineKernel.hh"

#include "Field/FieldList.hh"
#include "Field/FieldView.hh"
#include "Field/FieldListView.hh"
#include <Utilities/Logger.hh>

using DIM3 = Spheral::Dim<3>;
using Vector = Spheral::GeomVector<3>;
using FieldBase = Spheral::FieldBase<DIM3>;
using FieldDouble = Spheral::Field<DIM3, double>;
using FieldViewDouble = Spheral::FieldView<DIM3, double>;
using FieldListDouble = Spheral::FieldList<DIM3, double>;
using NodeList_t = Spheral::NodeList<DIM3>;
using FluidNodeList_t = Spheral::FluidNodeList<DIM3>;
using TreeNeighbor_t = Spheral::TreeNeighbor<DIM3>;
using GammaLawGas_t = Spheral::GammaLawGas<DIM3>;
using TableKernel_t = Spheral::TableKernel<DIM3>;
using BSplineKernel_t = Spheral::BSplineKernel<DIM3>;
using SymTensor = Spheral::GeomSymmetricTensor<3>;
using NPIT = Spheral::NodePairIdxType;
using NPLVec = std::vector<Spheral::NodePairIdxType>;
using NPL = Spheral::NodePairList;
using NPLV = Spheral::NodePairListView;

// Default Testing Size.
static constexpr size_t Nx = 20;
static constexpr size_t N = Nx*Nx*Nx;
static constexpr size_t N_table = 100;
static constexpr double hmin = 1.0e-20;
static constexpr double hmax = 1.0e20;
static constexpr double hminratio = 0.1;
static constexpr double nPerh = 2.01;
static constexpr size_t maxNumNeighbors = 500;
static constexpr double rhoMin = 1.0e-10;
static constexpr double rhoMax = 1.0e10;
static constexpr double kernelExtents = 2.0;
static Spheral::PhysicalConstants units(1.0, 1.0, 1.0);
static constexpr double dx = 0.1;
static constexpr double L = dx*(Nx - 1);

class LoopTest : public ::testing::Test {
public:
  // Constructor
  LoopTest(Vector xmin = Vector(0.0, 0.0, 0.0),
           Vector xmax = Vector(L, L, L),
           std::string name = "FNL", size_t count = N) :
    WT(BSplineKernel_t(), N_table),
    eos(2.0, 2.0, units, -1E200, 1E200,
        Spheral::MaterialPressureMinType::PressureFloor, 0.0),
    fluid_node_list(name, eos, count, 0,
                    hmin, hmax, hminratio,
                    nPerh, maxNumNeighbors,
                    rhoMin, rhoMax),
    tree_neighbor(fluid_node_list,
                  Spheral::NeighborSearchType::GatherScatter,
                  kernelExtents, xmin, xmax) {
    fluid_node_list.registerNeighbor(tree_neighbor);
  }

  NPL findNeighbors(double r2) {
    NPLVec nplvec;
    auto pos = fluid_node_list.positions();
    auto Ns = pos.size();
    for (auto i = 0u; i < Ns; ++i) {
      for (auto j = i + 1u; j < Ns; ++j) {
        auto dvec = pos[i] - pos[j];
        if (dvec.magnitude2() < r2) {
          NPIT nit(i, 0, j, 0);
          nplvec.push_back(nit);
        }
      }
    }
    return NPL(std::move(nplvec));
  }

  TableKernel_t WT;
  GammaLawGas_t eos;
  FluidNodeList_t fluid_node_list;
  TreeNeighbor_t tree_neighbor;
};

// Setting up G Test for Loops
TYPED_TEST_SUITE_P(LoopTypedTest);
template <typename T> class LoopTypedTest : public LoopTest {};

//------------------------------------------------------------------------------
// Start for utilizing hardware
//------------------------------------------------------------------------------
GPU_TYPED_TEST_P(LoopTypedTest, Start) {
  EXEC_IN_SPACE_BEGIN(TypeParam)
    Vector vec(1., 2., 3.);
  EXEC_IN_SPACE_END()
}

//------------------------------------------------------------------------------
// Basic forall
//------------------------------------------------------------------------------
GPU_TYPED_TEST_P(LoopTypedTest, BasicForAll) {
  auto& pos = gpu_this->fluid_node_list.positions();
  auto pos_v = pos.view();
  auto& H = gpu_this->fluid_node_list.Hfield();
  auto H_v = H.view();
  FieldDouble field("test", gpu_this->fluid_node_list, 0.0);
  FieldListDouble field_list;
  field_list.appendField(field);
  const int mx = Nx;
  RAJA::forall<TypeParam>(TRS_UINT(0, pos.numElements()),
     [=] SPHERAL_HOST_DEVICE (int i) {
       int idx = i;
       double x = int(idx%mx)*dx;
       idx /= mx;
       double y = int(idx%mx)*dx;
       idx /= mx;
       double z= int(idx)*dx;
       pos_v[i] = Vector(x, y, z);
       H_v[i] = SymTensor::one();
     });
  pos.move(chai::CPU);
  H.move(chai::CPU);
  gpu_this->fluid_node_list.neighbor().updateNodes();
  // Make a radius large enough to encompass diagonal nodes
  const double r2 = 3.*std::pow(1.01*dx, 2);
  auto field_v = field.view();
  auto pairs = gpu_this->findNeighbors(r2);
  auto pairs_v = pairs.view();
  auto fl_v = field_list.view();
  RAJA::forall<TypeParam>(TRS_UINT(0, pairs.size()),
     [=] SPHERAL_HOST_DEVICE (int kk) {
       auto i = pairs_v[kk].i_node;
       auto j = pairs_v[kk].j_node;
       auto nli = pairs_v[kk].i_list;
       auto nlj = pairs_v[kk].j_list;
       RAJA::atomicAdd<RAJA::auto_atomic>(&fl_v(nli, i), 1.);
       RAJA::atomicAdd<RAJA::auto_atomic>(&fl_v(nlj, j), 1.);
     });
  fl_v.move(chai::CPU);
  std::vector<int> ref_count(N, 0);
  for (auto kk = 0u; kk < pairs.size(); ++kk) {
    auto i = pairs[kk].i_node;
    auto j = pairs[kk].j_node;
    ref_count[i] += 1;
    ref_count[j] += 1;
  }
  for (auto kk = 0u; kk < N; ++kk) {
    // Make sure the node pair list is reasonable
    SPHERAL_ASSERT_TRUE(ref_count[kk] >= 7);
    SPHERAL_ASSERT_EQ(ref_count[kk], field_list(0, kk));
  }
}

REGISTER_TYPED_TEST_SUITE_P(LoopTypedTest, Start, BasicForAll);

INSTANTIATE_TYPED_TEST_SUITE_P(ForAllTests, LoopTypedTest,
                               typename Spheral::Test<EXEC_TYPES>::Types, );
