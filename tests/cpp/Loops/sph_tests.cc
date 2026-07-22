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
#include "ArtificialViscosity/LimitedMonaghanGingoldViscosity.hh"
#include "chai/managed_ptr.hpp"

#include "Field/FieldList.hh"
#include "Field/FieldView.hh"
#include "Field/FieldListView.hh"
#include <Utilities/Logger.hh>

using DIM3 = Spheral::Dim<3>;
using Vector = Spheral::GeomVector<3>;
using Scalar = double;
using Tensor = Spheral::GeomTensor<3>;
using FieldScalar = Spheral::Field<DIM3, Scalar>;
using FieldListScalarView = Spheral::FieldListView<DIM3, Scalar>;
using FieldListScalar = Spheral::FieldList<DIM3, Scalar>;
using FieldVector = Spheral::Field<DIM3, Vector>;
using FieldListVectorView = Spheral::FieldListView<DIM3, Vector>;
using FieldListVector = Spheral::FieldList<DIM3, Vector>;
using FieldTensor = Spheral::Field<DIM3, Tensor>;
using FieldListTensorView = Spheral::FieldListView<DIM3, Tensor>;
using FieldListTensor = Spheral::FieldList<DIM3, Tensor>;
using NodeList_t = Spheral::NodeList<DIM3>;
using FluidNodeList_t = Spheral::FluidNodeList<DIM3>;
using TreeNeighbor_t = Spheral::TreeNeighbor<DIM3>;
using GammaLawGas_t = Spheral::GammaLawGas<DIM3>;
using TableKernel_t = Spheral::TableKernel<DIM3>;
using BSplineKernel_t = Spheral::BSplineKernel<DIM3>;
using SymTensor = Spheral::GeomSymmetricTensor<3>;
using LimMonGVisc = Spheral::LimitedMonaghanGingoldViscosity<DIM3>;
using LimMonGView = Spheral::LimitedMonaghanGingoldViscosityView<DIM3>;
using MonGVisc = Spheral::MonaghanGingoldViscosity<DIM3>;
using MonGView = Spheral::MonaghanGingoldViscosityView<DIM3>;
using ArtVisc3D = Spheral::ArtificialViscosity<DIM3>;
using ArtViscView = Spheral::ArtificialViscosityView<DIM3, Scalar>;
using QPiType = LimMonGView::ReturnType;
using NPIT = Spheral::NodePairIdxType;
using NPLVec = std::vector<Spheral::NodePairIdxType>;
using NPL = Spheral::NodePairList;
using NPLV = Spheral::NodePairListView;

// Default Testing Size.
static constexpr size_t Nx = 40;
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

enum ArtViscType { MonG, LimMonG };

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
    Spheral::GPUUtils::initGPUs();
    fluid_node_list.registerNeighbor(tree_neighbor);
  }

  template<typename TypeParam>
  void initialize(ArtViscType a_visc) {
    av_type = a_visc;
    // Initialize the positions and H values
    auto& pos = fluid_node_list.positions();
    auto& H = fluid_node_list.Hfield();
    auto pos_v = pos.view();
    auto H_v = H.view();
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
    fluid_node_list.neighbor().updateNodes();
    if (a_visc == MonG) {
      monG = chai::make_managed<MonGView>(1.0, 1.0, false, false);
    } else if (a_visc == LimMonG) {
      limMonG = chai::make_managed<LimMonGView>(1.0, 1.0, false, false, 1.0, 0.2);
    }
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

  chai::managed_ptr<ArtViscView> getArtVisc() {
    if (av_type == MonG) {
      return chai::dynamic_pointer_cast<ArtViscView, MonGView>(monG);
    } else {
      return chai::dynamic_pointer_cast<ArtViscView, LimMonGView>(limMonG);
    }
  }
      
  ~LoopTest() { monG.free(); limMonG.free(); }

  TableKernel_t WT;
  GammaLawGas_t eos;
  FluidNodeList_t fluid_node_list;
  TreeNeighbor_t tree_neighbor;
  ArtViscType av_type;
  chai::managed_ptr<LimMonGView> limMonG;
  chai::managed_ptr<MonGView> monG;
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
GPU_TYPED_TEST_P(LoopTypedTest, SPHTest) {
  if (typeid(GPU_TEST_TYPE) != typeid(TypeParam)) {
    return;
  }
  volatile ArtViscType av = LimMonG;
  gpu_this->template initialize<RAJA::seq_exec>(av);
  FieldScalar fcl_f("fcl", gpu_this->fluid_node_list, 0.0);
  FieldListScalar fcl_list;
  fcl_list.appendField(fcl_f);
  FieldScalar fcq_f("fcq", gpu_this->fluid_node_list, 0.0);
  FieldListScalar fcq_list;
  fcq_list.appendField(fcq_f);
  FieldTensor dvdxq_f("dvdxq", gpu_this->fluid_node_list, Tensor::zero());
  FieldListTensor dvdxq_list;
  dvdxq_list.appendField(dvdxq_f);
  FieldListScalarView fcl = fcl_list.view();
  FieldListScalarView fcq = fcq_list.view();
  FieldListTensorView dvdxq = dvdxq_list.view();
  FieldTensor out_f("out", gpu_this->fluid_node_list, Tensor::zero());
  FieldListTensor out_list;
  out_list.appendField(out_f);
  FieldListTensorView out = out_list.view();
  // Make a radius large enough to encompass diagonal nodes
  const double r2 = 3.*std::pow(1.01*dx, 2);
  auto pairs = gpu_this->findNeighbors(r2);
  auto pairs_v = pairs.view();
  auto Q = gpu_this->getArtVisc();
  size_t npairs = pairs.size();
  RAJA::forall<TypeParam>(TRS_UINT(0u, npairs),
     [=] SPHERAL_HOST_DEVICE (size_t kk) {
       auto i = pairs_v[kk].i_node;
       auto j = pairs_v[kk].j_node;
       auto nodeListi = pairs_v[kk].i_list;
       auto nodeListj = pairs_v[kk].j_list;
       Vector xi(0.0);
       SymTensor Hi(0.0);
       Vector etai(0.0);
       Vector vi(0.0);
       Scalar rhoi(0.0);
       Scalar ci(0.0);
       Vector xj(0.0);
       SymTensor Hj(0.0);
       Vector etaj(0.0);
       Vector vj(0.0);
       Scalar rhoj(0.0);
       Scalar cj(0.0);
       QPiType QPiij(0.0);
       QPiType QPiji(0.0);
       Scalar Qi = 0.0;
       Scalar Qj = 0.0;
       Q->QPiij(QPiij, QPiji, Qi, Qj,
                nodeListi, i, nodeListj, j,
                xi, Hi, etai, vi, rhoi, ci,
                xj, Hj, etaj, vj, rhoj, cj,
                fcl, fcq, dvdxq);
       Tensor localQi(Qi);
       Tensor localQj(Qj);
       Scalar testval = 1.2;
       Tensor& outi = out(nodeListi, i);
       Tensor& outj = out(nodeListj, j);
       outi.atomicSub(testval*localQi);
       outj.atomicSub(testval*localQj);
     });
}

REGISTER_TYPED_TEST_SUITE_P(LoopTypedTest, Start, SPHTest);

INSTANTIATE_TYPED_TEST_SUITE_P(ForAllTests, LoopTypedTest,
                               typename Spheral::Test<EXEC_TYPES>::Types, );
