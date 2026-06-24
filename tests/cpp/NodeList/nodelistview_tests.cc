// Debug log printing can be quickly enabled for this unit test by uncommenting the
// definition below even if Spheral was not configured w/ SPHERAL_ENABLE_LOGGER=On.
// #define SPHERAL_ENABLE_LOGGER

#include "test-basic-exec-policies.hh"
#include "test-utilities.hh"

#include "NodeList/NodeList.hh"

#include <type_traits>

using DIM3 = Spheral::Dim<3>;
using Vector = typename DIM3::Vector;
using SymTensor = typename DIM3::SymTensor;
using NodeList_t = Spheral::NodeList<DIM3>;
using NodeListBase_t = Spheral::NodeListBase<DIM3>;
using NodeListView_t = Spheral::NodeListView<DIM3>;

static constexpr size_t Ninternal = 12;
static constexpr size_t Nghost = 3;
static constexpr size_t Ntotal = Ninternal + Nghost;

class NodeListViewTest : public ::testing::Test {
public:
  void initialize(NodeList_t& nodes) {
    auto& positions = nodes.positions();
    auto& H = nodes.Hfield();
    for (auto i = 0u; i < nodes.numNodes(); ++i) {
      positions[i] = Vector(double(i), double(i + 1u), double(i + 2u));
      H[i] = SymTensor::one();
    }
  }
};

TYPED_TEST_SUITE_P(NodeListViewTypedTest);
template <typename T> class NodeListViewTypedTest : public NodeListViewTest {};

//------------------------------------------------------------------------------
// Capture NodeListView in each configured execution space.
//------------------------------------------------------------------------------
GPU_TYPED_TEST_P(NodeListViewTypedTest, ViewCapture) {
  NodeList_t nodes("NodeListViewCapture", Ninternal, Nghost);
  gpu_this->initialize(nodes);

  auto nodesView = nodes.view();
  static_assert(std::is_same<decltype(nodesView), NodeListView_t>::value,
                "NodeList::view must return NodeListView");

  RAJA::forall<TypeParam>(TRS_UINT(0, nodesView.numNodes()),
    [=] SPHERAL_HOST_DEVICE (size_t i) {
      SPHERAL_ASSERT_EQ(nodesView.numNodes(), Ntotal);
      SPHERAL_ASSERT_EQ(nodesView.numInternalNodes(), Ninternal);
      SPHERAL_ASSERT_EQ(nodesView.firstGhostNode(), Ninternal);
      SPHERAL_ASSERT_EQ(nodesView.positions().numElements(), Ntotal);
      SPHERAL_ASSERT_EQ(nodesView.Hfield().numElements(), Ntotal);

      const auto& ri = nodesView.position(i);
      SPHERAL_ASSERT_EQ(ri.x(), double(i));
      SPHERAL_ASSERT_EQ(ri.y(), double(i + 1u));
      SPHERAL_ASSERT_EQ(ri.z(), double(i + 2u));

      const auto& Hi = nodesView.H(i);
      SPHERAL_ASSERT_EQ(Hi.xx(), 1.0);
      SPHERAL_ASSERT_EQ(Hi.yy(), 1.0);
      SPHERAL_ASSERT_EQ(Hi.zz(), 1.0);
    });
}

//------------------------------------------------------------------------------
// Reacquiring a view after resize should expose the new NodeList state.
//------------------------------------------------------------------------------
GPU_TYPED_TEST_P(NodeListViewTypedTest, ResizeAndReacquire) {
  NodeList_t nodes("NodeListViewResize", Ninternal, Nghost);
  gpu_this->initialize(nodes);

  auto oldView = nodes.view();
  SPHERAL_ASSERT_EQ(oldView.numNodes(), Ntotal);
  SPHERAL_ASSERT_EQ(oldView.numInternalNodes(), Ninternal);

  constexpr size_t newInternal = 2u*Ninternal;
  nodes.numInternalNodes(newInternal);
  gpu_this->initialize(nodes);

  auto newView = nodes.view();
  SPHERAL_ASSERT_EQ(oldView.numNodes(), Ntotal);
  SPHERAL_ASSERT_EQ(oldView.numInternalNodes(), Ninternal);
  SPHERAL_ASSERT_EQ(newView.numNodes(), newInternal + Nghost);
  SPHERAL_ASSERT_EQ(newView.numInternalNodes(), newInternal);

  RAJA::forall<TypeParam>(TRS_UINT(0, newView.numNodes()),
    [=] SPHERAL_HOST_DEVICE (size_t i) {
      SPHERAL_ASSERT_EQ(newView.numNodes(), newInternal + Nghost);
      SPHERAL_ASSERT_EQ(newView.firstGhostNode(), newInternal);
      SPHERAL_ASSERT_EQ(newView.positions().numElements(), newInternal + Nghost);
      SPHERAL_ASSERT_EQ(newView.position(i).x(), double(i));
    });
}

//------------------------------------------------------------------------------
// NodeList should remain usable through the new host base abstraction.
//------------------------------------------------------------------------------
GPU_TYPED_TEST_P(NodeListViewTypedTest, BaseContract) {
  static_assert(std::is_base_of<NodeListBase_t, NodeList_t>::value,
                "NodeList must inherit from NodeListBase");
  static_assert(std::is_base_of<NodeListView_t, NodeList_t>::value,
                "NodeList must inherit from NodeListView");

  NodeList_t nodes("NodeListBaseContract", Ninternal, Nghost);
  NodeListBase_t& base = nodes;

  SPHERAL_ASSERT_EQ(base.name(), "NodeListBaseContract");
  SPHERAL_ASSERT_EQ(base.numNodes(), Ntotal);
  SPHERAL_ASSERT_EQ(base.numInternalNodes(), Ninternal);
  SPHERAL_ASSERT_EQ(base.numGhostNodes(), Nghost);
  SPHERAL_ASSERT_EQ(base.firstGhostNode(), Ninternal);

  base.numGhostNodes(2u*Nghost);
  SPHERAL_ASSERT_EQ(base.numNodes(), Ninternal + 2u*Nghost);
  SPHERAL_ASSERT_EQ(nodes.view().numNodes(), Ninternal + 2u*Nghost);

  base.nodesPerSmoothingScale(3.25);
  base.maxNumNeighbors(123u);
  base.hmin(1.0e-8);
  base.hmax(1.0e8);
  base.hminratio(0.25);
  SPHERAL_ASSERT_EQ(base.nodesPerSmoothingScale(), 3.25);
  SPHERAL_ASSERT_EQ(base.maxNumNeighbors(), 123u);
  SPHERAL_ASSERT_EQ(base.hmin(), 1.0e-8);
  SPHERAL_ASSERT_EQ(base.hmax(), 1.0e8);
  SPHERAL_ASSERT_EQ(base.hminratio(), 0.25);
}

REGISTER_TYPED_TEST_SUITE_P(NodeListViewTypedTest,
                            ViewCapture,
                            ResizeAndReacquire,
                            BaseContract);

INSTANTIATE_TYPED_TEST_SUITE_P(NodeList,
                               NodeListViewTypedTest,
                               typename Spheral::Test<EXEC_TYPES>::Types, );
