// Debug log printing can be quickly enabled for this unit test by uncommenting the
// definition below even if Spheral was not configured w/ SPHERAL_ENABLE_LOGGER=On.
// #define SPHERAL_ENABLE_LOGGER

#include "test-basic-exec-policies.hh"
#include "test-utilities.hh"

#include "NodeList/NodeList.hh"

#include <type_traits>
#include <typeinfo>
#include <utility>
#include <vector>

using DIM3 = Spheral::Dim<3>;
using Vector = typename DIM3::Vector;
using SymTensor = typename DIM3::SymTensor;
using NodeList_t = Spheral::NodeList<DIM3>;
using NodeListView_t = Spheral::NodeListView<DIM3>;

static constexpr size_t Ninternal = 12;
static constexpr size_t Nghost = 3;
static constexpr size_t Ntotal = Ninternal + Nghost;

class NodeListViewTest : public ::testing::Test {
public:
  GPUCounters fieldCount;

  void initialize(NodeList_t& nodes) {
    auto& positions = nodes.positions();
    auto& H = nodes.Hfield();
    for (auto i = 0u; i < nodes.numNodes(); ++i) {
      positions[i] = Vector(double(i), double(i + 1u), double(i + 2u));
      H[i] = SymTensor::one();
    }
  }

  auto fieldCallback() {
    return [&](const chai::PointerRecord *, chai::Action action,
               chai::ExecutionSpace space) {
      if (action == chai::ACTION_MOVE)
        (space == chai::CPU) ? fieldCount.DToHCopies++ : fieldCount.HToDCopies++;
      if (action == chai::ACTION_ALLOC)
        (space == chai::CPU) ? fieldCount.HNumAlloc++ : fieldCount.DNumAlloc++;
      if (action == chai::ACTION_FREE)
        (space == chai::CPU) ? fieldCount.HNumFree++ : fieldCount.DNumFree++;
    };
  }
};

TYPED_TEST_SUITE_P(NodeListViewTypedTest);
template <typename T> class NodeListViewTypedTest : public NodeListViewTest {};

//------------------------------------------------------------------------------
// Default construction and explicit construction from FieldViews.
//------------------------------------------------------------------------------
GPU_TYPED_TEST_P(NodeListViewTypedTest, DefaultAndExplicitConstruction) {
  EXEC_IN_SPACE_BEGIN(TypeParam)
    NodeListView_t emptyView;
    SPHERAL_ASSERT_EQ(emptyView.numNodes(), 0u);
    SPHERAL_ASSERT_EQ(emptyView.numInternalNodes(), 0u);
    SPHERAL_ASSERT_EQ(emptyView.firstGhostNode(), 0u);
    SPHERAL_ASSERT_EQ(emptyView.mass().numElements(), 0u);
    SPHERAL_ASSERT_EQ(emptyView.positions().numElements(), 0u);
    SPHERAL_ASSERT_EQ(emptyView.velocity().numElements(), 0u);
    SPHERAL_ASSERT_EQ(emptyView.Hfield().numElements(), 0u);
    SPHERAL_ASSERT_EQ(emptyView.work().numElements(), 0u);

    NodeListView_t sizedView(Ntotal, Ninternal);
    SPHERAL_ASSERT_EQ(sizedView.numNodes(), Ntotal);
    SPHERAL_ASSERT_EQ(sizedView.numInternalNodes(), Ninternal);
    SPHERAL_ASSERT_EQ(sizedView.firstGhostNode(), Ninternal);
    SPHERAL_ASSERT_EQ(sizedView.numGhostNodes(), Nghost);
    SPHERAL_ASSERT_EQ(sizedView.nodeType(0u), Spheral::NodeType::InternalNode);
    SPHERAL_ASSERT_EQ(sizedView.nodeType(Ninternal), Spheral::NodeType::GhostNode);
    SPHERAL_ASSERT_EQ(sizedView.positions().numElements(), 0u);
    SPHERAL_ASSERT_EQ(sizedView.Hfield().numElements(), 0u);
    SPHERAL_ASSERT_EQ(sizedView.nodesPerSmoothingScale(), 2.01);
    SPHERAL_ASSERT_EQ(sizedView.maxNumNeighbors(), 500u);
  EXEC_IN_SPACE_END()

  NodeList_t nodes("NodeListViewExplicit", Ninternal, Nghost);
  gpu_this->initialize(nodes);

  NodeListView_t nodesView(nodes.numNodes(),
                           nodes.firstGhostNode(),
                           nodes.positions().view(),
                           nodes.Hfield().view());

  const NodeListView_t& constNodesView = nodesView;
  SPHERAL_ASSERT_EQ(constNodesView.positions().numElements(), Ntotal);
  SPHERAL_ASSERT_EQ(constNodesView.Hfield().numElements(), Ntotal);

  RAJA::forall<TypeParam>(TRS_UINT(0, nodesView.numNodes()),
    [=] SPHERAL_HOST_DEVICE (size_t i) {
      SPHERAL_ASSERT_EQ(nodesView.numNodes(), Ntotal);
      SPHERAL_ASSERT_EQ(nodesView.numInternalNodes(), Ninternal);
      SPHERAL_ASSERT_EQ(nodesView.firstGhostNode(), Ninternal);
      SPHERAL_ASSERT_EQ(nodesView.position(i).x(), double(i));
      SPHERAL_ASSERT_EQ(nodesView.H(i).xx(), 1.0);
    });
}

//------------------------------------------------------------------------------
// NodeListView copies and moves should preserve the value-view state.
//------------------------------------------------------------------------------
GPU_TYPED_TEST_P(NodeListViewTypedTest, ValueSemantics) {
  NodeList_t nodes("NodeListViewValueSemantics", Ninternal, Nghost);
  gpu_this->initialize(nodes);

  auto nodesView = nodes.view();

  NodeListView_t copiedView(nodesView);
  SPHERAL_ASSERT_EQ(copiedView.numNodes(), Ntotal);
  SPHERAL_ASSERT_EQ(copiedView.numInternalNodes(), Ninternal);
  SPHERAL_ASSERT_EQ(copiedView.firstGhostNode(), Ninternal);
  SPHERAL_ASSERT_EQ(copiedView.mass().numElements(), Ntotal);
  SPHERAL_ASSERT_EQ(copiedView.positions().numElements(), Ntotal);
  SPHERAL_ASSERT_EQ(copiedView.velocity().numElements(), Ntotal);
  SPHERAL_ASSERT_EQ(copiedView.Hfield().numElements(), Ntotal);
  SPHERAL_ASSERT_EQ(copiedView.work().numElements(), Ntotal);
  SPHERAL_ASSERT_EQ(copiedView.position(0).x(), 0.0);
  SPHERAL_ASSERT_EQ(copiedView.H(0).xx(), 1.0);

  NodeListView_t copyAssignedView;
  copyAssignedView = nodesView;
  SPHERAL_ASSERT_EQ(copyAssignedView.numNodes(), Ntotal);
  SPHERAL_ASSERT_EQ(copyAssignedView.numInternalNodes(), Ninternal);
  SPHERAL_ASSERT_EQ(copyAssignedView.firstGhostNode(), Ninternal);
  SPHERAL_ASSERT_EQ(copyAssignedView.positions().numElements(), Ntotal);
  SPHERAL_ASSERT_EQ(copyAssignedView.Hfield().numElements(), Ntotal);
  SPHERAL_ASSERT_EQ(copyAssignedView.position(1).x(), 1.0);
  SPHERAL_ASSERT_EQ(copyAssignedView.H(1).xx(), 1.0);

  NodeListView_t moveConstructedView(std::move(copiedView));
  SPHERAL_ASSERT_EQ(moveConstructedView.numNodes(), Ntotal);
  SPHERAL_ASSERT_EQ(moveConstructedView.numInternalNodes(), Ninternal);
  SPHERAL_ASSERT_EQ(moveConstructedView.firstGhostNode(), Ninternal);
  SPHERAL_ASSERT_EQ(moveConstructedView.positions().numElements(), Ntotal);
  SPHERAL_ASSERT_EQ(moveConstructedView.Hfield().numElements(), Ntotal);
  SPHERAL_ASSERT_EQ(moveConstructedView.position(2).x(), 2.0);
  SPHERAL_ASSERT_EQ(moveConstructedView.H(2).xx(), 1.0);

  NodeListView_t moveAssignedView;
  moveAssignedView = std::move(copyAssignedView);
  SPHERAL_ASSERT_EQ(moveAssignedView.numNodes(), Ntotal);
  SPHERAL_ASSERT_EQ(moveAssignedView.numInternalNodes(), Ninternal);
  SPHERAL_ASSERT_EQ(moveAssignedView.firstGhostNode(), Ninternal);
  SPHERAL_ASSERT_EQ(moveAssignedView.positions().numElements(), Ntotal);
  SPHERAL_ASSERT_EQ(moveAssignedView.Hfield().numElements(), Ntotal);
  SPHERAL_ASSERT_EQ(moveAssignedView.position(3).x(), 3.0);
  SPHERAL_ASSERT_EQ(moveAssignedView.H(3).xx(), 1.0);
}

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
// NodeListView::move should migrate device writes for all contained FieldViews.
//------------------------------------------------------------------------------
GPU_TYPED_TEST_P(NodeListViewTypedTest, MoveCopiesDeviceWritesToHost) {
  {
    NodeList_t nodes("NodeListViewMove", Ninternal, Nghost);
    gpu_this->initialize(nodes);
    nodes.positions().setCallback(gpu_this->fieldCallback());
    nodes.Hfield().setCallback(gpu_this->fieldCallback());

    auto nodesView = nodes.view();
    RAJA::forall<TypeParam>(TRS_UINT(0, nodesView.numNodes()),
      [=] SPHERAL_HOST_DEVICE (size_t i) {
        nodesView.position(i) = Vector(10.0 + double(i),
                                       20.0 + double(i),
                                       30.0 + double(i));
        nodesView.H(i).xx(40.0 + double(i));
        nodesView.H(i).yy(50.0 + double(i));
        nodesView.H(i).zz(60.0 + double(i));
      });

    nodesView.move(chai::CPU);

    for (auto i = 0u; i < nodes.numNodes(); ++i) {
      SPHERAL_ASSERT_EQ(nodes.positions()[i].x(), 10.0 + double(i));
      SPHERAL_ASSERT_EQ(nodes.positions()[i].y(), 20.0 + double(i));
      SPHERAL_ASSERT_EQ(nodes.positions()[i].z(), 30.0 + double(i));
      SPHERAL_ASSERT_EQ(nodes.Hfield()[i].xx(), 40.0 + double(i));
      SPHERAL_ASSERT_EQ(nodes.Hfield()[i].yy(), 50.0 + double(i));
      SPHERAL_ASSERT_EQ(nodes.Hfield()[i].zz(), 60.0 + double(i));
    }
  }

  GPUCounters refCount;
  if (typeid(GPU_TEST_TYPE) == typeid(TypeParam)) {
    refCount = {2, 2, 0, 2, 0, 2};
  }
  COMP_COUNTERS(gpu_this->fieldCount, refCount);
}

//------------------------------------------------------------------------------
// NodeListView::touch should mark both contained FieldViews as CPU-modified.
//------------------------------------------------------------------------------
GPU_TYPED_TEST_P(NodeListViewTypedTest, TouchForwardsToContainedViews) {
  {
    NodeList_t nodes("NodeListViewTouch", Ninternal, Nghost);
    gpu_this->initialize(nodes);
    nodes.positions().setCallback(gpu_this->fieldCallback());
    nodes.Hfield().setCallback(gpu_this->fieldCallback());

    auto nodesView = nodes.view();
    RAJA::forall<TypeParam>(TRS_UINT(0, nodesView.numNodes()),
      [=] SPHERAL_HOST_DEVICE (size_t i) {
        SPHERAL_ASSERT_EQ(nodesView.position(i).x(), double(i));
        SPHERAL_ASSERT_EQ(nodesView.H(i).xx(), 1.0);
      });

    nodesView.touch(chai::CPU);

    nodesView = nodes.view();
    RAJA::forall<TypeParam>(TRS_UINT(0, nodesView.numNodes()),
      [=] SPHERAL_HOST_DEVICE (size_t i) {
        SPHERAL_ASSERT_EQ(nodesView.position(i).x(), double(i));
        SPHERAL_ASSERT_EQ(nodesView.H(i).xx(), 1.0);
      });
    nodesView.touch(chai::CPU);
  }

  GPUCounters refCount;
  if (typeid(GPU_TEST_TYPE) == typeid(TypeParam)) {
    refCount = {4, 0, 0, 2, 0, 2};
  }
  COMP_COUNTERS(gpu_this->fieldCount, refCount);
}

//------------------------------------------------------------------------------
// NodeListView::touch is the usage hook after host-side Field mutation.
//------------------------------------------------------------------------------
GPU_TYPED_TEST_P(NodeListViewTypedTest, TouchAfterHostMutation) {
  NodeList_t nodes("NodeListViewHostTouch", Ninternal, Nghost);
  gpu_this->initialize(nodes);

  auto nodesView = nodes.view();
  RAJA::forall<TypeParam>(TRS_UINT(0, nodesView.numNodes()),
    [=] SPHERAL_HOST_DEVICE (size_t i) {
      SPHERAL_ASSERT_EQ(nodesView.position(i).x(), double(i));
      SPHERAL_ASSERT_EQ(nodesView.H(i).xx(), 1.0);
    });

  nodes.positions()[0] = Vector(100.0, 200.0, 300.0);
  nodes.Hfield()[0].xx(400.0);
  nodes.Hfield()[0].yy(500.0);
  nodes.Hfield()[0].zz(600.0);
  nodesView.touch(chai::CPU);

  nodesView = nodes.view();
  RAJA::forall<TypeParam>(TRS_UINT(0, nodesView.numNodes()),
    [=] SPHERAL_HOST_DEVICE (size_t i) {
      if (i == 0u) {
        SPHERAL_ASSERT_EQ(nodesView.position(i).x(), 100.0);
        SPHERAL_ASSERT_EQ(nodesView.position(i).y(), 200.0);
        SPHERAL_ASSERT_EQ(nodesView.position(i).z(), 300.0);
        SPHERAL_ASSERT_EQ(nodesView.H(i).xx(), 400.0);
        SPHERAL_ASSERT_EQ(nodesView.H(i).yy(), 500.0);
        SPHERAL_ASSERT_EQ(nodesView.H(i).zz(), 600.0);
      }
    });
  nodesView.touch(chai::CPU);
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
// Appending internal nodes should rebind the inherited view to the resized Fields.
//------------------------------------------------------------------------------
GPU_TYPED_TEST_P(NodeListViewTypedTest, AppendInternalNodesRefreshesInheritedView) {
  NodeList_t nodes("NodeListViewAppend", Ninternal, Nghost);
  gpu_this->initialize(nodes);

  NodeListView_t& inheritedView = nodes;
  auto oldView = nodes.view();

  const std::vector<size_t> sourceIDs = {1u, Ninternal - 1u};
  const auto packedValues = nodes.packNodeFieldValues(sourceIDs);
  nodes.appendInternalNodes(sourceIDs.size(), packedValues);

  const auto expectedInternal = Ninternal + sourceIDs.size();
  const auto expectedTotal = Ntotal + sourceIDs.size();
  SPHERAL_ASSERT_EQ(oldView.numNodes(), Ntotal);
  SPHERAL_ASSERT_EQ(nodes.numNodes(), expectedTotal);
  SPHERAL_ASSERT_EQ(nodes.numInternalNodes(), expectedInternal);
  SPHERAL_ASSERT_EQ(inheritedView.numNodes(), expectedTotal);
  SPHERAL_ASSERT_EQ(inheritedView.firstGhostNode(), expectedInternal);
  SPHERAL_ASSERT_EQ(inheritedView.mass().numElements(), expectedTotal);
  SPHERAL_ASSERT_EQ(inheritedView.positions().numElements(), expectedTotal);
  SPHERAL_ASSERT_EQ(inheritedView.velocity().numElements(), expectedTotal);
  SPHERAL_ASSERT_EQ(inheritedView.Hfield().numElements(), expectedTotal);
  SPHERAL_ASSERT_EQ(inheritedView.work().numElements(), expectedTotal);
  SPHERAL_ASSERT_EQ(inheritedView.position(Ninternal).x(), double(sourceIDs[0]));
  SPHERAL_ASSERT_EQ(inheritedView.position(Ninternal + 1u).x(), double(sourceIDs[1]));
  SPHERAL_ASSERT_EQ(inheritedView.H(Ninternal).xx(), 1.0);
  SPHERAL_ASSERT_EQ(inheritedView.H(Ninternal + 1u).xx(), 1.0);

  auto nodesView = nodes.view();
  SPHERAL_ASSERT_EQ(nodesView.position(Ninternal).x(), double(sourceIDs[0]));
  SPHERAL_ASSERT_EQ(nodesView.position(Ninternal + 1u).x(), double(sourceIDs[1]));
  SPHERAL_ASSERT_EQ(nodesView.H(Ninternal).xx(), 1.0);
  SPHERAL_ASSERT_EQ(nodesView.H(Ninternal + 1u).xx(), 1.0);

  // The first old ghost node should have moved after the appended internals.
  SPHERAL_ASSERT_EQ(nodesView.position(expectedInternal).x(), double(Ninternal));
}

//------------------------------------------------------------------------------
// Reordering should rebind the inherited view to the reordered internal Fields.
//------------------------------------------------------------------------------
GPU_TYPED_TEST_P(NodeListViewTypedTest, ReorderNodesRefreshesInheritedView) {
  NodeList_t nodes("NodeListViewReorder", Ninternal, Nghost);
  gpu_this->initialize(nodes);

  NodeListView_t& inheritedView = nodes;
  std::vector<size_t> newOrdering(Ninternal);
  for (auto i = 0u; i < Ninternal; ++i) newOrdering[i] = Ninternal - i - 1u;
  nodes.reorderNodes(newOrdering);

  SPHERAL_ASSERT_EQ(nodes.numNodes(), Ninternal);
  SPHERAL_ASSERT_EQ(nodes.numInternalNodes(), Ninternal);
  SPHERAL_ASSERT_EQ(nodes.numGhostNodes(), 0u);
  SPHERAL_ASSERT_EQ(inheritedView.numNodes(), Ninternal);
  SPHERAL_ASSERT_EQ(inheritedView.firstGhostNode(), Ninternal);
  SPHERAL_ASSERT_EQ(inheritedView.mass().numElements(), Ninternal);
  SPHERAL_ASSERT_EQ(inheritedView.positions().numElements(), Ninternal);
  SPHERAL_ASSERT_EQ(inheritedView.velocity().numElements(), Ninternal);
  SPHERAL_ASSERT_EQ(inheritedView.Hfield().numElements(), Ninternal);
  SPHERAL_ASSERT_EQ(inheritedView.work().numElements(), Ninternal);
  SPHERAL_ASSERT_EQ(inheritedView.position(0).x(), double(Ninternal - 1u));
  SPHERAL_ASSERT_EQ(inheritedView.position(Ninternal - 1u).x(), 0.0);

  auto nodesView = nodes.view();
  RAJA::forall<TypeParam>(TRS_UINT(0, nodesView.numNodes()),
    [=] SPHERAL_HOST_DEVICE (size_t i) {
      SPHERAL_ASSERT_EQ(nodesView.position(i).x(), double(Ninternal - i - 1u));
      SPHERAL_ASSERT_EQ(nodesView.H(i).xx(), 1.0);
    });
}

//------------------------------------------------------------------------------
// NodeList owns host state while inheriting allocation-free view state.
//------------------------------------------------------------------------------
GPU_TYPED_TEST_P(NodeListViewTypedTest, HostStateAndInheritedView) {
  static_assert(std::is_base_of<NodeListView_t, NodeList_t>::value,
                "NodeList should inherit NodeListView");

  NodeList_t nodes("NodeListHostState", Ninternal, Nghost);

  SPHERAL_ASSERT_EQ(nodes.name(), "NodeListHostState");
  SPHERAL_ASSERT_EQ(nodes.numNodes(), Ntotal);
  SPHERAL_ASSERT_EQ(nodes.numInternalNodes(), Ninternal);
  SPHERAL_ASSERT_EQ(nodes.numGhostNodes(), Nghost);
  SPHERAL_ASSERT_EQ(nodes.firstGhostNode(), Ninternal);

  NodeListView_t& inheritedView = nodes;
  SPHERAL_ASSERT_EQ(inheritedView.numNodes(), Ntotal);
  SPHERAL_ASSERT_EQ(inheritedView.mass().numElements(), Ntotal);
  SPHERAL_ASSERT_EQ(inheritedView.positions().numElements(), Ntotal);
  SPHERAL_ASSERT_EQ(inheritedView.velocity().numElements(), Ntotal);
  SPHERAL_ASSERT_EQ(inheritedView.Hfield().numElements(), Ntotal);
  SPHERAL_ASSERT_EQ(inheritedView.work().numElements(), Ntotal);

  nodes.numGhostNodes(2u*Nghost);
  SPHERAL_ASSERT_EQ(nodes.numNodes(), Ninternal + 2u*Nghost);
  SPHERAL_ASSERT_EQ(nodes.view().numNodes(), Ninternal + 2u*Nghost);
  SPHERAL_ASSERT_EQ(inheritedView.numNodes(), Ninternal + 2u*Nghost);
  SPHERAL_ASSERT_EQ(inheritedView.positions().numElements(), Ninternal + 2u*Nghost);

  nodes.nodesPerSmoothingScale(3.25);
  nodes.maxNumNeighbors(123u);
  nodes.hmin(1.0e-8);
  nodes.hmax(1.0e8);
  nodes.hminratio(0.25);
  SPHERAL_ASSERT_EQ(nodes.nodesPerSmoothingScale(), 3.25);
  SPHERAL_ASSERT_EQ(nodes.maxNumNeighbors(), 123u);
  SPHERAL_ASSERT_EQ(nodes.hmin(), 1.0e-8);
  SPHERAL_ASSERT_EQ(nodes.hmax(), 1.0e8);
  SPHERAL_ASSERT_EQ(nodes.hminratio(), 0.25);

  auto nodesView = nodes.view();
  SPHERAL_ASSERT_EQ(nodesView.nodesPerSmoothingScale(), 3.25);
  SPHERAL_ASSERT_EQ(nodesView.maxNumNeighbors(), 123u);
  SPHERAL_ASSERT_EQ(nodesView.hmin(), 1.0e-8);
  SPHERAL_ASSERT_EQ(nodesView.hmax(), 1.0e8);
  SPHERAL_ASSERT_EQ(nodesView.hminratio(), 0.25);
}

REGISTER_TYPED_TEST_SUITE_P(NodeListViewTypedTest,
                            DefaultAndExplicitConstruction,
                            ValueSemantics,
                            ViewCapture,
                            MoveCopiesDeviceWritesToHost,
                            TouchForwardsToContainedViews,
                            TouchAfterHostMutation,
                            ResizeAndReacquire,
                            AppendInternalNodesRefreshesInheritedView,
                            ReorderNodesRefreshesInheritedView,
                            HostStateAndInheritedView);

INSTANTIATE_TYPED_TEST_SUITE_P(NodeList,
                               NodeListViewTypedTest,
                               typename Spheral::Test<EXEC_TYPES>::Types, );
