#include "test-basic-exec-policies.hh"
#include "test-utilities.hh"

#include "Neighbor/TreeNeighbor.hh"
#include "NodeList/NodeList.hh"
#include "Threading/GPUUtils.hh"

#include <algorithm>
#include <type_traits>

using Dimension = Spheral::Dim<1>;
using Vector = Dimension::Vector;
using SymTensor = Dimension::SymTensor;
using NodeList = Spheral::NodeList<Dimension>;
using TreeNeighbor = Spheral::TreeNeighbor<Dimension>;

namespace {

TreeNeighbor
makeNeighbor(NodeList& nodes) {
  return TreeNeighbor(nodes,
                      Spheral::NeighborSearchType::GatherScatter,
                      2.0,
                      Vector(0.0),
                      Vector(1.0));
}

std::vector<int>
viewCandidates(const TreeNeighbor::ViewType& view,
               const Vector& position,
               const SymTensor& H) {
  std::vector<int> result(view.numMembers());
  size_t count = 0u;
  auto* resultPtr = result.data();
  auto* countPtr = &count;
  view.forEachCandidate(position, H,
    [=] SPHERAL_HOST_DEVICE (const int nodeID) {
      resultPtr[(*countPtr)++] = nodeID;
    });
  result.resize(count);
  return result;
}

// Exercise variable-length candidate generation using the same two-pass
// structure needed by a future ConnectivityMap GPU implementation.  The first
// pass counts each query's candidates.  The host converts those counts to
// disjoint output ranges, and the second pass repeats the traversal to fill
// those ranges without atomics.  ConnectivityMap can use the same structure
// with NodePairIdxType output and its exact gather/scatter test in each pass.
template<typename ExecPolicy>
std::vector<std::vector<int>>
twoPassViewCandidates(const TreeNeighbor::ViewType& view,
                      const std::vector<Vector>& positions,
                      const std::vector<SymTensor>& Hs) {
  using NodeID = TreeNeighbor::ViewType::NodeID;

  const auto numQueries = positions.size();
  EXPECT_EQ(Hs.size(), numQueries);
  if (Hs.size() != numQueries) return {};
  const auto executionSpace = (std::is_same<ExecPolicy, GPU_TEST_TYPE>::value ?
                               chai::GPU : chai::CPU);

  chai::ManagedArray<Vector> queryPositions(numQueries);
  chai::ManagedArray<SymTensor> queryHs(numQueries);
  chai::ManagedArray<size_t> candidateCounts(numQueries);
  for (auto query = 0u; query < numQueries; ++query) {
    queryPositions[query] = positions[query];
    queryHs[query] = Hs[query];
  }
  queryPositions.registerTouch(chai::CPU);
  queryHs.registerTouch(chai::CPU);

  // Count candidates on the selected execution space.
  auto activeView = view;
  activeView.move(executionSpace);
  queryPositions.move(executionSpace);
  queryHs.move(executionSpace);
  candidateCounts.move(executionSpace);
  const auto* activePositions = queryPositions.data(executionSpace, false);
  const auto* activeHs = queryHs.data(executionSpace, false);
  auto* activeCounts = candidateCounts.data(executionSpace, false);
  RAJA::forall<ExecPolicy>(TRS_UINT(0u, numQueries),
    [=] SPHERAL_HOST_DEVICE (const size_t query) {
      size_t count = 0u;
      activeView.forEachCandidate(activePositions[query], activeHs[query],
        [=, &count] SPHERAL_HOST_DEVICE (const NodeID) {
          ++count;
        });
      activeCounts[query] = count;
    });
  GPU_ERROR_CHECK;

  // Build an exclusive prefix sum on the host.  These offsets give each query
  // an exact-size, non-overlapping range in the device output buffer.
  candidateCounts.move(chai::CPU);
  const auto* hostCounts = candidateCounts.data(chai::CPU, false);
  std::vector<size_t> hostOffsets(numQueries + 1u, 0u);
  for (auto query = 0u; query < numQueries; ++query) {
    hostOffsets[query + 1u] = hostOffsets[query] + hostCounts[query];
  }

  chai::ManagedArray<size_t> candidateOffsets(numQueries + 1u);
  for (auto query = 0u; query <= numQueries; ++query) {
    candidateOffsets[query] = hostOffsets[query];
  }
  candidateOffsets.registerTouch(chai::CPU);

  // Repeat the traversal and fill the ranges determined by the count pass.
  const auto totalCandidates = hostOffsets.back();
  chai::ManagedArray<NodeID> candidateIDs(totalCandidates);
  if (totalCandidates > 0u) {
    candidateOffsets.move(executionSpace);
    candidateIDs.move(executionSpace);
    const auto* activeOffsets = candidateOffsets.data(executionSpace, false);
    auto* activeCandidates = candidateIDs.data(executionSpace, false);
    RAJA::forall<ExecPolicy>(TRS_UINT(0u, numQueries),
      [=] SPHERAL_HOST_DEVICE (const size_t query) {
        auto output = activeOffsets[query];
        activeView.forEachCandidate(activePositions[query], activeHs[query],
          [=, &output] SPHERAL_HOST_DEVICE (const NodeID nodeID) {
            activeCandidates[output++] = nodeID;
          });
      });
    GPU_ERROR_CHECK;
    candidateIDs.move(chai::CPU);
  }

  const auto* hostCandidates = (totalCandidates > 0u ?
                                candidateIDs.data(chai::CPU, false) : nullptr);
  std::vector<std::vector<int>> result(numQueries);
  for (auto query = 0u; query < numQueries; ++query) {
    if (hostOffsets[query] != hostOffsets[query + 1u]) {
      result[query].assign(hostCandidates + hostOffsets[query],
                           hostCandidates + hostOffsets[query + 1u]);
    }
  }

  queryPositions.free();
  queryHs.free();
  candidateCounts.free();
  candidateOffsets.free();
  candidateIDs.free();
  activeView.move(chai::CPU);
  return result;
}

} // anonymous namespace

//------------------------------------------------------------------------------
// An empty host tree produces an empty projection.
//------------------------------------------------------------------------------
TEST(TreeNeighborView, EmptyTree) {
  NodeList nodes("empty tree view", 0u, 0u);
  auto neighbor = makeNeighbor(nodes);
  const auto view = neighbor.view();

  EXPECT_TRUE(view.empty());
  EXPECT_EQ(view.numCells(), 0u);
  EXPECT_EQ(view.numMembers(), 0u);
  EXPECT_EQ(view.numDaughters(), 0u);
  EXPECT_EQ(view.numLevels(), 0u);
}

//------------------------------------------------------------------------------
// Cells are flattened by level and key.  Members and daughters use contiguous
// offset/count ranges, and daughter values are global cell indices.
//------------------------------------------------------------------------------
TEST(TreeNeighborView, DeterministicFlattening) {
  NodeList nodes("populated tree view", 3u, 0u);
  nodes.positions()[0] = Vector(0.10);
  nodes.positions()[1] = Vector(0.30);
  nodes.positions()[2] = Vector(0.90);
  for (auto i = 0u; i < nodes.numNodes(); ++i) {
    nodes.Hfield()[i] = 8.0*SymTensor::one();
  }

  auto neighbor = makeNeighbor(nodes);
  neighbor.updateNodes();
  const auto view = neighbor.view();

  ASSERT_EQ(view.numLevels(), 3u);
  ASSERT_EQ(view.numCells(), 6u);
  EXPECT_EQ(view.cell(0u).key, 0u);
  EXPECT_EQ(view.cell(1u).key, 0u);
  EXPECT_EQ(view.cell(2u).key, 1u);
  EXPECT_EQ(view.cell(3u).key, 0u);
  EXPECT_EQ(view.cell(4u).key, 1u);
  EXPECT_EQ(view.cell(5u).key, 3u);

  ASSERT_EQ(view.numDaughters(), 5u);
  EXPECT_EQ(view.cell(0u).daughterOffset, 0u);
  EXPECT_EQ(view.cell(0u).daughterCount, 2u);
  EXPECT_EQ(view.daughter(0u), 1u);
  EXPECT_EQ(view.daughter(1u), 2u);
  EXPECT_EQ(view.cell(1u).daughterOffset, 2u);
  EXPECT_EQ(view.cell(1u).daughterCount, 2u);
  EXPECT_EQ(view.daughter(2u), 3u);
  EXPECT_EQ(view.daughter(3u), 4u);
  EXPECT_EQ(view.cell(2u).daughterOffset, 4u);
  EXPECT_EQ(view.cell(2u).daughterCount, 1u);
  EXPECT_EQ(view.daughter(4u), 5u);

  ASSERT_EQ(view.numMembers(), 3u);
  EXPECT_EQ(view.cell(3u).memberOffset, 0u);
  EXPECT_EQ(view.cell(3u).memberCount, 1u);
  EXPECT_EQ(view.cell(4u).memberOffset, 1u);
  EXPECT_EQ(view.cell(4u).memberCount, 1u);
  EXPECT_EQ(view.cell(5u).memberOffset, 2u);
  EXPECT_EQ(view.cell(5u).memberCount, 1u);
  EXPECT_EQ(view.member(0u), 0);
  EXPECT_EQ(view.member(1u), 1);
  EXPECT_EQ(view.member(2u), 2);
}

//------------------------------------------------------------------------------
// Host-tree mutation invalidates the cached projection.
//------------------------------------------------------------------------------
TEST(TreeNeighborView, InvalidationLifecycle) {
  NodeList nodes("tree view lifecycle", 1u, 0u);
  nodes.positions()[0] = Vector(0.25);
  nodes.Hfield()[0] = 8.0*SymTensor::one();

  auto neighbor = makeNeighbor(nodes);
  neighbor.updateNodes();
  const auto first = neighbor.view();
  const auto reused = neighbor.view();
  EXPECT_EQ(reused.numCells(), first.numCells());
  EXPECT_EQ(reused.cell(2u).key, first.cell(2u).key);

  nodes.positions()[0] = Vector(0.75);
  neighbor.updateNodes();
  const auto rebuilt = neighbor.view();
  EXPECT_EQ(rebuilt.cell(2u).key, 3u);

  neighbor.reinitialize();
  const auto cleared = neighbor.view();
  EXPECT_TRUE(cleared.empty());
}

//------------------------------------------------------------------------------
// The flattened CPU traversal must produce the same scalar-H candidate sets
// as the established pointer traversal.  Ordering is deliberately ignored:
// the projection sorts daughter keys for deterministic flattening.
//------------------------------------------------------------------------------
TEST(TreeNeighborView, ScalarCandidateParity) {
  NodeList nodes("scalar candidate parity", 5u, 1u);
  nodes.positions()[0] = Vector(0.05);
  nodes.positions()[1] = Vector(0.20);
  nodes.positions()[2] = Vector(0.45);
  nodes.positions()[3] = Vector(0.70);
  nodes.positions()[4] = Vector(0.95);
  for (auto i = 0u; i < nodes.numNodes(); ++i) {
    nodes.Hfield()[i] = 8.0*SymTensor::one();
  }

  auto neighbor = makeNeighbor(nodes);
  neighbor.updateNodes();
  for (const auto position: {Vector(0.05), Vector(0.45), Vector(0.95)}) {
    for (const auto H: {4.0, 8.0, 16.0}) {
      for (const auto ghostConnectivity: {false, true}) {
        std::vector<int> hostMaster, hostCoarse;
        neighbor.setMasterList(position, H, hostMaster, hostCoarse,
                               ghostConnectivity);
        auto viewCoarse = viewCandidates(neighbor.view(), position,
                                         SymTensor::one()*H);
        std::sort(hostCoarse.begin(), hostCoarse.end());
        std::sort(viewCoarse.begin(), viewCoarse.end());
        EXPECT_EQ(viewCoarse, hostCoarse);
      }
    }
  }
}

//------------------------------------------------------------------------------
// Tensor-H uses the same effective scalar extent as the host implementation.
// Check it independently so device traversal cannot silently diverge.
//------------------------------------------------------------------------------
TEST(TreeNeighborView, TensorCandidateParity) {
  NodeList nodes("tensor candidate parity", 4u, 0u);
  nodes.positions()[0] = Vector(0.10);
  nodes.positions()[1] = Vector(0.30);
  nodes.positions()[2] = Vector(0.60);
  nodes.positions()[3] = Vector(0.90);
  for (auto i = 0u; i < nodes.numNodes(); ++i) {
    nodes.Hfield()[i] = 8.0*SymTensor::one();
  }

  auto neighbor = makeNeighbor(nodes);
  neighbor.updateNodes();
  for (const auto position: {Vector(0.10), Vector(0.60), Vector(0.90)}) {
    for (const auto H: {SymTensor::one()*4.0,
                        SymTensor::one()*8.0,
                        SymTensor::one()*16.0}) {
      std::vector<int> hostMaster, hostCoarse;
      neighbor.setMasterList(position, H, hostMaster, hostCoarse, false);
      auto viewCoarse = viewCandidates(neighbor.view(), position, H);
      std::sort(hostCoarse.begin(), hostCoarse.end());
      std::sort(viewCoarse.begin(), viewCoarse.end());
      EXPECT_EQ(viewCoarse, hostCoarse);
    }
  }
}

//------------------------------------------------------------------------------
// The read-only view can be captured and accessed in every configured
// execution space.  This checks transport and access, not neighbor traversal.
//------------------------------------------------------------------------------
template<typename ExecPolicy>
class TreeNeighborViewCapture: public ::testing::Test {
public:
  TreeNeighborViewCapture() {
    Spheral::GPUUtils::initGPUs();
  }
};

TYPED_TEST_SUITE(TreeNeighborViewCapture,
                 typename Spheral::Test<EXEC_TYPES>::Types);

GPU_TYPED_TEST(TreeNeighborViewCapture, ReadOnlyCapture) {
  NodeList nodes("captured tree view", 1u, 0u);
  nodes.positions()[0] = Vector(0.75);
  nodes.Hfield()[0] = 8.0*SymTensor::one();

  auto neighbor = makeNeighbor(nodes);
  neighbor.updateNodes();
  const auto view = neighbor.view();

  RAJA::forall<TypeParam>(TRS_UINT(0u, 1u),
    [=] SPHERAL_HOST_DEVICE (size_t) {
      SPHERAL_ASSERT_EQ(view.numLevels(), 3u);
      SPHERAL_ASSERT_EQ(view.numCells(), 3u);
      SPHERAL_ASSERT_EQ(view.numMembers(), 1u);
      SPHERAL_ASSERT_EQ(view.numDaughters(), 2u);
      SPHERAL_ASSERT_EQ(view.cell(2u).key, 3u);
      SPHERAL_ASSERT_EQ(view.member(0u), 0);
      SPHERAL_ASSERT_EQ(view.daughter(0u), 1u);
      SPHERAL_ASSERT_EQ(view.daughter(1u), 2u);
    });
}

//------------------------------------------------------------------------------
// Run complete candidate traversals in every configured execution space and
// compare them with the established pointer-based TreeNeighbor traversal.
//------------------------------------------------------------------------------
GPU_TYPED_TEST(TreeNeighborViewCapture, CandidateTraversalParity) {
  NodeList nodes("device candidate parity", 5u, 1u);
  nodes.positions()[0] = Vector(0.05);
  nodes.positions()[1] = Vector(0.20);
  nodes.positions()[2] = Vector(0.45);
  nodes.positions()[3] = Vector(0.70);
  nodes.positions()[4] = Vector(0.95);
  for (auto i = 0u; i < nodes.numNodes(); ++i) {
    nodes.Hfield()[i] = 8.0*SymTensor::one();
  }

  auto neighbor = makeNeighbor(nodes);
  neighbor.updateNodes();

  const std::vector<Vector> positions{Vector(0.05), Vector(0.45), Vector(0.95)};
  const std::vector<SymTensor> Hs{4.0*SymTensor::one(),
                                  8.0*SymTensor::one(),
                                  16.0*SymTensor::one()};
  auto actual = twoPassViewCandidates<TypeParam>(neighbor.view(), positions, Hs);

  for (auto query = 0u; query < positions.size(); ++query) {
    std::vector<int> hostMaster, expected;
    neighbor.setMasterList(positions[query], Hs[query], hostMaster, expected, true);
    std::sort(expected.begin(), expected.end());
    std::sort(actual[query].begin(), actual[query].end());
    EXPECT_EQ(actual[query], expected);
  }
}
