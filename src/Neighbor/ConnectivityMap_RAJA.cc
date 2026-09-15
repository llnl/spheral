//---------------------------------Spheral++----------------------------------//
// ConnectivityMap_RAJA -- TreeNeighbor candidate generation and exact
// symmetric pair tests using the configured RAJA execution policy.
//----------------------------------------------------------------------------//
#include "config.hh"
#include "Neighbor/ConnectivityMap_RAJA.hh"
#include "Neighbor/TreeNeighbor.hh"
#include "NodeList/NodeList.hh"
#include "Threading/GPUUtils.hh"
#include "Utilities/Timer.hh"

#include <utility>
#include <vector>

namespace Spheral {

namespace {

// A small host-owned buffer with a CHAI view for non-unified memory builds.
// Unified-memory builds can expose the host vector directly.  Do not compile
// this helper for a host-only build: chai::GPU is not a valid execution-space
// pointer slot when CHAI was configured without an accelerator backend.
#if defined(SPHERAL_ENABLE_HIP) || defined(SPHERAL_ENABLE_CUDA)
template<typename Value>
class ConnectivityRAJABuffer {
public:
  explicit ConnectivityRAJABuffer(const size_t size): mHost(size, Value()) {}
  ~ConnectivityRAJABuffer() {
#ifndef SPHERAL_UNIFIED_MEMORY
    GPUUtils::freeMAView(mView);
#endif
  }

  ConnectivityRAJABuffer(const ConnectivityRAJABuffer&) = delete;
  ConnectivityRAJABuffer& operator=(const ConnectivityRAJABuffer&) = delete;

  void moveToGPU() {
#ifndef SPHERAL_UNIFIED_MEMORY
    GPUUtils::initMAView(mView, mHost);
    GPUUtils::move(mView, chai::GPU);
#endif
  }

  void moveToCPU() {
#ifndef SPHERAL_UNIFIED_MEMORY
    GPUUtils::move(mView, chai::CPU);
#endif
  }

  Value* deviceData() {
#ifdef SPHERAL_UNIFIED_MEMORY
    return mHost.data();
#else
    return mView.data(chai::GPU, false);
#endif
  }

  std::vector<Value>& hostData() { return mHost; }

private:
  std::vector<Value> mHost;
#ifndef SPHERAL_UNIFIED_MEMORY
  chai::ManagedArray<Value> mView;
#endif
};
#endif

}

//------------------------------------------------------------------------------
// Build TreeNeighbor node pairs on the GPU using the same gather/scatter
// acceptance test as the CPU implementation. If any NodeList does not use
// TreeNeighbor, return false and let ConnectivityMap use its CPU path.
//------------------------------------------------------------------------------
#if defined(SPHERAL_ENABLE_HIP) || defined(SPHERAL_ENABLE_CUDA)
template<typename Dimension>
bool
ConnectivityMap_RAJA<Dimension>::
tryBuildNodePairs(const double kernelExtent,
                  const bool ghostConnectivity,
                  std::vector<NodePairIdxType>& nodePairs) {
  if (GPUUtils::deviceCount() == 0) return false;

  const auto& nodeLists = this->nodeLists();
  const auto numNodeLists = nodeLists.size();
  if (numNodeLists == 0u) return false;

  using NodeListViewType = typename NodeList<Dimension>::ViewType;
  using TreeViewType = typename TreeNeighbor<Dimension>::ViewType;

  // Use the GPU path only when every NodeList has a TreeNeighbor.
  std::vector<const TreeNeighbor<Dimension>*> treeNeighbors;
  treeNeighbors.reserve(numNodeLists);
  for (const auto* nodeListPtr: nodeLists) {
    const auto* treeNeighbor = dynamic_cast<const TreeNeighbor<Dimension>*>(&nodeListPtr->neighbor());
    if (treeNeighbor == nullptr) return false;
    treeNeighbors.push_back(treeNeighbor);
  }

  // Capture the existing NodeList fields and the lazily-created tree views.
  // The view objects are moved once for the whole connectivity operation and
  // returned to the CPU before host finalization and optional connectivity
  // construction.
  std::vector<NodeListViewType> nodeViews;
  std::vector<TreeViewType> treeViews;
  nodeViews.reserve(numNodeLists);
  treeViews.reserve(numNodeLists);
  for (auto k = 0u; k < numNodeLists; ++k) {
    auto* nodeListPtr = const_cast<NodeList<Dimension>*>(nodeLists[k]);
    nodeViews.push_back(nodeListPtr->view());
    treeViews.push_back(treeNeighbors[k]->view());

    TIME_BEGIN("ConnectivityMap_nodeViewCopyToGPU");
    nodeViews.back().move(chai::GPU);
    TIME_END("ConnectivityMap_nodeViewCopyToGPU");

    TIME_BEGIN("ConnectivityMap_treeViewCopyToGPU");
    treeViews.back().move(chai::GPU);
    TIME_END("ConnectivityMap_treeViewCopyToGPU");
  }

  const auto kernelExtent2 = kernelExtent*kernelExtent;
  const auto numBlocks = numNodeLists*numNodeLists;
  std::vector<std::vector<size_t>> blockCounts(numBlocks);

  // Every source node is present in its own TreeNeighbor master cell.  The
  // direct per-tree traversal therefore covers the same source set as the
  // host master-group loop; omitting host preculling may leave a superset of
  // coarse candidates, but the exact symmetric test below preserves results.
  // Only compact pair counts cross back to the host; candidate node IDs remain
  // in the device traversal and are consumed immediately by the exact tests.
  TIME_BEGIN("ConnectivityMap_gpuCandidatePairCount");
  for (auto iNodeList = 0u; iNodeList < numNodeLists; ++iNodeList) {
    const auto sourceNodes = nodeViews[iNodeList];
    const auto nsource = (ghostConnectivity ? sourceNodes.numNodes() : sourceNodes.numInternalNodes());
    for (auto jNodeList = 0u; jNodeList < numNodeLists; ++jNodeList) {
      if (nsource == 0u) continue;

      const auto targetNodes = nodeViews[jNodeList];
      const auto targetTree = treeViews[jNodeList];
      ConnectivityRAJABuffer<size_t> countBuffer(nsource);
      countBuffer.moveToGPU();
      auto* deviceCounts = countBuffer.deviceData();
      const auto firstGhostNodej = targetNodes.firstGhostNode();

      RAJA::forall<EXEC_POLICY>(TRS_UINT(0u, nsource),
        [=] SPHERAL_HOST_DEVICE (const size_t localIndex) {
          const auto i = localIndex;
          const auto& ri = sourceNodes.position(i);
          const auto& Hi = sourceNodes.H(i);
          size_t count = 0u;
          targetTree.forEachCandidate(ri, Hi,
            [=, &count] SPHERAL_HOST_DEVICE (const int j) {
              const auto& rj = targetNodes.position(j);
              const auto& Hj = targetNodes.H(j);
              const auto rij = ri - rj;
              if ((Hi*rij).magnitude2() <= kernelExtent2 or
                  (Hj*rij).magnitude2() <= kernelExtent2) {
                const bool notSelf = (iNodeList != jNodeList or i != size_t(j));
                const bool symmetric = (jNodeList > iNodeList or
                                        (jNodeList == iNodeList and size_t(j) > i) or
                                        (jNodeList < iNodeList and size_t(j) >= firstGhostNodej));
                if (notSelf and symmetric) ++count;
              }
            });
          deviceCounts[localIndex] = count;
        });

      countBuffer.moveToCPU();
      blockCounts[iNodeList*numNodeLists + jNodeList] = std::move(countBuffer.hostData());
    }
  }
  TIME_END("ConnectivityMap_gpuCandidatePairCount");

  std::vector<size_t> blockBases(numBlocks, 0u);
  size_t totalPairs = 0u;
  for (auto block = 0u; block < numBlocks; ++block) {
    blockBases[block] = totalPairs;
    for (const auto count: blockCounts[block]) totalPairs += count;
  }

  ConnectivityRAJABuffer<NodePairIdxType> pairBuffer(totalPairs);
  TIME_BEGIN("ConnectivityMap_pairBufferCopyToGPU");
  pairBuffer.moveToGPU();
  TIME_END("ConnectivityMap_pairBufferCopyToGPU");

  // Second pass: repeat the same device traversal and write each source
  // node's accepted pairs into its pre-sized, non-overlapping range.
  TIME_BEGIN("ConnectivityMap_gpuCandidatePairFill");
  for (auto iNodeList = 0u; iNodeList < numNodeLists; ++iNodeList) {
    const auto sourceNodes = nodeViews[iNodeList];
    const auto nsource = (ghostConnectivity ? sourceNodes.numNodes() : sourceNodes.numInternalNodes());
    for (auto jNodeList = 0u; jNodeList < numNodeLists; ++jNodeList) {
      const auto block = iNodeList*numNodeLists + jNodeList;
      if (nsource == 0u or blockCounts[block].empty()) continue;

      const auto targetNodes = nodeViews[jNodeList];
      const auto targetTree = treeViews[jNodeList];
      const auto& counts = blockCounts[block];
      ConnectivityRAJABuffer<size_t> offsetBuffer(nsource);
      auto& queryOffsets = offsetBuffer.hostData();
      size_t queryOffset = 0u;
      for (auto i = 0u; i < nsource; ++i) {
        queryOffsets[i] = queryOffset;
        queryOffset += counts[i];
      }

      offsetBuffer.moveToGPU();
      const auto* deviceOffsets = offsetBuffer.deviceData();
      auto* devicePairs = pairBuffer.deviceData();
      const auto firstGhostNodej = targetNodes.firstGhostNode();
      const auto pairBase = blockBases[block];

      RAJA::forall<EXEC_POLICY>(TRS_UINT(0u, nsource),
        [=] SPHERAL_HOST_DEVICE (const size_t i) {
          const auto& ri = sourceNodes.position(i);
          const auto& Hi = sourceNodes.H(i);
          size_t writeOffset = deviceOffsets[i];
          targetTree.forEachCandidate(ri, Hi,
            [=, &writeOffset] SPHERAL_HOST_DEVICE (const int j) {
              const auto& rj = targetNodes.position(j);
              const auto& Hj = targetNodes.H(j);
              const auto rij = ri - rj;
              if ((Hi*rij).magnitude2() <= kernelExtent2 or
                  (Hj*rij).magnitude2() <= kernelExtent2) {
                const bool notSelf = (iNodeList != jNodeList or i != size_t(j));
                const bool symmetric = (jNodeList > iNodeList or
                                        (jNodeList == iNodeList and size_t(j) > i) or
                                        (jNodeList < iNodeList and size_t(j) >= firstGhostNodej));
                if (notSelf and symmetric) {
                  devicePairs[pairBase + writeOffset++] = NodePairIdxType(i, iNodeList, j, jNodeList);
                }
              }
            });
        });

    }
  }
  TIME_END("ConnectivityMap_gpuCandidatePairFill");

  TIME_BEGIN("ConnectivityMap_nodeViewCopyToCPU");
  for (auto& nodeView: nodeViews) nodeView.move(chai::CPU);
  TIME_END("ConnectivityMap_nodeViewCopyToCPU");

  TIME_BEGIN("ConnectivityMap_treeViewCopyToCPU");
  for (auto& treeView: treeViews) treeView.move(chai::CPU);
  TIME_END("ConnectivityMap_treeViewCopyToCPU");

  if (totalPairs > 0u) {
    TIME_BEGIN("ConnectivityMap_pairBufferCopyToCPU");
    pairBuffer.moveToCPU();
    TIME_END("ConnectivityMap_pairBufferCopyToCPU");
    nodePairs = std::move(pairBuffer.hostData());
  } else {
    nodePairs.clear();
  }
  return true;
}
#else
template<typename Dimension>
bool
ConnectivityMap_RAJA<Dimension>::
tryBuildNodePairs(const double,
                  const bool,
                  std::vector<NodePairIdxType>&) {
  return false;
}
#endif

}
