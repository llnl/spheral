//---------------------------------Spheral++----------------------------------//
// ConnectivityMap
//
// Stores the full set of significant neighbors for a set of NodeLists.
//
// Created by J. Michael Owen, Sun Oct 30 15:36:33 PST 2005
//----------------------------------------------------------------------------//
#include "ConnectivityMap.hh"
#include "NodeList/NodeList.hh"
#include "Neighbor/Neighbor.hh"
#include "Neighbor/NestedGridNeighbor.hh"
#include "Neighbor/TreeNeighbor.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"
#include "Boundary/Boundary.hh"
#include "Utilities/globalNodeIDs.hh"
#include "Utilities/mortonOrderIndices.hh"
#include "Utilities/PairComparisons.hh"
#include "Utilities/pointDistances.hh"
#include "Utilities/Timer.hh"

#include <algorithm>
#include <ctime>
using std::vector;
using std::map;
using std::string;
using std::pair;

namespace Spheral {

namespace {
//------------------------------------------------------------------------------
// Append v2 to the end of v1
//------------------------------------------------------------------------------
template<typename T>
inline
void
appendSTLvectors(std::vector<T>& v1, std::vector<T>& v2) {
  if (not v2.empty()) {
    v1.reserve(v1.size() + v2.size());
    v1.insert(v1.end(), v2.begin(), v2.end());
  }
}

//------------------------------------------------------------------------------
// How should we compare pairs for sorting?
//------------------------------------------------------------------------------
inline
KeyTraits::Key
hashKeys(const KeyTraits::Key& a, const KeyTraits::Key& b) {
  // Szudzik's function
  // return (a >= b ?
  //         a * a + a + b :
  //         a + b * b);          // where a, b >= 0
  // return a + b;
  return ((KeyTraits::Key(a) << 16) | KeyTraits::Key(b));
}

//------------------------------------------------------------------------------
// Flatten the legacy ragged connectivity into CSR-like storage over
// (global node index, neighbor node list) entries.
//------------------------------------------------------------------------------
inline
void
flattenConnectivity(const std::vector<std::vector<std::vector<int>>>& source,
                    const size_t numNodeLists,
                    std::vector<int>& flatOffsets,
                    std::vector<int>& flatNeighbors) {
  const auto numEntries = source.size()*numNodeLists;
  flatOffsets.resize(numEntries + 1u);
  flatOffsets[0] = 0;
  auto offset = 0;
  auto entry = 0u;
  for (const auto& neighborsPerNode: source) {
    CHECK(neighborsPerNode.size() == numNodeLists);
    for (const auto& neighborsPerNodeList: neighborsPerNode) {
      offset += neighborsPerNodeList.size();
      ++entry;
      flatOffsets[entry] = offset;
    }
  }
  CHECK(entry == numEntries);

  flatNeighbors.clear();
  flatNeighbors.reserve(offset);
  for (const auto& neighborsPerNode: source) {
    for (const auto& neighborsPerNodeList: neighborsPerNode) {
      flatNeighbors.insert(flatNeighbors.end(),
                           neighborsPerNodeList.begin(),
                           neighborsPerNodeList.end());
    }
  }
  CHECK(flatOffsets.back() == int(flatNeighbors.size()));
}

//------------------------------------------------------------------------------
// Convert counts to CSR offsets.
//------------------------------------------------------------------------------
inline
void
countsToOffsets(const std::vector<int>& counts,
                std::vector<int>& offsets) {
  offsets.resize(counts.size() + 1u);
  if (counts.empty()) {
    offsets[0] = 0;
  } else {
    RAJA::exclusive_scan<RAJA::seq_exec>(RAJA::make_span(counts.data(), counts.size()),
                                         RAJA::make_span(offsets.data(), counts.size()),
                                         RAJA::operators::plus<int>{});
    offsets.back() = offsets[counts.size() - 1u] + counts.back();
  }
}

//------------------------------------------------------------------------------
// View-based coarse-neighbor helpers used by ConnectivityMap preprocessing.
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
void
extractTreeCellIndices(const typename TreeNeighbor<Dimension>::CellKey key,
                       typename TreeNeighbor<Dimension>::CellKey& ix,
                       typename TreeNeighbor<Dimension>::CellKey& iy,
                       typename TreeNeighbor<Dimension>::CellKey& iz) {
  using TreeType = TreeNeighbor<Dimension>;
  ix = key & TreeType::xKeyMask();
  iy = (key & TreeType::yKeyMask()) >> TreeType::num1DBits();
  iz = (key & TreeType::zKeyMask()) >> (2*TreeType::num1DBits());
}

template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
typename TreeNeighbor<Dimension>::CellKey
shiftTreeKeyLevel(const typename TreeNeighbor<Dimension>::CellKey ix,
                  const typename TreeNeighbor<Dimension>::LevelKey level0,
                  const typename TreeNeighbor<Dimension>::LevelKey level1) {
  return (level0 > level1 ? (ix >> (level0 - level1)) : (ix << (level1 - level0)));
}

template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
bool
treeKeyInRange(const typename TreeNeighbor<Dimension>::CellKey key,
               const typename TreeNeighbor<Dimension>::CellKey ix_min,
               const typename TreeNeighbor<Dimension>::CellKey iy_min,
               const typename TreeNeighbor<Dimension>::CellKey iz_min,
               const typename TreeNeighbor<Dimension>::CellKey ix_max,
               const typename TreeNeighbor<Dimension>::CellKey iy_max,
               const typename TreeNeighbor<Dimension>::CellKey iz_max) {
  typename TreeNeighbor<Dimension>::CellKey ix, iy, iz;
  extractTreeCellIndices<Dimension>(key, ix, iy, iz);
  return (ix >= ix_min and ix <= ix_max and
          iy >= iy_min and iy <= iy_max and
          iz >= iz_min and iz <= iz_max);
}

template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
size_t
countOrFillTreeGroupCoarseNeighbors(const TreeNeighborView<Dimension>& view,
                                    const int levelID,
                                    const uint64_t cellID,
                                    int* result) {
  using TreeType = TreeNeighbor<Dimension>;
  using CellKey = typename TreeType::CellKey;
  using LevelKey = typename TreeType::LevelKey;
  if (view.empty()) return 0u;
  CHECK(levelID >= 0);
  const auto masterLevel = LevelKey(levelID);
  const auto masterCell = CellKey(cellID);
  CellKey ix_master, iy_master, iz_master;
  extractTreeCellIndices<Dimension>(masterCell, ix_master, iy_master, iz_master);

  auto count = 0u;
  CHECK2(view.levelSize(0) > 0u, "TreeNeighbor root level is empty.");
  const auto rootCellIndex = view.levelBegin(0);
  CHECK2(view.memberSize(rootCellIndex) == 0u,
         "TreeNeighbor root cell occupied!  Will miss neighbors... " << view.memberSize(rootCellIndex));
  for (LevelKey ilevel = 1; ilevel < view.numLevels(); ++ilevel) {
    CellKey ix_min, iy_min, iz_min, ix_max, iy_max, iz_max;
    const auto delta = (ilevel <= masterLevel ? 1U : (1U << (ilevel - masterLevel)));
    const auto ix = shiftTreeKeyLevel<Dimension>(ix_master, masterLevel, ilevel);
    const auto iy = shiftTreeKeyLevel<Dimension>(iy_master, masterLevel, ilevel);
    const auto iz = shiftTreeKeyLevel<Dimension>(iz_master, masterLevel, ilevel);
    ix_min = (ix > delta                    ? ix - delta : 0U);
    iy_min = (iy > delta                    ? iy - delta : 0U);
    iz_min = (iz > delta                    ? iz - delta : 0U);
    ix_max = ((TreeType::max1dKeyValue() - ix) > delta ? ix + 2*delta - 1U : TreeType::max1dKeyValue());
    iy_max = ((TreeType::max1dKeyValue() - iy) > delta ? iy + 2*delta - 1U : TreeType::max1dKeyValue());
    iz_max = ((TreeType::max1dKeyValue() - iz) > delta ? iz + 2*delta - 1U : TreeType::max1dKeyValue());
    CHECK(ix_min <= ix_max and ix_max <= TreeType::max1dKeyValue());
    CHECK(iy_min <= iy_max and iy_max <= TreeType::max1dKeyValue());
    CHECK(iz_min <= iz_max and iz_max <= TreeType::max1dKeyValue());

    for (auto cellIndex = view.levelBegin(ilevel); cellIndex < view.levelEnd(ilevel); ++cellIndex) {
      if (treeKeyInRange<Dimension>(view.cellKey(cellIndex),
                                    ix_min, iy_min, iz_min,
                                    ix_max, iy_max, iz_max)) {
        const auto numMembers = view.memberSize(cellIndex);
        for (auto k = 0u; k < numMembers; ++k) {
          if (result != nullptr) result[count] = view.member(cellIndex, k);
          ++count;
        }
      }
    }
  }
  return count;
}

SPHERAL_HOST_DEVICE
inline
void
unpackNestedGridCellKey(const uint64_t key,
                        int& ix,
                        int& iy,
                        int& iz) {
  constexpr uint64_t mask = (1ull << 21) - 1ull;
  ix = int(key & mask) - (1 << 20);
  iy = int((key >> 21) & mask) - (1 << 20);
  iz = int((key >> 42) & mask) - (1 << 20);
}

SPHERAL_HOST_DEVICE
inline
bool
nestedKeyInRange(const uint64_t key,
                 const int ixMin,
                 const int iyMin,
                 const int izMin,
                 const int ixMax,
                 const int iyMax,
                 const int izMax,
                 const int ndim) {
  int ix, iy, iz;
  unpackNestedGridCellKey(key, ix, iy, iz);
  return (ix >= ixMin and ix <= ixMax and
          (ndim < 2 or (iy >= iyMin and iy <= iyMax)) and
          (ndim < 3 or (iz >= izMin and iz <= izMax)));
}

template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
size_t
countOrFillNestedGroupCoarseNeighbors(const NestedGridNeighborView<Dimension>& view,
                                      const int gridLevel,
                                      const uint64_t packedCellKey,
                                      const int gridCellInfluenceRadius,
                                      const int searchType,
                                      int* result) {
  if (view.empty()) return 0u;

  int cellX, cellY, cellZ;
  unpackNestedGridCellKey(packedCellKey, cellX, cellY, cellZ);
  const auto baseMinX = cellX - gridCellInfluenceRadius;
  const auto baseMinY = cellY - gridCellInfluenceRadius;
  const auto baseMinZ = cellZ - gridCellInfluenceRadius;
  const auto baseMaxX = cellX + gridCellInfluenceRadius;
  const auto baseMaxY = cellY + gridCellInfluenceRadius;
  const auto baseMaxZ = cellZ + gridCellInfluenceRadius;

  auto count = 0u;
  for (int currentGridLevel = 0; currentGridLevel != int(view.numLevels()); ++currentGridLevel) {
    if (view.levelSize(currentGridLevel) > 0u) {
      int ixMin, iyMin, izMin, ixMax, iyMax, izMax;
      if (currentGridLevel > gridLevel) {
        const int glfactor = 1 << (currentGridLevel - gridLevel);
        ixMin = baseMinX * glfactor;
        iyMin = baseMinY * glfactor;
        izMin = baseMinZ * glfactor;
        ixMax = (baseMaxX + 1) * glfactor - 1;
        iyMax = (baseMaxY + 1) * glfactor - 1;
        izMax = (baseMaxZ + 1) * glfactor - 1;
      } else {
        const int glfactor = 1 << (gridLevel - currentGridLevel);
        ixMin = baseMinX / glfactor;
        iyMin = baseMinY / glfactor;
        izMin = baseMinZ / glfactor;
        ixMax = baseMaxX / glfactor;
        iyMax = baseMaxY / glfactor;
        izMax = baseMaxZ / glfactor;
      }
      const int radius = (searchType == int(NeighborSearchType::GatherScatter) ?
                          gridCellInfluenceRadius * (1 << (currentGridLevel > gridLevel ? currentGridLevel - gridLevel : 0)) :
                          searchType == int(NeighborSearchType::Gather) ?
                          gridCellInfluenceRadius / (1 << (gridLevel > currentGridLevel ? gridLevel - currentGridLevel : 0)) + 1 :
                          gridCellInfluenceRadius);
      ixMin -= radius;
      iyMin -= radius;
      izMin -= radius;
      ixMax += radius;
      iyMax += radius;
      izMax += radius;
      for (auto cellIndex = view.levelBegin(currentGridLevel);
           cellIndex < view.levelEnd(currentGridLevel);
           ++cellIndex) {
        if (nestedKeyInRange(view.cellKey(cellIndex),
                             ixMin, iyMin, izMin,
                             ixMax, iyMax, izMax,
                             Dimension::nDim)) {
          const auto numMembers = view.memberSize(cellIndex);
          for (auto k = 0u; k < numMembers; ++k) {
            if (result != nullptr) result[count] = view.member(cellIndex, k);
            ++count;
          }
        }
      }
    }
  }
  return count;
}

//------------------------------------------------------------------------------
// Expand flat CSR-like connectivity back to the legacy ragged storage.
//------------------------------------------------------------------------------
inline
void
unflattenConnectivity(const size_t numEntries,
                      const size_t numNodeLists,
                      const std::vector<int>& flatOffsets,
                      const std::vector<int>& flatNeighbors,
                      std::vector<std::vector<std::vector<int>>>& connectivity) {
  connectivity = std::vector<std::vector<std::vector<int>>>(numEntries,
                                                            std::vector<std::vector<int>>(numNodeLists));
  for (auto globalNodeIndex = 0u; globalNodeIndex < numEntries; ++globalNodeIndex) {
    for (auto neighborNodeList = 0u; neighborNodeList < numNodeLists; ++neighborNodeList) {
      const auto entryIndex = globalNodeIndex*numNodeLists + neighborNodeList;
      connectivity[globalNodeIndex][neighborNodeList].assign(flatNeighbors.begin() + flatOffsets[entryIndex],
                                                             flatNeighbors.begin() + flatOffsets[entryIndex + 1u]);
    }
  }
}

//------------------------------------------------------------------------------
// Copy the flat connectivity slices for a given global node back into a
// temporary ragged representation.
//------------------------------------------------------------------------------
inline
std::vector<std::vector<int>>
flatConnectivityForNode(const ConnectivityMapFlatView& connectivity,
                        const size_t globalNodeIndex,
                        const size_t numNodeLists) {
  std::vector<std::vector<int>> result(numNodeLists);
  for (auto nodeList = 0u; nodeList < numNodeLists; ++nodeList) {
    const auto entryIndex = connectivity.entryIndex(globalNodeIndex, nodeList);
    const auto count = connectivity.size(entryIndex);
    result[nodeList].resize(count);
    for (auto k = 0u; k < count; ++k) result[nodeList][k] = connectivity(entryIndex, k);
  }
  return result;
}

//------------------------------------------------------------------------------
// Shared pair-selection rule for node-pair construction.
//------------------------------------------------------------------------------
SPHERAL_HOST_DEVICE
inline
bool
shouldCalculatePairInteraction(const int nodeListi, const int i,
                               const int nodeListj, const int j,
                               const int firstGhostNodej) {
  return ((nodeListj > nodeListi) or
          (nodeListj == nodeListi and j > i) or
          (nodeListj < nodeListi and j >= firstGhostNodej));
}

//------------------------------------------------------------------------------
// Check if a flat connectivity slice contains a target neighbor.
//------------------------------------------------------------------------------
inline
bool
flatConnectivityContains(const ConnectivityMapFlatView& connectivity,
                         const size_t globalNodeIndex,
                         const size_t neighborNodeList,
                         const int target) {
  const auto entryIndex = connectivity.entryIndex(globalNodeIndex, neighborNodeList);
  const auto count = connectivity.size(entryIndex);
  for (auto k = 0u; k < count; ++k) {
    if (connectivity(entryIndex, k) == target) return true;
  }
  return false;
}

//------------------------------------------------------------------------------
// Count or fill the raw overlap-neighbor candidates for one flat entry.
//------------------------------------------------------------------------------
template<typename PositionViewType, typename HViewType, typename OffsetViewType>
SPHERAL_HOST_DEVICE
inline
size_t
fillRawOverlapNeighborsForEntry(const ConnectivityMapFlatView& connectivity,
                                const OffsetViewType& offsets,
                                const PositionViewType& position,
                                const HViewType& H,
                                const double kernelExtent2,
                                const size_t numNodeLists,
                                const size_t globalNodeIndex,
                                const size_t iNodeList,
                                const int i,
                                const size_t targetNodeList,
                                int* result) {
  size_t count = 0u;

  // All direct gather/scatter neighbors are overlap neighbors.
  const auto baseEntryIndex = connectivity.entryIndex(globalNodeIndex, targetNodeList);
  const auto baseCount = connectivity.size(baseEntryIndex);
  for (auto k = 0u; k < baseCount; ++k) {
    if (result != nullptr) result[count] = connectivity(baseEntryIndex, k);
    ++count;
  }

  const auto& ri = position(iNodeList, i);
  const auto& Hi = H(iNodeList, i);
  for (auto jNodeList = 0u; jNodeList < numNodeLists; ++jNodeList) {
    const auto entryIndexjNodeList = connectivity.entryIndex(globalNodeIndex, jNodeList);
    const auto countjNodeList = connectivity.size(entryIndexjNodeList);
    for (auto neighborIndex = 0u; neighborIndex < countjNodeList; ++neighborIndex) {
      const auto j1 = connectivity(entryIndexjNodeList, neighborIndex);
      const auto& rj1 = position(jNodeList, j1);
      if ((Hi*(rj1 - ri)).magnitude2() <= kernelExtent2) {   // j1 is a gather neighbor of i.
        const auto globalj1 = size_t(offsets[jNodeList] + j1);
        const auto targetEntryIndex = connectivity.entryIndex(globalj1, targetNodeList);
        const auto targetCount = connectivity.size(targetEntryIndex);
        for (auto targetNeighborIndex = 0u; targetNeighborIndex < targetCount; ++targetNeighborIndex) {
          const auto j2 = connectivity(targetEntryIndex, targetNeighborIndex);
          if (targetNodeList != iNodeList or j2 != i) {
            const auto& rj2 = position(targetNodeList, j2);
            const auto& Hj2 = H(targetNodeList, j2);
            if ((Hj2*(rj2 - rj1)).magnitude2() <= kernelExtent2) {
              if (result != nullptr) result[count] = j2;
              ++count;
            }
          }
        }
      }
    }
  }

  return count;
}

//------------------------------------------------------------------------------
// Compare overlap neighbors in the desired final ordering.
//------------------------------------------------------------------------------
template<typename KeyViewType>
SPHERAL_HOST_DEVICE
inline
bool
overlapNeighborLess(const int a,
                    const int b,
                    const size_t targetNodeList,
                    const KeyViewType& keys,
                    const bool useKeys) {
  return (useKeys ? keys(targetNodeList, a) < keys(targetNodeList, b) : a < b);
}

//------------------------------------------------------------------------------
// Insertion sort a small overlap-neighbor slice in place.
//------------------------------------------------------------------------------
template<typename KeyViewType>
SPHERAL_HOST_DEVICE
inline
void
sortOverlapNeighbors(int* values,
                     const size_t count,
                     const size_t targetNodeList,
                     const KeyViewType& keys,
                     const bool useKeys) {
  for (auto i = 1u; i < count; ++i) {
    const auto value = values[i];
    auto j = i;
    while (j > 0u and overlapNeighborLess(value, values[j - 1u], targetNodeList, keys, useKeys)) {
      values[j] = values[j - 1u];
      --j;
    }
    values[j] = value;
  }
}

//------------------------------------------------------------------------------
// Remove duplicates from a sorted overlap-neighbor slice in place.
//------------------------------------------------------------------------------
SPHERAL_HOST_DEVICE
inline
size_t
uniqueOverlapNeighbors(int* values,
                       const size_t count) {
  if (count == 0u) return 0u;
  auto result = 1u;
  for (auto i = 1u; i < count; ++i) {
    if (values[i] != values[result - 1u]) values[result++] = values[i];
  }
  return result;
}

//------------------------------------------------------------------------------
// Local device-safe box predicates for preprocessing culling.
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
bool
pointInBox(const typename Dimension::Vector& point,
           const typename Dimension::Vector& xmin,
           const typename Dimension::Vector& xmax,
           const double tol = 1.0e-10) {
  for (auto k = 0; k < Dimension::nDim; ++k) {
    if (point(k) < xmin(k) - tol) return false;
    if (point(k) > xmax(k) + tol) return false;
  }
  return true;
}

template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
bool
boxesIntersect(const typename Dimension::Vector& xmin1,
               const typename Dimension::Vector& xmax1,
               const typename Dimension::Vector& xmin2,
               const typename Dimension::Vector& xmax2,
               const double tol = 1.0e-10) {
  for (auto k = 0; k < Dimension::nDim; ++k) {
    if (xmax1(k) < xmin2(k) - tol) return false;
    if (xmax2(k) < xmin1(k) - tol) return false;
  }
  return true;
}

template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
bool
keepCoarseNeighborForGroup(const int neighborSearchType,
                           const typename Dimension::Vector& nodePosition,
                           const typename Dimension::Vector& nodeExtent,
                           const typename Dimension::Vector& minMasterPosition,
                           const typename Dimension::Vector& maxMasterPosition,
                           const typename Dimension::Vector& minMasterExtent,
                           const typename Dimension::Vector& maxMasterExtent) {
  if (neighborSearchType == int(NeighborSearchType::GatherScatter)) {
    const auto minNodeExtent = nodePosition - nodeExtent;
    const auto maxNodeExtent = nodePosition + nodeExtent;
    return (pointInBox<Dimension>(nodePosition, minMasterExtent, maxMasterExtent) or
            boxesIntersect<Dimension>(minMasterPosition, maxMasterPosition, minNodeExtent, maxNodeExtent));

  } else if (neighborSearchType == int(NeighborSearchType::Gather)) {
    return pointInBox<Dimension>(nodePosition, minMasterExtent, maxMasterExtent);

  } else {
    const auto minNodeExtent = nodePosition - nodeExtent;
    const auto maxNodeExtent = nodePosition + nodeExtent;
    return boxesIntersect<Dimension>(minMasterPosition, maxMasterPosition, minNodeExtent, maxNodeExtent);
  }
}

//------------------------------------------------------------------------------
// Sort node pairs for domain independent ordering.
//------------------------------------------------------------------------------
template<typename KeyContainer>
inline
void
sortPairs(NodePairList& pairs,
          const KeyContainer& keys) {
  // Start by making sure the pairs themselves are deterministically arranged.
  for (auto& p: pairs) {
    if (keys(p.i_list, p.i_node) > keys(p.j_list, p.j_node)) {
      std::swap(p.i_list, p.j_list);
      std::swap(p.i_node, p.j_node);
    }
  }

  // Now sort the list as a whole.
  std::sort(pairs.begin(), pairs.end(),
            [&](const NodePairIdxType& a, const NodePairIdxType& b) {
              return hashKeys(keys(a.i_list, a.i_node), keys(a.j_list, a.j_node)) <
                hashKeys(keys(b.i_list, b.i_node), keys(b.j_list, b.j_node));
            });
}

}

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
ConnectivityMap<Dimension>::
ConnectivityMap():
  mNodeLists(),
  mBuildGhostConnectivity(false),
  mBuildOverlapConnectivity(false),
  mBuildIntersectionConnectivity(false),
  mOffsets(),
  mConnectivity(),
  mConnectivityFlatOffsets(),
  mConnectivityFlatNeighbors(),
  mConnectivityFlatOffsetsSpan(),
  mConnectivityFlatNeighborsSpan(),
  mNodePairListPtr(),
  mOverlapConnectivity(),
  mOverlapConnectivityFlatOffsets(),
  mOverlapConnectivityFlatNeighbors(),
  mOverlapConnectivityFlatOffsetsSpan(),
  mOverlapConnectivityFlatNeighborsSpan(),
  mNodeTraversalIndices(),
  mKeys(FieldStorageType::CopyFields),
  mCouplingPtr(std::make_shared<NodeCoupling>()),
  mIntersectionConnectivity() {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
ConnectivityMap<Dimension>::
~ConnectivityMap() {
  GPUUtils::freeMAView(mConnectivityFlatOffsetsSpan);
  GPUUtils::freeMAView(mConnectivityFlatNeighborsSpan);
  GPUUtils::freeMAView(mOverlapConnectivityFlatOffsetsSpan);
  GPUUtils::freeMAView(mOverlapConnectivityFlatNeighborsSpan);
}

//------------------------------------------------------------------------------
// Rebuild the flat connectivity views from the legacy ragged storage.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConnectivityMap<Dimension>::
rebuildFlatConnectivityViews() {
  const auto numNodeLists = mNodeLists.size();

  flattenConnectivity(mConnectivity,
                      numNodeLists,
                      mConnectivityFlatOffsets,
                      mConnectivityFlatNeighbors);
  GPUUtils::initMAView(mConnectivityFlatOffsetsSpan, mConnectivityFlatOffsets);
  GPUUtils::initMAView(mConnectivityFlatNeighborsSpan, mConnectivityFlatNeighbors);
  GPUUtils::touch(mConnectivityFlatOffsetsSpan, chai::CPU);
  GPUUtils::touch(mConnectivityFlatNeighborsSpan, chai::CPU);

  if (mBuildOverlapConnectivity) {
    flattenConnectivity(mOverlapConnectivity,
                        numNodeLists,
                        mOverlapConnectivityFlatOffsets,
                        mOverlapConnectivityFlatNeighbors);
    GPUUtils::initMAView(mOverlapConnectivityFlatOffsetsSpan, mOverlapConnectivityFlatOffsets);
    GPUUtils::initMAView(mOverlapConnectivityFlatNeighborsSpan, mOverlapConnectivityFlatNeighbors);
    GPUUtils::touch(mOverlapConnectivityFlatOffsetsSpan, chai::CPU);
    GPUUtils::touch(mOverlapConnectivityFlatNeighborsSpan, chai::CPU);
  } else {
    mOverlapConnectivityFlatOffsets.clear();
    mOverlapConnectivityFlatNeighbors.clear();
    GPUUtils::freeMAView(mOverlapConnectivityFlatOffsetsSpan);
    GPUUtils::freeMAView(mOverlapConnectivityFlatNeighborsSpan);
  }
}

//------------------------------------------------------------------------------
// Internal method to build the connectivity for the requested set of NodeLists.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConnectivityMap<Dimension>::
patchConnectivity(const FieldList<Dimension, size_t>& flags,
                  const FieldList<Dimension, size_t>& old2new) {
  TIME_BEGIN("ConnectivityMap_patch");

  const auto domainDecompIndependent = NodeListRegistrar<Dimension>::instance().domainDecompositionIndependent();

  // We have to recompute the keys to sort nodes by excluding the 
  // nodes that are being removed.
  const auto numNodeLists = mNodeLists.size();
  if (domainDecompIndependent) {
// #pragma omp parallel for collapse(2)
    for (auto iNodeList = 0u; iNodeList < numNodeLists; ++iNodeList) {
      for (auto i = 0u; i < mNodeLists[iNodeList]->numNodes(); ++i) {
        if (flags(iNodeList, i) == 0) mKeys(iNodeList, i) = KeyTraits::maxKey;
      }
    }
  }

  // Iterate over the Connectivity (NodeList).
  for (auto iNodeList = 0u; iNodeList != numNodeLists; ++iNodeList) {
    const auto ioff = mOffsets[iNodeList];
    const auto numNodes = ((domainDecompIndependent or mBuildGhostConnectivity or mBuildOverlapConnectivity) ? 
                           mNodeLists[iNodeList]->numNodes() :
                           mNodeLists[iNodeList]->numInternalNodes());

    vector<size_t> iNodesToKill;
    vector<pair<int, Key>> keys;
#pragma omp parallel
    {
      vector<size_t> iNodesToKill_thread;
      vector<pair<int, Key>> keys_thread;

      // Patch the traversal ordering and connectivity for this NodeList.
#pragma omp for schedule(dynamic)
      for (auto i = 0u; i < numNodes; ++i) {

        // Should we patch this set of neighbors?
        if (flags(iNodeList, i) == 0) {
          iNodesToKill_thread.push_back(i);
        } else {
          if (domainDecompIndependent) keys_thread.push_back(std::make_pair(old2new(iNodeList, i), mKeys(iNodeList, i)));
          mNodeTraversalIndices[iNodeList][i] = old2new(iNodeList, i);
          auto& neighbors = mConnectivity[ioff + i];
          CHECK(neighbors.size() == numNodeLists);
          for (auto jNodeList = 0u; jNodeList < numNodeLists; ++jNodeList) {
            vector<pair<int, Key>> nkeys;
            vector<size_t> jNodesToKill;
            for (auto k = 0u; k < neighbors[jNodeList].size(); ++k) {
              const auto j = neighbors[jNodeList][k];
              if (flags(jNodeList, j) == 0) {
                jNodesToKill.push_back(k);
              } else {
                if (domainDecompIndependent) nkeys.push_back(std::make_pair(old2new(jNodeList, j), mKeys(jNodeList, j)));
                neighbors[jNodeList][k] = old2new(jNodeList, j);
              }
            }
            removeElements(neighbors[jNodeList], jNodesToKill);

            // Recompute the ordering of the neighbors.
            if (domainDecompIndependent) {
              sort(nkeys.begin(), nkeys.end(), ComparePairsBySecondElement<pair<int, Key> >());
              for (size_t k = 0; k != neighbors[jNodeList].size(); ++k) {
                CHECK2(k == 0 or nkeys[k].second > nkeys[k-1].second,
                       "Incorrect neighbor ordering:  "
                       << i << " "
                       << k << " "
                       << nkeys[k-1].second << " "
                       << nkeys[k].second);
                neighbors[jNodeList][k] = nkeys[k].first;
              }
            } else {
              sort(neighbors[jNodeList].begin(), neighbors[jNodeList].end());
            }
          }
        }
      }

#pragma omp critical
      appendSTLvectors(iNodesToKill, iNodesToKill_thread);
      appendSTLvectors(keys, keys_thread);
    }
    removeElements(mNodeTraversalIndices[iNodeList], iNodesToKill);

    // Recompute the ordering for traversing the nodes.
    {
      const auto numNodes = mNodeTraversalIndices[iNodeList].size();
      if (domainDecompIndependent) {
        // keys = vector<pair<int, Key> >();
        // for (size_t k = 0; k != numNodes; ++k) {
        //   const int i = mNodeTraversalIndices[iNodeList][k];
        //   keys.push_back(std::make_pair(i, mKeys(iNodeList, i)));
        // }
        sort(keys.begin(), keys.end(), ComparePairsBySecondElement<pair<int, Key> >());
#pragma omp parallel for
        for (auto k = 0u; k < numNodes; ++k) {
          mNodeTraversalIndices[iNodeList][k] = keys[k].first;
        }
      } else {
#pragma omp parallel for
        for (auto i = 0u; i < numNodes; ++i) {
          mNodeTraversalIndices[iNodeList][i] = i;
        }
      }
    }
  }
  
  // We also need to patch the node pair structure
  // Note here we deliberately reallocate the NodePairList, which will invalidate any
  // PairFields pointing at the original pairs.
  REQUIRE(mNodePairListPtr);
  NodePairList& currentPairs = *mNodePairListPtr;
  std::vector<NodePairIdxType> culledPairs;
  culledPairs.reserve(currentPairs.size());
#pragma omp parallel
  {
    std::vector<NodePairIdxType> culledPairs_thread;
    const auto npairs = currentPairs.size();
    culledPairs_thread.reserve(npairs);
#pragma omp for
    for (auto k = 0u; k < npairs; ++k) {
      const auto iNodeList = currentPairs[k].i_list;
      const auto jNodeList = currentPairs[k].j_list;
      const auto i = currentPairs[k].i_node;
      const auto j = currentPairs[k].j_node;
      if (flags(iNodeList, i) != 0 and flags(jNodeList, j) != 0) {
        culledPairs_thread.push_back(NodePairIdxType(old2new(iNodeList, i), iNodeList,
                                                     old2new(jNodeList, j), jNodeList));
      }
    }
#pragma omp critical
    {
      culledPairs.insert(culledPairs.end(), culledPairs_thread.begin(), culledPairs_thread.end());
    }
  }
  mNodePairListPtr = std::make_shared<NodePairList>(std::move(culledPairs));

  // Sort the NodePairList in order to enforce domain decomposition independence.
  {
    auto& pairs = *mNodePairListPtr;
    if (domainDecompIndependent) {
      // sort(pairs.begin(), pairs.end(), [this](const NodePairIdxType& a, const NodePairIdxType& b) { return (mKeys(a.i_list, a.i_node) + mKeys(a.j_list, a.j_node)) < (mKeys(b.i_list, b.i_node) + mKeys(b.j_list, b.j_node)); });
      // sort(pairs.begin(), pairs.end(), [this](const NodePairIdxType& a, const NodePairIdxType& b) { return hashKeys(mKeys(a.i_list, a.i_node), mKeys(a.j_list, a.j_node)) < hashKeys(mKeys(b.i_list, b.i_node), mKeys(b.j_list, b.j_node)); });
      sortPairs(pairs, mKeys);
    } else {
      std::sort(pairs.begin(), pairs.end());
    }
  }
  // mNodePairListPtr->computeLookup();

  // Patch the intersection lists if we're maintaining them
  if (mBuildIntersectionConnectivity) {
    IntersectionConnectivityContainer intersection;
    for (const auto& element: mIntersectionConnectivity) {
      auto pair = element.first;
      const auto& oldintersect = element.second;
      if (flags(pair.i_list, pair.i_node) != 0 and flags(pair.j_list, pair.j_node) != 0) {
        pair.i_node = old2new(pair.i_list, pair.i_node);
        pair.j_node = old2new(pair.j_list, pair.j_node);
        vector<vector<int>> newintersect(numNodeLists);
        for (auto klist = 0u; klist < numNodeLists; ++klist) {
          for (const auto k: oldintersect[klist]) {
            if (flags(klist, k) != 0) newintersect[klist].push_back(old2new(klist, k));
          }
        }
        intersection[pair] = newintersect;
      }
    }
    mIntersectionConnectivity = std::move(intersection);
  }

  this->rebuildFlatConnectivityViews();

  // You can't check valid yet 'cause the NodeLists have not been resized
  // when we call patch!  The valid method should be checked by whoever called
  // this method after that point.
  //ENSURE(valid());
  TIME_END("ConnectivityMap_patch");
}

//------------------------------------------------------------------------------
// Compute the common neighbors for a pair of nodes.  Note this method 
// returns by value since this information is not stored by ConnectivityMap.
//------------------------------------------------------------------------------
template<typename Dimension>
vector<vector<int> >
ConnectivityMap<Dimension>::
connectivityIntersectionForNodes(const int nodeListi, const int i,
                                 const int nodeListj, const int j,
                                 const FieldList<Dimension, typename Dimension::Vector>& position) const {

  // Pre-conditions.
  TIME_BEGIN("ConnectivityMap_computeIntersectionConnectivity");
  const auto numNodeLists = mNodeLists.size();
  const auto domainDecompIndependent = NodeListRegistrar<Dimension>::instance().domainDecompositionIndependent();
  const auto ghostConnectivity = (mBuildGhostConnectivity or
                                  mBuildOverlapConnectivity or
                                  domainDecompIndependent);
  const auto firstGhostNodei = mNodeLists[nodeListi]->firstGhostNode();
  const auto firstGhostNodej = mNodeLists[nodeListj]->firstGhostNode();
  const bool usePosition = (position.numFields() == numNodeLists);
  REQUIRE(nodeListi < (int)numNodeLists and
          nodeListj < (int)numNodeLists);
  REQUIRE(ghostConnectivity or i < (int)firstGhostNodei or j < (int)firstGhostNodej);
  REQUIRE(position.numFields() == numNodeLists or position.numFields() == 0);

  // Prepare the result.
  vector<vector<int>> result(numNodeLists);
  const auto connectivity = this->connectivityFlatView();
  const auto globalNodei = size_t(mOffsets[nodeListi] + i);
  const auto globalNodej = size_t(mOffsets[nodeListj] + j);

  // If both nodes are internal, we simply intersect their neighbor lists.
  if (ghostConnectivity or (i < (int)firstGhostNodei and j < (int)firstGhostNodej)) {
    const auto neighborsi = flatConnectivityForNode(connectivity, globalNodei, numNodeLists);
    const auto neighborsj = flatConnectivityForNode(connectivity, globalNodej, numNodeLists);
    CHECK(neighborsi.size() == numNodeLists);
    CHECK(neighborsj.size() == numNodeLists);
    vector<int> neighborsijk;
    Vector posi, posj, b;
    if (usePosition) {
      posi = position(nodeListi, i);
      posj = position(nodeListj, j);
    }
    for (auto klist = 0u; klist < numNodeLists; ++klist) {
      neighborsijk.clear();
      if (domainDecompIndependent) {
        std::set_intersection(neighborsi[klist].begin(), neighborsi[klist].end(),
                              neighborsj[klist].begin(), neighborsj[klist].end(),
                              std::back_inserter(neighborsijk),
                              [&](const int a, const int b) { return mKeys(klist, a) < mKeys(klist, b); });
      } else {
        std::set_intersection(neighborsi[klist].begin(), neighborsi[klist].end(),
                              neighborsj[klist].begin(), neighborsj[klist].end(),
                              std::back_inserter(neighborsijk));
      }
      if (usePosition) {
        std::copy_if(neighborsijk.begin(), neighborsijk.end(), std::back_inserter(result[klist]),
                     [&](int k) { return (closestPointOnSegment(position(klist, k), posi, posj, b)); });
      } else {
        result[klist] = neighborsijk;
      }
    }
  } else if (i < (int)firstGhostNodei) {
    result = flatConnectivityForNode(connectivity, globalNodei, numNodeLists);
  } else {
    result = flatConnectivityForNode(connectivity, globalNodej, numNodeLists);
  }
  result[nodeListi].push_back(i);
  result[nodeListj].push_back(j);

  // That's it.
  TIME_END("ConnectivityMap_computeIntersectionConnectivity");
  return result;
}

//------------------------------------------------------------------------------
// Remove connectivity between neighbors.
// NOTE: this method assumes you are passing the indices of the neighbors to
//       remove!
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConnectivityMap<Dimension>::
removeConnectivity(const FieldList<Dimension, vector<vector<int>>>& neighborsToCut) {
  TIME_BEGIN("ConnectivityMap_cutConnectivity");

  const auto numNodeLists = mNodeLists.size();
  REQUIRE(neighborsToCut.numFields() == numNodeLists);

  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto n = mNodeLists[nodeListi]->numNodes();
    for (auto i = 0u; i < n; ++i) {
      const auto& allneighbors = neighborsToCut(nodeListi, i);
      CHECK(allneighbors.size() == 0 or allneighbors.size() == numNodeLists);
      for (auto nodeListj = 0u; nodeListj < allneighbors.size(); ++nodeListj) {
        auto& neighborsi = mConnectivity[mOffsets[nodeListi] + i][nodeListj];
        removeElements(neighborsi, allneighbors[nodeListj]);
      }
    }
  }

  this->rebuildFlatConnectivityViews();

  TIME_END("ConnectivityMap_cutConnectivity");
}

//------------------------------------------------------------------------------
// Compute the union of neighbors for a pair of nodes.  Note this method 
// returns by value since this information is not stored by ConnectivityMap.
//------------------------------------------------------------------------------
template<typename Dimension>
vector<vector<int> >
ConnectivityMap<Dimension>::
connectivityUnionForNodes(const int nodeListi, const int i,
                          const int nodeListj, const int j) const {

  // Pre-conditions.
  const unsigned numNodeLists = mNodeLists.size();
  const auto domainDecompIndependent = NodeListRegistrar<Dimension>::instance().domainDecompositionIndependent();
  const auto ghostConnectivity = (mBuildGhostConnectivity or
                                  mBuildOverlapConnectivity or
                                  domainDecompIndependent);
  const unsigned firstGhostNodei = mNodeLists[nodeListi]->firstGhostNode();
  const unsigned firstGhostNodej = mNodeLists[nodeListj]->firstGhostNode();
  CONTRACT_VAR(ghostConnectivity);
  CONTRACT_VAR(firstGhostNodei);
  CONTRACT_VAR(firstGhostNodej);
  REQUIRE(nodeListi < (int)numNodeLists and
          nodeListj < (int)numNodeLists);
  REQUIRE(ghostConnectivity or i < (int)firstGhostNodei or j < (int)firstGhostNodej);

  // Do the deed.
  vector<vector<int> > result(numNodeLists);
  const auto connectivity = this->connectivityFlatView();
  vector<vector<int> > neighborsi = flatConnectivityForNode(connectivity, mOffsets[nodeListi] + i, numNodeLists);
  vector<vector<int> > neighborsj = flatConnectivityForNode(connectivity, mOffsets[nodeListj] + j, numNodeLists);
  CHECK(neighborsi.size() == numNodeLists);
  CHECK(neighborsj.size() == numNodeLists);
  for (unsigned k = 0; k != numNodeLists; ++k) {
    sort(neighborsi[k].begin(), neighborsi[k].end());
    sort(neighborsj[k].begin(), neighborsj[k].end());
    set_union(neighborsi[k].begin(), neighborsi[k].end(),
              neighborsj[k].begin(), neighborsj[k].end(),
              back_inserter(result[k]));
  }

  // That's it.
  return result;
}

//------------------------------------------------------------------------------
// Return the connectivity in terms of global node IDs.
//------------------------------------------------------------------------------
template<typename Dimension>
map<int, vector<int> > 
ConnectivityMap<Dimension>::
globalConnectivity(vector<Boundary<Dimension>*>& boundaries) const {

  // Get the set of global node IDs.
  auto globalIDs = globalNodeIDs<Dimension, typename vector<const NodeList<Dimension>*>::const_iterator>
    (mNodeLists.begin(), mNodeLists.end());

  // Make sure all ghost nodes have the appropriate global IDs.
  for (typename vector<Boundary<Dimension>*>::iterator boundItr = boundaries.begin();
       boundItr != boundaries.end();
       ++boundItr) (*boundItr)->applyFieldListGhostBoundary(globalIDs);
  for (typename vector<Boundary<Dimension>*>::iterator boundItr = boundaries.begin();
       boundItr != boundaries.end();
       ++boundItr) (*boundItr)->finalizeGhostBoundary();

  // Now convert our connectivity to global IDs.
  map<int, vector<int> > result;
  const size_t numNodeLists = mNodeLists.size();
  const auto connectivity = this->connectivityFlatView();
  for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {

    const NodeList<Dimension>* nodeListPtr = mNodeLists[nodeListi];
    for (auto i = 0u; i != nodeListPtr->numInternalNodes(); ++i) {
      const int gid = globalIDs(nodeListi, i);
      result[gid] = vector<int>();
      const auto globalNodeIndex = size_t(mOffsets[nodeListi] + i);
      for (size_t nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
        const auto entryIndex = connectivity.entryIndex(globalNodeIndex, nodeListj);
        const auto count = connectivity.size(entryIndex);
        for (auto k = 0u; k < count; ++k) {
          result[gid].push_back(globalIDs(nodeListj, connectivity(entryIndex, k)));
        }
      }
      ENSURE(result[gid].size() == numNeighborsForNode(nodeListPtr, i));
    }
  }

  // That's it.
  return result;
}

//------------------------------------------------------------------------------
// Compute the index for the given NodeList in our known set.
//------------------------------------------------------------------------------
template<typename Dimension>
unsigned
ConnectivityMap<Dimension>::
nodeListIndex(const NodeList<Dimension>* nodeListPtr) const {
  return distance(mNodeLists.begin(), 
                  find(mNodeLists.begin(), mNodeLists.end(), nodeListPtr));
}

//------------------------------------------------------------------------------
// Valid test.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
ConnectivityMap<Dimension>::
valid() const {
  TIME_BEGIN("ConnectivityMap_valid");

  const auto domainDecompIndependent = NodeListRegistrar<Dimension>::instance().domainDecompositionIndependent();
  const auto ghostConnectivity = (mBuildGhostConnectivity or
                                  mBuildOverlapConnectivity or
                                  domainDecompIndependent);
  const auto connectivity = this->connectivityFlatView();

  // Check the offsets.
  const auto numNodeLists = mNodeLists.size();
  if (mOffsets.size() != numNodeLists) {
    cerr << "ConnectivityMap::valid: Failed mOffsets.size() == numNodeLists" << endl;
    return false;
  }
  {
    const auto numNodes = (ghostConnectivity ? 
                           mNodeLists.back()->numNodes() : 
                           mNodeLists.back()->numInternalNodes());
    if (mConnectivity.size() != mOffsets.back() + numNodes) {
      cerr << "ConnectivityMap::valid: Failed offset bounding: " << mConnectivity.size() << " != " << mOffsets.back() << " + " << numNodes << endl;
    }
    if (connectivity.numNodeLists() != numNodeLists or
        connectivity.numNodes() != mOffsets.back() + numNodes) {
      cerr << "ConnectivityMap::valid: Flat connectivity dimensions are inconsistent" << endl;
      return false;
    }
  }

  // Make sure that the NodeLists are listed in the correct sequence, and are
  // FluidNodeLists.
  {
    const NodeListRegistrar<Dimension>& registrar = NodeListRegistrar<Dimension>::instance();
    const vector<string> names = registrar.registeredNames();
    int lastPosition = -1;
    for (typename vector<const NodeList<Dimension>*>::const_iterator itr = mNodeLists.begin();
         itr != mNodeLists.end();
         ++itr) {
      const int newPosition = distance(names.begin(),
                                       find(names.begin(), names.end(), (*itr)->name()));
      if (newPosition <= lastPosition) {
        cerr << "ConnectivityMap::valid: Failed ordering of NodeLists" << endl;
        return false;
      }
      lastPosition = newPosition;
    }
  }

  // Iterate over each NodeList entered.
  for (auto nodeListIDi = 0u; nodeListIDi < numNodeLists; ++nodeListIDi) {

    // Are all internal nodes represented?
    const NodeList<Dimension>* nodeListPtri = mNodeLists[nodeListIDi];
    const auto numNodes = (ghostConnectivity ?
                           nodeListPtri->numNodes() : 
                           nodeListPtri->numInternalNodes());
    //const int firstGhostNodei = nodeListPtri->firstGhostNode();
    if (((nodeListIDi < numNodeLists - 1u) and ((mOffsets[nodeListIDi + 1] - mOffsets[nodeListIDi]) != (int)numNodes)) or
        ((nodeListIDi == numNodeLists - 1u) and ((mConnectivity.size() - (size_t)mOffsets[nodeListIDi]) != numNodes))) {
      cerr << "ConnectivityMap::valid: Failed test that all nodes set for NodeList "
           << mNodeLists[nodeListIDi]->name()
           << endl;
      return false;
    }

    // Iterate over the nodes for this NodeList.
    const int ioff = mOffsets[nodeListIDi];
    for (auto i = 0u; i < numNodes; ++i) {
      const auto globalNodeIndex = size_t(ioff + i);

      // Iterate over the sets of NodeList neighbors for this node.
      for (auto nodeListIDj = 0u; nodeListIDj < numNodeLists; ++nodeListIDj) {
        const NodeList<Dimension>* nodeListPtrj = mNodeLists[nodeListIDj];
        //const int firstGhostNodej = nodeListPtrj->firstGhostNode();
        const auto entryIndex = connectivity.entryIndex(globalNodeIndex, nodeListIDj);
        const auto numNeighbors = connectivity.size(entryIndex);

        // We require that the node IDs be sorted, unique, and of course in a valid range.
        if (numNeighbors > 0u) {
          auto minNeighbor = connectivity(entryIndex, 0u);
          auto maxNeighbor = minNeighbor;
          for (auto k = 1u; k < numNeighbors; ++k) {
            const auto neighbor = connectivity(entryIndex, k);
            minNeighbor = std::min(minNeighbor, neighbor);
            maxNeighbor = std::max(maxNeighbor, neighbor);
          }

          if (minNeighbor < 0 or (size_t)maxNeighbor >= nodeListPtrj->numNodes()) {
            cerr << "ConnectivityMap::valid: Failed test that neighbors must be valid IDs: " << minNeighbor << " " << maxNeighbor << " " << nodeListPtrj->numNodes() << endl;
            return false;
          }

          // // When enforcing domain independence the ith node may be a ghost, but all of it's neighbors should
          // // be internal.
          // if (domainDecompIndependent and (i >= firstGhostNodei) and (maxNeighbor > firstGhostNodej)) {
          //   cerr << "ConnectivityMap::valid: Failed test that all neighbors of a ghost node should be internal." << endl;
          //   return false;
          // }

          for (auto k = 1u; k < numNeighbors; ++k) {
            if (domainDecompIndependent) {
              // In the case of domain decomposition reproducibility, neighbors are sorted
              // by hashed IDs.
              if (mKeys(nodeListIDj, connectivity(entryIndex, k)) <
                  mKeys(nodeListIDj, connectivity(entryIndex, k - 1u))) {
                cerr << "ConnectivityMap::valid: Failed test that neighbors must be sorted for node "
                     << i << endl;
                for (auto kk = 0u; kk < numNeighbors; ++kk) {
                  const auto neighbor = connectivity(entryIndex, kk);
                  cerr << "(" << neighbor << " " << mKeys(nodeListIDj, neighbor) << ") ";
                }
                cerr << endl;
                return false;
              }

            } else {
              // Otherwise they should be sorted by local ID.
              if (connectivity(entryIndex, k) <= connectivity(entryIndex, k - 1u)) {
                cerr << "ConnectivityMap::valid: Failed test that neighbors must be sorted" << endl;
                for (auto kk = 0u; kk < numNeighbors; ++kk) cerr << " " << connectivity(entryIndex, kk);
                cerr << endl;
                return false;
              }
            }
          }
        }

        // Check that the connectivity is symmetric.
        for (auto k = 0u; k < numNeighbors; ++k) {
          const auto j = connectivity(entryIndex, k);
          if (ghostConnectivity or ((size_t)j < nodeListPtrj->numInternalNodes())) {
            if (not flatConnectivityContains(connectivity, size_t(mOffsets[nodeListIDj] + j), nodeListIDi, i)) {
              const auto otherEntryIndex = connectivity.entryIndex(size_t(mOffsets[nodeListIDj] + j), nodeListIDi);
              cerr << "ConnectivityMap::valid: Failed test that neighbors must be symmetric: " 
                   << i << " <> " << j 
                   << "  numneigbors(i)=" << numNeighbors
                   << "  numneigbors(j)=" << connectivity.size(otherEntryIndex)
                   << endl;
              cerr << "   " << i << " : ";
              for (auto kk = 0u; kk < numNeighbors; ++kk) cerr << connectivity(entryIndex, kk) << " ";
              cerr << endl
                   << "   " << j << " : ";
              for (auto kk = 0u; kk < connectivity.size(otherEntryIndex); ++kk) {
                cerr << connectivity(otherEntryIndex, kk) << " ";
              }
              cerr << endl;
              return false;
            }
          }
        }

      }
    }
  }

  // Check that if we are using domain decompostion independence then the keys 
  // have been calculated.
  if (domainDecompIndependent) {
    for (typename vector<const NodeList<Dimension>*>::const_iterator itr = mNodeLists.begin();
         itr != mNodeLists.end();
         ++itr) {
      if (not mKeys.haveNodeList(**itr)) {
        cerr << "ConnectivityMap::valid: missing information from Keys." << endl;
        return false;
      }
    }
  }

  // Make sure all nodes are listed in the node index traversal stuff.
  if (mNodeTraversalIndices.size() != mNodeLists.size()) {
    cerr << "ConnectivityMap::valid: mNodeTraversalIndices wrong size!" << endl;
    return false;
  }
  for (auto nodeList = 0u; nodeList < numNodeLists; ++nodeList) {
    const auto numExpected = domainDecompIndependent ? mNodeLists[nodeList]->numNodes() : mNodeLists[nodeList]->numInternalNodes();
    bool ok = mNodeTraversalIndices[nodeList].size() == numExpected;
    for (auto i = 0u; i < numExpected; ++i) {
      ok = ok and (count(mNodeTraversalIndices[nodeList].begin(),
                         mNodeTraversalIndices[nodeList].end(),
                         i) == 1);
    }
    if (not ok) {
      cerr << "ConnectivityMap::valid: mNodeTraversalIndices elements messed up!" << endl;
      return false;
    }
  }

  // // Check that the node traversal is ordered correctly.
  // for (int nodeList = 0; nodeList != numNodeLists; ++nodeList) {
  //   if ((domainDecompIndependent and mNodeLists[nodeList]->numNodes() > 0) or
  //       (not domainDecompIndependent and mNodeLists[nodeList]->numInternalNodes() > 0)) {
  //     const int firstGhostNode = mNodeLists[nodeList]->firstGhostNode();
  //     for (const_iterator itr = begin(nodeList); itr < end(nodeList) - 1; ++itr) {
  //       if (not calculatePairInteraction(nodeList, *itr,
  //                                        nodeList, *(itr + 1), 
  //                                        firstGhostNode)) {
  //         cerr << "ConnectivityMap::valid: mNodeTraversalIndices ordered incorrectly." << endl;
  //         cerr << *itr << " "
  //              << *(itr + 1) << " "
  //              << mKeys(nodeList, *itr) << " "
  //              << mKeys(nodeList, *(itr + 1)) << " "
  //              << mNodeLists[nodeList]->positions()(*itr) << " "
  //              << mNodeLists[nodeList]->positions()(*(itr + 1)) << " "
  //              << endl;
  //         for (int i = 0; i != 100; ++i) cerr << mKeys(nodeList, i) << " " << mNodeLists[nodeList]->positions()(i) << " ";
  //         cerr << endl;
  //         return false;
  //       }
  //     }
  //   }
  // }

  // Everything must be OK.
  TIME_END("ConnectivityMap_valid");
  return true;
}

//------------------------------------------------------------------------------
// Internal method to build the connectivity for the requested set of NodeLists.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConnectivityMap<Dimension>::
computeConnectivity() {
  TIME_BEGIN("ConnectivityMap_computeConnectivity");

  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;

  // Pre-conditions.
  BEGIN_CONTRACT_SCOPE
  {
    for (typename vector<const NodeList<Dimension>*>::const_iterator itr = mNodeLists.begin();
         itr != mNodeLists.end();
         ++itr) {
      REQUIRE((**itr).neighbor().valid());
    }
    REQUIRE(mOffsets.size() == mNodeLists.size());
  }
  END_CONTRACT_SCOPE

  const bool domainDecompIndependent = NodeListRegistrar<Dimension>::instance().domainDecompositionIndependent();
  // std::clock_t tpre = std::clock();

  // Do we need to build the ghost connectivity as well?
  const auto ghostConnectivity = (mBuildGhostConnectivity or
                                  mBuildOverlapConnectivity or
                                  domainDecompIndependent);

  // Build ourselves a temporary DataBase with the set of NodeLists.
  // Simultaneously find the maximum kernel extent.
  DataBase<Dimension> dataBase;
  double kernelExtent = 0.0;
  for (typename vector<const NodeList<Dimension>*>::const_iterator itr = mNodeLists.begin();
       itr != mNodeLists.end();
       ++itr) {
    dataBase.appendNodeList(const_cast<NodeList<Dimension>&>(**itr));
    kernelExtent = max(kernelExtent, (**itr).neighbor().kernelExtent());
  }
  const double kernelExtent2 = kernelExtent*kernelExtent;

  // Erase any prior information.
  const unsigned numNodeLists = dataBase.numNodeLists(),
             connectivitySize = mOffsets.back() + (ghostConnectivity ?
                                                   mNodeLists.back()->numNodes() :
                                                   mNodeLists.back()->numInternalNodes());
  const auto numConnectivityEntries = connectivitySize*numNodeLists;
  mNodeTraversalIndices = vector<vector<int> >(numNodeLists);
  mIntersectionConnectivity.clear();

  // If we're trying to be domain decomposition independent, we need a key to sort
  // by that will give us a unique ordering regardless of position.  The Morton ordered
  // hash fills the bill.
  using Key = typename KeyTraits::Key;
  if (domainDecompIndependent) mKeys = mortonOrderIndices(dataBase);

  // Fill in the ordering for walking the nodes.
  CHECK(mNodeTraversalIndices.size() == numNodeLists);
  if (domainDecompIndependent) {
    for (auto iNodeList = 0u; iNodeList != numNodeLists; ++iNodeList) {
      const NodeList<Dimension>& nodeList = *mNodeLists[iNodeList];
      mNodeTraversalIndices[iNodeList].resize(nodeList.numNodes());
      vector<pair<int, Key> > keys;
      keys.reserve(nodeList.numNodes());
      for (auto i = 0u; i != nodeList.numNodes(); ++i) keys.push_back(pair<int, Key>(i, mKeys(iNodeList, i)));
      sort(keys.begin(), keys.end(), ComparePairsBySecondElement<pair<int, Key> >());
      for (auto i = 0u; i != nodeList.numNodes(); ++i) mNodeTraversalIndices[iNodeList][i] = keys[i].first;
      CHECK(mNodeTraversalIndices[iNodeList].size() == nodeList.numNodes());
      // std::cerr << "Traversal: ";
      // std::copy(mNodeTraversalIndices[iNodeList].begin(), mNodeTraversalIndices[iNodeList].end(), std::ostream_iterator<int>(std::cerr, " "));
      // std::cerr << std::endl;
    }
  } else {
    for (auto iNodeList = 0u; iNodeList != numNodeLists; ++iNodeList) {
      const NodeList<Dimension>& nodeList = *mNodeLists[iNodeList];
      mNodeTraversalIndices[iNodeList].resize(nodeList.numInternalNodes());
      for (auto i = 0u; i != nodeList.numInternalNodes(); ++i) mNodeTraversalIndices[iNodeList][i] = i;
    }
  }

  // Get the position and H fields.
  FieldList<Dimension, Vector> position = dataBase.globalPosition();
  FieldList<Dimension, SymTensor> H = dataBase.globalHfield();
  FieldList<Dimension, Vector> extent = dataBase.globalNodeExtent();
  auto positionView = position.view();
  auto HView = H.view();
  auto extentView = extent.view();
  auto keysView = mKeys.view();
  using FlatIntSpan = typename ConnectivityMap<Dimension>::FlatConnectivitySpan;
#ifdef SPHERAL_UNIFIED_MEMORY
  using FlatVectorSpan = SPHERAL_SPAN_TYPE<Vector>;
  using FlatDescriptorSpan = SPHERAL_SPAN_TYPE<NeighborGroupDescriptor>;
  using FlatTreeNeighborViewSpan = SPHERAL_SPAN_TYPE<TreeNeighborView<Dimension>>;
  using FlatNestedNeighborViewSpan = SPHERAL_SPAN_TYPE<NestedGridNeighborView<Dimension>>;
#else
  using FlatVectorSpan = chai::ManagedArray<Vector>;
  using FlatDescriptorSpan = chai::ManagedArray<NeighborGroupDescriptor>;
  using FlatTreeNeighborViewSpan = chai::ManagedArray<TreeNeighborView<Dimension>>;
  using FlatNestedNeighborViewSpan = chai::ManagedArray<NestedGridNeighborView<Dimension>>;
#endif
  FlatIntSpan offsetsSpan;
  GPUUtils::initMAView(offsetsSpan, mOffsets);
  GPUUtils::touch(offsetsSpan, chai::CPU);
  std::vector<int> firstGhostNodes(numNodeLists);
  for (auto nodeList = 0u; nodeList < numNodeLists; ++nodeList) {
    firstGhostNodes[nodeList] = mNodeLists[nodeList]->firstGhostNode();
  }
  FlatIntSpan firstGhostNodesSpan;
  GPUUtils::initMAView(firstGhostNodesSpan, firstGhostNodes);
  GPUUtils::touch(firstGhostNodesSpan, chai::CPU);
  std::vector<int> neighborSearchTypes(numNodeLists);
  for (auto nodeList = 0u; nodeList < numNodeLists; ++nodeList) {
    neighborSearchTypes[nodeList] = int(mNodeLists[nodeList]->neighbor().neighborSearchType());
  }
  FlatIntSpan neighborSearchTypesSpan;
  GPUUtils::initMAView(neighborSearchTypesSpan, neighborSearchTypes);
  GPUUtils::touch(neighborSearchTypesSpan, chai::CPU);
  std::vector<int> neighborGroupKinds(numNodeLists, FallbackNeighborGroupKind);
  std::vector<int> nestedGridCellInfluenceRadii(numNodeLists, 0);
  std::vector<TreeNeighborView<Dimension>> treeNeighborViews(numNodeLists);
  std::vector<NestedGridNeighborView<Dimension>> nestedNeighborViews(numNodeLists);
  for (auto nodeList = 0u; nodeList < numNodeLists; ++nodeList) {
    const auto& neighbor = mNodeLists[nodeList]->neighbor();
    neighborGroupKinds[nodeList] = neighbor.groupKind();
    if (neighborGroupKinds[nodeList] == TreeNeighborGroupKind) {
      const auto& treeNeighbor = static_cast<const TreeNeighbor<Dimension>&>(neighbor);
      treeNeighborViews[nodeList] = treeNeighbor.view();
    } else if (neighborGroupKinds[nodeList] == NestedNeighborGroupKind) {
      const auto& nestedNeighbor = static_cast<const NestedGridNeighbor<Dimension>&>(neighbor);
      nestedNeighborViews[nodeList] = nestedNeighbor.view();
      nestedGridCellInfluenceRadii[nodeList] = nestedNeighbor.gridCellInfluenceRadius();
    }
  }
  FlatIntSpan neighborGroupKindsSpan;
  GPUUtils::initMAView(neighborGroupKindsSpan, neighborGroupKinds);
  GPUUtils::touch(neighborGroupKindsSpan, chai::CPU);
  FlatIntSpan nestedGridCellInfluenceRadiiSpan;
  GPUUtils::initMAView(nestedGridCellInfluenceRadiiSpan, nestedGridCellInfluenceRadii);
  GPUUtils::touch(nestedGridCellInfluenceRadiiSpan, chai::CPU);
  FlatTreeNeighborViewSpan treeNeighborViewsSpan;
  GPUUtils::initMAView(treeNeighborViewsSpan, treeNeighborViews);
  GPUUtils::touch(treeNeighborViewsSpan, chai::CPU);
  FlatNestedNeighborViewSpan nestedNeighborViewsSpan;
  GPUUtils::initMAView(nestedNeighborViewsSpan, nestedNeighborViews);
  GPUUtils::touch(nestedNeighborViewsSpan, chai::CPU);

  // Group nodes by identical per-NodeList descriptors of their own query
  // state. This replaces the old sequential "first seed claims a group"
  // ownership discovery with a canonical key-based grouping.
  std::vector<int> seedNodeListIDs;
  std::vector<int> seedNodeIDs;
  seedNodeListIDs.reserve(connectivitySize);
  seedNodeIDs.reserve(connectivitySize);
  for (auto iiNodeList = 0u; iiNodeList < numNodeLists; ++iiNodeList) {
    const auto nii = (ghostConnectivity ?
                      mNodeLists[iiNodeList]->numNodes() :
                      mNodeLists[iiNodeList]->numInternalNodes());
    for (auto ii = 0u; ii < nii; ++ii) {
      seedNodeListIDs.push_back(iiNodeList);
      seedNodeIDs.push_back(ii);
    }
  }
  CHECK(seedNodeIDs.size() == connectivitySize);
  FlatIntSpan seedNodeListIDsSpan;
  GPUUtils::initMAView(seedNodeListIDsSpan, seedNodeListIDs);
  GPUUtils::touch(seedNodeListIDsSpan, chai::CPU);
  FlatIntSpan seedNodeIDsSpan;
  GPUUtils::initMAView(seedNodeIDsSpan, seedNodeIDs);
  GPUUtils::touch(seedNodeIDsSpan, chai::CPU);

  std::vector<NeighborGroupDescriptor> seedDescriptors(connectivitySize*numNodeLists);
  FlatDescriptorSpan seedDescriptorsSpan;
  GPUUtils::initMAView(seedDescriptorsSpan, seedDescriptors);
  GPUUtils::touch(seedDescriptorsSpan, chai::CPU);
  RAJA::forall<RAJA::seq_exec>(RAJA::RangeSegment(0, int(connectivitySize*numNodeLists)),
  [&](int entryIndex) {
    const auto seedIndex = size_t(entryIndex) / numNodeLists;
    const auto nodeListi = size_t(entryIndex) % numNodeLists;
    const auto iiNodeList = size_t(seedNodeListIDsSpan[seedIndex]);
    const auto ii = seedNodeIDsSpan[seedIndex];
    seedDescriptorsSpan[entryIndex] =
      mNodeLists[nodeListi]->neighbor().groupDescriptor(position(iiNodeList, ii),
                                                        H(iiNodeList, ii),
                                                        seedIndex);
  });
  GPUUtils::touch(seedDescriptorsSpan, chai::CPU);

  std::vector<int> sortedSeedIndices(connectivitySize);
  for (auto i = 0u; i < connectivitySize; ++i) sortedSeedIndices[i] = i;
  std::sort(sortedSeedIndices.begin(), sortedSeedIndices.end(),
            [&](const int lhs, const int rhs) {
              for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
                const auto& ldesc = seedDescriptorsSpan[size_t(lhs)*numNodeLists + nodeListi];
                const auto& rdesc = seedDescriptorsSpan[size_t(rhs)*numNodeLists + nodeListi];
                if (ldesc < rdesc) return true;
                if (rdesc < ldesc) return false;
              }
              return lhs < rhs;
            });

  std::vector<int> descriptorRunOffsets(1, 0);
  for (auto sortedIndex = 1u; sortedIndex < connectivitySize; ++sortedIndex) {
    const auto seedIndex = sortedSeedIndices[sortedIndex];
    const auto repSeed = sortedSeedIndices[sortedIndex - 1u];
    auto same = true;
    for (auto nodeListi = 0u; nodeListi < numNodeLists and same; ++nodeListi) {
      same = (seedDescriptorsSpan[size_t(seedIndex)*numNodeLists + nodeListi] ==
              seedDescriptorsSpan[size_t(repSeed)*numNodeLists + nodeListi]);
    }
    if (not same) descriptorRunOffsets.push_back(sortedIndex);
  }
  descriptorRunOffsets.push_back(connectivitySize);

  const auto numNeighborGroups = descriptorRunOffsets.size() - 1u;
  std::vector<int> groupRepresentativeSeeds(numNeighborGroups);
  std::vector<int> groupSeedCounts(numNeighborGroups);
  for (auto groupIndex = 0u; groupIndex < numNeighborGroups; ++groupIndex) {
    const auto begin = descriptorRunOffsets[groupIndex];
    const auto end = descriptorRunOffsets[groupIndex + 1u];
    groupRepresentativeSeeds[groupIndex] = sortedSeedIndices[begin];
    groupSeedCounts[groupIndex] = end - begin;
  }

  std::vector<int> groupOrder(numNeighborGroups);
  for (auto i = 0u; i < numNeighborGroups; ++i) groupOrder[i] = i;
  std::sort(groupOrder.begin(), groupOrder.end(),
            [&](const int lhs, const int rhs) {
              return groupRepresentativeSeeds[lhs] < groupRepresentativeSeeds[rhs];
            });

  std::vector<int> orderedGroupRepresentativeSeeds(numNeighborGroups);
  std::vector<int> orderedGroupSeedCounts(numNeighborGroups);
  for (auto groupIndex = 0u; groupIndex < numNeighborGroups; ++groupIndex) {
    const auto sourceGroup = size_t(groupOrder[groupIndex]);
    orderedGroupRepresentativeSeeds[groupIndex] = groupRepresentativeSeeds[sourceGroup];
    orderedGroupSeedCounts[groupIndex] = groupSeedCounts[sourceGroup];
  }
  std::vector<int> groupedSeedOffsets;
  countsToOffsets(orderedGroupSeedCounts, groupedSeedOffsets);
  std::vector<int> groupedSeedIndices(groupedSeedOffsets.back());
  for (auto groupIndex = 0u; groupIndex < numNeighborGroups; ++groupIndex) {
    const auto sourceGroup = size_t(groupOrder[groupIndex]);
    const auto begin = descriptorRunOffsets[sourceGroup];
    const auto end = descriptorRunOffsets[sourceGroup + 1u];
    CHECK(groupedSeedOffsets[groupIndex + 1u] - groupedSeedOffsets[groupIndex] == end - begin);
    std::copy(sortedSeedIndices.begin() + begin,
              sortedSeedIndices.begin() + end,
              groupedSeedIndices.begin() + groupedSeedOffsets[groupIndex]);
  }

  std::vector<int> groupTaskOffsets = groupedSeedOffsets;
  std::vector<int> groupSeedNodeListIDs(numNeighborGroups);
  std::vector<int> groupSeedNodeIDs(numNeighborGroups);
  std::vector<int> taskGroupIDs(connectivitySize);
  std::vector<int> taskNodeListIDs(connectivitySize);
  std::vector<int> taskNodeIDs(connectivitySize);
  FlatIntSpan orderedGroupRepresentativeSeedsSpan;
  GPUUtils::initMAView(orderedGroupRepresentativeSeedsSpan, orderedGroupRepresentativeSeeds);
  GPUUtils::touch(orderedGroupRepresentativeSeedsSpan, chai::CPU);
  FlatIntSpan groupedSeedOffsetsSpan;
  GPUUtils::initMAView(groupedSeedOffsetsSpan, groupedSeedOffsets);
  GPUUtils::touch(groupedSeedOffsetsSpan, chai::CPU);
  FlatIntSpan groupedSeedIndicesSpan;
  GPUUtils::initMAView(groupedSeedIndicesSpan, groupedSeedIndices);
  GPUUtils::touch(groupedSeedIndicesSpan, chai::CPU);
  FlatIntSpan groupSeedNodeListIDsSpan;
  GPUUtils::initMAView(groupSeedNodeListIDsSpan, groupSeedNodeListIDs);
  GPUUtils::touch(groupSeedNodeListIDsSpan, chai::CPU);
  FlatIntSpan groupSeedNodeIDsSpan;
  GPUUtils::initMAView(groupSeedNodeIDsSpan, groupSeedNodeIDs);
  GPUUtils::touch(groupSeedNodeIDsSpan, chai::CPU);
  FlatIntSpan taskGroupIDsSpan;
  GPUUtils::initMAView(taskGroupIDsSpan, taskGroupIDs);
  GPUUtils::touch(taskGroupIDsSpan, chai::CPU);
  FlatIntSpan taskNodeListIDsSpan;
  GPUUtils::initMAView(taskNodeListIDsSpan, taskNodeListIDs);
  GPUUtils::touch(taskNodeListIDsSpan, chai::CPU);
  FlatIntSpan taskNodeIDsSpan;
  GPUUtils::initMAView(taskNodeIDsSpan, taskNodeIDs);
  GPUUtils::touch(taskNodeIDsSpan, chai::CPU);
  RAJA::forall<RAJA::seq_exec>(RAJA::RangeSegment(0, int(numNeighborGroups)),
  [&](int groupIndex) {
    const auto repSeed = size_t(orderedGroupRepresentativeSeedsSpan[groupIndex]);
    groupSeedNodeListIDsSpan[groupIndex] = seedNodeListIDsSpan[repSeed];
    groupSeedNodeIDsSpan[groupIndex] = seedNodeIDsSpan[repSeed];
    for (auto seedOffset = groupedSeedOffsetsSpan[groupIndex];
         seedOffset < groupedSeedOffsetsSpan[groupIndex + 1u];
         ++seedOffset) {
      const auto seedIndex = size_t(groupedSeedIndicesSpan[seedOffset]);
      taskGroupIDsSpan[seedOffset] = groupIndex;
      taskNodeListIDsSpan[seedOffset] = seedNodeListIDsSpan[seedIndex];
      taskNodeIDsSpan[seedOffset] = seedNodeIDsSpan[seedIndex];
    }
  });
  GPUUtils::touch(groupSeedNodeListIDsSpan, chai::CPU);
  GPUUtils::touch(groupSeedNodeIDsSpan, chai::CPU);
  GPUUtils::touch(taskGroupIDsSpan, chai::CPU);
  GPUUtils::touch(taskNodeListIDsSpan, chai::CPU);
  GPUUtils::touch(taskNodeIDsSpan, chai::CPU);

  auto countDescriptorCoarseNeighbors = [&](const size_t repSeed,
                                            const size_t nodeListi) -> int {
    const auto& descriptor = seedDescriptorsSpan[repSeed*numNodeLists + nodeListi];
    const auto neighborKind = neighborGroupKindsSpan[nodeListi];
    if (neighborKind == TreeNeighborGroupKind) {
      CHECK(descriptor.kind == TreeNeighborGroupKind);
      return countOrFillTreeGroupCoarseNeighbors<Dimension>(treeNeighborViewsSpan[nodeListi],
                                                            descriptor.level,
                                                            descriptor.key,
                                                            nullptr);
    } else if (neighborKind == NestedNeighborGroupKind) {
      CHECK(descriptor.kind == NestedNeighborGroupKind);
      return countOrFillNestedGroupCoarseNeighbors<Dimension>(nestedNeighborViewsSpan[nodeListi],
                                                              descriptor.level,
                                                              descriptor.key,
                                                              nestedGridCellInfluenceRadiiSpan[nodeListi],
                                                              neighborSearchTypesSpan[nodeListi],
                                                              nullptr);
    } else {
      const auto& neighbor = mNodeLists[nodeListi]->neighbor();
      const auto repNodeList = size_t(seedNodeListIDsSpan[repSeed]);
      const auto repNode = seedNodeIDsSpan[repSeed];
      return neighbor.countCoarseNeighbors(position(repNodeList, repNode),
                                           H(repNodeList, repNode),
                                           ghostConnectivity);
    }
  };

  auto fillDescriptorCoarseNeighbors = [&](const size_t repSeed,
                                           const size_t nodeListi,
                                           int* result) {
    const auto& descriptor = seedDescriptorsSpan[repSeed*numNodeLists + nodeListi];
    const auto neighborKind = neighborGroupKindsSpan[nodeListi];
    if (neighborKind == TreeNeighborGroupKind) {
      CHECK(descriptor.kind == TreeNeighborGroupKind);
      countOrFillTreeGroupCoarseNeighbors<Dimension>(treeNeighborViewsSpan[nodeListi],
                                                     descriptor.level,
                                                     descriptor.key,
                                                     result);
    } else if (neighborKind == NestedNeighborGroupKind) {
      CHECK(descriptor.kind == NestedNeighborGroupKind);
      countOrFillNestedGroupCoarseNeighbors<Dimension>(nestedNeighborViewsSpan[nodeListi],
                                                       descriptor.level,
                                                       descriptor.key,
                                                       nestedGridCellInfluenceRadiiSpan[nodeListi],
                                                       neighborSearchTypesSpan[nodeListi],
                                                       result);
    } else {
      const auto& neighbor = mNodeLists[nodeListi]->neighbor();
      const auto repNodeList = size_t(seedNodeListIDsSpan[repSeed]);
      const auto repNode = seedNodeIDsSpan[repSeed];
      neighbor.fillCoarseNeighbors(position(repNodeList, repNode),
                                   H(repNodeList, repNode),
                                   result,
                                   ghostConnectivity);
    }
  };

  std::vector<int> groupRawCoarseCounts(numNeighborGroups*numNodeLists);
  FlatIntSpan groupRawCoarseCountsSpan;
  GPUUtils::initMAView(groupRawCoarseCountsSpan, groupRawCoarseCounts);
  GPUUtils::touch(groupRawCoarseCountsSpan, chai::CPU);
  RAJA::forall<RAJA::seq_exec>(RAJA::RangeSegment(0, int(numNeighborGroups*numNodeLists)),
  [&](int entryIndex) {
    const auto groupIndex = size_t(entryIndex) / numNodeLists;
    const auto nodeListi = size_t(entryIndex) % numNodeLists;
    const auto repSeed = size_t(orderedGroupRepresentativeSeedsSpan[groupIndex]);
    groupRawCoarseCountsSpan[entryIndex] = countDescriptorCoarseNeighbors(repSeed, nodeListi);
  });
  GPUUtils::touch(groupRawCoarseCountsSpan, chai::CPU);
  std::vector<int> groupRawCoarseOffsets;
  countsToOffsets(groupRawCoarseCounts, groupRawCoarseOffsets);
  FlatIntSpan groupRawCoarseOffsetsSpan;
  GPUUtils::initMAView(groupRawCoarseOffsetsSpan, groupRawCoarseOffsets);
  GPUUtils::touch(groupRawCoarseOffsetsSpan, chai::CPU);
  std::vector<int> groupRawCoarseNeighbors(groupRawCoarseOffsets.back());
  FlatIntSpan groupRawCoarseNeighborsSpan;
  GPUUtils::initMAView(groupRawCoarseNeighborsSpan, groupRawCoarseNeighbors);
  GPUUtils::touch(groupRawCoarseNeighborsSpan, chai::CPU);
  RAJA::forall<RAJA::seq_exec>(RAJA::RangeSegment(0, int(numNeighborGroups*numNodeLists)),
  [&](int entryIndex) {
    const auto groupIndex = size_t(entryIndex) / numNodeLists;
    const auto nodeListi = size_t(entryIndex) % numNodeLists;
    const auto repSeed = size_t(orderedGroupRepresentativeSeedsSpan[groupIndex]);
    const auto count = groupRawCoarseCountsSpan[entryIndex];
    if (count > 0) {
      fillDescriptorCoarseNeighbors(repSeed,
                                    nodeListi,
                                    groupRawCoarseNeighborsSpan.data() + groupRawCoarseOffsetsSpan[entryIndex]);
    }
  });
  GPUUtils::touch(groupRawCoarseNeighborsSpan, chai::CPU);
  CHECK(groupTaskOffsets.size() == numNeighborGroups + 1u);
  CHECK(groupRawCoarseOffsets.size() == numNeighborGroups*numNodeLists + 1u);
  CHECK(taskNodeIDs.size() == connectivitySize);

  FlatIntSpan groupTaskOffsetsSpan;
  GPUUtils::initMAView(groupTaskOffsetsSpan, groupTaskOffsets);
  GPUUtils::touch(groupTaskOffsetsSpan, chai::CPU);

  std::vector<Vector> groupMinMasterPosition(numNeighborGroups);
  std::vector<Vector> groupMaxMasterPosition(numNeighborGroups);
  std::vector<Vector> groupMinMasterExtent(numNeighborGroups);
  std::vector<Vector> groupMaxMasterExtent(numNeighborGroups);
  FlatVectorSpan groupMinMasterPositionSpan;
  FlatVectorSpan groupMaxMasterPositionSpan;
  FlatVectorSpan groupMinMasterExtentSpan;
  FlatVectorSpan groupMaxMasterExtentSpan;
  GPUUtils::initMAView(groupMinMasterPositionSpan, groupMinMasterPosition);
  GPUUtils::initMAView(groupMaxMasterPositionSpan, groupMaxMasterPosition);
  GPUUtils::initMAView(groupMinMasterExtentSpan, groupMinMasterExtent);
  GPUUtils::initMAView(groupMaxMasterExtentSpan, groupMaxMasterExtent);
  GPUUtils::touch(groupMinMasterPositionSpan, chai::CPU);
  GPUUtils::touch(groupMaxMasterPositionSpan, chai::CPU);
  GPUUtils::touch(groupMinMasterExtentSpan, chai::CPU);
  GPUUtils::touch(groupMaxMasterExtentSpan, chai::CPU);

  RAJA::forall<EXEC_POLICY>(TRS_UINT(0u, numNeighborGroups),
  [=] SPHERAL_HOST_DEVICE (size_t groupIndex) {
    const auto seedNodeList = size_t(groupSeedNodeListIDsSpan[groupIndex]);
    const auto seedNode = groupSeedNodeIDsSpan[groupIndex];
    const auto& seedPosition = positionView(seedNodeList, seedNode);
    const auto& seedExtent = extentView(seedNodeList, seedNode);
    auto minMasterPosition = seedPosition;
    auto maxMasterPosition = seedPosition;
    auto minMasterExtent = seedPosition - seedExtent;
    auto maxMasterExtent = seedPosition + seedExtent;
    const auto taskBegin = groupTaskOffsetsSpan[groupIndex];
    const auto taskEnd = groupTaskOffsetsSpan[groupIndex + 1u];
    for (auto taskIndex = taskBegin; taskIndex < taskEnd; ++taskIndex) {
      const auto iNodeList = size_t(taskNodeListIDsSpan[taskIndex]);
      const auto i = taskNodeIDsSpan[taskIndex];
      const auto& ri = positionView(iNodeList, i);
      const auto& exti = extentView(iNodeList, i);
      const auto minExtenti = ri - exti;
      const auto maxExtenti = ri + exti;
      for (auto k = 0; k < Dimension::nDim; ++k) {
        if (ri(k) < minMasterPosition(k)) minMasterPosition(k) = ri(k);
        if (ri(k) > maxMasterPosition(k)) maxMasterPosition(k) = ri(k);
        if (minExtenti(k) < minMasterExtent(k)) minMasterExtent(k) = minExtenti(k);
        if (maxExtenti(k) > maxMasterExtent(k)) maxMasterExtent(k) = maxExtenti(k);
      }
    }
    groupMinMasterPositionSpan[groupIndex] = minMasterPosition;
    groupMaxMasterPositionSpan[groupIndex] = maxMasterPosition;
    groupMinMasterExtentSpan[groupIndex] = minMasterExtent;
    groupMaxMasterExtentSpan[groupIndex] = maxMasterExtent;
  });
  GPUUtils::touch(groupMinMasterPositionSpan, chai::CPU);
  GPUUtils::touch(groupMaxMasterPositionSpan, chai::CPU);
  GPUUtils::touch(groupMinMasterExtentSpan, chai::CPU);
  GPUUtils::touch(groupMaxMasterExtentSpan, chai::CPU);

  std::vector<int> groupCoarseCounts(numNeighborGroups*numNodeLists, 0);
  FlatIntSpan groupCoarseCountsSpan;
  GPUUtils::initMAView(groupCoarseCountsSpan, groupCoarseCounts);
  GPUUtils::touch(groupCoarseCountsSpan, chai::CPU);
  RAJA::forall<EXEC_POLICY>(TRS_UINT(0u, numNeighborGroups*numNodeLists),
  [=] SPHERAL_HOST_DEVICE (size_t entryIndex) {
    const auto groupIndex = entryIndex / numNodeLists;
    const auto nodeList = entryIndex % numNodeLists;
    const auto rawBegin = groupRawCoarseOffsetsSpan[entryIndex];
    const auto rawEnd = groupRawCoarseOffsetsSpan[entryIndex + 1u];
    auto count = 0;
    for (auto rawIndex = rawBegin; rawIndex < rawEnd; ++rawIndex) {
      const auto j = groupRawCoarseNeighborsSpan[rawIndex];
      if (keepCoarseNeighborForGroup<Dimension>(neighborSearchTypesSpan[nodeList],
                                                positionView(nodeList, j),
                                                extentView(nodeList, j),
                                                groupMinMasterPositionSpan[groupIndex],
                                                groupMaxMasterPositionSpan[groupIndex],
                                                groupMinMasterExtentSpan[groupIndex],
                                                groupMaxMasterExtentSpan[groupIndex])) {
        ++count;
      }
    }
    groupCoarseCountsSpan[entryIndex] = count;
  });
  GPUUtils::touch(groupCoarseCountsSpan, chai::CPU);

  std::vector<int> groupCoarseOffsets;
  countsToOffsets(groupCoarseCounts, groupCoarseOffsets);
  FlatIntSpan groupCoarseOffsetsSpan;
  GPUUtils::initMAView(groupCoarseOffsetsSpan, groupCoarseOffsets);
  GPUUtils::touch(groupCoarseOffsetsSpan, chai::CPU);
  std::vector<int> groupCoarseNeighbors(groupCoarseOffsets.back());
  FlatIntSpan groupCoarseNeighborsSpan;
  GPUUtils::initMAView(groupCoarseNeighborsSpan, groupCoarseNeighbors);
  GPUUtils::touch(groupCoarseNeighborsSpan, chai::CPU);
  RAJA::forall<EXEC_POLICY>(TRS_UINT(0u, numNeighborGroups*numNodeLists),
  [=] SPHERAL_HOST_DEVICE (size_t entryIndex) {
    const auto groupIndex = entryIndex / numNodeLists;
    const auto nodeList = entryIndex % numNodeLists;
    const auto rawBegin = groupRawCoarseOffsetsSpan[entryIndex];
    const auto rawEnd = groupRawCoarseOffsetsSpan[entryIndex + 1u];
    auto count = 0;
    for (auto rawIndex = rawBegin; rawIndex < rawEnd; ++rawIndex) {
      const auto j = groupRawCoarseNeighborsSpan[rawIndex];
      if (keepCoarseNeighborForGroup<Dimension>(neighborSearchTypesSpan[nodeList],
                                                positionView(nodeList, j),
                                                extentView(nodeList, j),
                                                groupMinMasterPositionSpan[groupIndex],
                                                groupMaxMasterPositionSpan[groupIndex],
                                                groupMinMasterExtentSpan[groupIndex],
                                                groupMaxMasterExtentSpan[groupIndex])) {
        groupCoarseNeighborsSpan[groupCoarseOffsetsSpan[entryIndex] + count] = j;
        ++count;
      }
    }
    CHECK(groupCoarseOffsetsSpan[entryIndex] + count == groupCoarseOffsetsSpan[entryIndex + 1u]);
  });
  GPUUtils::touch(groupCoarseNeighborsSpan, chai::CPU);
  for (auto groupIndex = 0u; groupIndex < numNeighborGroups; ++groupIndex) {
    for (auto nodeList = 0u; nodeList < numNodeLists; ++nodeList) {
      const auto entryIndex = groupIndex*numNodeLists + nodeList;
      const auto beginOffset = size_t(groupCoarseOffsets[entryIndex]);
      const auto endOffset = size_t(groupCoarseOffsets[entryIndex + 1u]);
      if (endOffset - beginOffset > 1u) {
        auto beginItr = groupCoarseNeighbors.begin() + beginOffset;
        auto endItr = groupCoarseNeighbors.begin() + endOffset;
        if (domainDecompIndependent) {
          std::sort(beginItr, endItr,
                    [&](const int a, const int b) { return mKeys(nodeList, a) < mKeys(nodeList, b); });
        } else {
          std::sort(beginItr, endItr);
        }
      }
    }
  }
  GPUUtils::touch(groupCoarseNeighborsSpan, chai::CPU);

  // Count/scan/fill the main connectivity in flat CSR-like storage first.
  std::vector<int> connectivityCounts(numConnectivityEntries, 0);
  std::vector<int> connectivityOffsets;
  std::vector<int> connectivityNeighbors;
  auto buildConnectivityPass = [&](std::vector<int>* flatCounts,
                                   const std::vector<int>* flatOffsets,
                                   std::vector<int>* flatNeighbors) {
    const auto fillingConnectivity = (flatNeighbors != nullptr);
    if (fillingConnectivity) CHECK(flatOffsets != nullptr);
    else CHECK(flatCounts != nullptr);

    FlatIntSpan flatCountsSpan;
    FlatIntSpan flatOffsetsSpan;
    FlatIntSpan flatNeighborsSpan;
    if (fillingConnectivity) {
      GPUUtils::initMAView(flatOffsetsSpan, const_cast<std::vector<int>&>(*flatOffsets));
      GPUUtils::touch(flatOffsetsSpan, chai::CPU);
      GPUUtils::initMAView(flatNeighborsSpan, *flatNeighbors);
      GPUUtils::touch(flatNeighborsSpan, chai::CPU);
    } else {
      GPUUtils::initMAView(flatCountsSpan, *flatCounts);
      GPUUtils::touch(flatCountsSpan, chai::CPU);
    }

    const auto numTasks = taskNodeIDs.size();
    RAJA::forall<EXEC_POLICY>(TRS_UINT(0u, numTasks),
    [=] SPHERAL_HOST_DEVICE (size_t taskIndex) {
      const auto groupIndex = size_t(taskGroupIDsSpan[taskIndex]);
      const auto iNodeList = size_t(taskNodeListIDsSpan[taskIndex]);
      const auto i = taskNodeIDsSpan[taskIndex];
      const auto globalNodeIndex = size_t(offsetsSpan[iNodeList] + i);
      CHECK(globalNodeIndex < connectivitySize);

      const auto& ri = positionView(iNodeList, i);
      const auto& Hi = HView(iNodeList, i);
      for (auto jNodeList = 0u; jNodeList != numNodeLists; ++jNodeList) {
        const auto entryIndex = globalNodeIndex*numNodeLists + jNodeList;
        const auto beginOffset = (fillingConnectivity ? flatOffsetsSpan[entryIndex] : 0);
        const auto coarseEntryIndex = groupIndex*numNodeLists + jNodeList;
        const auto coarseBegin = groupCoarseOffsetsSpan[coarseEntryIndex];
        const auto coarseEnd = groupCoarseOffsetsSpan[coarseEntryIndex + 1u];
        auto count = 0;

        for (auto coarseIndex = coarseBegin; coarseIndex < coarseEnd; ++coarseIndex) {
          const auto j = groupCoarseNeighborsSpan[coarseIndex];
          const auto& rj = positionView(jNodeList, j);
          const auto& Hj = HView(jNodeList, j);
          const auto rij = ri - rj;
          const auto eta2i = (Hi*rij).magnitude2();
          const auto eta2j = (Hj*rij).magnitude2();
          if ((eta2i <= kernelExtent2 or eta2j <= kernelExtent2) and
              ((iNodeList != jNodeList) or (i != j))) {
            if (fillingConnectivity) {
              flatNeighborsSpan[beginOffset + count] = j;
            }
            ++count;
          }
        }

        if (fillingConnectivity) {
          CHECK(beginOffset + count == flatOffsetsSpan[entryIndex + 1u]);
        } else {
          flatCountsSpan[entryIndex] = count;
        }
      }
    });

    if (fillingConnectivity) {
      GPUUtils::touch(flatNeighborsSpan, chai::CPU);
    } else {
      GPUUtils::touch(flatCountsSpan, chai::CPU);
    }
  };

  buildConnectivityPass(&connectivityCounts, nullptr, nullptr);
  countsToOffsets(connectivityCounts, connectivityOffsets);
  connectivityNeighbors.resize(connectivityOffsets.back());
  buildConnectivityPass(nullptr, &connectivityOffsets, &connectivityNeighbors);

  mConnectivityFlatOffsets = connectivityOffsets;
  mConnectivityFlatNeighbors = connectivityNeighbors;
  GPUUtils::initMAView(mConnectivityFlatOffsetsSpan, mConnectivityFlatOffsets);
  GPUUtils::initMAView(mConnectivityFlatNeighborsSpan, mConnectivityFlatNeighbors);
  GPUUtils::touch(mConnectivityFlatOffsetsSpan, chai::CPU);
  GPUUtils::touch(mConnectivityFlatNeighborsSpan, chai::CPU);
  const auto connectivityView = this->connectivityFlatView();

  unflattenConnectivity(connectivitySize,
                        numNodeLists,
                        mConnectivityFlatOffsets,
                        mConnectivityFlatNeighbors,
                        mConnectivity);

  std::vector<int> globalNodeListIDs(connectivitySize);
  std::vector<int> globalNodeIDs(connectivitySize);
  for (auto iNodeList = 0u; iNodeList != numNodeLists; ++iNodeList) {
    const auto ni = (ghostConnectivity ?
                     mNodeLists[iNodeList]->numNodes() :
                     mNodeLists[iNodeList]->numInternalNodes());
    for (auto i = 0u; i < ni; ++i) {
      const auto globalNodeIndex = size_t(mOffsets[iNodeList] + i);
      globalNodeListIDs[globalNodeIndex] = iNodeList;
      globalNodeIDs[globalNodeIndex] = i;
    }
  }
  FlatIntSpan globalNodeListIDsSpan;
  GPUUtils::initMAView(globalNodeListIDsSpan, globalNodeListIDs);
  GPUUtils::touch(globalNodeListIDsSpan, chai::CPU);
  FlatIntSpan globalNodeIDsSpan;
  GPUUtils::initMAView(globalNodeIDsSpan, globalNodeIDs);
  GPUUtils::touch(globalNodeIDsSpan, chai::CPU);

  std::vector<int> nodePairCounts(numConnectivityEntries, 0);
  auto buildNodePairPass = [&](std::vector<int>* pairCounts,
                               const std::vector<int>* pairOffsets,
                               NodePairList* nodePairs) {
    const auto fillingNodePairs = (nodePairs != nullptr);
    if (fillingNodePairs) CHECK(pairOffsets != nullptr);
    else CHECK(pairCounts != nullptr);

    FlatIntSpan pairCountsSpan;
    FlatIntSpan pairOffsetsSpan;
    NodePairListView nodePairsView;
    if (fillingNodePairs) {
      GPUUtils::initMAView(pairOffsetsSpan, const_cast<std::vector<int>&>(*pairOffsets));
      GPUUtils::touch(pairOffsetsSpan, chai::CPU);
      nodePairsView = nodePairs->view();
      nodePairsView.touch(chai::CPU);
    } else {
      GPUUtils::initMAView(pairCountsSpan, *pairCounts);
      GPUUtils::touch(pairCountsSpan, chai::CPU);
    }

    RAJA::forall<EXEC_POLICY>(TRS_UINT(0u, numConnectivityEntries),
    [=] SPHERAL_HOST_DEVICE (size_t entryIndex) {
      const auto globalNodeIndex = entryIndex / numNodeLists;
      const auto iNodeList = size_t(globalNodeListIDsSpan[globalNodeIndex]);
      const auto i = globalNodeIDsSpan[globalNodeIndex];
      const auto jNodeList = entryIndex % numNodeLists;
      const auto firstGhostNodej = firstGhostNodesSpan[jNodeList];
      const auto neighborCount = connectivityView.size(entryIndex);
      auto pairCount = 0;
      for (auto k = 0u; k < neighborCount; ++k) {
        const auto j = connectivityView(entryIndex, k);
        if (shouldCalculatePairInteraction(iNodeList, i, jNodeList, j, firstGhostNodej)) {
          if (fillingNodePairs) {
            nodePairsView[pairOffsetsSpan[entryIndex] + pairCount] =
              NodePairIdxType(i, iNodeList, j, jNodeList);
          }
          ++pairCount;
        }
      }
      if (fillingNodePairs) {
        CHECK(pairOffsetsSpan[entryIndex] + pairCount == pairOffsetsSpan[entryIndex + 1u]);
      } else {
        pairCountsSpan[entryIndex] = pairCount;
      }
    });

    if (fillingNodePairs) {
      nodePairsView.touch(chai::CPU);
    } else {
      GPUUtils::touch(pairCountsSpan, chai::CPU);
    }
  };

  buildNodePairPass(&nodePairCounts, nullptr, nullptr);
  std::vector<int> nodePairOffsets;
  countsToOffsets(nodePairCounts, nodePairOffsets);
  mNodePairListPtr = std::make_shared<NodePairList>(NodePairList::ContainerType(nodePairOffsets.back()));
  buildNodePairPass(nullptr, &nodePairOffsets, mNodePairListPtr.get());

  // // If necessary add ghost->internal connectivity.
  // if (ghostConnectivity) {
  //   for (auto iNodeList = 0; iNodeList < numNodeLists; ++iNodeList) {
  //     for (auto i = 0; i < mNodeLists[iNodeList]->numInternalNodes(); ++i) {
  //       const auto& neighborsi = mConnectivity[mOffsets[iNodeList] + i];
  //       CHECK(neighborsi.size() == numNodeLists);
  //       for (auto jNodeList = 0; jNodeList < numNodeLists; ++jNodeList) {
  //         const auto firstGhostNodej = mNodeLists[jNodeList]->firstGhostNode();
  //         for (auto jItr = neighborsi[jNodeList].begin(); jItr < neighborsi[jNodeList].end(); ++jItr) {
  //           const auto j = *jItr;
  //           if (j >= firstGhostNodej) {
  //             auto& neighborsj = mConnectivity[mOffsets[jNodeList] + j];
  //             CHECK(neighborsj.size() == numNodeLists);
  //             neighborsj[iNodeList].push_back(i);
  //             // mNodePairListPtr->push_back(NodePairIdxType(i, iNodeList, j, jNodeList));
  //           }
  //         }
  //       }
  //     }
  //   }

  //   // Flag ghost nodes as done if at least one neighbor was found
  //   for (auto iNodeList = 0; iNodeList < numNodeLists; ++iNodeList) {
  //     for (auto i = mNodeLists[iNodeList]->numInternalNodes(); i < mNodeLists[iNodeList]->numNodes(); ++i) {
  //       const auto& neighborsi = mConnectivity[mOffsets[iNodeList] + i];
  //       if (neighborsi.size() > 0) {
  //         flagNodeDone(iNodeList, i) = 1;
  //       }
  //     }
  //   }
  // }

  // Sort the NodePairList in order to enforce domain decomposition independence.
  if (domainDecompIndependent) {
    // sort(mNodePairListPtr->begin(), mNodePairListPtr->end(), [this](const NodePairIdxType& a, const NodePairIdxType& b) { return (mKeys(a.i_list, a.i_node) + mKeys(a.j_list, a.j_node)) < (mKeys(b.i_list, b.i_node) + mKeys(b.j_list, b.j_node)); });
    // sort(mNodePairListPtr->begin(), mNodePairListPtr->end(), [this](const NodePairIdxType& a, const NodePairIdxType& b) { return hashKeys(mKeys(a.i_list, a.i_node), mKeys(a.j_list, a.j_node)) < hashKeys(mKeys(b.i_list, b.i_node), mKeys(b.j_list, b.j_node)); });
    sortPairs(*mNodePairListPtr, mKeys);
  } else {
    std::sort(mNodePairListPtr->begin(), mNodePairListPtr->end());
  }
  // mNodePairListPtr->computeLookup();

  // Do we need overlap connectivity?
  if (mBuildOverlapConnectivity) {
    // VERIFY2(ghostConnectivity, "ghost connectivity is required for overlap connectivity");
    TIME_BEGIN("ConnectivityMap_computeOverlapConnectivity");

    const auto connectivity = this->connectivityFlatView();
    std::vector<int> overlapRawCounts(numConnectivityEntries, 0);
    FlatIntSpan overlapRawCountsSpan;
    GPUUtils::initMAView(overlapRawCountsSpan, overlapRawCounts);
    GPUUtils::touch(overlapRawCountsSpan, chai::CPU);

    for (auto iNodeList = 0u; iNodeList < numNodeLists; ++iNodeList) {
      const auto ni = mNodeLists[iNodeList]->numNodes();
      RAJA::forall<EXEC_POLICY>(TRS_UINT(0u, ni),
      [=] SPHERAL_HOST_DEVICE (size_t ii) {
        const auto i = int(ii);
        const auto globalNodeIndex = size_t(offsetsSpan[iNodeList] + i);
        for (auto jNodeList = 0u; jNodeList < numNodeLists; ++jNodeList) {
          const auto entryIndex = connectivity.entryIndex(globalNodeIndex, jNodeList);
          overlapRawCountsSpan[entryIndex] =
            fillRawOverlapNeighborsForEntry(connectivity,
                                           offsetsSpan,
                                           positionView,
                                           HView,
                                           kernelExtent2,
                                           numNodeLists,
                                           globalNodeIndex,
                                           iNodeList,
                                           i,
                                           jNodeList,
                                           nullptr);
        }
      });
    }
    GPUUtils::touch(overlapRawCountsSpan, chai::CPU);

    std::vector<int> overlapRawOffsets;
    countsToOffsets(overlapRawCounts, overlapRawOffsets);
    FlatIntSpan overlapRawOffsetsSpan;
    GPUUtils::initMAView(overlapRawOffsetsSpan, overlapRawOffsets);
    GPUUtils::touch(overlapRawOffsetsSpan, chai::CPU);

    std::vector<int> overlapRawNeighbors(overlapRawOffsets.back());
    FlatIntSpan overlapRawNeighborsSpan;
    GPUUtils::initMAView(overlapRawNeighborsSpan, overlapRawNeighbors);
    GPUUtils::touch(overlapRawNeighborsSpan, chai::CPU);

    for (auto iNodeList = 0u; iNodeList < numNodeLists; ++iNodeList) {
      const auto ni = mNodeLists[iNodeList]->numNodes();
      RAJA::forall<EXEC_POLICY>(TRS_UINT(0u, ni),
      [=] SPHERAL_HOST_DEVICE (size_t ii) {
        const auto i = int(ii);
        const auto globalNodeIndex = size_t(offsetsSpan[iNodeList] + i);
        for (auto jNodeList = 0u; jNodeList < numNodeLists; ++jNodeList) {
          const auto entryIndex = connectivity.entryIndex(globalNodeIndex, jNodeList);
          auto* entryPtr = (overlapRawNeighborsSpan.size() > 0u ?
                            overlapRawNeighborsSpan.data() + overlapRawOffsetsSpan[entryIndex] :
                            nullptr);
          fillRawOverlapNeighborsForEntry(connectivity,
                                         offsetsSpan,
                                         positionView,
                                         HView,
                                         kernelExtent2,
                                         numNodeLists,
                                         globalNodeIndex,
                                         iNodeList,
                                         i,
                                         jNodeList,
                                         entryPtr);
        }
      });
    }
    GPUUtils::touch(overlapRawNeighborsSpan, chai::CPU);

    std::vector<int> overlapCounts(numConnectivityEntries, 0);
    FlatIntSpan overlapCountsSpan;
    GPUUtils::initMAView(overlapCountsSpan, overlapCounts);
    GPUUtils::touch(overlapCountsSpan, chai::CPU);

    RAJA::forall<EXEC_POLICY>(TRS_UINT(0u, numConnectivityEntries),
    [=] SPHERAL_HOST_DEVICE (size_t entryIndex) {
      const auto beginOffset = overlapRawOffsetsSpan[entryIndex];
      const auto count = size_t(overlapRawOffsetsSpan[entryIndex + 1] - beginOffset);
      auto* values = (count > 0u ? overlapRawNeighborsSpan.data() + beginOffset : nullptr);
      const auto targetNodeList = size_t(entryIndex) % numNodeLists;
      if (count > 1u) {
        sortOverlapNeighbors(values, count, targetNodeList, keysView, domainDecompIndependent);
      }
      overlapCountsSpan[entryIndex] = uniqueOverlapNeighbors(values, count);
    });
    GPUUtils::touch(overlapCountsSpan, chai::CPU);

    std::vector<int> overlapOffsets;
    countsToOffsets(overlapCounts, overlapOffsets);
    FlatIntSpan overlapOffsetsSpan;
    GPUUtils::initMAView(overlapOffsetsSpan, overlapOffsets);
    GPUUtils::touch(overlapOffsetsSpan, chai::CPU);

    std::vector<int> overlapNeighbors(overlapOffsets.back());
    FlatIntSpan overlapNeighborsSpan;
    GPUUtils::initMAView(overlapNeighborsSpan, overlapNeighbors);
    GPUUtils::touch(overlapNeighborsSpan, chai::CPU);

    RAJA::forall<EXEC_POLICY>(TRS_UINT(0u, numConnectivityEntries),
    [=] SPHERAL_HOST_DEVICE (size_t entryIndex) {
      const auto beginSrc = overlapRawOffsetsSpan[entryIndex];
      const auto beginDst = overlapOffsetsSpan[entryIndex];
      const auto count = size_t(overlapCountsSpan[entryIndex]);
      for (auto k = 0u; k < count; ++k) {
        overlapNeighborsSpan[beginDst + k] = overlapRawNeighborsSpan[beginSrc + k];
      }
    });
    GPUUtils::touch(overlapNeighborsSpan, chai::CPU);

    unflattenConnectivity(connectivitySize,
                          numNodeLists,
                          overlapOffsets,
                          overlapNeighbors,
                          mOverlapConnectivity);
    TIME_END("ConnectivityMap_computeOverlapConnectivity");
  } else {
    mOverlapConnectivity.clear();
  }

  // Are we building intersection connectivity?
  if (mBuildIntersectionConnectivity) {
    TIME_BEGIN("ConnectivityMap_precomputeIntersectionConnectivity");
    auto& pairs = *mNodePairListPtr;
    const auto npairs = pairs.size();
#pragma omp parallel
    {
      IntersectionConnectivityContainer intersection_private;
#pragma omp for
      for (auto k = 0u; k < npairs; ++k) {
        const auto& pair = pairs[k];
        intersection_private[pair] = this->connectivityIntersectionForNodes(pair.i_list, pair.i_node,
                                                                            pair.j_list, pair.j_node,
                                                                            position);
      }
#pragma omp critical
      {
        for (const auto& element: intersection_private) {
          mIntersectionConnectivity[element.first] = element.second;
        }
      } // omp critical
    }   // omp parallel
    TIME_END("ConnectivityMap_precomputeIntersectionConnectivity");
  }

  // {
  //   tpre = allReduce(unsigned(tpre), SPHERAL_OP_SUM) / Process::getTotalNumberOfProcesses() / CLOCKS_PER_SEC;
  //   tmaster = allReduce(unsigned(tmaster), SPHERAL_OP_SUM) / Process::getTotalNumberOfProcesses() / CLOCKS_PER_SEC;
  //   trefine = allReduce(unsigned(trefine), SPHERAL_OP_SUM) / Process::getTotalNumberOfProcesses() / CLOCKS_PER_SEC;
  //   twalk = allReduce(unsigned(twalk), SPHERAL_OP_SUM) / Process::getTotalNumberOfProcesses() / CLOCKS_PER_SEC;
  //   if (Process::getRank() == 0) {
  //     std::cerr << "ConnectivityMap timings (pre, master, refine, walk) = " << tpre << " " << tmaster << " " << trefine << " " << twalk << std::endl;
  //   }
  // }

  // Post conditions.
  BEGIN_CONTRACT_SCOPE
  ENSURE2(taskNodeIDs.size() == connectivitySize,
          "Missed connectivity tasks: " << taskNodeIDs.size() << " " << connectivitySize);
  // Make sure we're ready to be used.
  ENSURE(valid());
  END_CONTRACT_SCOPE

  this->rebuildFlatConnectivityViews();

  TIME_END("ConnectivityMap_computeConnectivity");
}

}
