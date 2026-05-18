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
#include "care/KeyValueSorter_impl.h"

#include <algorithm>
#include <cmath>
#include <ctime>
#include <tuple>
#include <type_traits>
using std::vector;
using std::map;
using std::string;
using std::pair;

namespace Spheral {

namespace {
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
// Device-safe helper computing the max geometric neighbor extent from H.
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
double
maxNeighborExtent(const typename Dimension::SymTensor& H,
                  const double kernelExtent);

template<>
SPHERAL_HOST_DEVICE
inline
double
maxNeighborExtent<Dim<1>>(const Dim<1>::SymTensor& H,
                          const double kernelExtent) {
  return kernelExtent/H.xx();
}

template<>
SPHERAL_HOST_DEVICE
inline
double
maxNeighborExtent<Dim<2>>(const Dim<2>::SymTensor& H,
                          const double kernelExtent) {
  const auto Hdet = H.Determinant();
  const auto M = H.square();
  const auto ex = kernelExtent/Hdet*sqrt(M.yy());
  const auto ey = kernelExtent/Hdet*sqrt(M.xx());
  return (ex > ey ? ex : ey);
}

template<>
SPHERAL_HOST_DEVICE
inline
double
maxNeighborExtent<Dim<3>>(const Dim<3>::SymTensor& H,
                          const double kernelExtent) {
  const auto Hdet = H.Determinant();
  const auto M = H.square();
  const auto ex = kernelExtent/Hdet*sqrt(M.yy()*M.zz() - M.yz()*M.zy());
  const auto ey = kernelExtent/Hdet*sqrt(M.xx()*M.zz() - M.xz()*M.zx());
  const auto ez = kernelExtent/Hdet*sqrt(M.xx()*M.yy() - M.xy()*M.yx());
  return std::max(ex, std::max(ey, ez));
}

SPHERAL_HOST_DEVICE
inline
uint64_t
packGridCellKey(const int ix) {
  return uint64_t(ix + (1 << 20));
}

SPHERAL_HOST_DEVICE
inline
uint64_t
packGridCellKey(const int ix,
                const int iy) {
  const auto x = uint64_t(ix + (1 << 20));
  const auto y = uint64_t(iy + (1 << 20));
  return x | (y << 21);
}

SPHERAL_HOST_DEVICE
inline
uint64_t
packGridCellKey(const int ix,
                const int iy,
                const int iz) {
  const auto x = uint64_t(ix + (1 << 20));
  const auto y = uint64_t(iy + (1 << 20));
  const auto z = uint64_t(iz + (1 << 20));
  return x | (y << 21) | (z << 42);
}

template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
NeighborGroupDescriptor
makeTreeNeighborGroupDescriptor(const typename Dimension::Vector& position,
                                const typename Dimension::SymTensor& H,
                                const double boxLength,
                                const double gridLevelConst0,
                                const typename Dimension::Vector& xmin) {
  NeighborGroupDescriptor result;
  result.kind = TreeNeighborGroupKind;
  const auto h = 1.0/H.eigenValues().minElement();
  const auto maxLevel = int(TreeNeighbor<Dimension>::num1DBits()) - 1;
  const auto rawLevel = int(gridLevelConst0 - log(h)/log(2.0));
  result.level = std::max(0, std::min(maxLevel, rawLevel));
  const auto ncell = uint64_t(1U << result.level);
  const auto maxcell = ncell - 1U;
  const auto tx = std::max(0.0, std::min(1.0, (position.x() - xmin.x())/boxLength));
  const auto ty = std::max(0.0, std::min(1.0, (position.y() - xmin.y())/boxLength));
  const auto tz = std::max(0.0, std::min(1.0, (position.z() - xmin.z())/boxLength));
  const auto ix = std::min(maxcell, uint64_t(tx * ncell));
  const auto iy = std::min(maxcell, uint64_t(ty * ncell));
  const auto iz = std::min(maxcell, uint64_t(tz * ncell));
  result.key = ((std::max(uint64_t(0), std::min(TreeNeighbor<Dimension>::max1dKeyValue(), iz)) << 2*TreeNeighbor<Dimension>::num1DBits()) +
                (std::max(uint64_t(0), std::min(TreeNeighbor<Dimension>::max1dKeyValue(), iy)) <<   TreeNeighbor<Dimension>::num1DBits()) +
                (std::max(uint64_t(0), std::min(TreeNeighbor<Dimension>::max1dKeyValue(), ix))));
  return result;
}

template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
NeighborGroupDescriptor
makeNestedNeighborGroupDescriptor(const typename Dimension::Vector& position,
                                  const typename Dimension::SymTensor& H,
                                  const double kernelExtent,
                                  const int maxGridLevels,
                                  const double gridLevelConst0,
                                  const double topGridCellSizeInv,
                                  const typename Dimension::Vector& gridOrigin);

template<>
SPHERAL_HOST_DEVICE
inline
NeighborGroupDescriptor
makeNestedNeighborGroupDescriptor<Dim<1>>(const Dim<1>::Vector& position,
                                          const Dim<1>::SymTensor& H,
                                          const double kernelExtent,
                                          const int maxGridLevels,
                                          const double gridLevelConst0,
                                          const double topGridCellSizeInv,
                                          const Dim<1>::Vector& gridOrigin) {
  NeighborGroupDescriptor result;
  result.kind = NestedNeighborGroupKind;
  const auto h = maxNeighborExtent<Dim<1>>(H, kernelExtent);
  const auto rawLevel = int(gridLevelConst0 - log(h)/log(2.0));
  result.level = std::max(0, std::min(maxGridLevels - 1, rawLevel));
  const auto cellSizeInv = double(1u << result.level) * topGridCellSizeInv;
  const auto ix = int((position.x() - gridOrigin.x())*cellSizeInv) - (position.x() < gridOrigin.x() ? 1 : 0);
  result.key = packGridCellKey(ix);
  return result;
}

template<>
SPHERAL_HOST_DEVICE
inline
NeighborGroupDescriptor
makeNestedNeighborGroupDescriptor<Dim<2>>(const Dim<2>::Vector& position,
                                          const Dim<2>::SymTensor& H,
                                          const double kernelExtent,
                                          const int maxGridLevels,
                                          const double gridLevelConst0,
                                          const double topGridCellSizeInv,
                                          const Dim<2>::Vector& gridOrigin) {
  NeighborGroupDescriptor result;
  result.kind = NestedNeighborGroupKind;
  const auto h = maxNeighborExtent<Dim<2>>(H, kernelExtent);
  const auto rawLevel = int(gridLevelConst0 - log(h)/log(2.0));
  result.level = std::max(0, std::min(maxGridLevels - 1, rawLevel));
  const auto cellSizeInv = double(1u << result.level) * topGridCellSizeInv;
  const auto ix = int((position.x() - gridOrigin.x())*cellSizeInv) - (position.x() < gridOrigin.x() ? 1 : 0);
  const auto iy = int((position.y() - gridOrigin.y())*cellSizeInv) - (position.y() < gridOrigin.y() ? 1 : 0);
  result.key = packGridCellKey(ix, iy);
  return result;
}

template<>
SPHERAL_HOST_DEVICE
inline
NeighborGroupDescriptor
makeNestedNeighborGroupDescriptor<Dim<3>>(const Dim<3>::Vector& position,
                                          const Dim<3>::SymTensor& H,
                                          const double kernelExtent,
                                          const int maxGridLevels,
                                          const double gridLevelConst0,
                                          const double topGridCellSizeInv,
                                          const Dim<3>::Vector& gridOrigin) {
  NeighborGroupDescriptor result;
  result.kind = NestedNeighborGroupKind;
  const auto h = maxNeighborExtent<Dim<3>>(H, kernelExtent);
  const auto rawLevel = int(gridLevelConst0 - log(h)/log(2.0));
  result.level = std::max(0, std::min(maxGridLevels - 1, rawLevel));
  const auto cellSizeInv = double(1u << result.level) * topGridCellSizeInv;
  const auto ix = int((position.x() - gridOrigin.x())*cellSizeInv) - (position.x() < gridOrigin.x() ? 1 : 0);
  const auto iy = int((position.y() - gridOrigin.y())*cellSizeInv) - (position.y() < gridOrigin.y() ? 1 : 0);
  const auto iz = int((position.z() - gridOrigin.z())*cellSizeInv) - (position.z() < gridOrigin.z() ? 1 : 0);
  result.key = packGridCellKey(ix, iy, iz);
  return result;
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
// Number of global-node entries stored for the given NodeList in a flat layout.
//------------------------------------------------------------------------------
inline
size_t
nodeListEntryCount(const std::vector<int>& offsets,
                   const size_t totalEntries,
                   const size_t nodeList) {
  CHECK(nodeList < offsets.size());
  return size_t((nodeList + 1u < offsets.size() ? offsets[nodeList + 1u] : totalEntries) - offsets[nodeList]);
}

//------------------------------------------------------------------------------
// Build the canonical flat traversal ordering directly.
//------------------------------------------------------------------------------
template<typename Dimension, typename Key>
inline
void
buildTraversalOrdering(const std::vector<const NodeList<Dimension>*>& nodeLists,
                       const FieldList<Dimension, Key>& keys,
                       const bool domainDecompIndependent,
                       std::vector<int>& flatOffsets,
                       std::vector<int>& flatIndices) {
  const auto numNodeLists = nodeLists.size();
  flatOffsets.resize(numNodeLists + 1u);
  flatOffsets[0] = 0;
  for (auto iNodeList = 0u; iNodeList != numNodeLists; ++iNodeList) {
    flatOffsets[iNodeList + 1u] = flatOffsets[iNodeList] +
                                  (domainDecompIndependent ?
                                   nodeLists[iNodeList]->numNodes() :
                                   nodeLists[iNodeList]->numInternalNodes());
  }

  flatIndices.resize(flatOffsets.back());
  for (auto iNodeList = 0u; iNodeList != numNodeLists; ++iNodeList) {
    const auto beginOffset = size_t(flatOffsets[iNodeList]);
    const auto count = size_t(flatOffsets[iNodeList + 1u] - flatOffsets[iNodeList]);
    if (domainDecompIndependent) {
      std::vector<std::pair<int, Key>> orderedKeys;
      orderedKeys.reserve(count);
      for (auto i = 0u; i != count; ++i) orderedKeys.emplace_back(i, keys(iNodeList, i));
      sort(orderedKeys.begin(), orderedKeys.end(), ComparePairsBySecondElement<std::pair<int, Key>>());
      for (auto i = 0u; i != count; ++i) flatIndices[beginOffset + i] = orderedKeys[i].first;
    } else {
      for (auto i = 0u; i != count; ++i) flatIndices[beginOffset + i] = i;
    }
  }
}

//------------------------------------------------------------------------------
// Patch canonical flat connectivity after node deletion/reindexing.
//------------------------------------------------------------------------------
template<typename Dimension, typename Key>
inline
void
patchFlatConnectivity(const size_t oldNumEntries,
                      const std::vector<int>& sourceFlatOffsets,
                      const std::vector<int>& sourceFlatNeighbors,
                      const std::vector<int>& oldOffsets,
                      const std::vector<int>& newOffsets,
                      const FieldList<Dimension, size_t>& flags,
                      const FieldList<Dimension, size_t>& old2new,
                      const FieldList<Dimension, Key>& keys,
                      const bool domainDecompIndependent,
                      const size_t numNodeLists,
                      std::vector<int>& flatOffsets,
                      std::vector<int>& flatNeighbors) {
  const auto newNumEntries = size_t(newOffsets.back());
  std::vector<int> counts(newNumEntries*numNodeLists, 0);

  for (auto iNodeList = 0u; iNodeList != numNodeLists; ++iNodeList) {
    const auto oldBegin = size_t(oldOffsets[iNodeList]);
    const auto ni = nodeListEntryCount(oldOffsets, oldNumEntries, iNodeList);
    for (auto i = 0u; i != ni; ++i) {
      if (flags(iNodeList, i) == 0) continue;
      const auto srcNodeIndex = oldBegin + i;
      const auto dstNodeIndex = size_t(newOffsets[iNodeList] + old2new(iNodeList, i));
      for (auto jNodeList = 0u; jNodeList != numNodeLists; ++jNodeList) {
        const auto entryIndex = srcNodeIndex*numNodeLists + jNodeList;
        auto count = 0;
        for (auto k = sourceFlatOffsets[entryIndex]; k != sourceFlatOffsets[entryIndex + 1u]; ++k) {
          const auto j = sourceFlatNeighbors[k];
          if (flags(jNodeList, j) != 0) ++count;
        }
        counts[dstNodeIndex*numNodeLists + jNodeList] = count;
      }
    }
  }

  countsToOffsets(counts, flatOffsets);
  flatNeighbors.resize(flatOffsets.back());
  for (auto iNodeList = 0u; iNodeList != numNodeLists; ++iNodeList) {
    const auto oldBegin = size_t(oldOffsets[iNodeList]);
    const auto ni = nodeListEntryCount(oldOffsets, oldNumEntries, iNodeList);
    for (auto i = 0u; i != ni; ++i) {
      if (flags(iNodeList, i) == 0) continue;
      const auto srcNodeIndex = oldBegin + i;
      const auto dstNodeIndex = size_t(newOffsets[iNodeList] + old2new(iNodeList, i));
      for (auto jNodeList = 0u; jNodeList != numNodeLists; ++jNodeList) {
        const auto srcEntryIndex = srcNodeIndex*numNodeLists + jNodeList;
        const auto dstEntryIndex = dstNodeIndex*numNodeLists + jNodeList;
        auto dstOffset = size_t(flatOffsets[dstEntryIndex]);
        if (domainDecompIndependent) {
          std::vector<std::pair<int, Key>> neighborKeys;
          neighborKeys.reserve(sourceFlatOffsets[srcEntryIndex + 1u] - sourceFlatOffsets[srcEntryIndex]);
          for (auto k = sourceFlatOffsets[srcEntryIndex]; k != sourceFlatOffsets[srcEntryIndex + 1u]; ++k) {
            const auto j = sourceFlatNeighbors[k];
            if (flags(jNodeList, j) != 0) neighborKeys.emplace_back(old2new(jNodeList, j), keys(jNodeList, j));
          }
          sort(neighborKeys.begin(), neighborKeys.end(), ComparePairsBySecondElement<std::pair<int, Key>>());
          for (const auto& entry: neighborKeys) flatNeighbors[dstOffset++] = entry.first;
        } else {
          std::vector<int> patchedNeighbors;
          patchedNeighbors.reserve(sourceFlatOffsets[srcEntryIndex + 1u] - sourceFlatOffsets[srcEntryIndex]);
          for (auto k = sourceFlatOffsets[srcEntryIndex]; k != sourceFlatOffsets[srcEntryIndex + 1u]; ++k) {
            const auto j = sourceFlatNeighbors[k];
            if (flags(jNodeList, j) != 0) patchedNeighbors.push_back(old2new(jNodeList, j));
          }
          sort(patchedNeighbors.begin(), patchedNeighbors.end());
          for (const auto j: patchedNeighbors) flatNeighbors[dstOffset++] = j;
        }
        CHECK(dstOffset == size_t(flatOffsets[dstEntryIndex + 1u]));
      }
    }
  }
}

//------------------------------------------------------------------------------
// Patch canonical flat traversal ordering after node deletion/reindexing.
//------------------------------------------------------------------------------
template<typename Dimension, typename Key>
inline
void
patchTraversalOrdering(const std::vector<int>& currentOffsets,
                       const FieldList<Dimension, size_t>& flags,
                       const FieldList<Dimension, size_t>& old2new,
                       const FieldList<Dimension, Key>& keys,
                       const bool domainDecompIndependent,
                       std::vector<int>& flatOffsets,
                       std::vector<int>& flatIndices) {
  const auto numNodeLists = currentOffsets.size() - 1u;
  flatOffsets.resize(numNodeLists + 1u);
  flatOffsets[0] = 0;
  for (auto iNodeList = 0u; iNodeList != numNodeLists; ++iNodeList) {
    const auto n = size_t(currentOffsets[iNodeList + 1u] - currentOffsets[iNodeList]);
    auto count = 0;
    for (auto i = 0u; i != n; ++i) count += (flags(iNodeList, i) != 0 ? 1 : 0);
    flatOffsets[iNodeList + 1u] = flatOffsets[iNodeList] + count;
  }

  flatIndices.resize(flatOffsets.back());
  for (auto iNodeList = 0u; iNodeList != numNodeLists; ++iNodeList) {
    const auto beginOffset = size_t(flatOffsets[iNodeList]);
    const auto n = size_t(currentOffsets[iNodeList + 1u] - currentOffsets[iNodeList]);
    const auto count = size_t(flatOffsets[iNodeList + 1u] - flatOffsets[iNodeList]);
    if (domainDecompIndependent) {
      std::vector<std::pair<int, Key>> orderedKeys;
      orderedKeys.reserve(count);
      for (auto i = 0u; i != n; ++i) {
        if (flags(iNodeList, i) != 0) orderedKeys.emplace_back(old2new(iNodeList, i), keys(iNodeList, i));
      }
      sort(orderedKeys.begin(), orderedKeys.end(), ComparePairsBySecondElement<std::pair<int, Key>>());
      for (auto i = 0u; i != count; ++i) flatIndices[beginOffset + i] = orderedKeys[i].first;
    } else {
      for (auto i = 0u; i != count; ++i) flatIndices[beginOffset + i] = i;
    }
  }
}

//------------------------------------------------------------------------------
// Patch canonical flat intersection connectivity after node deletion/reindexing.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
patchIntersectionConnectivity(const std::vector<size_t>& sourcePairIndices,
                              const size_t numNodeLists,
                              const std::vector<int>& oldFlatOffsets,
                              const std::vector<int>& oldFlatNeighbors,
                              const FieldList<Dimension, size_t>& flags,
                              const FieldList<Dimension, size_t>& old2new,
                              std::vector<int>& flatOffsets,
                              std::vector<int>& flatNeighbors) {
  std::vector<int> counts(sourcePairIndices.size()*numNodeLists, 0);
  for (auto pairIndex = 0u; pairIndex != sourcePairIndices.size(); ++pairIndex) {
    const auto sourcePairIndex = sourcePairIndices[pairIndex];
    for (auto nodeList = 0u; nodeList != numNodeLists; ++nodeList) {
      const auto entryIndex = sourcePairIndex*numNodeLists + nodeList;
      auto count = 0;
      for (auto k = oldFlatOffsets[entryIndex]; k != oldFlatOffsets[entryIndex + 1u]; ++k) {
        if (flags(nodeList, oldFlatNeighbors[k]) != 0) ++count;
      }
      counts[pairIndex*numNodeLists + nodeList] = count;
    }
  }

  countsToOffsets(counts, flatOffsets);
  flatNeighbors.resize(flatOffsets.back());
  for (auto pairIndex = 0u; pairIndex != sourcePairIndices.size(); ++pairIndex) {
    const auto sourcePairIndex = sourcePairIndices[pairIndex];
    for (auto nodeList = 0u; nodeList != numNodeLists; ++nodeList) {
      const auto srcEntryIndex = sourcePairIndex*numNodeLists + nodeList;
      const auto dstEntryIndex = pairIndex*numNodeLists + nodeList;
      auto dstOffset = size_t(flatOffsets[dstEntryIndex]);
      for (auto k = oldFlatOffsets[srcEntryIndex]; k != oldFlatOffsets[srcEntryIndex + 1u]; ++k) {
        const auto oldNode = oldFlatNeighbors[k];
        if (flags(nodeList, oldNode) != 0) flatNeighbors[dstOffset++] = old2new(nodeList, oldNode);
      }
      CHECK(dstOffset == size_t(flatOffsets[dstEntryIndex + 1u]));
    }
  }
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
  mNodeTraversalOffsets(),
  mNodeTraversalIndicesFlat(),
  mNodeTraversalOffsetsSpan(),
  mNodeTraversalIndicesFlatSpan(),
  mNodeTraversalIndices(),
  mKeys(FieldStorageType::CopyFields),
  mCouplingPtr(std::make_shared<NodeCoupling>()),
  mIntersectionConnectivity(),
  mIntersectionConnectivityFlatOffsets(),
  mIntersectionConnectivityFlatNeighbors(),
  mIntersectionConnectivityFlatOffsetsSpan(),
  mIntersectionConnectivityFlatNeighborsSpan(),
  mConnectivityCacheValid(false),
  mOverlapConnectivityCacheValid(false),
  mNodeTraversalCacheValid(false),
  mIntersectionConnectivityCacheValid(false) {
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
  GPUUtils::freeMAView(mNodeTraversalOffsetsSpan);
  GPUUtils::freeMAView(mNodeTraversalIndicesFlatSpan);
  GPUUtils::freeMAView(mIntersectionConnectivityFlatOffsetsSpan);
  GPUUtils::freeMAView(mIntersectionConnectivityFlatNeighborsSpan);
}

//------------------------------------------------------------------------------
// Refresh the flat connectivity views from the canonical flat storage.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConnectivityMap<Dimension>::
rebuildFlatConnectivityViews() {
  GPUUtils::initMAView(mConnectivityFlatOffsetsSpan, mConnectivityFlatOffsets);
  GPUUtils::initMAView(mConnectivityFlatNeighborsSpan, mConnectivityFlatNeighbors);
  GPUUtils::touch(mConnectivityFlatOffsetsSpan, chai::CPU);
  GPUUtils::touch(mConnectivityFlatNeighborsSpan, chai::CPU);
  mConnectivityCacheValid = false;
  mConnectivity.clear();

  if (mBuildOverlapConnectivity) {
    GPUUtils::initMAView(mOverlapConnectivityFlatOffsetsSpan, mOverlapConnectivityFlatOffsets);
    GPUUtils::initMAView(mOverlapConnectivityFlatNeighborsSpan, mOverlapConnectivityFlatNeighbors);
    GPUUtils::touch(mOverlapConnectivityFlatOffsetsSpan, chai::CPU);
    GPUUtils::touch(mOverlapConnectivityFlatNeighborsSpan, chai::CPU);
    mOverlapConnectivityCacheValid = false;
    mOverlapConnectivity.clear();
  } else {
    mOverlapConnectivityFlatOffsets.clear();
    mOverlapConnectivityFlatNeighbors.clear();
    GPUUtils::freeMAView(mOverlapConnectivityFlatOffsetsSpan);
    GPUUtils::freeMAView(mOverlapConnectivityFlatNeighborsSpan);
    mOverlapConnectivity.clear();
    mOverlapConnectivityCacheValid = true;
  }
}

//------------------------------------------------------------------------------
// Refresh the flat traversal view from the canonical flat storage.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConnectivityMap<Dimension>::
rebuildFlatTraversalViews() {
  GPUUtils::initMAView(mNodeTraversalOffsetsSpan, mNodeTraversalOffsets);
  GPUUtils::initMAView(mNodeTraversalIndicesFlatSpan, mNodeTraversalIndicesFlat);
  GPUUtils::touch(mNodeTraversalOffsetsSpan, chai::CPU);
  GPUUtils::touch(mNodeTraversalIndicesFlatSpan, chai::CPU);
  mNodeTraversalIndices.clear();
  mNodeTraversalCacheValid = false;
}

//------------------------------------------------------------------------------
// Refresh the flat intersection view from the canonical flat storage.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConnectivityMap<Dimension>::
rebuildFlatIntersectionConnectivityViews() {
  if (mBuildIntersectionConnectivity) {
    GPUUtils::initMAView(mIntersectionConnectivityFlatOffsetsSpan, mIntersectionConnectivityFlatOffsets);
    GPUUtils::initMAView(mIntersectionConnectivityFlatNeighborsSpan, mIntersectionConnectivityFlatNeighbors);
    GPUUtils::touch(mIntersectionConnectivityFlatOffsetsSpan, chai::CPU);
    GPUUtils::touch(mIntersectionConnectivityFlatNeighborsSpan, chai::CPU);
    mIntersectionConnectivity.clear();
    mIntersectionConnectivityCacheValid = false;
  } else {
    mIntersectionConnectivityFlatOffsets.clear();
    mIntersectionConnectivityFlatNeighbors.clear();
    GPUUtils::freeMAView(mIntersectionConnectivityFlatOffsetsSpan);
    GPUUtils::freeMAView(mIntersectionConnectivityFlatNeighborsSpan);
    mIntersectionConnectivity.clear();
    mIntersectionConnectivityCacheValid = true;
  }
}

//------------------------------------------------------------------------------
// Rebuild the legacy ragged compatibility caches from the canonical flat data.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConnectivityMap<Dimension>::
ensureConnectivityCache() const {
  if (mConnectivityCacheValid) return;
  unflattenConnectivity(this->connectivityFlatView().numNodes(),
                        mNodeLists.size(),
                        mConnectivityFlatOffsets,
                        mConnectivityFlatNeighbors,
                        mConnectivity);
  mConnectivityCacheValid = true;
}

template<typename Dimension>
void
ConnectivityMap<Dimension>::
ensureOverlapConnectivityCache() const {
  if (mOverlapConnectivityCacheValid) return;
  if (not mBuildOverlapConnectivity) {
    mOverlapConnectivity.clear();
    mOverlapConnectivityCacheValid = true;
    return;
  }
  const auto numNodeLists = mNodeLists.size();
  const auto numEntries = (numNodeLists > 0u and mOverlapConnectivityFlatOffsets.size() > 0u ?
                           (mOverlapConnectivityFlatOffsets.size() - 1u)/numNodeLists :
                           0u);
  unflattenConnectivity(numEntries,
                        numNodeLists,
                        mOverlapConnectivityFlatOffsets,
                        mOverlapConnectivityFlatNeighbors,
                        mOverlapConnectivity);
  mOverlapConnectivityCacheValid = true;
}

template<typename Dimension>
void
ConnectivityMap<Dimension>::
ensureTraversalCache() const {
  if (mNodeTraversalCacheValid) return;
  const auto numNodeLists = mNodeLists.size();
  mNodeTraversalIndices.assign(numNodeLists, std::vector<int>());
  CHECK(mNodeTraversalOffsets.size() == numNodeLists + 1u);
  for (auto nodeList = 0u; nodeList < numNodeLists; ++nodeList) {
    const auto beginOffset = size_t(mNodeTraversalOffsets[nodeList]);
    const auto endOffset = size_t(mNodeTraversalOffsets[nodeList + 1u]);
    mNodeTraversalIndices[nodeList].assign(mNodeTraversalIndicesFlat.begin() + beginOffset,
                                           mNodeTraversalIndicesFlat.begin() + endOffset);
  }
  mNodeTraversalCacheValid = true;
}

template<typename Dimension>
void
ConnectivityMap<Dimension>::
ensureIntersectionConnectivityCache() const {
  if (mIntersectionConnectivityCacheValid) return;
  mIntersectionConnectivity.clear();
  if (not mBuildIntersectionConnectivity) {
    mIntersectionConnectivityCacheValid = true;
    return;
  }
  REQUIRE(mNodePairListPtr);
  const auto& pairs = *mNodePairListPtr;
  const auto numNodeLists = mNodeLists.size();
  const auto expectedOffsets = pairs.size()*numNodeLists + 1u;
  if (mIntersectionConnectivityFlatOffsets.size() != expectedOffsets) {
    mIntersectionConnectivityCacheValid = true;
    return;
  }
  for (auto pairIndex = 0u; pairIndex < pairs.size(); ++pairIndex) {
    std::vector<std::vector<int>> intersection(numNodeLists);
    for (auto nodeList = 0u; nodeList < numNodeLists; ++nodeList) {
      const auto entryIndex = pairIndex*numNodeLists + nodeList;
      const auto beginOffset = size_t(mIntersectionConnectivityFlatOffsets[entryIndex]);
      const auto endOffset = size_t(mIntersectionConnectivityFlatOffsets[entryIndex + 1u]);
      intersection[nodeList].assign(mIntersectionConnectivityFlatNeighbors.begin() + beginOffset,
                                    mIntersectionConnectivityFlatNeighbors.begin() + endOffset);
    }
    mIntersectionConnectivity[pairs[pairIndex]] = std::move(intersection);
  }
  mIntersectionConnectivityCacheValid = true;
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
  const auto oldOffsets = mOffsets;
  const auto oldConnectivityFlatOffsets = mConnectivityFlatOffsets;
  const auto oldConnectivityFlatNeighbors = mConnectivityFlatNeighbors;
  const auto oldConnectivityNumEntries = (numNodeLists > 0u and oldConnectivityFlatOffsets.size() > 0u ?
                                          (oldConnectivityFlatOffsets.size() - 1u)/numNodeLists :
                                          0u);
  const auto oldOverlapConnectivityFlatOffsets = mOverlapConnectivityFlatOffsets;
  const auto oldOverlapConnectivityFlatNeighbors = mOverlapConnectivityFlatNeighbors;
  const auto oldTraversalOffsets = mNodeTraversalOffsets;
  const auto oldIntersectionOffsets = mIntersectionConnectivityFlatOffsets;
  const auto oldIntersectionNeighbors = mIntersectionConnectivityFlatNeighbors;
  if (domainDecompIndependent) {
// #pragma omp parallel for collapse(2)
    for (auto iNodeList = 0u; iNodeList < numNodeLists; ++iNodeList) {
      for (auto i = 0u; i < mNodeLists[iNodeList]->numNodes(); ++i) {
        if (flags(iNodeList, i) == 0) mKeys(iNodeList, i) = KeyTraits::maxKey;
      }
    }
  }

  std::vector<int> newConnectivityOffsets(numNodeLists + 1u, 0);
  for (auto iNodeList = 0u; iNodeList != numNodeLists; ++iNodeList) {
    const auto numNodes = nodeListEntryCount(oldOffsets, oldConnectivityNumEntries, iNodeList);
    auto count = 0;
    for (auto i = 0u; i != numNodes; ++i) count += (flags(iNodeList, i) != 0 ? 1 : 0);
    newConnectivityOffsets[iNodeList + 1u] = newConnectivityOffsets[iNodeList] + count;
  }

  patchFlatConnectivity(oldConnectivityNumEntries,
                        oldConnectivityFlatOffsets,
                        oldConnectivityFlatNeighbors,
                        oldOffsets,
                        newConnectivityOffsets,
                        flags,
                        old2new,
                        mKeys,
                        domainDecompIndependent,
                        numNodeLists,
                        mConnectivityFlatOffsets,
                        mConnectivityFlatNeighbors);
  if (mBuildOverlapConnectivity) {
    patchFlatConnectivity(oldConnectivityNumEntries,
                          oldOverlapConnectivityFlatOffsets,
                          oldOverlapConnectivityFlatNeighbors,
                          oldOffsets,
                          newConnectivityOffsets,
                          flags,
                          old2new,
                          mKeys,
                          domainDecompIndependent,
                          numNodeLists,
                          mOverlapConnectivityFlatOffsets,
                          mOverlapConnectivityFlatNeighbors);
  }

  mOffsets.resize(numNodeLists);
  for (auto iNodeList = 0u; iNodeList != numNodeLists; ++iNodeList) {
    mOffsets[iNodeList] = newConnectivityOffsets[iNodeList];
  }

  patchTraversalOrdering(oldTraversalOffsets,
                         flags,
                         old2new,
                         mKeys,
                         domainDecompIndependent,
                         mNodeTraversalOffsets,
                         mNodeTraversalIndicesFlat);

  REQUIRE(mNodePairListPtr);
  const auto& currentPairs = *mNodePairListPtr;
  const auto currentPairCount = currentPairs.size();
  struct PairWithSource {
    NodePairIdxType pair;
    size_t sourceIndex;
  };
  std::vector<PairWithSource> culledPairs;
  culledPairs.reserve(currentPairs.size());
  for (auto k = 0u; k < currentPairs.size(); ++k) {
    const auto iNodeList = currentPairs[k].i_list;
    const auto jNodeList = currentPairs[k].j_list;
    const auto i = currentPairs[k].i_node;
    const auto j = currentPairs[k].j_node;
    if (flags(iNodeList, i) != 0 and flags(jNodeList, j) != 0) {
      culledPairs.push_back({NodePairIdxType(old2new(iNodeList, i), iNodeList,
                                             old2new(jNodeList, j), jNodeList),
                             k});
    }
  }

  if (domainDecompIndependent) {
    std::sort(culledPairs.begin(), culledPairs.end(),
              [this, &currentPairs](const PairWithSource& a, const PairWithSource& b) {
                const auto& apair = currentPairs[a.sourceIndex];
                const auto& bpair = currentPairs[b.sourceIndex];
                return hashKeys(mKeys(apair.i_list, apair.i_node), mKeys(apair.j_list, apair.j_node)) <
                       hashKeys(mKeys(bpair.i_list, bpair.i_node), mKeys(bpair.j_list, bpair.j_node));
              });
  } else {
    std::sort(culledPairs.begin(), culledPairs.end(),
              [](const PairWithSource& a, const PairWithSource& b) {
                return std::tie(a.pair.i_list, a.pair.i_node, a.pair.j_list, a.pair.j_node) <
                       std::tie(b.pair.i_list, b.pair.i_node, b.pair.j_list, b.pair.j_node);
              });
  }

  std::vector<NodePairIdxType> sortedPairs;
  std::vector<size_t> sourcePairIndices;
  sortedPairs.reserve(culledPairs.size());
  sourcePairIndices.reserve(culledPairs.size());
  for (const auto& entry: culledPairs) {
    sortedPairs.push_back(entry.pair);
    sourcePairIndices.push_back(entry.sourceIndex);
  }
  mNodePairListPtr = std::make_shared<NodePairList>(std::move(sortedPairs));

  if (mBuildIntersectionConnectivity) {
    if (oldIntersectionOffsets.size() == currentPairCount*numNodeLists + 1u) {
      patchIntersectionConnectivity(sourcePairIndices,
                                    numNodeLists,
                                    oldIntersectionOffsets,
                                    oldIntersectionNeighbors,
                                    flags,
                                    old2new,
                                    mIntersectionConnectivityFlatOffsets,
                                    mIntersectionConnectivityFlatNeighbors);
    } else {
      mIntersectionConnectivityFlatOffsets.clear();
      mIntersectionConnectivityFlatNeighbors.clear();
    }
  } else {
    mIntersectionConnectivityFlatOffsets.clear();
    mIntersectionConnectivityFlatNeighbors.clear();
  }

  this->rebuildFlatConnectivityViews();
  this->rebuildFlatTraversalViews();
  this->rebuildFlatIntersectionConnectivityViews();
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
  const auto oldConnectivityFlatOffsets = mConnectivityFlatOffsets;
  const auto oldConnectivityFlatNeighbors = mConnectivityFlatNeighbors;
  const auto numEntries = (numNodeLists > 0u and oldConnectivityFlatOffsets.size() > 0u ?
                           (oldConnectivityFlatOffsets.size() - 1u)/numNodeLists :
                           0u);
  std::vector<int> counts(numEntries*numNodeLists, 0);
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto n = nodeListEntryCount(mOffsets, numEntries, nodeListi);
    for (auto i = 0u; i < n; ++i) {
      const auto& allneighbors = neighborsToCut(nodeListi, i);
      CHECK(allneighbors.size() == 0 or allneighbors.size() == numNodeLists);
      for (auto nodeListj = 0u; nodeListj < numNodeLists; ++nodeListj) {
        const auto entryIndex = (mOffsets[nodeListi] + i)*numNodeLists + nodeListj;
        const auto oldCount = size_t(oldConnectivityFlatOffsets[entryIndex + 1u] - oldConnectivityFlatOffsets[entryIndex]);
        if (nodeListj >= allneighbors.size() or allneighbors[nodeListj].empty()) {
          counts[entryIndex] = oldCount;
        } else {
          auto cutIndices = allneighbors[nodeListj];
          sort(cutIndices.begin(), cutIndices.end());
          cutIndices.erase(unique(cutIndices.begin(), cutIndices.end()), cutIndices.end());
          auto keepCount = int(oldCount);
          for (const auto k: cutIndices) {
            if (k >= 0 and size_t(k) < oldCount) --keepCount;
          }
          counts[entryIndex] = keepCount;
        }
      }
    }
  }

  countsToOffsets(counts, mConnectivityFlatOffsets);
  mConnectivityFlatNeighbors.resize(mConnectivityFlatOffsets.back());
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto n = nodeListEntryCount(mOffsets, numEntries, nodeListi);
    for (auto i = 0u; i < n; ++i) {
      const auto& allneighbors = neighborsToCut(nodeListi, i);
      for (auto nodeListj = 0u; nodeListj < numNodeLists; ++nodeListj) {
        const auto entryIndex = (mOffsets[nodeListi] + i)*numNodeLists + nodeListj;
        auto dstOffset = size_t(mConnectivityFlatOffsets[entryIndex]);
        const auto srcBegin = oldConnectivityFlatOffsets[entryIndex];
        const auto srcEnd = oldConnectivityFlatOffsets[entryIndex + 1u];
        if (nodeListj >= allneighbors.size() or allneighbors[nodeListj].empty()) {
          for (auto k = srcBegin; k != srcEnd; ++k) {
            mConnectivityFlatNeighbors[dstOffset++] = oldConnectivityFlatNeighbors[k];
          }
        } else {
          auto cutIndices = allneighbors[nodeListj];
          sort(cutIndices.begin(), cutIndices.end());
          cutIndices.erase(unique(cutIndices.begin(), cutIndices.end()), cutIndices.end());
          auto cutItr = cutIndices.begin();
          auto k = 0u;
          for (auto srcOffset = srcBegin; srcOffset != srcEnd; ++srcOffset) {
            while (cutItr != cutIndices.end() and *cutItr < int(k)) ++cutItr;
            if (cutItr == cutIndices.end() or *cutItr != int(k)) {
              mConnectivityFlatNeighbors[dstOffset++] = oldConnectivityFlatNeighbors[srcOffset];
            }
            ++k;
          }
        }
        CHECK(dstOffset == size_t(mConnectivityFlatOffsets[entryIndex + 1u]));
      }
    }
  }

  this->rebuildFlatConnectivityViews();
  mIntersectionConnectivity.clear();
  mIntersectionConnectivityFlatOffsets.clear();
  mIntersectionConnectivityFlatNeighbors.clear();
  this->rebuildFlatIntersectionConnectivityViews();

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
        ((nodeListIDi == numNodeLists - 1u) and ((connectivity.numNodes() - (size_t)mOffsets[nodeListIDi]) != numNodes))) {
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
  if (mNodeTraversalOffsets.size() != mNodeLists.size() + 1u) {
    cerr << "ConnectivityMap::valid: mNodeTraversalOffsets wrong size!" << endl;
    return false;
  }
  for (auto nodeList = 0u; nodeList < numNodeLists; ++nodeList) {
    const auto numExpected = domainDecompIndependent ? mNodeLists[nodeList]->numNodes() : mNodeLists[nodeList]->numInternalNodes();
    const auto beginOffset = size_t(mNodeTraversalOffsets[nodeList]);
    const auto endOffset = size_t(mNodeTraversalOffsets[nodeList + 1u]);
    bool ok = endOffset - beginOffset == numExpected;
    for (auto i = 0u; i < numExpected; ++i) {
      ok = ok and (count(mNodeTraversalIndicesFlat.begin() + beginOffset,
                         mNodeTraversalIndicesFlat.begin() + endOffset,
                         i) == 1);
    }
    if (not ok) {
      cerr << "ConnectivityMap::valid: mNodeTraversalIndicesFlat elements messed up!" << endl;
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
  mIntersectionConnectivity.clear();
  mIntersectionConnectivityCacheValid = false;

  // If we're trying to be domain decomposition independent, we need a key to sort
  // by that will give us a unique ordering regardless of position.  The Morton ordered
  // hash fills the bill.
  using Key = typename KeyTraits::Key;
  if (domainDecompIndependent) mKeys = mortonOrderIndices(dataBase);

  // Fill in the ordering for walking the nodes.
  buildTraversalOrdering(mNodeLists,
                         mKeys,
                         domainDecompIndependent,
                         mNodeTraversalOffsets,
                         mNodeTraversalIndicesFlat);

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
  using FlatDoubleSpan = SPHERAL_SPAN_TYPE<double>;
  using FlatDescriptorSpan = SPHERAL_SPAN_TYPE<NeighborGroupDescriptor>;
  using FlatUInt64Span = SPHERAL_SPAN_TYPE<uint64_t>;
  using FlatTreeNeighborViewSpan = SPHERAL_SPAN_TYPE<TreeNeighborView<Dimension>>;
  using FlatNestedNeighborViewSpan = SPHERAL_SPAN_TYPE<NestedGridNeighborView<Dimension>>;
#else
  using FlatVectorSpan = chai::ManagedArray<Vector>;
  using FlatDoubleSpan = chai::ManagedArray<double>;
  using FlatDescriptorSpan = chai::ManagedArray<NeighborGroupDescriptor>;
  using FlatUInt64Span = chai::ManagedArray<uint64_t>;
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
  std::vector<int> nestedMaxGridLevels(numNodeLists, 0);
  std::vector<double> neighborKernelExtents(numNodeLists, 0.0);
  std::vector<double> treeBoxLengths(numNodeLists, 1.0);
  std::vector<double> treeGridLevelConst0(numNodeLists, 0.0);
  std::vector<Vector> treeXmins(numNodeLists, Vector::zero());
  std::vector<double> nestedGridLevelConst0(numNodeLists, 0.0);
  std::vector<double> nestedTopGridCellSizeInv(numNodeLists, 0.0);
  std::vector<Vector> nestedGridOrigins(numNodeLists, Vector::zero());
  std::vector<TreeNeighborView<Dimension>> treeNeighborViews(numNodeLists);
  std::vector<NestedGridNeighborView<Dimension>> nestedNeighborViews(numNodeLists);
  for (auto nodeList = 0u; nodeList < numNodeLists; ++nodeList) {
    const auto& neighbor = mNodeLists[nodeList]->neighbor();
    neighborKernelExtents[nodeList] = neighbor.kernelExtent();
    neighborGroupKinds[nodeList] = neighbor.groupKind();
    if (neighborGroupKinds[nodeList] == TreeNeighborGroupKind) {
      const auto& treeNeighbor = static_cast<const TreeNeighbor<Dimension>&>(neighbor);
      treeNeighborViews[nodeList] = treeNeighbor.view();
      treeBoxLengths[nodeList] = treeNeighbor.boxLength();
      treeGridLevelConst0[nodeList] = log(treeNeighbor.boxLength()/neighbor.kernelExtent())/log(2.0);
      treeXmins[nodeList] = treeNeighbor.xmin();
    } else if (neighborGroupKinds[nodeList] == NestedNeighborGroupKind) {
      const auto& nestedNeighbor = static_cast<const NestedGridNeighbor<Dimension>&>(neighbor);
      nestedNeighborViews[nodeList] = nestedNeighbor.view();
      nestedGridCellInfluenceRadii[nodeList] = nestedNeighbor.gridCellInfluenceRadius();
      nestedMaxGridLevels[nodeList] = nestedNeighbor.numGridLevels();
      nestedGridLevelConst0[nodeList] = log(nestedNeighbor.topGridSize()*nestedNeighbor.gridCellInfluenceRadius())/log(2.0);
      nestedTopGridCellSizeInv[nodeList] = 1.0/nestedNeighbor.topGridSize();
      nestedGridOrigins[nodeList] = nestedNeighbor.origin();
    }
  }
  FlatIntSpan neighborGroupKindsSpan;
  GPUUtils::initMAView(neighborGroupKindsSpan, neighborGroupKinds);
  GPUUtils::touch(neighborGroupKindsSpan, chai::CPU);
  FlatIntSpan nestedGridCellInfluenceRadiiSpan;
  GPUUtils::initMAView(nestedGridCellInfluenceRadiiSpan, nestedGridCellInfluenceRadii);
  GPUUtils::touch(nestedGridCellInfluenceRadiiSpan, chai::CPU);
  FlatIntSpan nestedMaxGridLevelsSpan;
  GPUUtils::initMAView(nestedMaxGridLevelsSpan, nestedMaxGridLevels);
  GPUUtils::touch(nestedMaxGridLevelsSpan, chai::CPU);
  FlatDoubleSpan neighborKernelExtentsSpan;
  GPUUtils::initMAView(neighborKernelExtentsSpan, neighborKernelExtents);
  GPUUtils::touch(neighborKernelExtentsSpan, chai::CPU);
  FlatDoubleSpan treeBoxLengthsSpan;
  GPUUtils::initMAView(treeBoxLengthsSpan, treeBoxLengths);
  GPUUtils::touch(treeBoxLengthsSpan, chai::CPU);
  FlatDoubleSpan treeGridLevelConst0Span;
  GPUUtils::initMAView(treeGridLevelConst0Span, treeGridLevelConst0);
  GPUUtils::touch(treeGridLevelConst0Span, chai::CPU);
  FlatVectorSpan treeXminsSpan;
  GPUUtils::initMAView(treeXminsSpan, treeXmins);
  GPUUtils::touch(treeXminsSpan, chai::CPU);
  FlatDoubleSpan nestedGridLevelConst0Span;
  GPUUtils::initMAView(nestedGridLevelConst0Span, nestedGridLevelConst0);
  GPUUtils::touch(nestedGridLevelConst0Span, chai::CPU);
  FlatDoubleSpan nestedTopGridCellSizeInvSpan;
  GPUUtils::initMAView(nestedTopGridCellSizeInvSpan, nestedTopGridCellSizeInv);
  GPUUtils::touch(nestedTopGridCellSizeInvSpan, chai::CPU);
  FlatVectorSpan nestedGridOriginsSpan;
  GPUUtils::initMAView(nestedGridOriginsSpan, nestedGridOrigins);
  GPUUtils::touch(nestedGridOriginsSpan, chai::CPU);
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

  using SortExecPolicy = std::conditional_t<std::is_same_v<EXEC_POLICY, RAJA::seq_exec>,
                                            RAJA::seq_exec,
                                            EXEC_POLICY>;

  std::vector<NeighborGroupDescriptor> seedDescriptors(connectivitySize*numNodeLists);
  FlatDescriptorSpan seedDescriptorsSpan;
  GPUUtils::initMAView(seedDescriptorsSpan, seedDescriptors);
  GPUUtils::touch(seedDescriptorsSpan, chai::CPU);
  RAJA::forall<EXEC_POLICY>(TRS_UINT(0u, connectivitySize*numNodeLists),
  [=] SPHERAL_HOST_DEVICE (size_t entryIndex) {
    const auto seedIndex = entryIndex / numNodeLists;
    const auto nodeListi = entryIndex % numNodeLists;
    const auto iiNodeList = size_t(seedNodeListIDsSpan[seedIndex]);
    const auto ii = seedNodeIDsSpan[seedIndex];
    const auto& ri = positionView(iiNodeList, ii);
    const auto& Hi = HView(iiNodeList, ii);
    const auto neighborKind = neighborGroupKindsSpan[nodeListi];
    NeighborGroupDescriptor descriptor;
    if (neighborKind == TreeNeighborGroupKind) {
      descriptor = makeTreeNeighborGroupDescriptor<Dimension>(ri,
                                                              Hi,
                                                              treeBoxLengthsSpan[nodeListi],
                                                              treeGridLevelConst0Span[nodeListi],
                                                              treeXminsSpan[nodeListi]);
    } else if (neighborKind == NestedNeighborGroupKind) {
      descriptor = makeNestedNeighborGroupDescriptor<Dimension>(ri,
                                                                Hi,
                                                                neighborKernelExtentsSpan[nodeListi],
                                                                nestedMaxGridLevelsSpan[nodeListi],
                                                                nestedGridLevelConst0Span[nodeListi],
                                                                nestedTopGridCellSizeInvSpan[nodeListi],
                                                                nestedGridOriginsSpan[nodeListi]);
    } else {
      descriptor.kind = FallbackNeighborGroupKind;
      descriptor.level = 0;
      descriptor.key = seedIndex;
    }
    seedDescriptorsSpan[entryIndex] = descriptor;
  });
  GPUUtils::touch(seedDescriptorsSpan, chai::CPU);

  care::host_device_ptr<int> sortedSeedIndices(connectivitySize, "sortedSeedIndices");
  care::host_device_ptr<int> descriptorSortKeysInt(connectivitySize, "descriptorSortKeysInt");
  care::host_device_ptr<uint64_t> descriptorSortKeysUInt64(connectivitySize, "descriptorSortKeysUInt64");
  RAJA::forall<EXEC_POLICY>(TRS_UINT(0u, connectivitySize),
  [=] SPHERAL_HOST_DEVICE (size_t seedIndex) {
    sortedSeedIndices[seedIndex] = seedIndex;
  });
  for (int nodeListi = int(numNodeLists) - 1; nodeListi >= 0; --nodeListi) {
    RAJA::forall<EXEC_POLICY>(TRS_UINT(0u, connectivitySize),
    [=] SPHERAL_HOST_DEVICE (size_t sortedIndex) {
      const auto seedIndex = size_t(sortedSeedIndices[sortedIndex]);
      descriptorSortKeysUInt64[sortedIndex] = seedDescriptorsSpan[seedIndex*numNodeLists + size_t(nodeListi)].key;
    });
    care::sortKeyValueArrays<SortExecPolicy, uint64_t, int>(descriptorSortKeysUInt64, sortedSeedIndices, 0, connectivitySize);

    RAJA::forall<EXEC_POLICY>(TRS_UINT(0u, connectivitySize),
    [=] SPHERAL_HOST_DEVICE (size_t sortedIndex) {
      const auto seedIndex = size_t(sortedSeedIndices[sortedIndex]);
      descriptorSortKeysInt[sortedIndex] = seedDescriptorsSpan[seedIndex*numNodeLists + size_t(nodeListi)].level;
    });
    care::sortKeyValueArrays<SortExecPolicy, int, int>(descriptorSortKeysInt, sortedSeedIndices, 0, connectivitySize);

    RAJA::forall<EXEC_POLICY>(TRS_UINT(0u, connectivitySize),
    [=] SPHERAL_HOST_DEVICE (size_t sortedIndex) {
      const auto seedIndex = size_t(sortedSeedIndices[sortedIndex]);
      descriptorSortKeysInt[sortedIndex] = seedDescriptorsSpan[seedIndex*numNodeLists + size_t(nodeListi)].kind;
    });
    care::sortKeyValueArrays<SortExecPolicy, int, int>(descriptorSortKeysInt, sortedSeedIndices, 0, connectivitySize);
  }
  sortedSeedIndices.registerTouch(care::CPU);

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

  std::vector<int> orderedGroupRepresentativeSeeds(numNeighborGroups);
  std::vector<int> orderedGroupSeedCounts(numNeighborGroups);
  std::vector<int> representativeSeedToGroup(connectivitySize, -1);
  for (auto groupIndex = 0u; groupIndex < numNeighborGroups; ++groupIndex) {
    representativeSeedToGroup[groupRepresentativeSeeds[groupIndex]] = groupIndex;
  }
  for (auto repSeed = 0u, orderedGroupIndex = 0u; repSeed < connectivitySize; ++repSeed) {
    const auto sourceGroup = representativeSeedToGroup[repSeed];
    if (sourceGroup >= 0) {
      orderedGroupRepresentativeSeeds[orderedGroupIndex] = repSeed;
      orderedGroupSeedCounts[orderedGroupIndex] = groupSeedCounts[sourceGroup];
      ++orderedGroupIndex;
    }
  }
  CHECK(std::count_if(representativeSeedToGroup.begin(),
                      representativeSeedToGroup.end(),
                      [](const int groupIndex) { return groupIndex >= 0; }) == int(numNeighborGroups));
  std::vector<int> groupedSeedOffsets;
  countsToOffsets(orderedGroupSeedCounts, groupedSeedOffsets);
  std::vector<int> groupedSeedIndices(groupedSeedOffsets.back());
  for (auto groupIndex = 0u; groupIndex < numNeighborGroups; ++groupIndex) {
    const auto sourceGroup = size_t(representativeSeedToGroup[orderedGroupRepresentativeSeeds[groupIndex]]);
    const auto begin = descriptorRunOffsets[sourceGroup];
    const auto end = descriptorRunOffsets[sourceGroup + 1u];
    CHECK(groupedSeedOffsets[groupIndex + 1u] - groupedSeedOffsets[groupIndex] == end - begin);
    for (auto offset = 0u; offset < end - begin; ++offset) {
      groupedSeedIndices[groupedSeedOffsets[groupIndex] + offset] = sortedSeedIndices[begin + offset];
    }
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
  RAJA::forall<EXEC_POLICY>(TRS_UINT(0u, numNeighborGroups),
  [=] SPHERAL_HOST_DEVICE (size_t groupIndex) {
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

  const auto hasFallbackNeighborKinds =
    std::find(neighborGroupKinds.begin(), neighborGroupKinds.end(), FallbackNeighborGroupKind) != neighborGroupKinds.end();

  std::vector<int> groupRawCoarseCounts(numNeighborGroups*numNodeLists);
  FlatIntSpan groupRawCoarseCountsSpan;
  GPUUtils::initMAView(groupRawCoarseCountsSpan, groupRawCoarseCounts);
  GPUUtils::touch(groupRawCoarseCountsSpan, chai::CPU);
  RAJA::forall<EXEC_POLICY>(TRS_UINT(0u, numNeighborGroups*numNodeLists),
  [=] SPHERAL_HOST_DEVICE (size_t entryIndex) {
    const auto groupIndex = entryIndex / numNodeLists;
    const auto nodeListi = entryIndex % numNodeLists;
    const auto repSeed = size_t(orderedGroupRepresentativeSeedsSpan[groupIndex]);
    const auto& descriptor = seedDescriptorsSpan[repSeed*numNodeLists + nodeListi];
    const auto neighborKind = neighborGroupKindsSpan[nodeListi];
    if (neighborKind == TreeNeighborGroupKind) {
      CHECK(descriptor.kind == TreeNeighborGroupKind);
      groupRawCoarseCountsSpan[entryIndex] = countOrFillTreeGroupCoarseNeighbors<Dimension>(treeNeighborViewsSpan[nodeListi],
                                                                                             descriptor.level,
                                                                                             descriptor.key,
                                                                                             nullptr);
    } else if (neighborKind == NestedNeighborGroupKind) {
      CHECK(descriptor.kind == NestedNeighborGroupKind);
      groupRawCoarseCountsSpan[entryIndex] = countOrFillNestedGroupCoarseNeighbors<Dimension>(nestedNeighborViewsSpan[nodeListi],
                                                                                               descriptor.level,
                                                                                               descriptor.key,
                                                                                               nestedGridCellInfluenceRadiiSpan[nodeListi],
                                                                                               neighborSearchTypesSpan[nodeListi],
                                                                                               nullptr);
    } else {
      groupRawCoarseCountsSpan[entryIndex] = 0;
    }
  });
  GPUUtils::touch(groupRawCoarseCountsSpan, chai::CPU);
  if (hasFallbackNeighborKinds) {
    for (auto groupIndex = 0u; groupIndex < numNeighborGroups; ++groupIndex) {
      const auto repSeed = size_t(orderedGroupRepresentativeSeeds[groupIndex]);
      const auto repNodeList = size_t(seedNodeListIDs[repSeed]);
      const auto repNode = seedNodeIDs[repSeed];
      for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
        if (neighborGroupKinds[nodeListi] == FallbackNeighborGroupKind) {
          const auto& neighbor = mNodeLists[nodeListi]->neighbor();
          groupRawCoarseCounts[groupIndex*numNodeLists + nodeListi] =
            neighbor.countCoarseNeighbors(position(repNodeList, repNode),
                                          H(repNodeList, repNode),
                                          ghostConnectivity);
        }
      }
    }
  }
  std::vector<int> groupRawCoarseOffsets;
  countsToOffsets(groupRawCoarseCounts, groupRawCoarseOffsets);
  FlatIntSpan groupRawCoarseOffsetsSpan;
  GPUUtils::initMAView(groupRawCoarseOffsetsSpan, groupRawCoarseOffsets);
  GPUUtils::touch(groupRawCoarseOffsetsSpan, chai::CPU);
  std::vector<int> groupRawCoarseNeighbors(groupRawCoarseOffsets.back());
  FlatIntSpan groupRawCoarseNeighborsSpan;
  GPUUtils::initMAView(groupRawCoarseNeighborsSpan, groupRawCoarseNeighbors);
  GPUUtils::touch(groupRawCoarseNeighborsSpan, chai::CPU);
  RAJA::forall<EXEC_POLICY>(TRS_UINT(0u, numNeighborGroups*numNodeLists),
  [=] SPHERAL_HOST_DEVICE (size_t entryIndex) {
    const auto groupIndex = entryIndex / numNodeLists;
    const auto nodeListi = entryIndex % numNodeLists;
    const auto repSeed = size_t(orderedGroupRepresentativeSeedsSpan[groupIndex]);
    const auto count = groupRawCoarseCountsSpan[entryIndex];
    if (count > 0) {
      const auto& descriptor = seedDescriptorsSpan[repSeed*numNodeLists + nodeListi];
      const auto neighborKind = neighborGroupKindsSpan[nodeListi];
      if (neighborKind == TreeNeighborGroupKind) {
        CHECK(descriptor.kind == TreeNeighborGroupKind);
        countOrFillTreeGroupCoarseNeighbors<Dimension>(treeNeighborViewsSpan[nodeListi],
                                                       descriptor.level,
                                                       descriptor.key,
                                                       groupRawCoarseNeighborsSpan.data() + groupRawCoarseOffsetsSpan[entryIndex]);
      } else if (neighborKind == NestedNeighborGroupKind) {
        CHECK(descriptor.kind == NestedNeighborGroupKind);
        countOrFillNestedGroupCoarseNeighbors<Dimension>(nestedNeighborViewsSpan[nodeListi],
                                                         descriptor.level,
                                                         descriptor.key,
                                                         nestedGridCellInfluenceRadiiSpan[nodeListi],
                                                         neighborSearchTypesSpan[nodeListi],
                                                         groupRawCoarseNeighborsSpan.data() + groupRawCoarseOffsetsSpan[entryIndex]);
      }
    }
  });
  GPUUtils::touch(groupRawCoarseNeighborsSpan, chai::CPU);
  if (hasFallbackNeighborKinds) {
    for (auto groupIndex = 0u; groupIndex < numNeighborGroups; ++groupIndex) {
      const auto repSeed = size_t(orderedGroupRepresentativeSeeds[groupIndex]);
      const auto repNodeList = size_t(seedNodeListIDs[repSeed]);
      const auto repNode = seedNodeIDs[repSeed];
      for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
        if (neighborGroupKinds[nodeListi] == FallbackNeighborGroupKind) {
          const auto entryIndex = groupIndex*numNodeLists + nodeListi;
          const auto count = groupRawCoarseCounts[entryIndex];
          if (count > 0) {
            mNodeLists[nodeListi]->neighbor().fillCoarseNeighbors(position(repNodeList, repNode),
                                                                  H(repNodeList, repNode),
                                                                  groupRawCoarseNeighbors.data() + groupRawCoarseOffsets[entryIndex],
                                                                  ghostConnectivity);
          }
        }
      }
    }
  }
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
  care::host_device_ptr<Key> groupCoarseSortKeys(groupCoarseNeighbors.size(), "groupCoarseSortKeys");
#ifndef SPHERAL_UNIFIED_MEMORY
  care::host_device_ptr<int> groupCoarseSortValues(groupCoarseNeighborsSpan);
#else
  care::host_device_ptr<int> groupCoarseSortValues(groupCoarseNeighbors.size(), "groupCoarseSortValues");
#endif
  RAJA::forall<EXEC_POLICY>(TRS_UINT(0u, numNeighborGroups*numNodeLists),
  [=] SPHERAL_HOST_DEVICE (size_t entryIndex) {
    const auto nodeList = entryIndex % numNodeLists;
    const auto beginOffset = size_t(groupCoarseOffsetsSpan[entryIndex]);
    const auto endOffset = size_t(groupCoarseOffsetsSpan[entryIndex + 1u]);
    for (auto coarseIndex = beginOffset; coarseIndex < endOffset; ++coarseIndex) {
      const auto j = groupCoarseNeighborsSpan[coarseIndex];
#ifdef SPHERAL_UNIFIED_MEMORY
      groupCoarseSortValues[coarseIndex] = j;
#endif
      groupCoarseSortKeys[coarseIndex] = (domainDecompIndependent ? keysView(nodeList, j) : Key(j));
    }
  });
  for (auto entryIndex = 0u; entryIndex < numNeighborGroups*numNodeLists; ++entryIndex) {
    const auto beginOffset = size_t(groupCoarseOffsets[entryIndex]);
    const auto sortCount = size_t(groupCoarseOffsets[entryIndex + 1u] - groupCoarseOffsets[entryIndex]);
    if (sortCount > 1u) {
      care::sortKeyValueArrays<SortExecPolicy, Key, int>(groupCoarseSortKeys,
                                                         groupCoarseSortValues,
                                                         beginOffset,
                                                         sortCount);
    }
  }
#ifdef SPHERAL_UNIFIED_MEMORY
  RAJA::forall<EXEC_POLICY>(TRS_UINT(0u, groupCoarseNeighbors.size()),
  [=] SPHERAL_HOST_DEVICE (size_t coarseIndex) {
    groupCoarseNeighborsSpan[coarseIndex] = groupCoarseSortValues[coarseIndex];
  });
  GPUUtils::touch(groupCoarseNeighborsSpan, chai::CPU);
#endif

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
  this->rebuildFlatConnectivityViews();
  const auto connectivityView = this->connectivityFlatView();

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

    RAJA::forall<EXEC_POLICY>(TRS_UINT(0u, numConnectivityEntries),
    [=] SPHERAL_HOST_DEVICE (size_t entryIndex) {
      const auto globalNodeIndex = entryIndex / numNodeLists;
      const auto iNodeList = size_t(globalNodeListIDsSpan[globalNodeIndex]);
      const auto i = globalNodeIDsSpan[globalNodeIndex];
      const auto jNodeList = entryIndex % numNodeLists;
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
    });
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

    RAJA::forall<EXEC_POLICY>(TRS_UINT(0u, numConnectivityEntries),
    [=] SPHERAL_HOST_DEVICE (size_t entryIndex) {
      const auto globalNodeIndex = entryIndex / numNodeLists;
      const auto iNodeList = size_t(globalNodeListIDsSpan[globalNodeIndex]);
      const auto i = globalNodeIDsSpan[globalNodeIndex];
      const auto jNodeList = entryIndex % numNodeLists;
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
    });
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

    mOverlapConnectivityFlatOffsets = std::move(overlapOffsets);
    mOverlapConnectivityFlatNeighbors = std::move(overlapNeighbors);
    mOverlapConnectivity.clear();
    TIME_END("ConnectivityMap_computeOverlapConnectivity");
  } else {
    mOverlapConnectivity.clear();
  }

  this->rebuildFlatConnectivityViews();

  // Are we building intersection connectivity?
  if (mBuildIntersectionConnectivity) {
    TIME_BEGIN("ConnectivityMap_precomputeIntersectionConnectivity");
    auto& pairs = *mNodePairListPtr;
    const auto npairs = pairs.size();
    std::vector<int> intersectionCounts(npairs*numNodeLists, 0);
#pragma omp parallel for
    for (auto k = 0u; k < npairs; ++k) {
      const auto& pair = pairs[k];
      const auto intersection = this->connectivityIntersectionForNodes(pair.i_list, pair.i_node,
                                                                       pair.j_list, pair.j_node,
                                                                       position);
      CHECK(intersection.size() == numNodeLists);
      for (auto nodeList = 0u; nodeList < numNodeLists; ++nodeList) {
        intersectionCounts[k*numNodeLists + nodeList] = intersection[nodeList].size();
      }
    }
    countsToOffsets(intersectionCounts, mIntersectionConnectivityFlatOffsets);
    mIntersectionConnectivityFlatNeighbors.resize(mIntersectionConnectivityFlatOffsets.back());
#pragma omp parallel for
    for (auto k = 0u; k < npairs; ++k) {
      const auto& pair = pairs[k];
      const auto intersection = this->connectivityIntersectionForNodes(pair.i_list, pair.i_node,
                                                                       pair.j_list, pair.j_node,
                                                                       position);
      for (auto nodeList = 0u; nodeList < numNodeLists; ++nodeList) {
        auto dstOffset = size_t(mIntersectionConnectivityFlatOffsets[k*numNodeLists + nodeList]);
        for (const auto neighbor: intersection[nodeList]) {
          mIntersectionConnectivityFlatNeighbors[dstOffset++] = neighbor;
        }
        CHECK(dstOffset == size_t(mIntersectionConnectivityFlatOffsets[k*numNodeLists + nodeList + 1u]));
      }
    }
    TIME_END("ConnectivityMap_precomputeIntersectionConnectivity");
  } else {
    mIntersectionConnectivityFlatOffsets.clear();
    mIntersectionConnectivityFlatNeighbors.clear();
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
  this->rebuildFlatTraversalViews();
  this->rebuildFlatIntersectionConnectivityViews();

  BEGIN_CONTRACT_SCOPE
  ENSURE2(taskNodeIDs.size() == connectivitySize,
          "Missed connectivity tasks: " << taskNodeIDs.size() << " " << connectivitySize);
  // Make sure we're ready to be used.
  ENSURE(valid());
  END_CONTRACT_SCOPE

  mConnectivityCacheValid = false;
  mOverlapConnectivityCacheValid = false;
  mIntersectionConnectivityCacheValid = false;

  TIME_END("ConnectivityMap_computeConnectivity");
}

}
