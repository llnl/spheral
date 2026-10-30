//---------------------------------Spheral++----------------------------------//
// TreeNeighborView
//
// An SPH neighbor finder based on an oct-tree like structure with specialized
// cell membership criteria.
//
// Based on the algorithm originally described in 
// Owen, Villumsen, Shapiro, & Martel 1998, ApJS, 116, 155
//
// Created by J. Michael Owen, Fri Mar 23 15:50:35 PDT 2012
// Split into TreeNeighborView/TreeNeighbor by JMO Tue Sep 29 10:26:17 PDT 2026
//----------------------------------------------------------------------------//
#include "Neighbor/TreeNeighborView.hh"
// #include "NodeList/NodeListView.hh"
// #include "Utilities/globalBoundingVolumes.hh"
#include "Utilities/FastMath.hh"
#include "Boundary/mapPositionThroughPlanes.hh"
#include "Geometry/Dimension.hh"
#include "Geometry/GeomPlane.hh"
#include "Utilities/DBC.hh"

#include <array>
#include <algorithm>
#include <cmath>

namespace Spheral {

namespace { // anonymous

//------------------------------------------------------------------------------
// Compute the vertex coordinates for a cell.
//------------------------------------------------------------------------------
// 1D
inline
SPHERAL_HOST_DEVICE 
std::array<Dim<1>::Vector, 2>
findCellVertices(const Dim<1>::Vector& xmin,
                 const double boxLength,
                 const uint32_t ilevel,
                 const uint64_t ix,
                 const uint64_t /*iy*/,
                 const uint64_t /*iz*/) {
  using Vector = Dim<1>::Vector;
  const auto cellSize = boxLength/(1U << ilevel);
  return std::array<Vector, 2>
  {{
      xmin + Vector(ix       *cellSize),
      xmin + Vector((ix + 1U)*cellSize)
  }};
}

// 2D
inline
SPHERAL_HOST_DEVICE
std::array<Dim<2>::Vector, 4>
findCellVertices(const Dim<2>::Vector& xmin,
                 const double boxLength,
                 const uint32_t ilevel,
                 const uint64_t ix,
                 const uint64_t iy,
                 const uint64_t /*iz*/) {
  using Vector = Dim<2>::Vector;
  const auto cellSize = boxLength/(1U << ilevel);
  return std::array<Vector, 4>
  {{
      xmin + Vector(ix       *cellSize, iy       *cellSize),
      xmin + Vector((ix + 1U)*cellSize, iy       *cellSize),
      xmin + Vector((ix + 1U)*cellSize, (iy + 1U)*cellSize),
      xmin + Vector(ix       *cellSize, (iy + 1U)*cellSize)
  }};
}

// 3D
inline
SPHERAL_HOST_DEVICE 
std::array<Dim<3>::Vector, 8>
findCellVertices(const Dim<3>::Vector& xmin,
                 const double boxLength,
                 const uint32_t ilevel,
                 const uint64_t ix,
                 const uint64_t iy,
                 const uint64_t iz) {
  using Vector = Dim<3>::Vector;
  const auto cellSize = boxLength/(1U << ilevel);
  return std::array<Vector, 8>
  {{
      xmin + Vector(ix       *cellSize, iy       *cellSize, iz       *cellSize),
      xmin + Vector((ix + 1U)*cellSize, iy       *cellSize, iz       *cellSize),
      xmin + Vector((ix + 1U)*cellSize, (iy + 1U)*cellSize, iz       *cellSize),
      xmin + Vector(ix       *cellSize, (iy + 1U)*cellSize, iz       *cellSize),
      xmin + Vector(ix       *cellSize, iy       *cellSize, (iz + 1U)*cellSize),
      xmin + Vector((ix + 1U)*cellSize, iy       *cellSize, (iz + 1U)*cellSize),
      xmin + Vector((ix + 1U)*cellSize, (iy + 1U)*cellSize, (iz + 1U)*cellSize),
      xmin + Vector(ix       *cellSize, (iy + 1U)*cellSize, (iz + 1U)*cellSize)
  }};
}

}  // anonymous

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE 
TreeNeighborView<Dimension>::
TreeNeighborView(DataBaseView<Dimension> dataBase,
                 const NeighborSearchType searchType,
                 const double kernelExtent,
                 const Vector& xmin,
                 const Vector& xmax):
  NeighborView<Dimension, TreeNeighborView<Dimension>>(dataBase, searchType, kernelExtent),
  mBoxLength((xmax - xmin).maxElement()),
  mGridLevelConst0(std::log((xmax - xmin).maxElement()/kernelExtent)*log2inv),
  mXmin(xmin),
  mXmax(xmax),
  mFlatTreeView(),
  mFlatTreeDaughtersIndicesView(),
  mFlatTreeMembersView(),
  mFlatTreeLevelOffsetsView() {
}

//------------------------------------------------------------------------------
// Set the master list for the given position & H
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE 
std::tuple<bool, size_t, size_t>
TreeNeighborView<Dimension>::
setMasterList(const Vector& position,
              const SymTensor& H,
              Container<NodeID>& masterList,
              Container<NodeID>& coarseNeighbors,
              const bool ghostConnectivity) const {
  REQUIRE(H.Determinant() > 0.0);
  const auto hinvValues = H.eigenValues();
  CHECK(hinvValues.minElement() > 0.0);
  const auto h = 1.0/hinvValues.minElement();
  return this->setTreeMasterList(position, h, masterList, coarseNeighbors, ghostConnectivity);
}

template<typename Dimension>
SPHERAL_HOST_DEVICE 
std::tuple<bool, size_t, size_t>
TreeNeighborView<Dimension>::
setMasterList(const Vector& position,
              Container<NodeID>& masterList,
              Container<NodeID>& coarseNeighbors,
              const bool ghostConnectivity) const {
  return this->setTreeMasterList(position, 1.0e-30*mBoxLength, masterList, coarseNeighbors, ghostConnectivity);
}

//------------------------------------------------------------------------------
// Set the refine neighbor list for position & H
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE 
size_t
TreeNeighborView<Dimension>::
setRefineNeighborList(const Vector& position,
                      const SymTensor& H,
                      const size_t numCoarse,
                      const Container<NodeID>& coarseNeighbors,
                      Container<NodeID>& refineNeighbors) const {
  REQUIRE(H.Determinant() > 0.0);
  return this->setTreeRefineNeighborList(position, H, numCoarse, coarseNeighbors, refineNeighbors);
}

template<typename Dimension>
SPHERAL_HOST_DEVICE 
size_t
TreeNeighborView<Dimension>::
setRefineNeighborList(const Vector& position,
                      const size_t numCoarse,
                      const Container<NodeID>& coarseNeighbors,
                      Container<NodeID>& refineNeighbors) const {
  return this->setTreeRefineNeighborList(position, 1.0e30*mBoxLength*SymTensor::one(), numCoarse, coarseNeighbors, refineNeighbors);
}

//------------------------------------------------------------------------------
// Set the master list based on proximity to planes.
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE 
std::tuple<bool, size_t, size_t>
TreeNeighborView<Dimension>::
setMasterList(const GeomPlane<Dimension>& enterPlane,
              const GeomPlane<Dimension>& exitPlane,
              Container<NodeID>& masterList,
              Container<NodeID>& coarseNeighbors) const {

  // Maximum sizes
  const auto maxMaster = masterList.size();
  const auto maxCoarse = coarseNeighbors.size();
  constexpr size_t maxDaughters = Dimension::pownu(4u);  // Down 2 levels

  // Get the master and coarse lists.
  bool ok = true;
  size_t numMaster = 0u, numCoarse = 0u, numNewDaughters;
  std::array<size_t, maxDaughters> remainingDaughters, newDaughters;
  std::array<CellKey, maxDaughters> mappedKeys;
  if (mFlatTreeView.size() > 0) {
    CHECK(mFlatTreeLevelOffsetsView[0] == 0u and mFlatTreeLevelOffsetsView[1] == 1u);    // Only one root cell
    CHECK(mFlatTreeView[0].numMembers == 0u);                                            // No points are members of the root cell

    // Declare a bunch of variables we're going to need.
    LevelKey ilevel = 0u;
    CellKey ix, iy, iz;
    Scalar cellSize;
    size_t numRemainingDaughters = mFlatTreeView[0].numDaughters;
    CHECK(numRemainingDaughters <= maxDaughters);
    std::copy_n(mFlatTreeDaughtersIndicesView.begin() + mFlatTreeView[0].daughterOffset, numRemainingDaughters, remainingDaughters.begin());

    // Walk the tree, looking for any master cells that are in range of the entrance plane.
    while (numRemainingDaughters > 0u) {
      numNewDaughters = 0u;
      ++ilevel;
      cellSize = this->cellSize(ilevel);
    
      // Walk the candidates.
      for (auto candidateID = 0u; candidateID < numRemainingDaughters; ++candidateID) {
        const auto daughterIndex = remainingDaughters[candidateID];
        CHECK(daughterIndex < mFlatTreeView.size());
        const auto& cell = mFlatTreeView[daughterIndex];

        // Check if we're in range of either plane.
        const bool entranceCheck = (this->distanceToCell(ilevel, cell.key, enterPlane) <= cellSize);
        const bool exitCheck = (this->distanceToCell(ilevel, cell.key, exitPlane) <= cellSize);

        // If so, add the daughters to check on the next pass.
        if (entranceCheck or exitCheck) {
          CHECK(numNewDaughters + cell.numDaughters <= maxDaughters);
          std::copy_n(mFlatTreeDaughtersIndicesView.begin() + cell.daughterOffset, cell.numDaughters, newDaughters.begin() + numNewDaughters);
          numNewDaughters += cell.numDaughters;
        }

        // Does this cell have members in range of the entrance plane?
        if (entranceCheck and cell.numMembers > 0u) {
          ok &= numMaster + cell.numMembers <= maxMaster;
          if (ok) std::copy_n(mFlatTreeMembersView.begin() + cell.memberOffset, cell.numMembers, masterList.begin() + numMaster);
          numMaster += cell.numMembers;

          // Map the cell key through to the exit plane (which may result in more
          // than one equivalent key).  Then find all the neighbors for those cells,
          // and add them to the coarse set.
          auto numMapped = this->mapKey(ilevel, cell.key, enterPlane, exitPlane, mappedKeys);
          for (auto k = 0u; k < numMapped; ++k) {
            this->extractCellIndices(mappedKeys[k], ix, iy, iz);
            const auto coarseOffset = std::min(numCoarse, maxCoarse);
            auto [neighborsOK, numNeighbors] = this->findTreeNeighbors(ilevel, ix, iy, iz,
                                                                       std::span<NodeID>(coarseNeighbors.begin() + coarseOffset,
                                                                                         maxCoarse - coarseOffset));
            ok &= neighborsOK;
            numCoarse += numNeighbors;
            CHECK(numCoarse <= maxCoarse or not ok);
          }
        }

        // Does this cell have members in range of the exit plane?
        if (exitCheck and cell.numMembers > 0u) {
          ok &= numCoarse + cell.numMembers <= maxCoarse;
          if (ok) std::copy_n(mFlatTreeMembersView.begin() + cell.memberOffset, cell.numMembers, coarseNeighbors.begin() + numCoarse);
          numCoarse += cell.numMembers;

          // Map the cell key through to the entrance plane.  Any nodes we interact
          // with on that side are potential masters.
          auto numMapped = this->mapKey(ilevel, cell.key, exitPlane, enterPlane, mappedKeys);
          for (auto k = 0u; k < numMapped; ++k) {
            this->extractCellIndices(mappedKeys[k], ix, iy, iz);
            const auto masterOffset = std::min(numMaster, maxMaster);
            auto [neighborsOK, numNeighbors] = this->findTreeNeighbors(ilevel, ix, iy, iz,
                                                                       std::span<NodeID>(masterList.begin() + masterOffset,
                                                                                         maxMaster - masterOffset));
            ok &= neighborsOK;
            numMaster += numNeighbors;
            CHECK(numMaster <= maxMaster or not ok);
          }
        }
      }

      // Update the set of daughters to check on the next pass.
      numRemainingDaughters = numNewDaughters;
      remainingDaughters = newDaughters;
    }

    // Remove duplicates from the master and coarse sets.
    if (ok) {
      std::sort(masterList.begin(), masterList.begin() + numMaster);
      numMaster = std::distance(masterList.begin(), std::unique(masterList.begin(), masterList.begin() + numMaster));
      std::sort(coarseNeighbors.begin(), coarseNeighbors.begin() + numCoarse);
      numCoarse = std::distance(coarseNeighbors.begin(), std::unique(coarseNeighbors.begin(), coarseNeighbors.begin() + numCoarse));
    }
  }
  return std::make_tuple(ok, numMaster, numCoarse);
}

//------------------------------------------------------------------------------
// Do any existing points occupy cells overlapping the given position and
// extent?
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE 
bool
TreeNeighborView<Dimension>::
occupied(const Vector& xi,
         const SymTensor& Hi) const {
  if (mFlatTreeView.size() == 0) return false;
  const auto gridLeveli = std::min(gridLevel(Hi), numGridLevels() - 1u);
  CellKey keyi, ix, iy, iz;
  for (auto ilevel = 0u; ilevel < gridLeveli; ++ilevel) {
    buildCellKey(ilevel, xi, keyi, ix, iy, iz);
    auto k = findCell(ilevel, keyi);
    if (k != mFlatTreeView.size() and
        mFlatTreeView[k].numMembers > 0u) return true;
  }
  buildCellKey(gridLeveli, xi, keyi, ix, iy, iz);
  auto k = findCell(gridLeveli, keyi);
  return k < mFlatTreeView.size();
}

template<typename Dimension>
SPHERAL_HOST_DEVICE 
bool
TreeNeighborView<Dimension>::
occupied(const Vector& xi,
         const Scalar& hi) const {
  return this->occupied(xi, SymTensor::one() * 1.0/hi);
}

//------------------------------------------------------------------------------
// Set the master & coarse neighbor sets for a given tree level and cell ID.
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE 
std::tuple<bool, size_t, size_t>
TreeNeighborView<Dimension>::
setTreeMasterList(const LevelKey levelID,
                  const CellKey cellID,
                  Container<NodeID>& masterList,
                  Container<NodeID>& coarseNeighbors,
                  const bool ghostConnectivity) const {
  REQUIRE(levelID < num1dbits);
  const auto maxMaster = masterList.size();

  // Get the per dimension cell indices.
  CellKey ix_master, iy_master, iz_master;
  extractCellIndices(cellID, ix_master, iy_master, iz_master);

  // Prep the result
  bool ok = true;
  size_t numMaster = 0u;
  size_t numCoarse = 0u;

  // Set the master list.
  if (levelID < numGridLevels()) {
    auto k = findCell(levelID, cellID);
    if (k < mFlatTreeView.size()) {
      const auto& cell = mFlatTreeView[k];
      numMaster = cell.numMembers;
      ok &= numMaster <= maxMaster;
      if (ok) std::copy_n(mFlatTreeMembersView.begin() + cell.memberOffset, cell.numMembers, masterList.begin());
    }
  }

  // Set the coarse list.
  if (mFlatTreeView.size() > 0) {
    bool ok_coarse;
    std::tie(ok_coarse, numCoarse) = findTreeNeighbors(levelID,
                                                       ix_master,
                                                       iy_master,
                                                       iz_master,
                                                       std::span<NodeID>(coarseNeighbors.begin(), coarseNeighbors.size()));
    ok &= ok_coarse;
  }

  // Remove all ghost nodes from the master list.
  if (numMaster <= maxMaster) {
    std::sort(masterList.begin(), masterList.begin() + numMaster);
  }
  if (numMaster <= maxMaster and not ghostConnectivity) {
    const auto& nodeLists = mDataBaseView.nodeListViews();
    size_t writeIndex = 0u;
    for (size_t readIndex = 0u; readIndex < numMaster; ++readIndex) {
      const auto& nodeID = masterList[readIndex];
      CHECK(nodeID.first < nodeLists.size());
      if (nodeID.second < nodeLists[nodeID.first].firstGhostNode()) {
        masterList[writeIndex++] = nodeID;
      }
    }
    numMaster = writeIndex;
  }

  // Post conditions.
  ENSURE2(numMaster <= numCoarse or not ok, "Strange sizes: (numMaster, numCoarse) = " << numMaster << " " << numCoarse);
  return std::make_tuple(ok, numMaster, numCoarse);
}

//------------------------------------------------------------------------------
// Find the index of a FlatCell in the tree given (level, key)
// Return mFlatTreeView.size() on failure
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE
size_t
TreeNeighborView<Dimension>::
findCell(const LevelKey levelID, const CellKey cellID) const {
  if (levelID + 1u >= mFlatTreeLevelOffsetsView.size()) return mFlatTreeView.size();
  const auto begin = mFlatTreeLevelOffsetsView[levelID];
  const auto end = mFlatTreeLevelOffsetsView[levelID + 1u];
  size_t low = begin, high = end;
  while (low < high) {
    const auto mid = low + (high - low)/2u;
    if (mFlatTreeView[mid].key < cellID) {
      low = mid + 1u;
    } else {
      high = mid;
    }
  }
  if (low < end and mFlatTreeView[low].key == cellID) {
    return low;
  } else {
    return mFlatTreeView.size();
  }
}

//------------------------------------------------------------------------------
// Set the master & coarse neighbor sets by walking the tree.
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE 
std::tuple<bool, size_t, size_t>
TreeNeighborView<Dimension>::
setTreeMasterList(const Vector& position,
                  const Scalar& h,
                  Container<NodeID>& masterList,
                  Container<NodeID>& coarseNeighbors,
                  const bool ghostConnectivity) const {

  // Set the working master grid level and cell.
  CellKey masterKey, ix_master, iy_master, iz_master;
  const LevelKey masterLevel = gridLevel(h);
  buildCellKey(masterLevel, position, masterKey, ix_master, iy_master, iz_master);
  CHECK(masterLevel >= 0 and masterLevel < num1dbits);

  // We can just use the method based on IDs at this point.
  return setTreeMasterList(masterLevel, masterKey, masterList, coarseNeighbors, ghostConnectivity);
}

//------------------------------------------------------------------------------
// Produce the refined list of potential neighbors for a single node.
//------------------------------------------------------------------------------
template<typename Dimension>
size_t
SPHERAL_HOST_DEVICE 
TreeNeighborView<Dimension>::
setTreeRefineNeighborList(const Vector& /*position*/,
                          const SymTensor& /*H*/,
                          const size_t numCoarse,
                          const Container<NodeID>& coarseNeighbors,
                          Container<NodeID>& refineNeighbors) const {

  // // Determine the maximum extent of this H tensor in each dimension.
  // const Vector extent = this->HExtent(H, this->kernelExtent());
  // const Vector minExtent = position - extent;
  // const Vector maxExtent = position + extent;

  // Use precull to set the refined neighbor list.
  REQUIRE(numCoarse <= coarseNeighbors.size());
  REQUIRE(numCoarse <= refineNeighbors.size());
  std::copy_n(coarseNeighbors.begin(), numCoarse, refineNeighbors.begin()); // this->precullList(position, position, minExtent, maxExtent, coarseNeighbors);
  return numCoarse;
}

//------------------------------------------------------------------------------
// Walk the tree to find the neighbors for the given (gridlevel, gridcell) 
// pair.
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE 
std::tuple<bool, size_t>
TreeNeighborView<Dimension>::
findTreeNeighbors(const LevelKey masterLevel,
                  const CellKey ix_master,
                  const CellKey iy_master,
                  const CellKey iz_master,
                  std::span<NodeID> neighbors) const {

  // Declare variables.
  const auto maxNeighbors = neighbors.size();
  constexpr size_t maxDaughters = Dimension::pownu(4u);  // Down 2 levels
  LevelKey ilevel = 0;
  CellKey ix, iy, iz, ix_min, iy_min, iz_min, ix_max, iy_max, iz_max, delta;
  size_t numRemainingDaughters, numNewDaughters;
  std::array<size_t, maxDaughters> remainingDaughters, newDaughters;

  // Walk the tree until we run out of daughters to check.
  bool ok = true;
  size_t numNeighbors = 0u;
  if (mFlatTreeView.size() > 0u) {
    CHECK2(mFlatTreeView[0].numMembers == 0u, "TreeNeighbor root cell occupied!  Will miss neighbors... " << mFlatTreeView[0].numMembers);
    numRemainingDaughters = mFlatTreeView[0].numDaughters;
    CHECK(numRemainingDaughters <= maxDaughters);
    std::copy_n(mFlatTreeDaughtersIndicesView.begin() + mFlatTreeView[0].daughterOffset, numRemainingDaughters, remainingDaughters.begin());
                
    while (numRemainingDaughters > 0u) {
      numNewDaughters = 0u;
      ++ilevel;
      delta = (ilevel <= masterLevel ? 1U : (1U << (ilevel - masterLevel)));

      // Find the target range of keys on this level.
      ix = shiftKeyLevel(ix_master, masterLevel, ilevel);
      iy = shiftKeyLevel(iy_master, masterLevel, ilevel);
      iz = shiftKeyLevel(iz_master, masterLevel, ilevel);
      ix_min = (ix > delta              ? ix - delta : 0U);
      iy_min = (iy > delta              ? iy - delta : 0U);
      iz_min = (iz > delta              ? iz - delta : 0U);
      ix_max = ((max1dKey - ix) > delta ? ix + 2*delta - 1U : max1dKey);
      iy_max = ((max1dKey - iy) > delta ? iy + 2*delta - 1U : max1dKey);
      iz_max = ((max1dKey - iz) > delta ? iz + 2*delta - 1U : max1dKey);
      CHECK(ix_min <= ix_max and ix_max <= max1dKey);
      CHECK(iy_min <= iy_max and iy_max <= max1dKey);
      CHECK(iz_min <= iz_max and iz_max <= max1dKey);
    
      // Walk the candidate daughters on this level.
      for (auto candidateID = 0u; candidateID < numRemainingDaughters; ++candidateID) {
        const auto cellIndex = remainingDaughters[candidateID];
        CHECK(cellIndex < mFlatTreeView.size());
        const auto& cell = mFlatTreeView[cellIndex];

        // Is this daughter in range?
        if (keyInRange(cell.key, ix_min, iy_min, iz_min, ix_max, iy_max, iz_max)) {
        
          // Copy this cells members to the result.
          ok &= numNeighbors + cell.numMembers <= maxNeighbors;
          if (ok) std::copy_n(mFlatTreeMembersView.begin() + cell.memberOffset, cell.numMembers, neighbors.begin() + numNeighbors);
          numNeighbors += cell.numMembers;

          // Add any daughters of this cell to our candidates to check on the next level.
          CHECK(numNewDaughters + cell.numDaughters <= maxDaughters);
          std::copy_n(mFlatTreeDaughtersIndicesView.begin() + cell.daughterOffset, cell.numDaughters, newDaughters.begin() + numNewDaughters);
          numNewDaughters += cell.numDaughters;
        }
      }

      // Update the daughters to check on the next pass.
      numRemainingDaughters = numNewDaughters;
      remainingDaughters = newDaughters;
    }
  }

  // That's it.
  return std::make_tuple(ok, numNeighbors);
}

//------------------------------------------------------------------------------
// Find the minimum distance from a cell to a plane.
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE 
typename Dimension::Scalar
TreeNeighborView<Dimension>::
distanceToCell(const LevelKey ilevel,
               const CellKey key,
               const GeomPlane<Dimension>& plane) const {
  CellKey ix, iy, iz;
  extractCellIndices(key, ix, iy, iz);
  const auto cellVertices = findCellVertices(mXmin, mBoxLength, 
                                             ilevel, ix, iy, iz);
  CHECK(cellVertices.size() == (1U << Dimension::nDim));
  double dist, minDist = plane.signedDistance(cellVertices[0]);
  for (auto i = 1u; i < cellVertices.size(); ++i) {
    dist = plane.signedDistance(cellVertices[i]);
    if (dist*minDist < 0.0) {
      minDist = 0.0;
    } else if (std::abs(dist) < std::abs(minDist)) {
      minDist = dist;
    }
  }
  return std::abs(minDist);
}

//------------------------------------------------------------------------------
// Map a cell key/level through an entrance/exit plane pair, returning the set 
// of cells that overlap the mapped cell.
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE 
size_t
TreeNeighborView<Dimension>::
mapKey(const LevelKey ilevel,
       const CellKey key,
       const GeomPlane<Dimension>& enterPlane,
       const GeomPlane<Dimension>& exitPlane,
       std::span<CellKey> cells) const {
  CellKey ix, iy, iz;
  extractCellIndices(key, ix, iy, iz);
  const auto vertices = findCellVertices(mXmin, mBoxLength,
                                         ilevel, ix, iy, iz);
  constexpr unsigned n = (1U << Dimension::nDim);
  CHECK(vertices.size() == n);

  // Find the range of (ix,iy,iz) the mapped vertices cover.
  CellKey ixmin = max1dKey, 
          iymin = max1dKey, 
          izmin = max1dKey,
          ixmax = CellKey(0),
          iymax = CellKey(0),
          izmax = CellKey(0),
          newKey;
  for (auto i = 0u; i < n; ++i) {
    buildCellKey(ilevel,
                 mapPositionThroughPlanes(vertices[i], enterPlane, exitPlane),
                 newKey, ix, iy, iz);
    ixmin = std::min(ixmin, ix);
    iymin = std::min(iymin, iy);
    izmin = std::min(izmin, iz);
    ixmax = std::max(ixmax, ix);
    iymax = std::max(iymax, iy);
    izmax = std::max(izmax, iz);
  }
  CHECK((ixmax - ixmin + 1)*(iymax - iymin + 1)*(izmax - izmin + 1) <= cells.size());

  // Now fill in the set of cells we map to.
  size_t result = 0u;
  for (ix = ixmin; ix <= ixmax; ++ix) {
    for (iy = iymin; iy <= iymax; ++iy) {
      for (iz = izmin; iz <= izmax; ++iz) {
        cells[result++] = ((std::max(CellKey(0), std::min(max1dKey, iz)) << 2*num1dbits) +
                           (std::max(CellKey(0), std::min(max1dKey, iy)) <<   num1dbits) +
                           (std::max(CellKey(0), std::min(max1dKey, ix))));
      }
    }
  }

  // That's it.
  return result;
}

}
