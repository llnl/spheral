#include "Utilities/copy_n_conditional.hh"

#include <array>

namespace Spheral {

//------------------------------------------------------------------------------
// Find the grid level for the given smoothing scale.
//------------------------------------------------------------------------------
// Argument units of length
template<typename Dimension>
SPHERAL_HOST_DEVICE 
inline
typename TreeNeighborView<Dimension>::LevelKey
TreeNeighborView<Dimension>::
gridLevel(const Scalar& h) const {   
  // REQUIRE2(this->kernelExtent()*h <= mBoxLength,
  //          "h larger than box size: " << this->kernelExtent()*h << " " << mBoxLength);
  const LevelKey result = std::clamp(int(mGridLevelConst0 - std::log(h)*log2inv), 0, int(num1dbits) - 1);
  ENSURE(result < num1dbits);
  ENSURE2((result == 0) or (result == num1dbits - 1) or
          (fuzzyLessThanOrEqual(h*this->kernelExtent(), mBoxLength/(1U << result), 1.0e-10) and
           fuzzyGreaterThanOrEqual(h*this->kernelExtent(), mBoxLength/(1U << (result + 1U)), 1.0e-10)),
          result << " " << h*this->kernelExtent() << " in? ["
          << mBoxLength/(1U << (result + 1U)) << " " 
          << mBoxLength/(1U << result) << "]\n"
          << mBoxLength << " " << mGridLevelConst0);
  return result;
}

// Argument units of 1/length
template<typename Dimension>
SPHERAL_HOST_DEVICE 
inline
typename TreeNeighborView<Dimension>::LevelKey
TreeNeighborView<Dimension>::
gridLevel(const SymTensor& H) const {   
  REQUIRE(H.Determinant() > 0.0);
  return this->gridLevel(1.0/H.eigenValues().minElement());
}

//------------------------------------------------------------------------------
// Find the nearest cell center on the given grid level.
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE 
inline
typename TreeNeighborView<Dimension>::Vector
TreeNeighborView<Dimension>::
nearestCellCenter(const Vector& xi,
                  const SymTensor& Hi) const {
  const auto leveli = gridLevel(Hi);
  CellKey keyi, ix, iy, iz;
  buildCellKey(leveli, xi, keyi, ix, iy, iz);
  return mXmin + cellSize(leveli)*Vector(double(ix) + 0.5,
                                         double(iy) + 0.5,
                                         double(iz) + 0.5);
}

template<typename Dimension>
SPHERAL_HOST_DEVICE 
inline
typename TreeNeighborView<Dimension>::Vector
TreeNeighborView<Dimension>::
nearestCellCenter(const Vector& xi,
                  const Scalar& hi) const {
  return this->nearestCellCenter(xi, SymTensor::one() * 1.0/hi);
}

//------------------------------------------------------------------------------
// Build a cell key from coordinate indices.
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE 
inline
void
TreeNeighborView<Dimension>::
buildCellKey(const LevelKey ilevel,
             const Vector& xi,
             CellKey& key,
             CellKey& ix,
             CellKey& iy,
             CellKey& iz) const {
  const CellKey ncell = (1U << ilevel);
  const CellKey maxcell = ncell - 1U;
  ix = std::min(maxcell, CellKey(std::clamp((xi.x() - mXmin.x())/mBoxLength, 0.0, 1.0) * ncell));
  iy = std::min(maxcell, CellKey(std::clamp((xi.y() - mXmin.y())/mBoxLength, 0.0, 1.0) * ncell));
  iz = std::min(maxcell, CellKey(std::clamp((xi.z() - mXmin.z())/mBoxLength, 0.0, 1.0) * ncell));
  key = ((std::max(CellKey(0), std::min(max1dKey, iz)) << 2*num1dbits) +
         (std::max(CellKey(0), std::min(max1dKey, iy)) <<   num1dbits) +
         (std::max(CellKey(0), std::min(max1dKey, ix))));
}

//------------------------------------------------------------------------------
// Extract the individual coordinate indices from a cell index.
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE 
inline
void
TreeNeighborView<Dimension>::
extractCellIndices(const CellKey& key,
                   CellKey& ix,
                   CellKey& iy,
                   CellKey& iz) const {
  ix = key & xkeymask;
  iy = (key & ykeymask) >> num1dbits;
  iz = (key & zkeymask) >> 2*num1dbits;
}

//------------------------------------------------------------------------------
// Shift a 1D key from one level to another.
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE 
inline
typename TreeNeighborView<Dimension>::CellKey
TreeNeighborView<Dimension>::
shiftKeyLevel(const CellKey ix,
              const LevelKey level0,
              const LevelKey level1) const {
  return (level1 <= level0 ? 
          ix >> (level0 - level1) :
          ix << (level1 - level0));
}

//------------------------------------------------------------------------------
// Check if the given key is in a specified range.
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE 
inline
bool
TreeNeighborView<Dimension>::
keyInRange(const CellKey key,
           const CellKey ix_min,
           const CellKey iy_min,
           const CellKey iz_min,
           const CellKey ix_max,
           const CellKey iy_max,
           const CellKey iz_max) const {
  CellKey ix, iy, iz;
  extractCellIndices(key, ix, iy, iz);
  return (ix >= ix_min and ix <= ix_max and
          iy >= iy_min and iy <= iy_max and
          iz >= iz_min and iz <= iz_max);
}

//------------------------------------------------------------------------------
// Set the master list for the given position & H
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename Container>
SPHERAL_HOST_DEVICE 
inline
std::tuple<bool, size_t, size_t>
TreeNeighborView<Dimension>::
setMasterList(const Vector& position,
              const SymTensor& H,
              Container& masterList,
              Container& coarseNeighbors,
              const bool ghostConnectivity) const {
  REQUIRE(H.Determinant() > 0.0);
  const auto hinvValues = H.eigenValues();
  CHECK(hinvValues.minElement() > 0.0);
  const auto h = 1.0/hinvValues.minElement();
  return this->setTreeMasterList(position, h, masterList, coarseNeighbors, ghostConnectivity);
}

template<typename Dimension>
template<typename Container>
SPHERAL_HOST_DEVICE 
inline
std::tuple<bool, size_t, size_t>
TreeNeighborView<Dimension>::
setMasterList(const Vector& position,
              Container& masterList,
              Container& coarseNeighbors,
              const bool ghostConnectivity) const {
  return this->setTreeMasterList(position, 1.0e-30*mBoxLength, masterList, coarseNeighbors, ghostConnectivity);
}

//------------------------------------------------------------------------------
// Set the master & coarse neighbor sets by walking the tree.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename Container>
SPHERAL_HOST_DEVICE 
inline
std::tuple<bool, size_t, size_t>
TreeNeighborView<Dimension>::
setTreeMasterList(const Vector& position,
                  const Scalar& h,
                  Container& masterList,
                  Container& coarseNeighbors,
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
// Set the master & coarse neighbor sets for a given tree level and cell ID.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename Container>
SPHERAL_HOST_DEVICE 
inline
std::tuple<bool, size_t, size_t>
TreeNeighborView<Dimension>::
setTreeMasterList(const LevelKey levelID,
                  const CellKey cellID,
                  Container& masterList,
                  Container& coarseNeighbors,
                  const bool ghostConnectivity) const {
  REQUIRE(levelID < num1dbits);

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
      ok &= copy_n_conditional(mFlatTreeMembersView, cell.numMembers, masterList, cell.memberOffset);
    }
  }

  // Set the coarse list.
  if (mFlatTreeView.size() > 0) {
    bool ok_coarse;
    std::tie(ok_coarse, numCoarse) = findTreeNeighbors(levelID,
                                                       ix_master,
                                                       iy_master,
                                                       iz_master,
                                                       coarseNeighbors);
    ok &= ok_coarse;
  }

  // Remove all ghost nodes from the master list.
  if (ok) {
    std::sort(masterList.begin(), masterList.begin() + numMaster);
    if (not ghostConnectivity) {
      size_t writeIndex = 0u;
      for (size_t readIndex = 0u; readIndex < numMaster; ++readIndex) {
        const auto& nodeID = masterList[readIndex];
        CHECK(nodeID.first < mFirstGhostNodeView.size());
        if (nodeID.second < mFirstGhostNodeView[nodeID.first]) {
          masterList[writeIndex++] = nodeID;
        }
      }
      numMaster = writeIndex;
    }
  }

  // Post conditions.
  ENSURE2(numMaster <= numCoarse or not ok, "Strange sizes: (numMaster, numCoarse) = " << numMaster << " " << numCoarse);
  return std::make_tuple(ok, numMaster, numCoarse);
}

//------------------------------------------------------------------------------
// Set the refine neighbor list for position & H
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename Container>
SPHERAL_HOST_DEVICE 
inline
size_t
TreeNeighborView<Dimension>::
setRefineNeighborList(const Vector& position,
                      const SymTensor& H,
                      const size_t numCoarse,
                      const Container& coarseNeighbors,
                      Container& refineNeighbors) const {
  REQUIRE(H.Determinant() > 0.0);
  return this->setTreeRefineNeighborList(position, H, numCoarse, coarseNeighbors, refineNeighbors);
}

template<typename Dimension>
template<typename Container>
SPHERAL_HOST_DEVICE 
inline
size_t
TreeNeighborView<Dimension>::
setRefineNeighborList(const Vector& position,
                      const size_t numCoarse,
                      const Container& coarseNeighbors,
                      Container& refineNeighbors) const {
  return this->setTreeRefineNeighborList(position, 1.0e30*mBoxLength*SymTensor::one(), numCoarse, coarseNeighbors, refineNeighbors);
}

//------------------------------------------------------------------------------
// Set the master list based on proximity to planes.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename Container>
SPHERAL_HOST_DEVICE 
inline
std::tuple<bool, size_t, size_t>
TreeNeighborView<Dimension>::
setMasterList(const GeomPlane<Dimension>& enterPlane,
              const GeomPlane<Dimension>& exitPlane,
              Container& masterList,
              Container& coarseNeighbors) const {

  // A depth-first walk has bounded candidate storage based only on the tree's
  // maximum branching factor and depth.  The mapped-cell buffer has a separate
  // bound: mapping one cell through the planes can overlap at most four cells
  // per dimension.
  constexpr size_t maxDaughters = Dimension::pownu(2u);
  constexpr size_t maxStackSize = 1u + (maxDaughters - 1u)*num1dbits;
  constexpr size_t maxMappedCells = Dimension::pownu(4u);

  // Get the master and coarse lists.
  bool ok = true;
  size_t numMaster = 0u, numCoarse = 0u, stackSize = 0u;
  std::array<size_t, maxStackSize> candidateCells;
  std::array<LevelKey, maxStackSize> candidateLevels;
  std::array<CellKey, maxMappedCells> mappedKeys;
  if (mFlatTreeView.size() > 0) {
    CHECK(mFlatTreeLevelOffsetsView[0] == 0u and mFlatTreeLevelOffsetsView[1] == 1u);    // Only one root cell
    CHECK(mFlatTreeView[0].numMembers == 0u);                                            // No points are members of the root cell

    // Seed the stack with the root daughters in reverse order so popping the
    // stack preserves the original daughter traversal order.
    const auto& root = mFlatTreeView[0];
    CHECK(root.numDaughters <= maxDaughters);
    for (auto daughterID = root.numDaughters; daughterID > 0u; --daughterID) {
      candidateCells[stackSize] = mFlatTreeDaughtersIndicesView[root.daughterOffset + daughterID - 1u];
      candidateLevels[stackSize] = 1u;
      ++stackSize;
    }

    LevelKey ilevel;
    CellKey ix, iy, iz;
    while (stackSize > 0u) {
      --stackSize;
      const auto cellIndex = candidateCells[stackSize];
      ilevel = candidateLevels[stackSize];
      CHECK(ilevel < num1dbits);
      CHECK(cellIndex < mFlatTreeView.size());
      const auto& cell = mFlatTreeView[cellIndex];
      const auto cellSize = this->cellSize(ilevel);

      // Check if we're in range of either plane.
      const bool entranceCheck = (this->distanceToCell(ilevel, cell.key, enterPlane) <= cellSize);
      const bool exitCheck = (this->distanceToCell(ilevel, cell.key, exitPlane) <= cellSize);

      // If so, add the daughters to the depth-first stack.
      if (entranceCheck or exitCheck) {
        CHECK(cell.numDaughters <= maxDaughters);
        CHECK(stackSize + cell.numDaughters <= maxStackSize);
        CHECK(cell.numDaughters == 0u or ilevel + 1u < num1dbits);
        for (auto daughterID = cell.numDaughters; daughterID > 0u; --daughterID) {
          candidateCells[stackSize] = mFlatTreeDaughtersIndicesView[cell.daughterOffset + daughterID - 1u];
          candidateLevels[stackSize] = ilevel + 1u;
          ++stackSize;
        }
      }

      // Does this cell have members in range of the entrance plane?
      if (entranceCheck and cell.numMembers > 0u) {
        ok &= copy_n_conditional(mFlatTreeMembersView, cell.numMembers, masterList, cell.memberOffset, numMaster);
        numMaster += cell.numMembers;

        // Map the cell key through to the exit plane (which may result in more
        // than one equivalent key).  Then find all the neighbors for those cells,
        // and add them to the coarse set.
        const auto numMapped = this->mapKey(ilevel, cell.key, enterPlane, exitPlane, mappedKeys);
        for (auto k = 0u; k < numMapped; ++k) {
          this->extractCellIndices(mappedKeys[k], ix, iy, iz);
          auto [neighborsOK, numNeighbors] = this->findTreeNeighbors(ilevel, ix, iy, iz,
                                                                     coarseNeighbors, numCoarse);
          ok &= neighborsOK;
          numCoarse = numNeighbors;
        }
      }

      // Does this cell have members in range of the exit plane?
      if (exitCheck and cell.numMembers > 0u) {
        ok &= copy_n_conditional(mFlatTreeMembersView, cell.numMembers, coarseNeighbors, cell.memberOffset, numCoarse);
        numCoarse += cell.numMembers;

        // Map the cell key through to the entrance plane.  Any nodes we interact
        // with on that side are potential masters.
        const auto numMapped = this->mapKey(ilevel, cell.key, exitPlane, enterPlane, mappedKeys);
        for (auto k = 0u; k < numMapped; ++k) {
          this->extractCellIndices(mappedKeys[k], ix, iy, iz);
          auto [neighborsOK, numNeighbors] = this->findTreeNeighbors(ilevel, ix, iy, iz,
                                                                     masterList, numMaster);
          ok &= neighborsOK;
          numMaster = numNeighbors;
        }
      }
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
// Produce the refined list of potential neighbors for a single node.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename Container>
SPHERAL_HOST_DEVICE 
size_t
inline
TreeNeighborView<Dimension>::
setTreeRefineNeighborList(const Vector& /*position*/,
                          const SymTensor& /*H*/,
                          const size_t numCoarse,
                          const Container& coarseNeighbors,
                          Container& refineNeighbors) const {

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
template<typename Container>
SPHERAL_HOST_DEVICE 
inline
std::tuple<bool, size_t>
TreeNeighborView<Dimension>::
findTreeNeighbors(const LevelKey masterLevel,
                  const CellKey ix_master,
                  const CellKey iy_master,
                  const CellKey iz_master,
                  Container& neighbors,
                  size_t numNeighbors) const {

  // Declare variables.  Walk depth-first so the candidate storage is bounded
  // by the maximum branching factor and tree depth, rather than by the number
  // of occupied cells on any one level.
  constexpr size_t maxDaughters = Dimension::pownu(2u);
  constexpr size_t maxStackSize = 1u + (maxDaughters - 1u)*num1dbits;
  LevelKey ilevel;
  CellKey ix, iy, iz, ix_min, iy_min, iz_min, ix_max, iy_max, iz_max, delta;
  size_t stackSize = 0u;
  std::array<size_t, maxStackSize> candidateCells;
  std::array<LevelKey, maxStackSize> candidateLevels;

  // Walk the tree until we run out of daughters to check.
  bool ok = true;
  if (mFlatTreeView.size() > 0u) {
    CHECK2(mFlatTreeView[0].numMembers == 0u, "TreeNeighbor root cell occupied!  Will miss neighbors... " << mFlatTreeView[0].numMembers);
    const auto& root = mFlatTreeView[0];
    CHECK(root.numDaughters <= maxDaughters);

    // Push in reverse order so popping the stack preserves the original
    // daughter traversal order.
    for (auto daughterID = root.numDaughters; daughterID > 0u; --daughterID) {
      candidateCells[stackSize] = mFlatTreeDaughtersIndicesView[root.daughterOffset + daughterID - 1u];
      candidateLevels[stackSize] = 1u;
      ++stackSize;
    }

    while (stackSize > 0u) {
      --stackSize;
      const auto cellIndex = candidateCells[stackSize];
      ilevel = candidateLevels[stackSize];
      CHECK(ilevel < num1dbits);
      CHECK(cellIndex < mFlatTreeView.size());
      const auto& cell = mFlatTreeView[cellIndex];

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

      // If this cell is in range, collect its members and descend into its
      // daughters.  A depth-first stack is bounded by
      // 1 + (branching factor - 1)*maximum depth.
      if (keyInRange(cell.key, ix_min, iy_min, iz_min, ix_max, iy_max, iz_max)) {
        ok &= copy_n_conditional(mFlatTreeMembersView, cell.numMembers, neighbors, cell.memberOffset, numNeighbors);
        numNeighbors += cell.numMembers;

        CHECK(cell.numDaughters <= maxDaughters);
        CHECK(stackSize + cell.numDaughters <= maxStackSize);
        CHECK(cell.numDaughters == 0u or ilevel + 1u < num1dbits);
        for (auto daughterID = cell.numDaughters; daughterID > 0u; --daughterID) {
          candidateCells[stackSize] = mFlatTreeDaughtersIndicesView[cell.daughterOffset + daughterID - 1u];
          candidateLevels[stackSize] = ilevel + 1u;
          ++stackSize;
        }
      }
    }
  }

  // That's it.
  return std::make_tuple(ok, numNeighbors);
}

}
