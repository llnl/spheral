//---------------------------------Spheral++----------------------------------//
// TreeNeighborView
//
// An SPH neighbor finder based on an oct-tree like structure with specialized
// cell membership criteria.  This View object can be used on CPU or GPU devices.
//
// Based on the algorithm originally described in 
// Owen, Villumsen, Shapiro, & Martel 1998, ApJS, 116, 155
//
// Created by J. Michael Owen, Fri Mar 23 15:50:35 PDT 2012
// Split into TreeNeighborView/TreeNeighbor by JMO Tue Sep 29 10:26:17 PDT 2026
//----------------------------------------------------------------------------//
#ifndef __Spheral_TreeNeighborView__
#define __Spheral_TreeNeighborView__

#include "Neighbor/NeighborView.hh"

#include <algorithm>
#include <cmath>
#include <span>
#include <cstdint>
#include <numbers>

namespace Spheral {

template<typename Dimension>
class TreeNeighborView: public NeighborView<Dimension, TreeNeighborView<Dimension>> {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  using DimensionType = Dimension;
  using BaseType = NeighborView<Dimension, TreeNeighborView<Dimension>>;
  using NodeID = typename BaseType::NodeID;
  template<typename T> using ContainerView = typename BaseType::template ContainerView<T>;

  // Data types we use to build the internal tree structure.
  using LevelKey = uint32_t;
  using CellKey = uint64_t;

  // Constructors and destructors
  SPHERAL_HOST_DEVICE TreeNeighborView(const NeighborSearchType searchType,
                                       const double kernelExtent,
                                       const Vector& xmin,
                                       const Vector& xmax);
  SPHERAL_HOST_DEVICE TreeNeighborView(const TreeNeighborView&) = default;
  SPHERAL_HOST_DEVICE TreeNeighborView(TreeNeighborView&&) = default;
  SPHERAL_HOST_DEVICE ~TreeNeighborView() = default;

  // -------------------- NeighborViews must provide these methods --------------------
  // Set or refine the neighbor lists for the given position and smoothing scale.
  template<typename Container>
  SPHERAL_HOST_DEVICE std::tuple<bool, size_t, size_t> setMasterList(const Vector& position,
                                                                     const SymTensor& H,
                                                                     Container& masterList,
                                                                     Container& coarseNeighbors,
                                                                     const bool ghostConnectivity = false) const;

  template<typename Container>
  SPHERAL_HOST_DEVICE size_t setRefineNeighborList(const Vector& position,
                                                   const SymTensor& H,
                                                   const size_t numCoarse,
                                                   const Container& coarseNeighbors,
                                                   Container& refineNeighbors) const;

  // Set Neighbors for the given position.
  template<typename Container>
  SPHERAL_HOST_DEVICE std::tuple<bool, size_t, size_t> setMasterList(const Vector& position,
                                                                     Container& masterList,
                                                                     Container& coarseNeighbors,
                                                                     const bool ghostConnectivity = false) const;
  template<typename Container>
  SPHERAL_HOST_DEVICE size_t setRefineNeighborList(const Vector& position,
                                                   const size_t numCoarse,
                                                   const Container& coarseNeighbors,
                                                   Container& refineNeighbors) const;

  // Set the neighbor lists based on proximity to planes.
  template<typename Container>
  SPHERAL_HOST_DEVICE std::tuple<bool, size_t, size_t> setMasterList(const GeomPlane<Dimension>& enterPlane,
                                                                     const GeomPlane<Dimension>& exitPlane,
                                                                     Container& masterList,
                                                                     Container& coarseNeighbors) const;
  // -----------------------------------------------------------------------------------

  // Compute the grid level appropriate for the given smoothing scale.
  SPHERAL_HOST_DEVICE LevelKey gridLevel(const Scalar& h) const;        // units of length
  SPHERAL_HOST_DEVICE LevelKey gridLevel(const SymTensor& H) const;     // units of 1/length

  // Access internal attributes.
  SPHERAL_HOST_DEVICE const Vector& xmin()                    const { return mXmin; }
  SPHERAL_HOST_DEVICE const Vector& xmax()                    const { return mXmax; }
  SPHERAL_HOST_DEVICE Scalar boxLength()                      const { return mBoxLength; }

  // Cell size on the given level
  SPHERAL_HOST_DEVICE Scalar cellSize(const LevelKey levelID) const { return mBoxLength/(1u << levelID); }

  // How many grid levels in the current tree?
  SPHERAL_HOST_DEVICE LevelKey numGridLevels()                const { return mFlatTreeLevelOffsetsView.size() > 0u ? LevelKey(mFlatTreeLevelOffsetsView.size() - 1u) : 0u; }

  // Find the nearest cell center on the given grid level.
  SPHERAL_HOST_DEVICE Vector nearestCellCenter(const Vector& xi, const Scalar& hi) const;
  SPHERAL_HOST_DEVICE Vector nearestCellCenter(const Vector& xi, const SymTensor& Hi) const;

  // Test if the given position and extent is occupied, in the sense that
  // any tree cell is occupied at any coarser level down to it's native
  // level.
  SPHERAL_HOST_DEVICE bool occupied(const Vector& xi, const Scalar& hi) const;
  SPHERAL_HOST_DEVICE bool occupied(const Vector& xi, const SymTensor& Hi) const;

  // For our parallel algorithm it is useful to be able to set the master/coarse
  // information based on the given (level, cell).
  template<typename Container>
  SPHERAL_HOST_DEVICE std::tuple<bool, size_t, size_t> setTreeMasterList(const LevelKey levelID,
                                                                         const CellKey cellID,
                                                                         Container& masterList,
                                                                         Container& coarseNeighbors,
                                                                         const bool ghostConnectivity) const;

  // Expose hidden Neighbor methods
  using BaseType::setMasterList;
  using BaseType::setRefineNeighborList;

protected:
  //--------------------------- Protected Interface ---------------------------//
  static constexpr unsigned num1dbits = 21U;                  // The number of bits we quantize 1D coordinates to.  We have to fit three of these in 64 bits.
  static constexpr CellKey max1dKey = 1U << num1dbits;        // The maximum number of cells this corresponds to in a direction.
  static constexpr CellKey xkeymask = (1U << num1dbits) - 1U; // Bit masks we can use to extract the coordinate specific indices from a cell key.
  static constexpr CellKey ykeymask = xkeymask << num1dbits;
  static constexpr CellKey zkeymask = ykeymask << num1dbits;                                 
  static constexpr double log2inv = std::numbers::log2e;      // 1.0/std::log(2.0);

  //----------------------------------------------------------------------------
  // FlatCell holds the properties of cells in the tree.
  // This is an immutable version of the TreeCell for use after the Tree is
  // fully constructed for use on the GPU when walking the Tree.
  //----------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  // The FlatCell value type, which we build the FlatTree from for reference
  // in the View object
  //--------------------------------------------------------------------------
  struct FlatCell {
    CellKey key;
    size_t daughterOffset, numDaughters;
    size_t memberOffset, numMembers;

    // Constructors
    SPHERAL_HOST_DEVICE FlatCell(): key(0), daughterOffset(0), numDaughters(0), memberOffset(0), numMembers(0) {}
    SPHERAL_HOST_DEVICE FlatCell(const CellKey key0,
                                 const size_t daughterOff,
                                 const size_t nDaughters,
                                 const size_t memberOff,
                                 const size_t nMembers):
      key(key0),
      daughterOffset(daughterOff),
      numDaughters(nDaughters),
      memberOffset(memberOff),
      numMembers(nMembers) {}
    SPHERAL_HOST_DEVICE FlatCell(const FlatCell&) = default;
    SPHERAL_HOST_DEVICE FlatCell(const FlatCell&&) = default;
    SPHERAL_HOST_DEVICE FlatCell& operator=(const FlatCell&) = default;

    // Throw in comparison operators for help sorting.
    SPHERAL_HOST_DEVICE bool operator==(const FlatCell& rhs) const { return key == rhs.key; }
    SPHERAL_HOST_DEVICE bool operator<(const FlatCell& rhs) const { return key < rhs.key; }
  };

  // Find the index of the given (level, cellKey) if present. Returns mFlatTreeView.size() if not found.
  SPHERAL_HOST_DEVICE size_t findCell(const LevelKey ilevel, const CellKey key) const;

  // Build a cell key based on a level and position.
  SPHERAL_HOST_DEVICE void buildCellKey(const LevelKey ilevel,
                                        const Vector& xi,
                                        CellKey& key,
                                        CellKey& ix,
                                        CellKey& iy,
                                        CellKey& iz) const;

  // Extract the individual coordinate indices from a cell key.
  SPHERAL_HOST_DEVICE void extractCellIndices(const CellKey& key,
                                              CellKey& ix,
                                              CellKey& iy,
                                              CellKey& iz) const;

  // Actual method for setting the master list.
  template<typename Container>
  SPHERAL_HOST_DEVICE std::tuple<bool, size_t, size_t> setTreeMasterList(const Vector& position,
                                                                         const Scalar& h,
                                                                         Container& masterList,
                                                                         Container& coarseNeighbors,
                                                                         const bool ghostConnectivity) const;

  // Set the refine neighbors
  template<typename Container>
  SPHERAL_HOST_DEVICE size_t setTreeRefineNeighborList(const Vector& position,
                                                       const SymTensor& H,
                                                       const size_t numCoarse,
                                                       const Container& coarseNeighbors,
                                                       Container& refineNeighbors) const;
  // Return all the members in range of the specified cell.
  template<typename Container>
  SPHERAL_HOST_DEVICE std::tuple<bool, size_t> findTreeNeighbors(const LevelKey masterLevel,
                                                                 const CellKey ix_master,
                                                                 const CellKey iy_master,
                                                                 const CellKey iz_master,
                                                                 Container& neighbors,
                                                                 size_t numNeighbors = 0u) const;

  // Shift a 1D key from one level to another.
  SPHERAL_HOST_DEVICE CellKey shiftKeyLevel(const CellKey ix, const LevelKey level0, const LevelKey level1) const;

  // Determine if a full cell key is in the set of 1D key ranges.
  SPHERAL_HOST_DEVICE bool keyInRange(const CellKey key, 
                                      const CellKey ix_min, const CellKey iy_min, const CellKey iz_min,
                                      const CellKey ix_max, const CellKey iy_max, const CellKey iz_max) const;

  // Find the minimum distance from a cell/level to a plane.
  SPHERAL_HOST_DEVICE Scalar distanceToCell(const LevelKey ilevel, const CellKey key, const GeomPlane<Dimension>& plane) const;

  // Map a cell key/level through an entrance/exit plane pair, returning the set of cells
  // that overlap the mapped cell.
  SPHERAL_HOST_DEVICE size_t mapKey(const LevelKey ilevel, const CellKey key, 
                                    const GeomPlane<Dimension>& enterPlane,
                                    const GeomPlane<Dimension>& exitPlane,
                                    std::span<CellKey> cells) const;

  // Protected data.
  Scalar mBoxLength, mGridLevelConst0;
  Vector mXmin, mXmax;
  
  // The FlatTreeView data
  ContainerView<FlatCell> mFlatTreeView;
  ContainerView<size_t> mFlatTreeDaughtersIndicesView;
  ContainerView<NodeID> mFlatTreeMembersView;
  ContainerView<size_t> mFlatTreeLevelOffsetsView;
  ContainerView<size_t> mFirstGhostNodeView;

  using BaseType::mNodeExtentView;
  using BaseType::mSearchType;
  using BaseType::mKernelExtent;
};

}

#include "Neighbor/TreeNeighborViewInline.hh"

#endif
