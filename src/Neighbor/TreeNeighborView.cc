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

// //------------------------------------------------------------------------------
// // Append elements to a container. Grow resizable containers when necessary;
// // fixed-capacity views report failure so their owning caller can resize and retry.
// //------------------------------------------------------------------------------
// template<typename InContainer,
//          typename OutContainer>
// SPHERAL_HOST_DEVICE
// inline
// bool
// copy_n_conditional(const InContainer& source,
//                    const size_t numSource,
//                    OutContainer& dest,
//                    const size_t sourceOffset = 0u,
//                    const size_t destOffset = 0u) {
//   REQUIRE(sourceOffset <= source.size());
//   REQUIRE(numSource <= source.size() - sourceOffset);
//   if (numSource == 0u) return true;
//   const auto sizeAvail = destOffset < dest.size() ? dest.size() - destOffset : 0u;
//   if (sizeAvail >= numSource) {
//     std::copy_n(source.begin() + sourceOffset, numSource, dest.begin() + destOffset);
//     return true;
//   }

//   if constexpr (requires { dest.resize(std::size_t{0}); }) {
//     dest.resize(destOffset + numSource);
//     std::copy_n(source.begin() + sourceOffset, numSource, dest.begin() + destOffset);
//     return true;
//   }

//   return false;
// }

}  // anonymous

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE 
TreeNeighborView<Dimension>::
TreeNeighborView(const NeighborSearchType searchType,
                 const double kernelExtent,
                 const Vector& xmin,
                 const Vector& xmax):
  NeighborView<Dimension, TreeNeighborView<Dimension>>(searchType, kernelExtent),
  mBoxLength((xmax - xmin).maxElement()),
  mGridLevelConst0(0.0),
  mXmin(xmin),
  mXmax(xmax),
  mFlatTreeView(),
  mFlatTreeDaughtersIndicesView(),
  mFlatTreeMembersView(),
  mFlatTreeLevelOffsetsView(),
  mFirstGhostNodeView() {
  REQUIRE(kernelExtent > 0.0);
  mGridLevelConst0 = std::log((xmax - xmin).maxElement()/kernelExtent)*log2inv;
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
