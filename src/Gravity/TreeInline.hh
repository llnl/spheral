//---------------------------------Spheral++----------------------------------//
// Tree -- a implementation of a quad/oct (2D/3D) tree structure.
//
// Extracted from our original implementation in TreeGravity.
// This class carries along a few data members from that heritage (like the mass
// and velocity) which we retain for convenience. Those attributes can be ignored
// for purely geometrical applications.
//
// Created by JMO, Tue Oct  4 10:17:41 PDT 2022
//----------------------------------------------------------------------------//

namespace Spheral {

//------------------------------------------------------------------------------
// xmin
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const typename Tree<Dimension>::Vector&
Tree<Dimension>::
xmin() const {
  return mXmin;
}

//------------------------------------------------------------------------------
// xmax
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const typename Tree<Dimension>::Vector&
Tree<Dimension>::
xmax() const {
  return mXmax;
}

//------------------------------------------------------------------------------
// Return the central position of a cell
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Tree<Dimension>::Vector
Tree<Dimension>::
cellCenter(const LevelKey& level, const CellKey& key) const {
  CellKey ix, iy, iz;
  extractCellIndices(key, ix, iy, iz);
  const auto dx = cellSize(level);
  return mXmin + dx*Vector(double(ix) + 0.5,
                           double(iy) + 0.5,
                           double(iz) + 0.5);
}

//------------------------------------------------------------------------------
// Return the lower bound coordinate of a cell
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Tree<Dimension>::Vector
Tree<Dimension>::
lowerBound(const LevelKey& level, const CellKey& key) const {
  CellKey ix, iy, iz;
  extractCellIndices(key, ix, iy, iz);
  const auto dx = cellSize(level);
  return mXmin + dx*Vector(double(ix),
                           double(iy),
                           double(iz));
}

//------------------------------------------------------------------------------
// Return the upper bound coordainte of a cell
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Tree<Dimension>::Vector
Tree<Dimension>::
upperBound(const LevelKey& level, const CellKey& key) const {
  CellKey ix, iy, iz;
  extractCellIndices(key, ix, iy, iz);
  const auto dx = cellSize(level);
  return mXmin + dx*Vector(double(ix) + 1.0,
                           double(iy) + 1.0,
                           double(iz) + 1.0);
}

//------------------------------------------------------------------------------
// Return the cell size on a level
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
Tree<Dimension>::
cellSize(const LevelKey& level) const {
  return mBoxLength/(1U << level);
}

//------------------------------------------------------------------------------
// Number of levels in the tree
//------------------------------------------------------------------------------
template<typename Dimension>
inline
size_t
Tree<Dimension>::
numLevels() const {
  return mLevels.size();
}

template<typename Dimension>
inline
void
Tree<Dimension>::
numLevels(const size_t n) {
  mLevels.resize(n);
}

//------------------------------------------------------------------------------
// Get the occupied cells on the given level
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Tree<Dimension>::TreeLevel&
Tree<Dimension>::
operator[](const size_t level) {
  return mLevels[level];
}

template<typename Dimension>
inline
const typename Tree<Dimension>::TreeLevel&
Tree<Dimension>::
operator[](const size_t level) const {
  return mLevels[level];
}

//------------------------------------------------------------------------------
// Build a cell key from coordinate indices.
//------------------------------------------------------------------------------
template<>
inline
void
Tree<Dim<2>>::
buildCellKey(const Tree<Dim<2>>::LevelKey ilevel,
             const Tree<Dim<2>>::Vector& xi,
             Tree<Dim<2>>::CellKey& key,
             Tree<Dim<2>>::CellKey& ix,
             Tree<Dim<2>>::CellKey& iy,
             Tree<Dim<2>>::CellKey& iz) const {
  REQUIRE(xi.x() >= mXmin.x() and xi.x() <= mXmax.x());
  REQUIRE(xi.y() >= mXmin.y() and xi.y() <= mXmax.y());
  const CellKey ncell = (1U << ilevel);
  const CellKey maxcell = ncell - 1U;
  ix = std::min(maxcell, CellKey((xi.x() - mXmin.x())/mBoxLength * ncell));
  iy = std::min(maxcell, CellKey((xi.y() - mXmin.y())/mBoxLength * ncell));
  iz = 0U;
  key = ((std::max(CellKey(0), std::min(max1dKey, iy)) <<   num1dbits) +
         (std::max(CellKey(0), std::min(max1dKey, ix))));
}

template<>
inline
void
Tree<Dim<3>>::
buildCellKey(const Tree<Dim<3>>::LevelKey ilevel,
             const Tree<Dim<3>>::Vector& xi,
             Tree<Dim<3>>::CellKey& key,
             Tree<Dim<3>>::CellKey& ix,
             Tree<Dim<3>>::CellKey& iy,
             Tree<Dim<3>>::CellKey& iz) const {
  REQUIRE(xi.x() >= mXmin.x() and xi.x() <= mXmax.x());
  REQUIRE(xi.y() >= mXmin.y() and xi.y() <= mXmax.y());
  REQUIRE(xi.z() >= mXmin.z() and xi.z() <= mXmax.z());
  const CellKey ncell = (1U << ilevel);
  const CellKey maxcell = ncell - 1U;
  ix = std::min(maxcell, CellKey((xi.x() - mXmin.x())/mBoxLength * ncell));
  iy = std::min(maxcell, CellKey((xi.y() - mXmin.y())/mBoxLength * ncell));
  iz = std::min(maxcell, CellKey((xi.z() - mXmin.z())/mBoxLength * ncell));
  key = ((std::max(CellKey(0), std::min(max1dKey, iz)) << 2*num1dbits) +
         (std::max(CellKey(0), std::min(max1dKey, iy)) <<   num1dbits) +
         (std::max(CellKey(0), std::min(max1dKey, ix))));
}

//------------------------------------------------------------------------------
// Extract the individual coordinate indices from a cell index.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
Tree<Dimension>::
extractCellIndices(const Tree<Dimension>::CellKey& key,
                   Tree<Dimension>::CellKey& ix,
                   Tree<Dimension>::CellKey& iy,
                   Tree<Dimension>::CellKey& iz) const {
  ix = key & xkeymask;
  iy = (key & ykeymask) >> num1dbits;
  iz = (key & zkeymask) >> 2*num1dbits;
}

//------------------------------------------------------------------------------
// Add a daughter to a cell if not present.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
Tree<Dimension>::
addDaughter(Tree<Dimension>::Cell& cell,
            const Tree<Dimension>::CellKey daughterKey) const {
  if (std::find(cell.daughters.begin(), cell.daughters.end(), daughterKey) == cell.daughters.end())
    cell.daughters.push_back(daughterKey);
  ENSURE(cell.daughters.size() <= Dimension::pownu(2));
}

//------------------------------------------------------------------------------
// Add a node to the internal Tree structure.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
Tree<Dimension>::
addNodeToTree(const double mi,
              const Tree<Dimension>::Vector& xi,
              const Tree<Dimension>::Vector& vi) {
  mLevels.reserve(num1dbits); // This is necessary to avoid memory errors!

  LevelKey ilevel = 0;
  bool terminated = false;
  CellKey key, parentKey, otherKey, ix, iy, iz;
  typename TreeLevel::iterator itr;
  while (ilevel < Tree<Dimension>::num1dbits and not terminated) {

    // Do we need to add another level to the tree?
    if (ilevel == mLevels.size()) mLevels.push_back(TreeLevel());

    // Create the key for the cell containing this particle on this level.
    buildCellKey(ilevel, xi, key, ix, iy, iz);
    itr = mLevels[ilevel].find(key);

    if (itr == mLevels[ilevel].end()) {
      // If this is an unregistered cell, add it with this node as the sole leaf
      // and we're done.
      terminated = true;
      mLevels[ilevel][key] = Cell(mi, xi, vi, key);

    } else {
      Cell& cell = itr->second;

      // Is this cell a single leaf already?
      if (cell.masses.size() > 0) {
        CHECK(cell.masses.size() == cell.positions.size() and
              cell.masses.size() == cell.velocities.size());
        CHECK(cell.daughters.size() == 0);

        // Yep, so we need to split it unless we're at the maximum refinement.
        if (ilevel < Tree<Dimension>::num1dbits - 1) {
          CHECK(cell.masses.size() == 1);
          const LevelKey ilevel1 = ilevel + 1;
          CHECK(ilevel1 < Tree<Dimension>::num1dbits);
          if (ilevel1 == mLevels.size()) mLevels.push_back(TreeLevel());
          buildCellKey(ilevel1, cell.xcm, otherKey, ix, iy, iz);
          mLevels[ilevel1][otherKey] = Cell(cell.M, cell.xcm, cell.vcm, otherKey);
          cell.daughters = std::vector<CellKey>(1, otherKey);
          cell.masses = std::vector<double>();
          cell.positions = std::vector<Vector>();
          cell.velocities = std::vector<Vector>();

        } else {
          // If we've maxed out the levels, then we just huck this node in 
          // the members of this cell.
          cell.masses.push_back(mi);
          cell.positions.push_back(xi);
          cell.velocities.push_back(vi);
        }
      }

      // Increment the cell moments.
      cell.xcm = (cell.M*cell.xcm + mi*xi)/(cell.M + mi);
      cell.vcm = (cell.M*cell.vcm + mi*vi)/(cell.M + mi);
      cell.M += mi;
      cell.Mglobal = cell.M;
    }

    // Link this cell as a daughter of its parent.
    if (ilevel > 0) {
      CHECK(mLevels[ilevel - 1].find(parentKey) != mLevels[ilevel - 1].end());
      addDaughter(mLevels[ilevel - 1][parentKey], key);
    }

    parentKey = key;
    ++ilevel;
  }
}

//------------------------------------------------------------------------------
// Geometry only
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
Tree<Dimension>::
addNodeToTree(const Tree<Dimension>::Vector& xi) {
  this->addNodeToTree(0.0, xi, Vector::zero());
}

//------------------------------------------------------------------------------
// Construct the daughter pointers in a tree.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
Tree<Dimension>::
constructDaughterPtrs() {
  const unsigned nlevels = mLevels.size();
  const unsigned n = nlevels > 0 ? nlevels - 1 : nlevels;
  unsigned ilevel, ilevel1;
  for (ilevel = 0; ilevel != n; ++ilevel) {
    ilevel1 = ilevel + 1;
    for (typename TreeLevel::iterator itr = mLevels[ilevel].begin();
         itr != mLevels[ilevel].end();
         ++itr) {
      Cell& cell = itr->second;
      cell.daughterPtrs = std::vector<Cell*>();
      for (std::vector<CellKey>::const_iterator ditr = cell.daughters.begin();
           ditr != cell.daughters.end();
           ++ditr) {
        cell.daughterPtrs.push_back(&(mLevels[ilevel1][*ditr]));
      }
      CHECK(cell.daughters.size() == cell.daughterPtrs.size());
    }
  }
}

}
