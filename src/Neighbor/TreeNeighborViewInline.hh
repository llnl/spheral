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
  const LevelKey result = std::clamp(int(mGridLevelConst0 - log(h)*log2inv), 0, int(num1dbits) - 1);
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

}
