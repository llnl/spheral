namespace Spheral {

//------------------------------------------------------------------------------
// Get specialized hash functions
//------------------------------------------------------------------------------
// For an array, we are hashing a normal, and the final element is either 0 or 1
template<>
inline
int
NormalHash<std::array<int, 1>>::
operator()(const std::array<int, 1>& v) const {
  return v[0];
}
template<>
inline
int
NormalHash<std::array<int, 2>>::
operator()(const std::array<int, 2>& v) const {
  return (v[0] << 1) ^ v[1];
}
template<>
inline
int
NormalHash<std::array<int, 3>>::
operator()(const std::array<int, 3>& v) const {
  return (v[0] << std::numeric_limits<int>::digits / 2) ^ (v[1] << 1) ^ v[2];
}

//------------------------------------------------------------------------------
// Return whether things have been initialized
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
FlatConnectivity<Dimension>::
indexingInitialized() const {
  return mIndexingInitialized;
}

template<typename Dimension>
inline
bool
FlatConnectivity<Dimension>::
overlapIndexingInitialized() const {
  return mOverlapIndexingInitialized;
}

template<typename Dimension>
inline
bool
FlatConnectivity<Dimension>::
globalIndexingInitialized() const {
  return mGlobalIndexingInitialized;
}

template<typename Dimension>
inline
bool
FlatConnectivity<Dimension>::
surfaceIndexingInitialized() const {
  return mSurfaceIndexingInitialized;
}

template<typename Dimension>
inline
bool
FlatConnectivity<Dimension>::
boundaryInformationInitialized() const {
  return mBoundaryInformationInitialized;
}

//------------------------------------------------------------------------------
// Get global indices
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
FlatConnectivity<Dimension>::
firstGlobalIndex() const {
  CHECK(mGlobalIndexingInitialized);
  return mFirstGlobalIndex;
}

template<typename Dimension>
inline
int
FlatConnectivity<Dimension>::
lastGlobalIndex() const {
  CHECK(mGlobalIndexingInitialized);
  return mLastGlobalIndex;
}

//------------------------------------------------------------------------------
// Get number of elements
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
FlatConnectivity<Dimension>::
numNodes() const {
  CHECK(mIndexingInitialized);
  return mNumLocalNodes;
}

template<typename Dimension>
inline
int
FlatConnectivity<Dimension>::
numInternalNodes() const {
  CHECK(mIndexingInitialized);
  return mNumInternalLocalNodes;
}

template<typename Dimension>
inline
int
FlatConnectivity<Dimension>::
numGlobalNodes() const {
  CHECK(mGlobalIndexingInitialized);
  return mNumGlobalNodes;
}

template<typename Dimension>
inline
int
FlatConnectivity<Dimension>::
numBoundaryNodes() const {
  CHECK(mBoundaryInformationInitialized);
  return mNumBoundaryNodes;
}

template<typename Dimension>
inline
int
FlatConnectivity<Dimension>::
numGlobalBoundaryNodes() const {
  CHECK(mBoundaryInformationInitialized);
  return mNumGlobalBoundaryNodes;
}

//------------------------------------------------------------------------------
// Convert between local and node indices
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
FlatConnectivity<Dimension>::
nodeToLocal(const int nodeListi, const int nodei) const {
  CHECK(mIndexingInitialized);
  CHECK(size_t(nodeListi) < mNodeToLocalIndex.size());
  CHECK(size_t(nodei) < mNodeToLocalIndex[nodeListi].size());
  return mNodeToLocalIndex[nodeListi][nodei];
}

template<typename Dimension>
inline
std::pair<int, int>
FlatConnectivity<Dimension>::
localToNode(const int locali) const {
  CHECK(mIndexingInitialized);
  CHECK(locali < mNumLocalNodes);
  CHECK(locali < int(mLocalToNodeIndex.size()));
  return mLocalToNodeIndex[locali];
}

//------------------------------------------------------------------------------
// Get the global index for this local index
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
FlatConnectivity<Dimension>::
localToGlobal(const int locali) const {
  CHECK(mGlobalIndexingInitialized);
  CHECK(locali < mNumLocalNodes);
  CHECK(size_t(locali) < mLocalToGlobalIndex.size());
  return mLocalToGlobalIndex[locali];
}

//------------------------------------------------------------------------------
// Get number of neighbors, including the point itself
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
FlatConnectivity<Dimension>::
numNeighbors(const int locali) const {
  CHECK(mIndexingInitialized);
  CHECK(locali < mNumConnectivityNodes);
  CHECK(locali < int(mNumNeighbors.size()));
  return mNumNeighbors[locali];
}

template<typename Dimension>
inline
int
FlatConnectivity<Dimension>::
numOverlapNeighbors(const int locali) const {
  CHECK(mOverlapIndexingInitialized);
  CHECK(locali < mNumConnectivityNodes);
  CHECK(size_t(locali) < mNumOverlapNeighbors.size());
  return mNumOverlapNeighbors[locali];
}

template<typename Dimension>
inline
int
FlatConnectivity<Dimension>::
numConstNeighbors(const int locali) const {
  CHECK(mIndexingInitialized);
  CHECK(mBoundaryInformationInitialized);
  CHECK(locali < mNumConnectivityNodes);
  CHECK(size_t(locali) < mNumConstantBoundaryNeighbors.size());
  return mNumConstantBoundaryNeighbors[locali];
}

template<typename Dimension>
inline
int
FlatConnectivity<Dimension>::
numConstOverlapNeighbors(const int locali) const {
  CHECK(mOverlapIndexingInitialized);
  CHECK(mBoundaryInformationInitialized);
  CHECK(locali < mNumConnectivityNodes);
  CHECK(size_t(locali) < mNumConstantBoundaryOverlapNeighbors.size());
  return mNumConstantBoundaryOverlapNeighbors[locali];
}

template<typename Dimension>
inline
int
FlatConnectivity<Dimension>::
numNonConstNeighbors(const int locali) const {
  CHECK(mIndexingInitialized);
  CHECK(mBoundaryInformationInitialized);
  CHECK(locali < mNumConnectivityNodes);
  CHECK(size_t(locali) < mNumNeighbors.size());
  CHECK(size_t(locali) < mNumConstantBoundaryNeighbors.size());
  return mNumNeighbors[locali] - mNumConstantBoundaryNeighbors[locali];
}

template<typename Dimension>
inline
int
FlatConnectivity<Dimension>::
numNonConstOverlapNeighbors(const int locali) const {
  CHECK(mOverlapIndexingInitialized);
  CHECK(mBoundaryInformationInitialized);
  CHECK(locali < mNumConnectivityNodes);
  CHECK(size_t(locali) < mNumOverlapNeighbors.size());
  CHECK(size_t(locali) < mNumConstantBoundaryOverlapNeighbors.size());
  return mNumOverlapNeighbors[locali] - mNumConstantBoundaryOverlapNeighbors[locali];
}

template<typename Dimension>
inline
int
FlatConnectivity<Dimension>::
totalNumNonConstNeighbors() const {
  CHECK(mIndexingInitialized);
  CHECK(mBoundaryInformationInitialized);
  auto total = 0;
  for (auto locali = 0; locali < mNumInternalLocalNodes; ++locali)
  {
    total += mNumNeighbors[locali] - mNumConstantBoundaryNeighbors[locali];
  }
  return total;
}

//------------------------------------------------------------------------------
// For the point i, for its neighbor j, get the flattened index for point j
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
FlatConnectivity<Dimension>::
localToFlat(const int locali,  const int localj) const {
  CHECK(mIndexingInitialized);
  CHECK(locali < mNumConnectivityNodes);
  CHECK(size_t(locali) < mLocalToFlatIndex.size());
  const auto& localMap = mLocalToFlatIndex[locali];
  const auto it = localMap.find(localj);
  if (it == localMap.end()) {
    return NoConnectivity;
  }
  else {
    return it->second;
  }
}

template<typename Dimension>
inline
int
FlatConnectivity<Dimension>::
localToFlatOverlap(const int locali,  const int localj) const {
  CHECK(mOverlapIndexingInitialized);
  CHECK(locali < mNumConnectivityNodes);
  CHECK(size_t(locali) < mLocalToFlatOverlapIndex.size());
  const auto& localMap = mLocalToFlatOverlapIndex[locali];
  const auto it = localMap.find(localj);
  if (it == localMap.end()) {
    return NoConnectivity;
  }
  else {
    return it->second;
  }
}

//------------------------------------------------------------------------------
// For the point i, get the node index for the flattened point j
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
FlatConnectivity<Dimension>::
flatToLocal(const int locali,  const int flatj) const {
  CHECK(mIndexingInitialized);
  CHECK(locali < mNumConnectivityNodes);
  CHECK(flatj < mNumNeighbors[locali]);
  CHECK(locali < int(mFlatToLocalIndex.size()));
  CHECK(flatj < int(mFlatToLocalIndex[locali].size()));
  return mFlatToLocalIndex[locali][flatj];
}

template<typename Dimension>
inline
int
FlatConnectivity<Dimension>::
flatOverlapToLocal(const int locali,  const int flatj) const {
  CHECK(mIndexingInitialized);
  CHECK(locali < mNumConnectivityNodes);
  CHECK(flatj < mNumOverlapNeighbors[locali]);
  CHECK(size_t(locali) < mFlatOverlapToLocalIndex.size());
  CHECK(size_t(flatj) < mFlatOverlapToLocalIndex[locali].size());
  return mFlatOverlapToLocalIndex[locali][flatj];
}

//------------------------------------------------------------------------------
// For the point i, get the node index for the flattened point j
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
FlatConnectivity<Dimension>::
isConstantBoundaryNode(const int locali) const {
  CHECK(mBoundaryInformationInitialized);
  CHECK(locali < mNumLocalNodes);
  CHECK(size_t(locali) < mIsConstantBoundaryNode.size());
  return mIsConstantBoundaryNode[locali];
}

//------------------------------------------------------------------------------
// Get the number of unique surfaces included in this point's neighbors
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
FlatConnectivity<Dimension>::
numSurfaces(const int locali) const {
  CHECK(mSurfaceIndexingInitialized);
  CHECK(locali < mNumLocalNodes);
  CHECK(locali < int(mSurfaceNormal.size()));
  return mSurfaceNormal[locali].size();
}

//------------------------------------------------------------------------------
// Get the unique surface index from the normal direction
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
FlatConnectivity<Dimension>::
surfaceIndex(const int locali,
             const Vector& normal) const {
  CHECK(fuzzyEqual(normal.magnitude(), 1.0));
  ArrayDim array;
  normalToArray(normal, array); 
  return surfaceIndex(locali, array);
}

template<typename Dimension>
inline
int
FlatConnectivity<Dimension>::
surfaceIndex(const int locali,
             const ArrayDim& values) const {
  CHECK(mSurfaceIndexingInitialized);
  CHECK(locali <  mNumLocalNodes);
  CHECK(size_t(locali) < mSurfaceFlatIndex.size());
  const auto& localMap = mSurfaceFlatIndex[locali];
  const auto it = localMap.find(values);
  if (it == localMap.end()) {
    return NoConnectivity;
  }
  else {
    return it->second;
  }
}

//------------------------------------------------------------------------------
// Get the unique surface index from the normal direction
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const typename FlatConnectivity<Dimension>::Vector&
FlatConnectivity<Dimension>::
normal(const int locali, const int flats) const {
  CHECK(mSurfaceIndexingInitialized);
  CHECK(locali < mNumLocalNodes);
  CHECK(locali < int(mSurfaceNormal.size()));
  CHECK(flats < int(mSurfaceNormal[locali].size()));
  return mSurfaceNormal[locali][flats];
}

//------------------------------------------------------------------------------
// Given the Voronoi cell for point i, which surfaces have voids on the other side?
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
FlatConnectivity<Dimension>::
numSurfacesForCell(const int locali) const {
  CHECK(mSurfaceIndexingInitialized);
  CHECK(locali <  mNumLocalNodes);
  CHECK(size_t(locali) < mVoidSurfaces.size());
  return mVoidSurfaces[locali].size();
}

//------------------------------------------------------------------------------
// Given the Voronoi cell for point i, which surfaces have voids on the other side?
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
FlatConnectivity<Dimension>::
surfaceIndexForCell(const int locali, const int flats) const {
  CHECK(mSurfaceIndexingInitialized);
  CHECK(locali <  mNumLocalNodes);
  CHECK(size_t(locali) < mVoidSurfaces.size());
  CHECK(size_t(flats) < mVoidSurfaces[locali].size());
  return mVoidSurfaces[locali][flats];
}

//------------------------------------------------------------------------------
// Convert from the normal to rounded integer values
// We can convert one dimension to only sign dependence due to normalization
//------------------------------------------------------------------------------
template<>
inline
void
FlatConnectivity<Dim<1>>::
normalToArray(const Vector& normal,
              ArrayDim& values) const {
  CHECK(fuzzyEqual(normal.magnitude(), 1.0));
  values[0] = normal[0] > 0;
}

template<>
inline
void
FlatConnectivity<Dim<2>>::
normalToArray(const Vector& normal,
              ArrayDim& values) const {
  CHECK(fuzzyEqual(normal.magnitude(), 1.0));
  values[0] = static_cast<int>(normal[0] * mRoundValue + 0.5);
  values[1] = static_cast<int>(normal[1] * mRoundValue + 0.5 > 0);
}

template<>
inline
void
FlatConnectivity<Dim<3>>::
normalToArray(const Vector& normal,
              ArrayDim& values) const {
  CHECK(fuzzyEqual(normal.magnitude(), 1.0));
  values[0] = static_cast<int>(normal[0] * mRoundValue + 0.5);
  values[1] = static_cast<int>(normal[1] * mRoundValue + 0.5);
  values[2] = static_cast<int>(normal[2] * mRoundValue + 0.5 > 0);
}

//------------------------------------------------------------------------------
// Convert from the rounded integer values to the normal
//------------------------------------------------------------------------------
template<>
inline
void
FlatConnectivity<Dim<1>>::
arrayToNormal(const ArrayDim& values,
              Vector& normal) const {
  CHECK(values[0] == 0 || values[0] == 1);
  normal[0] = values[0] ? 1 : -1;
  ENSURE(-1 <= normal[0] && normal[0] <= 1);
}

template<>
inline
void
FlatConnectivity<Dim<2>>::
arrayToNormal(const ArrayDim& values,
              Vector& normal) const {
  normal[0] = values[0] * mRoundValueInv;
  const auto normal1Mag = sqrt(1 - normal[0] * normal[0]);
  const auto normal1Sign = values[0] ? 1 : -1;
  normal[1] = normal1Mag * normal1Sign;
  ENSURE(fuzzyEqual(normal.magnitude(), 1.0));
}

template<>
inline
void
FlatConnectivity<Dim<3>>::
arrayToNormal(const ArrayDim& values,
              Vector& normal) const {
  normal[0] = values[0] * mRoundValueInv;
  normal[1] = values[1] * mRoundValueInv;
  const auto normal2Mag = sqrt(1 - normal[0] * normal[0] - normal[1] * normal[1]);
  const auto normal2Sign = values[0] ? 1 : -1;
  normal[2] = normal2Mag * normal2Sign;
  ENSURE(fuzzyEqual(normal.magnitude(), 1.0));
}

} // end namespace Spheral
