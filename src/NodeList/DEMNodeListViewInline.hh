#include "Utilities/DBC.hh"

namespace Spheral {

template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
DEMNodeListView<Dimension>::
DEMNodeListView(const size_t numNodes,
                const size_t firstGhostNode,
                const Scalar hmin,
                const Scalar hmax,
                const Scalar hminratio,
                const Scalar nPerh,
                const Scalar neighborSearchBuffer,
                const size_t maxNumNeighbors):
  NodeListViewType(numNodes, firstGhostNode, hmin, hmax, hminratio, nPerh, maxNumNeighbors),
  mNeighborSearchBuffer(neighborSearchBuffer) {
}

template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
DEMNodeListView<Dimension>::
DEMNodeListView(const size_t numNodes,
                const size_t firstGhostNode,
                const ScalarView& particleRadius,
                const IntView& compositeParticleIndex,
                const IndexView& uniqueIndex,
                const Scalar neighborSearchBuffer):
  NodeListViewType(numNodes, firstGhostNode),
  mNeighborSearchBuffer(neighborSearchBuffer),
  mParticleRadiusView(particleRadius),
  mCompositeParticleIndexView(compositeParticleIndex),
  mUniqueIndexView(uniqueIndex) {
  REQUIRE(mParticleRadiusView.numElements() == this->mNumNodes);
  REQUIRE(mCompositeParticleIndexView.numElements() == this->mNumNodes);
  REQUIRE(mUniqueIndexView.numElements() == this->mNumNodes);
}

template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
typename DEMNodeListView<Dimension>::Scalar&
DEMNodeListView<Dimension>::particleRadius(const size_t nodeID) const {
  REQUIRE(nodeID < mParticleRadiusView.numElements());
  return mParticleRadiusView[nodeID];
}

template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
int&
DEMNodeListView<Dimension>::compositeParticleIndex(const size_t nodeID) const {
  REQUIRE(nodeID < mCompositeParticleIndexView.numElements());
  return mCompositeParticleIndexView[nodeID];
}

template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
size_t&
DEMNodeListView<Dimension>::uniqueIndex(const size_t nodeID) const {
  REQUIRE(nodeID < mUniqueIndexView.numElements());
  return mUniqueIndexView[nodeID];
}

template<typename Dimension>
SPHERAL_HOST
inline
void
DEMNodeListView<Dimension>::move(chai::ExecutionSpace space) {
  NodeListViewType::move(space);
#ifndef SPHERAL_UNIFIED_MEMORY
  mParticleRadiusView.move(space);
  mCompositeParticleIndexView.move(space);
  mUniqueIndexView.move(space);
#endif
}

template<typename Dimension>
SPHERAL_HOST
inline
void
DEMNodeListView<Dimension>::touch(chai::ExecutionSpace space) {
  NodeListViewType::touch(space);
#ifndef SPHERAL_UNIFIED_MEMORY
  mParticleRadiusView.touch(space);
  mCompositeParticleIndexView.touch(space);
  mUniqueIndexView.touch(space);
#endif
}

}
