#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructors.
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
NodeListView<Dimension>::
NodeListView(const size_t numNodes,
             const size_t firstGhostNode,
             const Scalar hmin,
             const Scalar hmax,
             const Scalar hminratio,
             const Scalar nPerh,
             const size_t maxNumNeighbors):
  mNumNodes(numNodes),
  mFirstGhostNode(firstGhostNode),
  mhmin(hmin),
  mhmax(hmax),
  mhminratio(hminratio),
  mNodesPerSmoothingScale(nPerh),
  mMaxNumNeighbors(maxNumNeighbors),
  mMassView(),
  mPositionsView(),
  mVelocityView(),
  mHfieldView(),
  mWorkView() {
  REQUIRE(mFirstGhostNode <= mNumNodes);
}

template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
NodeListView<Dimension>::
NodeListView(const size_t numNodes,
             const size_t firstGhostNode,
             const PositionView& positions,
             const HfieldView& Hfield,
             const Scalar hmin,
             const Scalar hmax,
             const Scalar hminratio,
             const Scalar nPerh,
             const size_t maxNumNeighbors):
  mNumNodes(numNodes),
  mFirstGhostNode(firstGhostNode),
  mhmin(hmin),
  mhmax(hmax),
  mhminratio(hminratio),
  mNodesPerSmoothingScale(nPerh),
  mMaxNumNeighbors(maxNumNeighbors),
  mMassView(),
  mPositionsView(positions),
  mVelocityView(),
  mHfieldView(Hfield),
  mWorkView() {
  REQUIRE(mFirstGhostNode <= mNumNodes);
  REQUIRE(mPositionsView.numElements() == mNumNodes);
  REQUIRE(mHfieldView.numElements() == mNumNodes);
}

template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
NodeListView<Dimension>::
NodeListView(const size_t numNodes,
             const size_t firstGhostNode,
             const MassView& mass,
             const PositionView& positions,
             const VelocityView& velocity,
             const HfieldView& Hfield,
             const WorkView& work,
             const Scalar hmin,
             const Scalar hmax,
             const Scalar hminratio,
             const Scalar nPerh,
             const size_t maxNumNeighbors):
  mNumNodes(numNodes),
  mFirstGhostNode(firstGhostNode),
  mhmin(hmin),
  mhmax(hmax),
  mhminratio(hminratio),
  mNodesPerSmoothingScale(nPerh),
  mMaxNumNeighbors(maxNumNeighbors),
  mMassView(mass),
  mPositionsView(positions),
  mVelocityView(velocity),
  mHfieldView(Hfield),
  mWorkView(work) {
  REQUIRE(mFirstGhostNode <= mNumNodes);
  REQUIRE(mMassView.numElements() == mNumNodes);
  REQUIRE(mPositionsView.numElements() == mNumNodes);
  REQUIRE(mVelocityView.numElements() == mNumNodes);
  REQUIRE(mHfieldView.numElements() == mNumNodes);
  REQUIRE(mWorkView.numElements() == mNumNodes);
}

//------------------------------------------------------------------------------
// CHAI/ManagedArray specific operations.
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST
inline
void
NodeListView<Dimension>::
move(chai::ExecutionSpace space) {
#ifndef SPHERAL_UNIFIED_MEMORY
  mMassView.move(space);
  mPositionsView.move(space);
  mVelocityView.move(space);
  mHfieldView.move(space);
  mWorkView.move(space);
#endif
}

template<typename Dimension>
SPHERAL_HOST
inline
void
NodeListView<Dimension>::
touch(chai::ExecutionSpace space) {
#ifndef SPHERAL_UNIFIED_MEMORY
  mMassView.touch(space);
  mPositionsView.touch(space);
  mVelocityView.touch(space);
  mHfieldView.touch(space);
  mWorkView.touch(space);
#endif
}

}
