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
// Node count accessors.
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
size_t
NodeListView<Dimension>::
numNodes() const {
  return mNumNodes;
}

template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
size_t
NodeListView<Dimension>::
numInternalNodes() const {
  return mFirstGhostNode;
}

template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
size_t
NodeListView<Dimension>::
firstGhostNode() const {
  return mFirstGhostNode;
}

template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
size_t
NodeListView<Dimension>::
numGhostNodes() const {
  CHECK(mFirstGhostNode <= mNumNodes);
  return mNumNodes - mFirstGhostNode;
}

template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
NodeType
NodeListView<Dimension>::
nodeType(size_t i) const {
  REQUIRE(i < mNumNodes);
  return (i < mFirstGhostNode ?
          NodeType::InternalNode :
          NodeType::GhostNode);
}

//------------------------------------------------------------------------------
// Smoothing scale controls.
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
typename Dimension::Scalar
NodeListView<Dimension>::
nodesPerSmoothingScale() const {
  return mNodesPerSmoothingScale;
}

template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
void
NodeListView<Dimension>::
nodesPerSmoothingScale(Scalar val) {
  mNodesPerSmoothingScale = val;
}

template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
size_t
NodeListView<Dimension>::
maxNumNeighbors() const {
  return mMaxNumNeighbors;
}

template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
void
NodeListView<Dimension>::
maxNumNeighbors(size_t val) {
  mMaxNumNeighbors = val;
}

template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
typename Dimension::Scalar
NodeListView<Dimension>::
hmin() const {
  return mhmin;
}

template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
void
NodeListView<Dimension>::
hmin(Scalar val) {
  mhmin = val;
}

template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
typename Dimension::Scalar
NodeListView<Dimension>::
hmax() const {
  return mhmax;
}

template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
void
NodeListView<Dimension>::
hmax(Scalar val) {
  mhmax = val;
}

template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
typename Dimension::Scalar
NodeListView<Dimension>::
hminratio() const {
  return mhminratio;
}

template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
void
NodeListView<Dimension>::
hminratio(Scalar val) {
  mhminratio = val;
}

//------------------------------------------------------------------------------
// Field views.
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
typename NodeListView<Dimension>::MassView&
NodeListView<Dimension>::
mass() {
  return mMassView;
}

template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
const typename NodeListView<Dimension>::MassView&
NodeListView<Dimension>::
mass() const {
  return mMassView;
}

template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
typename NodeListView<Dimension>::PositionView&
NodeListView<Dimension>::
positions() {
  return mPositionsView;
}

template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
const typename NodeListView<Dimension>::PositionView&
NodeListView<Dimension>::
positions() const {
  return mPositionsView;
}

template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
typename NodeListView<Dimension>::VelocityView&
NodeListView<Dimension>::
velocity() {
  return mVelocityView;
}

template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
const typename NodeListView<Dimension>::VelocityView&
NodeListView<Dimension>::
velocity() const {
  return mVelocityView;
}

template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
typename NodeListView<Dimension>::HfieldView&
NodeListView<Dimension>::
Hfield() {
  return mHfieldView;
}

template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
const typename NodeListView<Dimension>::HfieldView&
NodeListView<Dimension>::
Hfield() const {
  return mHfieldView;
}

template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
typename NodeListView<Dimension>::WorkView&
NodeListView<Dimension>::
work() {
  return mWorkView;
}

template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
const typename NodeListView<Dimension>::WorkView&
NodeListView<Dimension>::
work() const {
  return mWorkView;
}

template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
typename Dimension::Scalar&
NodeListView<Dimension>::
mass(const size_t nodeID) const {
  REQUIRE(nodeID < mMassView.numElements());
  return mMassView[nodeID];
}

template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
typename Dimension::Vector&
NodeListView<Dimension>::
position(const size_t nodeID) const {
  REQUIRE(nodeID < mPositionsView.numElements());
  return mPositionsView[nodeID];
}

template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
typename Dimension::Vector&
NodeListView<Dimension>::
velocity(const size_t nodeID) const {
  REQUIRE(nodeID < mVelocityView.numElements());
  return mVelocityView[nodeID];
}

template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
typename Dimension::SymTensor&
NodeListView<Dimension>::
H(const size_t nodeID) const {
  REQUIRE(nodeID < mHfieldView.numElements());
  return mHfieldView[nodeID];
}

template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
typename Dimension::Scalar&
NodeListView<Dimension>::
work(const size_t nodeID) const {
  REQUIRE(nodeID < mWorkView.numElements());
  return mWorkView[nodeID];
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
