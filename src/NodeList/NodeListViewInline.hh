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
             const size_t firstGhostNode):
  mNumNodes(numNodes),
  mFirstGhostNode(firstGhostNode),
  mPositionsView(),
  mHfieldView() {
}

template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
NodeListView<Dimension>::
NodeListView(const size_t numNodes,
             const size_t firstGhostNode,
             const PositionView& positions,
             const HfieldView& Hfield):
  mNumNodes(numNodes),
  mFirstGhostNode(firstGhostNode),
  mPositionsView(positions),
  mHfieldView(Hfield) {
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

//------------------------------------------------------------------------------
// Field views.
//------------------------------------------------------------------------------
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
typename Dimension::Vector&
NodeListView<Dimension>::
position(const size_t nodeID) const {
  REQUIRE(nodeID < mPositionsView.numElements());
  return mPositionsView[nodeID];
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
  mPositionsView.move(space);
  mHfieldView.move(space);
#endif
}

template<typename Dimension>
SPHERAL_HOST
inline
void
NodeListView<Dimension>::
touch(chai::ExecutionSpace space) {
#ifndef SPHERAL_UNIFIED_MEMORY
  mPositionsView.touch(space);
  mHfieldView.touch(space);
#endif
}

}
