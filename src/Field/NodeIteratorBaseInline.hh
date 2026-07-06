#include "Utilities/DBC.hh"
#include "Neighbor/Neighbor.hh"
#include "NodeList/NodeListView.hh"
#include "NodeList/NodeList.hh"

#include <vector>

namespace Spheral {

//------------------------------------------------------------------------------
// State info we can extract from the NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
NodeIteratorBase<Dimension>::
mass() const {
  return nodeListPtr()->mass()(mNodeID);
}

template<typename Dimension>
inline
const typename Dimension::Vector&
NodeIteratorBase<Dimension>::
position() const {
  return nodeListPtr()->positions()(mNodeID);
}

template<typename Dimension>
inline
const typename Dimension::Vector&
NodeIteratorBase<Dimension>::
velocity() const {
  return nodeListPtr()->velocity()(mNodeID);
}

template<typename Dimension>
inline
const typename Dimension::SymTensor&
NodeIteratorBase<Dimension>::
H() const {
  return nodeListPtr()->Hfield()(mNodeID);
}

//------------------------------------------------------------------------------
// operator==
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
NodeIteratorBase<Dimension>::
operator==(const NodeIteratorBase<Dimension>& rhs) const {
  return (nodeListPtr() == rhs.nodeListPtr() &&
          nodeID() == rhs.nodeID());
}

//------------------------------------------------------------------------------
// operator<
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
NodeIteratorBase<Dimension>::
operator<(const NodeIteratorBase<Dimension>& rhs) const {

  // Do a pointer comparison between the NodeLists.
  if (nodeListIterator() < mNodeListEnd &&
      rhs.nodeListIterator() < rhs.mNodeListEnd) {

    const NodeList<Dimension>* lhsPtr = nodeListPtr();
    const NodeList<Dimension>* rhsPtr = rhs.nodeListPtr();
    if (lhsPtr < rhsPtr) {
      return true;
    } else if (lhsPtr == rhsPtr) {
      return nodeID() < rhs.nodeID();
    } else {
      return false;
    }

  } else {

    if (nodeListIterator() < mNodeListEnd &&
        rhs.nodeListIterator() == rhs.mNodeListEnd) {
      return true;
    } else {
      return false;
    }

  }
}

//------------------------------------------------------------------------------
// Return the NodeList iterator.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename std::vector<NodeList<Dimension>*>::const_iterator
NodeIteratorBase<Dimension>::
nodeListIterator() const {
  return mNodeListItr;
}

//------------------------------------------------------------------------------
// Return a pointer to the NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const NodeList<Dimension>*
NodeIteratorBase<Dimension>::
nodeListPtr() const {
  if (mNodeListItr < mNodeListEnd) {
    return *mNodeListItr;
  } else {
    return 0;
  }
}

//------------------------------------------------------------------------------
// Return the node ID.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
NodeIteratorBase<Dimension>::
nodeID() const {
  return mNodeID;
}

//------------------------------------------------------------------------------
// Return the index of the current Field.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
NodeIteratorBase<Dimension>::
fieldID() const {
  return mFieldID;
}

//------------------------------------------------------------------------------
// Test if this is an internal node.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
NodeIteratorBase<Dimension>::
internalNode() const {
  return nodeListPtr()->nodeType(nodeID()) == NodeType::InternalNode;
}

//------------------------------------------------------------------------------
// Test if this is a ghost node.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
NodeIteratorBase<Dimension>::
ghostNode() const {
  return nodeListPtr()->nodeType(nodeID()) == NodeType::GhostNode;
}

//------------------------------------------------------------------------------
// operator!=
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
NodeIteratorBase<Dimension>::
operator!=(const NodeIteratorBase<Dimension>& rhs) const {
  return !(operator==(rhs));
}

//------------------------------------------------------------------------------
// operator>
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
NodeIteratorBase<Dimension>::
operator>(const NodeIteratorBase<Dimension>& rhs) const {
  return !(operator<(rhs) || operator==(rhs));
}

//------------------------------------------------------------------------------
// operator<=
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
NodeIteratorBase<Dimension>::
operator<=(const NodeIteratorBase<Dimension>& rhs) const {
  return operator==(rhs) || operator<(rhs);
}

//------------------------------------------------------------------------------
// operator>=
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
NodeIteratorBase<Dimension>::
operator>=(const NodeIteratorBase<Dimension>& rhs) const {
  return operator==(rhs) || operator>(rhs);
}

}
