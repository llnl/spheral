//---------------------------------Spheral++----------------------------------//
// NodeIteratorBase -- This is a helper class for the DataBase.  This is the
// base class for the node iterators.  This class implements the basic methods,
// most importantly how to index by such iterators into Fields/FieldLists.
// Based on the old NodeIDIterator class.
//
// Created by J. Michael Owen, Mon Mar 17 13:11:19 PST 2003
//----------------------------------------------------------------------------//

#include "NodeIteratorBase.hh"

#include "NodeList/FluidNodeList.hh"

#include <algorithm>
#include <cstdlib>
using std::vector;

namespace Spheral {

//------------------------------------------------------------------------------
// Default constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
NodeIteratorBase<Dimension>::
NodeIteratorBase():
  mNodeID(0),
  mFieldID(0),
  mNodeListBegin(),
  mNodeListEnd(),
  mNodeListItr() {
}

//------------------------------------------------------------------------------
// Copy constructor
//------------------------------------------------------------------------------
template<typename Dimension>
NodeIteratorBase<Dimension>::
NodeIteratorBase(const NodeIteratorBase<Dimension>& rhs):
  mNodeID(rhs.mNodeID),
  mFieldID(rhs.mFieldID),
  mNodeListBegin(rhs.mNodeListBegin),
  mNodeListEnd(rhs.mNodeListEnd),
  mNodeListItr(rhs.mNodeListItr) {
  ENSURE(valid());
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
NodeIteratorBase<Dimension>::
~NodeIteratorBase() {
}

//------------------------------------------------------------------------------
// Assignment
//------------------------------------------------------------------------------
template<typename Dimension>
NodeIteratorBase<Dimension>&
NodeIteratorBase<Dimension>::
operator=(const NodeIteratorBase<Dimension>& rhs) {
  if (&rhs != this) {
    mNodeID = rhs.mNodeID;
    mFieldID = rhs.mFieldID;
    mNodeListBegin = rhs.mNodeListBegin;
    mNodeListEnd = rhs.mNodeListEnd;
    mNodeListItr = rhs.mNodeListItr;
  }
  ENSURE(valid());
  return *this;
}

//------------------------------------------------------------------------------
// Return a pointer to the NodeList cast as a FluidNodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
const FluidNodeList<Dimension>*
NodeIteratorBase<Dimension>::
fluidNodeListPtr() const {
  const FluidNodeList<Dimension>* result = dynamic_cast<const FluidNodeList<Dimension>*>(nodeListPtr());
  ENSURE(result != 0);
  return result;
}

//------------------------------------------------------------------------------
// Base valid test.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
NodeIteratorBase<Dimension>::
valid() const {

  // We want real NodeList iterators.
  // const bool validNodeListIterators = (mNodeListItr != typename vector<NodeList<Dimension>*>::const_iterator() &&
  //                                      mNodeListBegin != typename vector<NodeList<Dimension>*>::const_iterator() &&
  //                                      mNodeListEnd != typename vector<NodeList<Dimension>*>::const_iterator());

  // Test if the NodeList iterator is in the appropriate range.
  const bool nodeListRange = (mNodeListItr >= mNodeListBegin &&
                              mNodeListItr <= mNodeListEnd);

  // Are we on a valid NodeList?
  if (nodeListRange && mNodeListItr < mNodeListEnd) {

    const bool nodeIDTest = mNodeID >= 0 && mNodeID < (int)nodeListPtr()->numNodes();
    const bool fieldIDTest = mFieldID == distance(mNodeListBegin, mNodeListItr);
    return nodeListRange && nodeIDTest && fieldIDTest;

  } else {

    return nodeListRange;

  }

}

}

