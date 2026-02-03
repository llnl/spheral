//---------------------------------Spheral++----------------------------------//
// PairwiseField
//
// Stores a value per node pair in a NodePairList. Because connectivity is
// allowed to change step to step in our meshfree methods, PairwiseField is
// ephemeral and will be invalidated when topology is updated.
//
// Created by J. Michael Owen, Wed Nov 20 14:44:44 PST 2024
//----------------------------------------------------------------------------//

#include "Neighbor/ConnectivityMap.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor (connectivity)
//------------------------------------------------------------------------------
template<typename Dimension, typename Value, size_t numElements>
PairwiseField<Dimension, Value, numElements>::PairwiseField(const ConnectivityMap<Dimension>& connectivity):
  PairwiseField<Dimension, Value, numElements>(connectivity.nodePairListPtr()) {
}

//------------------------------------------------------------------------------
// Constructor (pairs)
//------------------------------------------------------------------------------
template<typename Dimension, typename Value, size_t numElements>
PairwiseField<Dimension, Value, numElements>::PairwiseField(std::shared_ptr<NodePairList> pairsPtr):
  PairwiseFieldView<Value, numElements>(),
  mPairsPtr(pairsPtr),
  mArray() {
  if (auto p = mPairsPtr.lock()) {
    mArray.resize(numElements * p->size());
  } else {
    VERIFY2(false, "PairwiseField constructed with invalid NodePairList");
  }
  assignDataSpan();
}

//------------------------------------------------------------------------------
// Copy Constructor.
// Note we deliberately do not use the View copy constructor here, but
// instead reassign its span view in our own assignDataSpan call.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value, size_t numElements>
inline
PairwiseField<Dimension, Value, numElements>::PairwiseField(const PairwiseField& rhs):
  PairwiseFieldView<Value, numElements>(),
  mPairsPtr(rhs.mPairsPtr),
  mArray(rhs.mArray) {
  assignDataSpan();
  DEBUG_LOG << "PairwiseField::copy : " << rhs.mArray.data() << " -> " << mArray.data() << " : " << rhs.mSpan.data() << " -> " << mSpan.data();
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension, typename Value, size_t numElements>
inline
PairwiseField<Dimension, Value, numElements>::
~PairwiseField() {
  DEBUG_LOG << " --> PairwiseField::~PairwiseField() " << this;
#ifndef SPHERAL_UNIFIED_MEMORY
  mSpan.free();
#endif
}

//------------------------------------------------------------------------------
// Assignment operator
//------------------------------------------------------------------------------
template<typename Dimension, typename Value, size_t numElements>
inline
PairwiseField<Dimension, Value, numElements>&
PairwiseField<Dimension, Value, numElements>::operator=(const PairwiseField& rhs) {
  if (this != &rhs) {
    mArray = rhs.mArray;
    assignDataSpan();
  }
  DEBUG_LOG << "PairwiseField::assign : " << rhs.mArray.data() << " -> " << mArray.data() << " : " << rhs.mSpan.data() << " -> " <<  mSpan.data();
  return *this;
}

//------------------------------------------------------------------------------
// Index by pair
//------------------------------------------------------------------------------
template<typename Dimension, typename Value, size_t numElements>
inline
typename PairwiseField<Dimension, Value, numElements>::const_reference
PairwiseField<Dimension, Value, numElements>::operator()(const NodePairIdxType& x) const {
  auto p = mPairsPtr.lock();
  if (!p) {
    VERIFY2(false, "PairwiseField ERROR: attempt to index with invalid pair " << x);
  }
  return Accessor::at(mArray, p->index(x));
}

template<typename Dimension, typename Value, size_t numElements>
inline
typename PairwiseField<Dimension, Value, numElements>::reference
PairwiseField<Dimension, Value, numElements>::operator()(const NodePairIdxType& x) {
  auto p = mPairsPtr.lock();
  if (!p) {
    VERIFY2(false, "PairwiseField ERROR: attempt to index with invalid pair " << x);
  }
  return Accessor::at(mArray, p->index(x));
}

//------------------------------------------------------------------------------
// NodePairList
//------------------------------------------------------------------------------
template<typename Dimension, typename Value, size_t numElements>
inline
const NodePairList&
PairwiseField<Dimension, Value, numElements>::pairs() const {
  auto p = mPairsPtr.lock();
  if (!p) {
    VERIFY2(false, "Orphaned PairwiseField without NodePairList");
  }
  return *p;
}

//------------------------------------------------------------------------------
// Assign the view span
//------------------------------------------------------------------------------
template<typename Dimension, typename Value, size_t numElements>
inline
void
PairwiseField<Dimension, Value, numElements>::assignDataSpan() {
#ifdef SPHERAL_UNIFIED_MEMORY
  mSpan = mArray;
#else
  if (mSpan.size() != mArray.size() or
      mSpan.data(chai::CPU, false) != mArray.data()) {
    DEBUG_LOG << "PairwiseField::assignDataSpan " << this;
    initMAView(mSpan, mArray);
  }
  mSpan.registerTouch(chai::CPU);
#endif
  ENSURE(mSpan.size() == mArray.size() and (mArray.empty() or mSpan.data() == &mArray[0]));
}

}
