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
// Constructor
//------------------------------------------------------------------------------
template<typename Dimension, typename Value, size_t numElements>
PairwiseField<Dimension, Value, numElements>::PairwiseField(const ConnectivityMap<Dimension>& connectivity):
  mPairsPtr(connectivity.nodePairListPtr()),
  mValues() {
  if (auto p = mPairsPtr.lock()) {
    mValues.resize(numElements * p->size());
  } else {
    VERIFY2(false, "PairwiseField constructed with invalid NodePairList");
  }
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
  return Accessor::at(mValues, p->index(x));
}

template<typename Dimension, typename Value, size_t numElements>
inline
typename PairwiseField<Dimension, Value, numElements>::reference
PairwiseField<Dimension, Value, numElements>::operator()(const NodePairIdxType& x) {
  auto p = mPairsPtr.lock();
  if (!p) {
    VERIFY2(false, "PairwiseField ERROR: attempt to index with invalid pair " << x);
  }
  return Accessor::at(mValues, p->index(x));
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

}
