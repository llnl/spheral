#include "Field/Field.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// buffer distance for contact detection
//------------------------------------------------------------------------------

template<typename Dimension>
inline
typename Dimension::Scalar
DEMNodeList<Dimension>::
neighborSearchBuffer() const {
  return mNeighborSearchBuffer;
}

template<typename Dimension>
inline
void
DEMNodeList<Dimension>::
neighborSearchBuffer(typename Dimension::Scalar x) {
  mNeighborSearchBuffer = x;
}

//------------------------------------------------------------------------------
// particle radius per node.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
Field<Dimension, typename Dimension::Scalar>&
DEMNodeList<Dimension>::particleRadius() {
  REQUIRE(mParticleRadius.nodeListPtr() == this);
  return mParticleRadius;
}

template<typename Dimension>
inline
const Field<Dimension, typename Dimension::Scalar>&
DEMNodeList<Dimension>::particleRadius() const {
  REQUIRE(mParticleRadius.nodeListPtr() == this);
  return mParticleRadius;
}


//------------------------------------------------------------------------------
// composite particle indices
//------------------------------------------------------------------------------
template<typename Dimension>
inline
Field<Dimension, int>&
DEMNodeList<Dimension>::compositeParticleIndex() {
  REQUIRE(mCompositeParticleIndex.nodeListPtr() == this);
  return mCompositeParticleIndex;
}

template<typename Dimension>
inline
const Field<Dimension, int>&
DEMNodeList<Dimension>::compositeParticleIndex() const {
  REQUIRE(mCompositeParticleIndex.nodeListPtr() == this);
  return mCompositeParticleIndex;
}

//------------------------------------------------------------------------------
// unique particle indices
//------------------------------------------------------------------------------
template<typename Dimension>
inline
Field<Dimension, size_t>&
DEMNodeList<Dimension>::uniqueIndex() {
  REQUIRE(mUniqueIndex.nodeListPtr() == this);
  return mUniqueIndex;
}

template<typename Dimension>
inline
const Field<Dimension, size_t>&
DEMNodeList<Dimension>::uniqueIndex() const {
  REQUIRE(mUniqueIndex.nodeListPtr() == this);
  return mUniqueIndex;
}
}
