#include "Utilities/DBC.hh"
#include "Field/Field.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Get the name of the NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
std::string
NodeList<Dimension>::name() const {
  return mName;
}

//------------------------------------------------------------------------------
// Return the number of nodes.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
size_t 
NodeList<Dimension>::numNodes() const {
  return mNumNodes;
}

//------------------------------------------------------------------------------
// Return the number of internal nodes.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
size_t 
NodeList<Dimension>::numInternalNodes() const {
  CHECK2(mFirstGhostNode <= numNodes(), "NodeList sizes: " << mFirstGhostNode << " " << numNodes());
  return mFirstGhostNode;
}

//------------------------------------------------------------------------------
// Return the number of ghost nodes.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
size_t 
NodeList<Dimension>::numGhostNodes() const {
  CHECK(mFirstGhostNode <= numNodes());
  return numNodes() - mFirstGhostNode;
}

//------------------------------------------------------------------------------
// Access the mass field.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
Field<Dimension, typename Dimension::Scalar>&
NodeList<Dimension>::mass() {
  CHECK(mMass.nodeListPtr() == this);
  return mMass;
}

template<typename Dimension>
inline
const Field<Dimension, typename Dimension::Scalar>&
NodeList<Dimension>::mass() const {
  CHECK(mMass.nodeListPtr() == this);
  return mMass;
}

//------------------------------------------------------------------------------
// Access the position field.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
Field<Dimension, typename Dimension::Vector>&
NodeList<Dimension>::positions() {
  CHECK(mPositions.nodeListPtr() == this);
  return mPositions;
}

template<typename Dimension>
inline
const Field<Dimension, typename Dimension::Vector>&
NodeList<Dimension>::positions() const {
  CHECK(mPositions.nodeListPtr() == this);
  return mPositions;
}

//------------------------------------------------------------------------------
// Access the velocity field.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
Field<Dimension, typename Dimension::Vector>&
NodeList<Dimension>::velocity() {
  CHECK(mVelocity.nodeListPtr() == this);
  return mVelocity;
}

template<typename Dimension>
inline
const Field<Dimension, typename Dimension::Vector>&
NodeList<Dimension>::velocity() const {
  CHECK(mVelocity.nodeListPtr() == this);
  return mVelocity;
}

//------------------------------------------------------------------------------
// Smoothing scale per node.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
Field<Dimension, typename Dimension::SymTensor>&
NodeList<Dimension>::Hfield() {
  CHECK(mH.nodeListPtr() == this);
  return mH;
}

template<typename Dimension>
inline
const Field<Dimension, typename Dimension::SymTensor>&
NodeList<Dimension>::Hfield() const {
  CHECK(mH.nodeListPtr() == this);
  return mH;
}

//------------------------------------------------------------------------------
// Work per node.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
Field<Dimension, typename Dimension::Scalar>&
NodeList<Dimension>::work() const {
  CHECK(mWork.nodeListPtr() == this);
  return mWork;
}

//------------------------------------------------------------------------------
// Iterators over the Fields registered on this NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename std::vector<FieldBase<Dimension>*>::iterator
NodeList<Dimension>::
registeredFieldsBegin() {
  return mFieldBaseList.begin();
}

template<typename Dimension>
inline
typename std::vector<FieldBase<Dimension>*>::iterator
NodeList<Dimension>::
registeredFieldsEnd() {
  return mFieldBaseList.end();
}

template<typename Dimension>
inline
typename std::vector<FieldBase<Dimension>*>::const_iterator
NodeList<Dimension>::
registeredFieldsBegin() const {
  return mFieldBaseList.begin();
}

template<typename Dimension>
inline
typename std::vector<FieldBase<Dimension>*>::const_iterator
NodeList<Dimension>::
registeredFieldsEnd() const {
  return mFieldBaseList.end();
}

//------------------------------------------------------------------------------
// Access the target number of nodes to maintain per smoothing scale.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
NodeList<Dimension>::nodesPerSmoothingScale() const {
  return mNodesPerSmoothingScale;
}

template<typename Dimension>
void
NodeList<Dimension>::nodesPerSmoothingScale(typename Dimension::Scalar val) {
  mNodesPerSmoothingScale = val;
}

//------------------------------------------------------------------------------
// Access the maximum number of neighbors to allow.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
size_t 
NodeList<Dimension>::maxNumNeighbors() const {
  return mMaxNumNeighbors;
}

template<typename Dimension>
void
NodeList<Dimension>::maxNumNeighbors(size_t val) {
  mMaxNumNeighbors = val;
}

//------------------------------------------------------------------------------
// The min allowed smoothing scale.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
NodeList<Dimension>::hmin() const {
  return mhmin;
}

template<typename Dimension>
inline
void
NodeList<Dimension>::hmin(typename Dimension::Scalar val) {
  mhmin = val;
}

//------------------------------------------------------------------------------
// The max allowed smoothing scale.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
NodeList<Dimension>::hmax() const {
  return mhmax;
}

template<typename Dimension>
inline
void
NodeList<Dimension>::hmax(typename Dimension::Scalar val) {
  mhmax = val;
}

//------------------------------------------------------------------------------
// The min allowed ratio of the eigen values of the H tensor.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
NodeList<Dimension>::hminratio() const {
  return mhminratio;
}

template<typename Dimension>
inline
void
NodeList<Dimension>::hminratio(typename Dimension::Scalar val) {
  mhminratio = val;
}

//------------------------------------------------------------------------------
// Comparison operators.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
NodeList<Dimension>::operator==(const NodeList<Dimension>& rhs) const {
  return this == &rhs;
}

template<typename Dimension>
inline
bool
NodeList<Dimension>::operator!=(const NodeList<Dimension>& rhs) const {
  return not (*this == rhs);
}

}
