#include "Field/FieldBase.hh"
#include "Utilities/DBC.hh"
#include "Utilities/Logger.hh"
#include "Utilities/SpheralMessage.hh"

#include <algorithm>

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
NodeListBase<Dimension>::
NodeListBase(std::string name,
             const Scalar hmin,
             const Scalar hmax,
             const Scalar hminratio,
             const Scalar nPerh,
             const size_t maxNumNeighbors):
  mName(name),
  mFieldBases(),
  mNeighborPtr(nullptr),
  mhmin(hmin),
  mhmax(hmax),
  mhminratio(hminratio),
  mNodesPerSmoothingScale(nPerh),
  mMaxNumNeighbors(maxNumNeighbors) {
}

//------------------------------------------------------------------------------
// Access the name.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
std::string
NodeListBase<Dimension>::
name() const {
  return mName;
}

template<typename Dimension>
inline
void
NodeListBase<Dimension>::
name(std::string name) {
  mName = name;
}

//------------------------------------------------------------------------------
// Registered Fields.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename NodeListBase<Dimension>::FieldBaseIterator
NodeListBase<Dimension>::
registeredFieldsBegin() {
  return mFieldBases.begin();
}

template<typename Dimension>
inline
typename NodeListBase<Dimension>::FieldBaseIterator
NodeListBase<Dimension>::
registeredFieldsEnd() {
  return mFieldBases.end();
}

template<typename Dimension>
inline
typename NodeListBase<Dimension>::const_FieldBaseIterator
NodeListBase<Dimension>::
registeredFieldsBegin() const {
  return mFieldBases.begin();
}

template<typename Dimension>
inline
typename NodeListBase<Dimension>::const_FieldBaseIterator
NodeListBase<Dimension>::
registeredFieldsEnd() const {
  return mFieldBases.end();
}

template<typename Dimension>
inline
typename NodeListBase<Dimension>::FieldBaseSpan
NodeListBase<Dimension>::
registeredFields() const {
  return FieldBaseSpan(mFieldBases.data(), mFieldBases.size());
}

template<typename Dimension>
inline
size_t
NodeListBase<Dimension>::
numFields() const {
  return mFieldBases.size();
}

//------------------------------------------------------------------------------
// Check if the given field is registered with this NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
NodeListBase<Dimension>::
haveField(const FieldBase<Dimension>& field) const {
  return (std::find_if(mFieldBases.begin(),
                       mFieldBases.end(),
                       [&](const std::reference_wrapper<FieldBase<Dimension>>& x) { return &(x.get()) == &field; }) != mFieldBases.end());
}

//------------------------------------------------------------------------------
// Access the neighbor object.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
Neighbor<Dimension>&
NodeListBase<Dimension>::
neighbor() const {
  CHECK(mNeighborPtr != nullptr);
  return *mNeighborPtr;
}

//------------------------------------------------------------------------------
// Smoothing scale controls.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename NodeListBase<Dimension>::Scalar
NodeListBase<Dimension>::
nodesPerSmoothingScale() const {
  return mNodesPerSmoothingScale;
}

template<typename Dimension>
inline
void
NodeListBase<Dimension>::
nodesPerSmoothingScale(Scalar val) {
  mNodesPerSmoothingScale = val;
}

template<typename Dimension>
inline
size_t
NodeListBase<Dimension>::
maxNumNeighbors() const {
  return mMaxNumNeighbors;
}

template<typename Dimension>
inline
void
NodeListBase<Dimension>::
maxNumNeighbors(size_t val) {
  mMaxNumNeighbors = val;
}

template<typename Dimension>
inline
typename NodeListBase<Dimension>::Scalar
NodeListBase<Dimension>::
hmin() const {
  return mhmin;
}

template<typename Dimension>
inline
void
NodeListBase<Dimension>::
hmin(Scalar val) {
  mhmin = val;
}

template<typename Dimension>
inline
typename NodeListBase<Dimension>::Scalar
NodeListBase<Dimension>::
hmax() const {
  return mhmax;
}

template<typename Dimension>
inline
void
NodeListBase<Dimension>::
hmax(Scalar val) {
  mhmax = val;
}

template<typename Dimension>
inline
typename NodeListBase<Dimension>::Scalar
NodeListBase<Dimension>::
hminratio() const {
  return mhminratio;
}

template<typename Dimension>
inline
void
NodeListBase<Dimension>::
hminratio(Scalar val) {
  mhminratio = val;
}

//------------------------------------------------------------------------------
// Some operators.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
NodeListBase<Dimension>::
operator==(const NodeListBase& rhs) const {
  return this == &rhs;
}

template<typename Dimension>
inline
bool
NodeListBase<Dimension>::
operator!=(const NodeListBase& rhs) const {
  return !(*this == rhs);
}

//------------------------------------------------------------------------------
// Register a field with this NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
NodeListBase<Dimension>::
registerField(FieldBase<Dimension>& field) const {
  DEBUG_LOG << "NodeListBase::registerField : " << mName << " " << this << " : " << field.name() << " " << &field;
  if (haveField(field)) {
    SpheralMessage("WARNING: Attempt to register field " << &field << " (" << field.name()
                   << ") with NodeList " << this << " (" << this->name() << ") that already has it.");
  } else {
    mFieldBases.push_back(std::ref(field));
  }
}

//------------------------------------------------------------------------------
// Unregister a field that is listed with this NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
NodeListBase<Dimension>::
unregisterField(FieldBase<Dimension>& field) const {
  DEBUG_LOG << "NodeListBase::unregisterField : " << mName << " " << this << " : " << field.name() << " " << &field;
#pragma omp critical
  {
    auto itr = std::find_if(mFieldBases.begin(),
                            mFieldBases.end(),
                            [&](const std::reference_wrapper<FieldBase<Dimension>>& x) { return &(x.get()) == &field; });
    if (itr != mFieldBases.end()) mFieldBases.erase(itr);
  }
}

//------------------------------------------------------------------------------
// Register the given neighbor object with this node list.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
NodeListBase<Dimension>::
registerNeighbor(Neighbor<Dimension>& neighbor) {
  DEBUG_LOG << "NodeListBase::registerNeighbor : " << mName << " " << this << " : " << &neighbor;
  mNeighborPtr = &neighbor;
}

//------------------------------------------------------------------------------
// Unregister the current neighbor object from this node list.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
NodeListBase<Dimension>::
unregisterNeighbor() {
  DEBUG_LOG << "NodeListBase::unregisterNeighbor : " << mName << " " << this;
  mNeighborPtr = nullptr;
}

}
