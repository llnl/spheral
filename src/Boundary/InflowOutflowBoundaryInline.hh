#include "Boundary/ConstantBoundaryUtilities.hh"

#include <typeinfo>

namespace Spheral {
//------------------------------------------------------------------------------
// beam parameters
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
InflowOutflowBoundary<Dimension>::
inflowRadius() const {
  return mInflowRadius;
}
template<typename Dimension>
inline
void
InflowOutflowBoundary<Dimension>::
inflowRadius(const typename Dimension::Scalar x) {
  mInflowRadius=x;
}


//------------------------------------------------------------------------------
// The effective timestep vote
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
InflowOutflowBoundary<Dimension>::
dtmin() const {
  return mDT;
}

//------------------------------------------------------------------------------
// Return the DataBase this Boundary is defined on.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const DataBase<Dimension>&
InflowOutflowBoundary<Dimension>::
dataBase() const {
  return mDataBase;
}

//------------------------------------------------------------------------------
// The plane for points to enter by.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const GeomPlane<Dimension>& 
InflowOutflowBoundary<Dimension>::
plane() const {
  return mPlane;
}

//------------------------------------------------------------------------------
// Return the number of nodes this Boundary is going to create.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
size_t
InflowOutflowBoundary<Dimension>::
numInflowNodes(const NodeList<Dimension>& nodeList) const {
  const auto itr = mNumInflowNodes.find(nodeList.name());
  VERIFY2(itr != mNumInflowNodes.end(), "InflowOutflowBoundary::numInflowNodes no entry for " << nodeList.name());
  return itr->second;
}

//------------------------------------------------------------------------------
// Access the stored template field values for ghost points.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename DataType>
inline
std::vector<DataType>
InflowOutflowBoundary<Dimension>::
storedValues(const KeyType key, const DataType&) {
  const auto itr = mBufferedValues.find(key);
  VERIFY2(itr != mBufferedValues.end(), "InflowOutflowBoundary ERROR: attempt to extract stored value for " << key << ", which was not found.");
  return extractBufferedValues<DataType>(itr->second);
}

template<typename Dimension>
template<typename DataType>
inline
std::vector<DataType>
InflowOutflowBoundary<Dimension>::
storedValues(const Field<Dimension, DataType>& field) {
  const auto key = StateBase<Dimension>::key(field);
  return storedValues(key, DataTypeTraits<DataType>::zero());
}

//------------------------------------------------------------------------------
// Set the stored template field values for ghost points.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename DataType>
inline
void
InflowOutflowBoundary<Dimension>::
setStoredValues(const KeyType key, const std::vector<DataType>& values) {
  auto itr = mBufferedValues.find(key);
  VERIFY2(itr != mBufferedValues.end(), "InflowOutflowBoundary ERROR: attempt to set stored value for " << key << ", which was not found.");
  auto currentVals = this->storedValues(key, DataType());
  const auto n = currentVals.size();
  VERIFY2(values.size() == n, "InflowOutflowBoundary ERROR: attempt to set stored value for " << key << " with vector of wrong size");
  auto& buffer = itr->second;
  buffer.clear();
  for (auto i = 0u; i < n; ++i) packElement(values[i], buffer);
}

template<typename Dimension>
template<typename DataType>
inline
void
InflowOutflowBoundary<Dimension>::
setStoredValues(const Field<Dimension, DataType>& field, const std::vector<DataType>& values) {
  const auto key = StateBase<Dimension>::key(field);
  this->setStoredValues<DataType>(key, values);
}

//------------------------------------------------------------------------------
// Set the stored template field values for ghost points (constant value).
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename DataType>
inline
void
InflowOutflowBoundary<Dimension>::
setStoredValues(const KeyType key, const DataType& value) {
  auto itr = mBufferedValues.find(key);
  VERIFY2(itr != mBufferedValues.end(), "InflowOutflowBoundary ERROR: attempt to extract stored value for " << key << ", which was not found.");
  auto currentVals = this->storedValues(key, DataType());
  const auto n = currentVals.size();
  auto& buffer = itr->second;
  buffer.clear();
  for (auto i = 0u; i < n; ++i) packElement(value, buffer);
}

template<typename Dimension>
template<typename DataType>
inline
void
InflowOutflowBoundary<Dimension>::
setStoredValues(const Field<Dimension, DataType>& field, const DataType& value) {
  const auto key = StateBase<Dimension>::key(field);
  this->setStoredValues<DataType>(key, value);
}

}
