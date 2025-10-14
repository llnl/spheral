#include "Field/Field.hh"
#include "Geometry/Dimension.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/DataTypeTraits.hh"
#include "Utilities/Logger.hh"

#include <cmath>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <sstream>
#include <limits>

// Inlined methods.
namespace Spheral {

//------------------------------------------------------------------------------
// Assignment operator with a constant value of DataType
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE 
inline
FieldView<Dimension, DataType>&
FieldView<Dimension, DataType>::
operator=(const DataType& rhs) {
  for (auto& x: mDataSpan) x = rhs;
  return *this;
}

//------------------------------------------------------------------------------
// Element access by integer index.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE 
inline
DataType&
FieldView<Dimension, DataType>::operator()(size_t index) {
  CHECK2(index < mDataSpan.size(), "FieldView index out of range: " << index << " " << mDataSpan.size());
  return mDataSpan[index];
}

template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE 
inline
DataType&
FieldView<Dimension, DataType>::operator()(size_t index) const {
  CHECK2(index < mDataSpan.size(), "FieldView index out of range: " << index << " " << mDataSpan.size());
  return mDataSpan[index];
}

//------------------------------------------------------------------------------
// at version, for consistency with STL interface.
// Since std::span doesn't support at() until C++26, we emulate it with VERIFY
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE 
inline
DataType&
FieldView<Dimension, DataType>::at(size_t index) {
  VERIFY2(index < mDataSpan.size(), "FieldView index out of range: " << index << " " << mDataSpan.size());
  return mDataSpan[index];
}

template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE 
inline
DataType&
FieldView<Dimension, DataType>::at(size_t index) const {
  VERIFY2(index < mDataSpan.size(), "FieldView index out of range: " << index << " " << mDataSpan.size());
  return mDataSpan[index];
}

//------------------------------------------------------------------------------
// Index operators.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE 
inline
DataType&
FieldView<Dimension, DataType>::
operator[](const size_t index) {
  CHECK2(index < mDataSpan.size(), "FieldView index out of range: " << index << " " << mDataSpan.size());
  return mDataSpan[index];
}

template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE 
inline
DataType&
FieldView<Dimension, DataType>::
operator[](const size_t index) const {
  CHECK2(index < mDataSpan.size(), "FieldView index out of range: " << index << " " << mDataSpan.size());
  return mDataSpan[index];
}

//------------------------------------------------------------------------------
// Apply a minimum value to the elements.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE 
inline
void
FieldView<Dimension, DataType>::applyMin(const DataType& dataMin) {
  for (auto& x: *this) x = std::max(x, dataMin);
}

//------------------------------------------------------------------------------
// Apply a maximum value to the elements.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE 
inline
void
FieldView<Dimension, DataType>::applyMax(const DataType& dataMax) {
  for (auto& x: *this) x = std::min(x, dataMax);
}

//------------------------------------------------------------------------------
// Apply a scalar minimum value  to the elements.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE 
inline
void
FieldView<Dimension, DataType>::applyScalarMin(const Scalar& dataMin) {
  for (auto& x: *this) x = std::max(x, dataMin);
}

//------------------------------------------------------------------------------
// Apply a scalar maximum value to the elements.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE 
inline
void
FieldView<Dimension, DataType>::applyScalarMax(const Scalar& dataMax) {
  for (auto& x: *this) x = std::min(x, dataMax);
}

//------------------------------------------------------------------------------
// Addition with another FieldView in place
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE 
inline
FieldView<Dimension, DataType>&
FieldView<Dimension, DataType>::
operator+=(const FieldView<Dimension, DataType>& rhs) {
  const auto n = mDataSpan.size();
  REQUIRE(rhs.numElements() == n);
  for (auto i = 0u; i < n; ++i) (*this)(i) += rhs(i);
  return *this;
}

//------------------------------------------------------------------------------
// Subtract another FieldView from this one in place.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE 
inline
FieldView<Dimension, DataType>&
FieldView<Dimension, DataType>::
operator-=(const FieldView<Dimension, DataType>& rhs) {
  const auto n = mDataSpan.size();
  REQUIRE(rhs.numElements() == n);
  for (auto i = 0u; i < n; ++i) (*this)(i) -= rhs(i);
  return *this;
}

//------------------------------------------------------------------------------
// Addition with a single value in place
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE 
inline
FieldView<Dimension, DataType>&
FieldView<Dimension, DataType>::
operator+=(const DataType& rhs) {
  for (auto& x: *this) x += rhs;
  return *this;
}

//------------------------------------------------------------------------------
// Subtract a single value in place
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE 
inline
FieldView<Dimension, DataType>&
FieldView<Dimension, DataType>::
operator-=(const DataType& rhs) {
  for (auto& x: *this) x -= rhs;
  return *this;
}

//------------------------------------------------------------------------------
// Multiplication by a Scalar FieldView in place.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE 
inline
FieldView<Dimension, DataType>&
FieldView<Dimension, DataType>::
operator*=(const FieldView<Dimension, Scalar>& rhs) {
  const auto n = mDataSpan.size();
  REQUIRE(rhs.numElements() == n);
  for (auto i = 0u; i < n; ++i) (*this)(i) *= rhs(i);
  return *this;
}

//------------------------------------------------------------------------------
// Division by a Scalar FieldView in place.
//-------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE 
inline
FieldView<Dimension, DataType>&
FieldView<Dimension, DataType>::
operator/=(const FieldView<Dimension, typename Dimension::Scalar>& rhs) {
  const auto n = mDataSpan.size();
  REQUIRE(rhs.numElements() == n);
  for (auto i = 0u; i < n; ++i) (*this)(i) *= safeInvVar(rhs(i), 1.0e-60);
  return *this;
}

//------------------------------------------------------------------------------
// Multiplication by a Scalar in place
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE 
inline
FieldView<Dimension, DataType>&
FieldView<Dimension, DataType>::
operator*=(const Scalar& rhs) {
  for (auto& x: *this) x *= rhs;
  return *this;
}

//------------------------------------------------------------------------------
// Division by a Scalar value in place
//-------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE 
inline
FieldView<Dimension, DataType>&
FieldView<Dimension, DataType>::
operator/=(const Scalar& rhs) {
  const auto rhsInv = safeInvVar(rhs, 1.0e-60);
  for (auto& x: *this) x *= rhsInv;
  return *this;
}

//------------------------------------------------------------------------------
// Sum the elements of the FieldView (assumes the DataType::operator+= is 
// available).  
// LOCAL to processor!
//-------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE 
inline
DataType
FieldView<Dimension, DataType>::
localSumElements(const bool includeGhosts) const {
  const auto* start = &this->front();
  const auto n = includeGhosts ? numElements() : numInternalElements();
  return std::accumulate(start, start + n, DataTypeTraits<DataType>::zero());
}

//------------------------------------------------------------------------------
// Minimum.
// LOCAL to processor!
//-------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE 
inline
DataType
FieldView<Dimension, DataType>::
localMin(const bool includeGhosts) const {
  const auto* start = &this->front();
  const auto n = includeGhosts ? numElements() : numInternalElements();
  return (this->empty() ?
          std::numeric_limits<DataType>::max() :
          *std::min_element(start, start + n));
}

//------------------------------------------------------------------------------
// Maximum.
// LOCAL to processor!
//-------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE 
inline
DataType
FieldView<Dimension, DataType>::
localMax(const bool includeGhosts) const {
  const auto* start = &this->front();
  const auto n = includeGhosts ? numElements() : numInternalElements();
  return (this->empty() ?
          std::numeric_limits<DataType>::lowest() :
          *std::max_element(start, start + n));
}

//------------------------------------------------------------------------------
// operator==(FieldView)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE 
inline
bool
FieldView<Dimension, DataType>::
operator==(const FieldView<Dimension, DataType>& rhs) const {
  const auto n = mDataSpan.size();
  if (rhs.numElements() != n) return false;
  auto result = true;
  size_t i = 0u;
  while (i < n and result) {
    result = (*this)[i] == rhs[i];
    ++i;
  }
  return result;
}

//------------------------------------------------------------------------------
// operator!=(FieldView)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE 
inline
bool
FieldView<Dimension, DataType>::
operator!=(const FieldView<Dimension, DataType>& rhs) const {
  return !((*this) == rhs);
}

//------------------------------------------------------------------------------
// operator==(value)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE 
inline
bool
FieldView<Dimension, DataType>::
operator==(const DataType& rhs) const {
  const auto n = mDataSpan.size();
  auto result = true;
  size_t i = 0u;
  while (i < n and result) {
    result = (*this)[i] == rhs;
    ++i;
  }
  return result;
}

//------------------------------------------------------------------------------
// operator!=(value)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE 
inline
bool
FieldView<Dimension, DataType>::
operator!=(const DataType& rhs) const {
  return !((*this) == rhs);
}

//------------------------------------------------------------------------------
// operator>(value)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE 
inline
bool
FieldView<Dimension, DataType>::
operator>(const DataType& rhs) const {
  const auto n = mDataSpan.size();
  auto result = true;
  size_t i = 0u;
  while (i < n and result) {
    result = (*this)[i] > rhs;
    ++i;
  }
  return result;
}

//------------------------------------------------------------------------------
// operator<(value)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE 
inline
bool
FieldView<Dimension, DataType>::
operator<(const DataType& rhs) const {
  const auto n = mDataSpan.size();
  auto result = true;
  size_t i = 0u;
  while (i < n and result) {
    result = (*this)[i] < rhs;
    ++i;
  }
  return result;
}

//------------------------------------------------------------------------------
// operator>=(value)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE 
inline
bool
FieldView<Dimension, DataType>::
operator>=(const DataType& rhs) const {
  const auto n = mDataSpan.size();
  auto result = true;
  size_t i = 0u;
  while (i < n and result) {
    result = (*this)[i] >= rhs;
    ++i;
  }
  return result;
}

//------------------------------------------------------------------------------
// operator<=(value)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE 
inline
bool
FieldView<Dimension, DataType>::
operator<=(const DataType& rhs) const {
  const auto n = mDataSpan.size();
  auto result = true;
  size_t i = 0u;
  while (i < n and result) {
    result = (*this)[i] <= rhs;
    ++i;
  }
  return result;
}

//------------------------------------------------------------------------------
// data pointer
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE 
inline
DataType*
FieldView<Dimension, DataType>::
data() const {
#ifdef SPHERAL_UNIFIED_MEMORY
  return mDataSpan.data();
#else
  return mDataSpan.getActivePointer();
#endif
}
  
//------------------------------------------------------------------------------
// data pointer
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST
inline
DataType*
FieldView<Dimension, DataType>::
data(chai::ExecutionSpace space,
     bool do_move) const {
#ifdef SPHERAL_UNIFIED_MEMORY
  return mDataSpan.data();
#else
  return mDataSpan.data(space, do_move);
#endif
}
  
//------------------------------------------------------------------------------
// move
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST
inline
void
FieldView<Dimension, DataType>::
move(chai::ExecutionSpace space) {
#ifndef SPHERAL_UNIFIED_MEMORY
  mDataSpan.move(space);
#endif
}
  
//------------------------------------------------------------------------------
// shallowCopy
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE 
inline
void
FieldView<Dimension, DataType>::
shallowCopy(FieldView const& other) const {
#ifdef SPHERAL_UNIFIED_MEMORY
  const_cast<ContainerType&>(mDataSpan) = other.mDataSpan;
#else
  mDataSpan.shallowCopy(other.mDataSpan);
#endif
}
  
//------------------------------------------------------------------------------
// touch
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST
inline
void
FieldView<Dimension, DataType>::
touch(chai::ExecutionSpace space) {
#ifndef SPHERAL_UNIFIED_MEMORY
  mDataSpan.registerTouch(space);
#endif
}

//****************************** Global Functions ******************************
//------------------------------------------------------------------------------
// Output (ostream) operator.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST 
std::ostream&
operator<<(std::ostream& os, const FieldView<Dimension, DataType>& fieldSpan) {

  // Write the number of internal elements.
  os << fieldSpan.numInternalElements() << " ";

  // Write the internal elements.
  for (auto itr = fieldSpan.internalBegin(); itr < fieldSpan.internalEnd(); ++itr) {
    os << *itr << " ";
  }
//   os << endl;
  return os;
}

} // namespace Spheral
