// Includes.
#include "Geometry/MathTraits.hh"

#include "Field/FieldSpan.hh"
#include "Distributed/allReduce.hh"

#include <algorithm>
#include <limits>

namespace Spheral {

//------------------------------------------------------------------------------
// Assignment with a constant.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
FieldSpanList<Dimension, DataType>&
FieldSpanList<Dimension, DataType>::
operator=(const DataType& rhs) {
  for (auto* fspan: mSpanFieldSpans) *fspan = rhs;
  return *this;
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST
inline
FieldSpanList<Dimension, DataType>::
~FieldSpanList() {
#ifndef SPHERAL_UNIFIED_MEMORY
  mSpanFieldSpans.free();
#endif
}

//------------------------------------------------------------------------------
// Index operator.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
typename FieldSpanList<Dimension, DataType>::value_type
FieldSpanList<Dimension, DataType>::
operator[](const size_t index) const {
  REQUIRE2(index < this->size(), "FieldSpanList index ERROR: out of bounds " << index << " !< " << this->size());
  return mSpanFieldSpans[index];
}

//------------------------------------------------------------------------------
// at version, for consistency with the STL.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
typename FieldSpanList<Dimension, DataType>::value_type
FieldSpanList<Dimension, DataType>::
at(const size_t index) const {
  return (*this)[index];
}

//------------------------------------------------------------------------------
// Provide direct access to FieldSpan elements
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
DataType&
FieldSpanList<Dimension, DataType>::
operator()(const size_t fieldIndex,
           const size_t nodeIndex) const {
  REQUIRE2(fieldIndex < mSpanFieldSpans.size(), "FieldSpanList index ERROR: out of bounds " << fieldIndex << " !< " << mSpanFieldSpans.size());
  REQUIRE2(nodeIndex < mSpanFieldSpans[fieldIndex]->numElements(), "FieldSpanList node index ERROR: out of bounds " << nodeIndex << " !< " << mSpanFieldSpans[fieldIndex]->numElements());
  return (*mSpanFieldSpans[fieldIndex])[nodeIndex];
}

//------------------------------------------------------------------------------
// Apply a minimum data value.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
void
FieldSpanList<Dimension, DataType>::applyMin(const DataType& dataMin) {
  for (auto* x: mSpanFieldSpans) x->applyMin(dataMin);
}

//------------------------------------------------------------------------------
// Apply a maximum data value.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
void
FieldSpanList<Dimension, DataType>::applyMax(const DataType& dataMax) {
  for (auto* x: mSpanFieldSpans) x->applyMax(dataMax);
}

//------------------------------------------------------------------------------
// Apply a (scalar) minimum data value.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
void
FieldSpanList<Dimension, DataType>::applyScalarMin(const Scalar dataMin) {
  for (auto* x: mSpanFieldSpans) x->applyScalarMin(dataMin);
}

//------------------------------------------------------------------------------
// Apply a (scalar) maximum data value.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
void
FieldSpanList<Dimension, DataType>::applyScalarMax(const Scalar dataMax) {
  for (auto x: mSpanFieldSpans) x->applyScalarMax(dataMax);
}

//------------------------------------------------------------------------------
// Add two FieldSpanLists in place.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
FieldSpanList<Dimension, DataType>&
FieldSpanList<Dimension, DataType>::operator+=(const FieldSpanList<Dimension, DataType>& rhs) {

  // Pre-conditions.
  const auto n = this->size();
  BEGIN_CONTRACT_SCOPE
  {
    REQUIRE(rhs.size() == n);
    for (size_t i = 0u; i < n; ++i) REQUIRE(mSpanFieldSpans[i]->numElements() == rhs[i]->numElements());
  }
  END_CONTRACT_SCOPE

  for (size_t i = 0u; i < n; ++i) *(*this)[i] += *rhs[i];
  return *this;
}

//------------------------------------------------------------------------------
// Subtract a FieldList from another in place.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
FieldSpanList<Dimension, DataType>&
FieldSpanList<Dimension, DataType>::operator-=(const FieldSpanList<Dimension, DataType>& rhs) {

  // Pre-conditions.
  const auto n = this->size();
  BEGIN_CONTRACT_SCOPE
  {
    REQUIRE(rhs.size() == n);
    for (size_t i = 0u; i < n; ++i) REQUIRE(mSpanFieldSpans[i]->numElements() == rhs[i]->numElements());
  }
  END_CONTRACT_SCOPE

  for (size_t i = 0u; i < n; ++i) *(*this)[i] -= *rhs[i];
  return *this;
}

//------------------------------------------------------------------------------
// Add a single value to the FieldSpanList in place.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
FieldSpanList<Dimension, DataType>&
FieldSpanList<Dimension, DataType>::operator+=(const DataType& rhs) {
  const auto n = this->size();
  for (size_t i = 0u; i < n; ++i) *(*this)[i] += rhs;
  return *this;
}

//------------------------------------------------------------------------------
// Subtract a single value from the FieldSpanList in place.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
FieldSpanList<Dimension, DataType>&
FieldSpanList<Dimension, DataType>::operator-=(const DataType& rhs) {
  const auto n = this->size();
  for (size_t i = 0u; i < n; ++i) *(*this)[i] -= rhs;
  return *this;
}

//------------------------------------------------------------------------------
// Multiply this FieldSpanList by a Scalar FieldSpanList in place.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
FieldSpanList<Dimension, DataType>&
FieldSpanList<Dimension, DataType>::
operator*=(const FieldSpanList<Dimension, typename Dimension::Scalar>& rhs) {

  // Pre-conditions.
  const auto n = this->size();
  BEGIN_CONTRACT_SCOPE
  {
    REQUIRE(rhs.size() == n);
    for (size_t i = 0u; i < n; ++i) REQUIRE(mSpanFieldSpans[i]->numElements() == rhs[i]->numElements());
  }
  END_CONTRACT_SCOPE

  for (size_t i = 0u; i < n; ++i) *(*this)[i] *= *rhs[i];
  return *this;
}

//------------------------------------------------------------------------------
// Multiply this FieldSpanList by a Scalar in place.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
FieldSpanList<Dimension, DataType>&
FieldSpanList<Dimension, DataType>::operator*=(const Scalar& rhs) {
  const auto n = this->size();
  for (size_t i = 0u; i < n; ++i) *(*this)[i] *= rhs;
  return *this;
}

//------------------------------------------------------------------------------
// Divide this FieldSpanList by a Scalar FieldSpanList in place.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
FieldSpanList<Dimension, DataType>&
FieldSpanList<Dimension, DataType>::
operator/=(const FieldSpanList<Dimension, typename Dimension::Scalar>& rhs) {

  // Pre-conditions.
  const auto n = this->size();
  BEGIN_CONTRACT_SCOPE
  {
    REQUIRE(rhs.size() == n);
    for (size_t i = 0u; i < n; ++i) REQUIRE(mSpanFieldSpans[i]->numElements() == rhs[i]->numElements());
  }
  END_CONTRACT_SCOPE

  for (size_t i = 0u; i < n; ++i) *(*this)[i] /= *rhs[i];
  return *this;
}

//------------------------------------------------------------------------------
// Divide this FieldSpanList by a Scalar in place.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
FieldSpanList<Dimension, DataType>&
FieldSpanList<Dimension, DataType>::operator/=(const typename Dimension::Scalar& rhs) {
  const auto n = this->size();
  for (size_t i = 0u; i < n; ++i) *(*this)[i] /= rhs;
  return *this;
}

//------------------------------------------------------------------------------
// Sum the field elements.
// LOCAL to processor!
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
DataType
FieldSpanList<Dimension, DataType>::
localSumElements() const {
  auto result = DataTypeTraits<DataType>::zero();
  for (auto* x: mSpanFieldSpans) result += x->localSumElements();
  return result;
}

//------------------------------------------------------------------------------
// Find the minimum.
// LOCAL to processor!
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
DataType
FieldSpanList<Dimension, DataType>::
localMin() const {
  auto result = std::numeric_limits<DataType>::max();
  for (auto* x: mSpanFieldSpans) result = std::min(result, x->localMin());
  return result;
}

//------------------------------------------------------------------------------
// Find the maximum.
// LOCAL to processor!
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
DataType
FieldSpanList<Dimension, DataType>::
localMax() const {
  auto result = std::numeric_limits<DataType>::lowest();
  for (auto* x: mSpanFieldSpans) result = std::max(result, x->localMax());
  return result;
}

//------------------------------------------------------------------------------
// operator==
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
bool
FieldSpanList<Dimension, DataType>::
operator==(const FieldSpanList<Dimension, DataType>& rhs) const {
  // Pre-conditions.
  const auto n = this->size();
  BEGIN_CONTRACT_SCOPE
  {
    REQUIRE(rhs.size() == n);
    for (size_t i = 0u; i < n; ++i) REQUIRE(mSpanFieldSpans[i]->numElements() == rhs[i]->numElements());
  }
  END_CONTRACT_SCOPE

  auto result = true;
  size_t i = 0u;
  while (result and i < n) {
    result = *mSpanFieldSpans[i] == *rhs[i];
    ++i;
  }
  return result;
}

//------------------------------------------------------------------------------
// operator!=
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
bool
FieldSpanList<Dimension, DataType>::
operator!=(const FieldSpanList<Dimension, DataType>& rhs) const {
  return !(operator==(rhs));
}

//------------------------------------------------------------------------------
// operator==
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
bool
FieldSpanList<Dimension, DataType>::
operator==(const DataType& rhs) const {
  const auto n = this->size();
  bool result = true;
  size_t i = 0u;
  while (result and i < n) {
    result = *mSpanFieldSpans[i] == rhs;
    ++i;
  }
  return result;
}

//------------------------------------------------------------------------------
// operator!=
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
bool
FieldSpanList<Dimension, DataType>::
operator!=(const DataType& rhs) const {
  return !(operator==(rhs));
}

//------------------------------------------------------------------------------
// operator>
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
bool
FieldSpanList<Dimension, DataType>::
operator>(const DataType& rhs) const {
  const auto n = this->size();
  bool result = true;
  size_t i = 0u;
  while (result and i < n) {
    result = *mSpanFieldSpans[i] > rhs;
    ++i;
  }
  return result;
}

//------------------------------------------------------------------------------
// operator<
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
bool
FieldSpanList<Dimension, DataType>::
operator<(const DataType& rhs) const {
  const auto n = this->size();
  bool result = true;
  size_t i = 0u;
  while (result and i < n) {
    result = *mSpanFieldSpans[i] < rhs;
    ++i;
  }
  return result;
}

//------------------------------------------------------------------------------
// operator>=
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
bool
FieldSpanList<Dimension, DataType>::
operator>=(const DataType& rhs) const {
  const auto n = this->size();
  bool result = true;
  size_t i = 0u;
  while (result and i < n) {
    result = *mSpanFieldSpans[i] >= rhs;
    ++i;
  }
  return result;
}

//------------------------------------------------------------------------------
// operator<=
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
bool
FieldSpanList<Dimension, DataType>::
operator<=(const DataType& rhs) const {
  const auto n = this->size();
  bool result = true;
  size_t i = 0u;
  while (result and i < n) {
    result = *mSpanFieldSpans[i] <= rhs;
    ++i;
  }
  return result;
}

//------------------------------------------------------------------------------
// numElements
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
size_t
FieldSpanList<Dimension, DataType>::
numElements() const {
  size_t result = 0u;
  for (auto* x: mSpanFieldSpans) result += x->numElements();
  return result;
}

//------------------------------------------------------------------------------
// numInternalElements
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
size_t
FieldSpanList<Dimension, DataType>::
numInternalElements() const {
  size_t result = 0u;
  for (auto* x: mSpanFieldSpans) result += x->numInternalElements();
  return result;
}

//------------------------------------------------------------------------------
// numGhostElements
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
size_t
FieldSpanList<Dimension, DataType>::
numGhostElements() const {
  size_t result = 0u;
  for (auto* x: mSpanFieldSpans) result += x->numGhostElements();
  return result;
}

#ifndef SPHERAL_UNIFIED_MEMORY

//------------------------------------------------------------------------------
// move
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
void
FieldSpanList<Dimension, DataType>::
move(chai::ExecutionSpace space, bool recursive) {
  mSpanFieldSpans.move(space);
  if (recursive) {
    for (auto& d: mSpanFieldSpans) {
      d.move(space);
    }
  }
}

//------------------------------------------------------------------------------
// touch
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
void
FieldSpanList<Dimension, DataType>::
touch(chai::ExecutionSpace space, bool recursive) {
  mSpanFieldSpans.registerTouch(space);
  if (recursive) {
    for (auto& d : mSpanFieldSpans) {
      d.touch(space);
    }
  }
}

#endif

}
