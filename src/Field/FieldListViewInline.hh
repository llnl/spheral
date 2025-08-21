// Includes.
#include "Geometry/MathTraits.hh"

#include "Field/FieldView.hh"
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
FieldListView<Dimension, DataType>&
FieldListView<Dimension, DataType>::
operator=(const DataType& rhs) {
  for (auto* fspan: mSpanFieldViews) *fspan = rhs;
  return *this;
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST
inline
FieldListView<Dimension, DataType>::
~FieldListView() {
#ifndef SPHERAL_UNIFIED_MEMORY
  mSpanFieldViews.free();
#endif
}

//------------------------------------------------------------------------------
// Index operator.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
typename FieldListView<Dimension, DataType>::value_type
FieldListView<Dimension, DataType>::
operator[](const size_t index) const {
  REQUIRE2(index < this->size(), "FieldListView index ERROR: out of bounds " << index << " !< " << this->size());
  return mSpanFieldViews[index];
}

//------------------------------------------------------------------------------
// at version, for consistency with the STL.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
typename FieldListView<Dimension, DataType>::value_type
FieldListView<Dimension, DataType>::
at(const size_t index) const {
  return (*this)[index];
}

//------------------------------------------------------------------------------
// Provide direct access to FieldView elements
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
DataType&
FieldListView<Dimension, DataType>::
operator()(const size_t fieldIndex,
           const size_t nodeIndex) const {
  REQUIRE2(fieldIndex < mSpanFieldViews.size(), "FieldListView index ERROR: out of bounds " << fieldIndex << " !< " << mSpanFieldViews.size());
  REQUIRE2(nodeIndex < mSpanFieldViews[fieldIndex]->numElements(), "FieldListView node index ERROR: out of bounds " << nodeIndex << " !< " << mSpanFieldViews[fieldIndex]->numElements());
  return (*mSpanFieldViews[fieldIndex])[nodeIndex];
}

//------------------------------------------------------------------------------
// Apply a minimum data value.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
void
FieldListView<Dimension, DataType>::applyMin(const DataType& dataMin) {
  for (auto* x: mSpanFieldViews) x->applyMin(dataMin);
}

//------------------------------------------------------------------------------
// Apply a maximum data value.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
void
FieldListView<Dimension, DataType>::applyMax(const DataType& dataMax) {
  for (auto* x: mSpanFieldViews) x->applyMax(dataMax);
}

//------------------------------------------------------------------------------
// Apply a (scalar) minimum data value.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
void
FieldListView<Dimension, DataType>::applyScalarMin(const Scalar dataMin) {
  for (auto* x: mSpanFieldViews) x->applyScalarMin(dataMin);
}

//------------------------------------------------------------------------------
// Apply a (scalar) maximum data value.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
void
FieldListView<Dimension, DataType>::applyScalarMax(const Scalar dataMax) {
  for (auto x: mSpanFieldViews) x->applyScalarMax(dataMax);
}

//------------------------------------------------------------------------------
// Add two FieldListViews in place.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
FieldListView<Dimension, DataType>&
FieldListView<Dimension, DataType>::operator+=(const FieldListView<Dimension, DataType>& rhs) {

  // Pre-conditions.
  const auto n = this->size();
  BEGIN_CONTRACT_SCOPE
  {
    REQUIRE(rhs.size() == n);
    for (size_t i = 0u; i < n; ++i) REQUIRE(mSpanFieldViews[i]->numElements() == rhs[i]->numElements());
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
FieldListView<Dimension, DataType>&
FieldListView<Dimension, DataType>::operator-=(const FieldListView<Dimension, DataType>& rhs) {

  // Pre-conditions.
  const auto n = this->size();
  BEGIN_CONTRACT_SCOPE
  {
    REQUIRE(rhs.size() == n);
    for (size_t i = 0u; i < n; ++i) REQUIRE(mSpanFieldViews[i]->numElements() == rhs[i]->numElements());
  }
  END_CONTRACT_SCOPE

  for (size_t i = 0u; i < n; ++i) *(*this)[i] -= *rhs[i];
  return *this;
}

//------------------------------------------------------------------------------
// Add a single value to the FieldListView in place.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
FieldListView<Dimension, DataType>&
FieldListView<Dimension, DataType>::operator+=(const DataType& rhs) {
  const auto n = this->size();
  for (size_t i = 0u; i < n; ++i) *(*this)[i] += rhs;
  return *this;
}

//------------------------------------------------------------------------------
// Subtract a single value from the FieldListView in place.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
FieldListView<Dimension, DataType>&
FieldListView<Dimension, DataType>::operator-=(const DataType& rhs) {
  const auto n = this->size();
  for (size_t i = 0u; i < n; ++i) *(*this)[i] -= rhs;
  return *this;
}

//------------------------------------------------------------------------------
// Multiply this FieldListView by a Scalar FieldListView in place.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
FieldListView<Dimension, DataType>&
FieldListView<Dimension, DataType>::
operator*=(const FieldListView<Dimension, typename Dimension::Scalar>& rhs) {

  // Pre-conditions.
  const auto n = this->size();
  BEGIN_CONTRACT_SCOPE
  {
    REQUIRE(rhs.size() == n);
    for (size_t i = 0u; i < n; ++i) REQUIRE(mSpanFieldViews[i]->numElements() == rhs[i]->numElements());
  }
  END_CONTRACT_SCOPE

  for (size_t i = 0u; i < n; ++i) *(*this)[i] *= *rhs[i];
  return *this;
}

//------------------------------------------------------------------------------
// Multiply this FieldListView by a Scalar in place.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
FieldListView<Dimension, DataType>&
FieldListView<Dimension, DataType>::operator*=(const Scalar& rhs) {
  const auto n = this->size();
  for (size_t i = 0u; i < n; ++i) *(*this)[i] *= rhs;
  return *this;
}

//------------------------------------------------------------------------------
// Divide this FieldListView by a Scalar FieldListView in place.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
FieldListView<Dimension, DataType>&
FieldListView<Dimension, DataType>::
operator/=(const FieldListView<Dimension, typename Dimension::Scalar>& rhs) {

  // Pre-conditions.
  const auto n = this->size();
  BEGIN_CONTRACT_SCOPE
  {
    REQUIRE(rhs.size() == n);
    for (size_t i = 0u; i < n; ++i) REQUIRE(mSpanFieldViews[i]->numElements() == rhs[i]->numElements());
  }
  END_CONTRACT_SCOPE

  for (size_t i = 0u; i < n; ++i) *(*this)[i] /= *rhs[i];
  return *this;
}

//------------------------------------------------------------------------------
// Divide this FieldListView by a Scalar in place.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
FieldListView<Dimension, DataType>&
FieldListView<Dimension, DataType>::operator/=(const typename Dimension::Scalar& rhs) {
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
FieldListView<Dimension, DataType>::
localSumElements() const {
  auto result = DataTypeTraits<DataType>::zero();
  for (auto* x: mSpanFieldViews) result += x->localSumElements();
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
FieldListView<Dimension, DataType>::
localMin() const {
  auto result = std::numeric_limits<DataType>::max();
  for (auto* x: mSpanFieldViews) result = std::min(result, x->localMin());
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
FieldListView<Dimension, DataType>::
localMax() const {
  auto result = std::numeric_limits<DataType>::lowest();
  for (auto* x: mSpanFieldViews) result = std::max(result, x->localMax());
  return result;
}

//------------------------------------------------------------------------------
// operator==
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
bool
FieldListView<Dimension, DataType>::
operator==(const FieldListView<Dimension, DataType>& rhs) const {
  // Pre-conditions.
  const auto n = this->size();
  BEGIN_CONTRACT_SCOPE
  {
    REQUIRE(rhs.size() == n);
    for (size_t i = 0u; i < n; ++i) REQUIRE(mSpanFieldViews[i]->numElements() == rhs[i]->numElements());
  }
  END_CONTRACT_SCOPE

  auto result = true;
  size_t i = 0u;
  while (result and i < n) {
    result = *mSpanFieldViews[i] == *rhs[i];
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
FieldListView<Dimension, DataType>::
operator!=(const FieldListView<Dimension, DataType>& rhs) const {
  return !(operator==(rhs));
}

//------------------------------------------------------------------------------
// operator==
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
bool
FieldListView<Dimension, DataType>::
operator==(const DataType& rhs) const {
  const auto n = this->size();
  bool result = true;
  size_t i = 0u;
  while (result and i < n) {
    result = *mSpanFieldViews[i] == rhs;
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
FieldListView<Dimension, DataType>::
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
FieldListView<Dimension, DataType>::
operator>(const DataType& rhs) const {
  const auto n = this->size();
  bool result = true;
  size_t i = 0u;
  while (result and i < n) {
    result = *mSpanFieldViews[i] > rhs;
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
FieldListView<Dimension, DataType>::
operator<(const DataType& rhs) const {
  const auto n = this->size();
  bool result = true;
  size_t i = 0u;
  while (result and i < n) {
    result = *mSpanFieldViews[i] < rhs;
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
FieldListView<Dimension, DataType>::
operator>=(const DataType& rhs) const {
  const auto n = this->size();
  bool result = true;
  size_t i = 0u;
  while (result and i < n) {
    result = *mSpanFieldViews[i] >= rhs;
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
FieldListView<Dimension, DataType>::
operator<=(const DataType& rhs) const {
  const auto n = this->size();
  bool result = true;
  size_t i = 0u;
  while (result and i < n) {
    result = *mSpanFieldViews[i] <= rhs;
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
FieldListView<Dimension, DataType>::
numElements() const {
  size_t result = 0u;
  for (auto* x: mSpanFieldViews) result += x->numElements();
  return result;
}

//------------------------------------------------------------------------------
// numInternalElements
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
size_t
FieldListView<Dimension, DataType>::
numInternalElements() const {
  size_t result = 0u;
  for (auto* x: mSpanFieldViews) result += x->numInternalElements();
  return result;
}

//------------------------------------------------------------------------------
// numGhostElements
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
size_t
FieldListView<Dimension, DataType>::
numGhostElements() const {
  size_t result = 0u;
  for (auto* x: mSpanFieldViews) result += x->numGhostElements();
  return result;
}

#ifndef SPHERAL_UNIFIED_MEMORY

//------------------------------------------------------------------------------
// move
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST
inline
void
FieldListView<Dimension, DataType>::
move(chai::ExecutionSpace space, bool recursive) {
  mSpanFieldViews.move(space);
  if (recursive) {
    for (auto& d: mSpanFieldViews) {
      d.move(space);
    }
  }
}

//------------------------------------------------------------------------------
// touch
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST
inline
void
FieldListView<Dimension, DataType>::
touch(chai::ExecutionSpace space, bool recursive) {
  mSpanFieldViews.registerTouch(space);
  if (recursive) {
    for (auto& d : mSpanFieldViews) {
      d.touch(space);
    }
  }
}

#endif

}
