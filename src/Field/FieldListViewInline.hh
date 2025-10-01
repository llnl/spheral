// Includes.
#include "Geometry/MathTraits.hh"

#include "Field/FieldView.hh"
#include "Field/FieldList.hh"
#include "Distributed/allReduce.hh"

#include <algorithm>
#include <limits>

namespace Spheral {

//------------------------------------------------------------------------------
// Empty constructor.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
FieldListView<Dimension, DataType>::
FieldListView():
  mFieldViews(),
  mChaiCallback([](const chai::PointerRecord*, chai::Action, chai::ExecutionSpace) {}) {
  mFieldViews.setUserCallback(getCallback());
  DEBUG_LOG << "FieldListView::FieldListView() : " << this;
}

//------------------------------------------------------------------------------
// Construct from a FieldList
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST
inline
FieldListView<Dimension, DataType>::
FieldListView(const FieldList<Dimension, DataType>& rhs):
  mFieldViews(),
  mChaiCallback([](const chai::PointerRecord*, chai::Action, chai::ExecutionSpace) {}) {
  const auto n = rhs.size();
  if (n > 0u) {
    mFieldViews.reallocate(n);
    for (auto i = 0u; i < n; ++i) mFieldViews[i] = rhs[i]->view();
  } else {
    mFieldViews.free();
  }
  mFieldViews.setUserCallback(getCallback());
  ENSURE(this->size() == rhs.size());
  DEBUG_LOG << "FieldListView::FieldListView(const FieldList& " << &rhs << ") : " << this;
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
FieldListView<Dimension, DataType>::
~FieldListView() {
  // mFieldViews.free();
  DEBUG_LOG << "FieldListView::~FieldListView() : " << this;
}

//------------------------------------------------------------------------------
// Assignment with a constant.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
FieldListView<Dimension, DataType>&
FieldListView<Dimension, DataType>::
operator=(const DataType& rhs) {
  for (auto& x: mFieldViews) x = rhs;
  return *this;
}

//------------------------------------------------------------------------------
// Index operator.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
typename FieldListView<Dimension, DataType>::value_type&
FieldListView<Dimension, DataType>::
operator[](const size_t index) const {
  REQUIRE2(index < this->size(), "FieldListView index ERROR: out of bounds " << index << " !< " << this->size());
  return mFieldViews[index];
}

//------------------------------------------------------------------------------
// at version, for consistency with the STL.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
typename FieldListView<Dimension, DataType>::value_type&
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
  REQUIRE2(fieldIndex < mFieldViews.size(), "FieldListView index ERROR: out of bounds " << fieldIndex << " !< " << mFieldViews.size());
  REQUIRE2(nodeIndex < mFieldViews[fieldIndex].numElements(), "FieldListView node index ERROR: out of bounds " << nodeIndex << " !< " << mFieldViews[fieldIndex].numElements());
  return mFieldViews[fieldIndex][nodeIndex];
}

//------------------------------------------------------------------------------
// Apply a minimum data value.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
void
FieldListView<Dimension, DataType>::
applyMin(const DataType& dataMin) {
  for (auto& x: mFieldViews) x.applyMin(dataMin);
}

//------------------------------------------------------------------------------
// Apply a maximum data value.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
void
FieldListView<Dimension, DataType>::
applyMax(const DataType& dataMax) {
  for (auto& x: mFieldViews) x.applyMax(dataMax);
}

//------------------------------------------------------------------------------
// Apply a (scalar) minimum data value.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
void
FieldListView<Dimension, DataType>::
applyScalarMin(const Scalar dataMin) {
  for (auto& x: mFieldViews) x.applyScalarMin(dataMin);
}

//------------------------------------------------------------------------------
// Apply a (scalar) maximum data value.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
void
FieldListView<Dimension, DataType>::
applyScalarMax(const Scalar dataMax) {
  for (auto& x: mFieldViews) x.applyScalarMax(dataMax);
}

//------------------------------------------------------------------------------
// Add two FieldListViews in place.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
FieldListView<Dimension, DataType>&
FieldListView<Dimension, DataType>::
operator+=(const FieldListView<Dimension, DataType>& rhs) {

  // Pre-conditions.
  const auto n = this->size();
  BEGIN_CONTRACT_SCOPE
  {
    REQUIRE(rhs.size() == n);
    for (size_t i = 0u; i < n; ++i) REQUIRE(mFieldViews[i].numElements() == rhs[i].numElements());
  }
  END_CONTRACT_SCOPE

  for (size_t i = 0u; i < n; ++i) (*this)[i] += rhs[i];
  return *this;
}

//------------------------------------------------------------------------------
// Subtract a FieldList from another in place.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
FieldListView<Dimension, DataType>&
FieldListView<Dimension, DataType>::
operator-=(const FieldListView<Dimension, DataType>& rhs) {

  // Pre-conditions.
  const auto n = this->size();
  BEGIN_CONTRACT_SCOPE
  {
    REQUIRE(rhs.size() == n);
    for (size_t i = 0u; i < n; ++i) REQUIRE(mFieldViews[i].numElements() == rhs[i].numElements());
  }
  END_CONTRACT_SCOPE

  for (size_t i = 0u; i < n; ++i) (*this)[i] -= rhs[i];
  return *this;
}

//------------------------------------------------------------------------------
// Add a single value to the FieldListView in place.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
FieldListView<Dimension, DataType>&
FieldListView<Dimension, DataType>::
operator+=(const DataType& rhs) {
  const auto n = this->size();
  for (size_t i = 0u; i < n; ++i) (*this)[i] += rhs;
  return *this;
}

//------------------------------------------------------------------------------
// Subtract a single value from the FieldListView in place.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
FieldListView<Dimension, DataType>&
FieldListView<Dimension, DataType>::
operator-=(const DataType& rhs) {
  const auto n = this->size();
  for (size_t i = 0u; i < n; ++i) (*this)[i] -= rhs;
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
    for (size_t i = 0u; i < n; ++i) REQUIRE(mFieldViews[i].numElements() == rhs[i].numElements());
  }
  END_CONTRACT_SCOPE

  for (size_t i = 0u; i < n; ++i) (*this)[i] *= rhs[i];
  return *this;
}

//------------------------------------------------------------------------------
// Multiply this FieldListView by a Scalar in place.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
FieldListView<Dimension, DataType>&
FieldListView<Dimension, DataType>::
operator*=(const Scalar& rhs) {
  const auto n = this->size();
  for (size_t i = 0u; i < n; ++i) (*this)[i] *= rhs;
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
    for (size_t i = 0u; i < n; ++i) REQUIRE(mFieldViews[i].numElements() == rhs[i].numElements());
  }
  END_CONTRACT_SCOPE

  for (size_t i = 0u; i < n; ++i) (*this)[i] /= rhs[i];
  return *this;
}

//------------------------------------------------------------------------------
// Divide this FieldListView by a Scalar in place.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
FieldListView<Dimension, DataType>&
FieldListView<Dimension, DataType>::
operator/=(const typename Dimension::Scalar& rhs) {
  const auto n = this->size();
  for (size_t i = 0u; i < n; ++i) (*this)[i] /= rhs;
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
localSumElements(const bool includeGhosts) const {
  auto result = DataTypeTraits<DataType>::zero();
  for (auto& x: mFieldViews) result += x.localSumElements(includeGhosts);
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
localMin(const bool includeGhosts) const {
  auto result = std::numeric_limits<DataType>::max();
  for (auto& x: mFieldViews) result = std::min(result, x.localMin(includeGhosts));
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
localMax(const bool includeGhosts) const {
  auto result = std::numeric_limits<DataType>::lowest();
  for (auto& x: mFieldViews) result = std::max(result, x.localMax(includeGhosts));
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
    for (size_t i = 0u; i < n; ++i) REQUIRE(mFieldViews[i].numElements() == rhs[i].numElements());
  }
  END_CONTRACT_SCOPE

  auto result = true;
  size_t i = 0u;
  while (result and i < n) {
    result = (mFieldViews[i] == rhs[i]);
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
    result = (mFieldViews[i] == rhs);
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
    result = (mFieldViews[i] > rhs);
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
    result = (mFieldViews[i] < rhs);
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
    result = (mFieldViews[i] >= rhs);
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
    result = (mFieldViews[i] <= rhs);
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
  for (auto& x: mFieldViews) result += x.numElements();
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
  for (auto& x: mFieldViews) result += x.numInternalElements();
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
  for (auto& x: mFieldViews) result += x.numGhostElements();
  return result;
}

//------------------------------------------------------------------------------
// move
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST
inline
void
FieldListView<Dimension, DataType>::
move(chai::ExecutionSpace space, bool recursive) {
  mFieldViews.move(space);
  if (recursive) {
    for (auto& d: mFieldViews) {
      d.move(space);
    }
  }
}

//------------------------------------------------------------------------------
// data
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST_DEVICE
inline
typename FieldListView<Dimension, DataType>::value_type*
FieldListView<Dimension, DataType>::
data() const {
  return mFieldViews.data();
}

//------------------------------------------------------------------------------
// data
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
SPHERAL_HOST
inline
typename FieldListView<Dimension, DataType>::value_type*
FieldListView<Dimension, DataType>::
data(chai::ExecutionSpace space, bool do_move) const {
  return mFieldViews.data(space, do_move);
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
  mFieldViews.registerTouch(space);
  if (recursive) {
    for (auto& d : mFieldViews) {
      d.touch(space);
    }
  }
}

//------------------------------------------------------------------------------
// Default callback action to be used with chai Managed containers. An
// additional calback can be passed to extend this functionality. Useful for
// debugging, testing and probing for performance counters / metrics.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
auto
FieldListView<Dimension, DataType>::
getCallback() {
  return [callback = mChaiCallback](
    const chai::PointerRecord * record,
    chai::Action action,
    chai::ExecutionSpace space) {
      if (action == chai::ACTION_MOVE) {
        if (space == chai::CPU)
          DEBUG_LOG << "FieldList : MOVED to the CPU";
        if (space == chai::GPU)
          DEBUG_LOG << "FieldList : MOVED to the GPU";
      }
      else if (action == chai::ACTION_ALLOC) {
        if (space == chai::CPU)
          DEBUG_LOG << "FieldList : ALLOC on the CPU";
        if (space == chai::GPU)
          DEBUG_LOG << "FieldList : ALLOC on the GPU";
      }
      else if (action == chai::ACTION_FREE) {
        if (space == chai::CPU)
          DEBUG_LOG << "FieldList : FREE on the CPU";
        if (space == chai::GPU)
          DEBUG_LOG << "FieldList : FREE on the GPU";
      }
      callback(record, action, space);
    };
}

}
