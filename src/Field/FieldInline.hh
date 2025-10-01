#include "Field/NodeIterators.hh"

#include "Geometry/Dimension.hh"
#include "NodeList/NodeList.hh"
#include "Utilities/packElement.hh"
#include "Utilities/removeElements.hh"
#include "Utilities/safeInv.hh"
#include "Distributed/allReduce.hh"
#include "Distributed/Communicator.hh"

#include <cmath>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <sstream>
#include <limits>

#include "config.hh"
#ifdef ENABLE_MPI
extern "C" {
#include <mpi.h>
}
#endif

// Inlined methods.
namespace Spheral {

//------------------------------------------------------------------------------
// Construct with name.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>::
Field(typename FieldBase<Dimension>::FieldName name):
  FieldBase<Dimension>(name),
  mDataArray() {
}

//------------------------------------------------------------------------------
// Construct with name and field values.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>::
Field(typename FieldBase<Dimension>::FieldName name,
      const Field<Dimension, DataType>& field):
  FieldBase<Dimension>(name, *field.nodeListPtr()),
  mDataArray(field.mDataArray) {}

//------------------------------------------------------------------------------
// Construct with the given name and NodeList.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>::
Field(typename FieldBase<Dimension>::FieldName name,
      const NodeList<Dimension>& nodeList):
  FieldBase<Dimension>(name, nodeList),
  mDataArray(nodeList.numNodes(), DataTypeTraits<DataType>::zero()) {
  REQUIRE(numElements() == nodeList.numNodes());
}

//------------------------------------------------------------------------------
// Construct with given name, NodeList, and value.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>::
Field(typename FieldBase<Dimension>::FieldName name,
      const NodeList<Dimension>& nodeList,
      DataType value):
  FieldBase<Dimension>(name, nodeList),
  mDataArray(nodeList.numNodes(), value) {
  REQUIRE(size() == nodeList.numNodes());
}

//------------------------------------------------------------------------------
// Construct for a given name, NodeList, and vector of values.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>::
Field(typename FieldBase<Dimension>::FieldName name, 
      const NodeList<Dimension>& nodeList,
      const std::vector<DataType,DataAllocator<DataType>>& array):
  FieldBase<Dimension>(name, nodeList),
  mDataArray(nodeList.numNodes()) {
  REQUIRE(size() == nodeList.numNodes());
  REQUIRE(size() == array.size());
  mDataArray = array;
}

//------------------------------------------------------------------------------
// Construct by copying the values of another Field, but using a different
// node list.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>::Field(const NodeList<Dimension>& nodeList,
                                  const Field<Dimension, DataType>& field):
  FieldBase<Dimension>(field.name(), nodeList),
  mDataArray(field.mDataArray) {
  ENSURE(size() == nodeList.numNodes());
}

//------------------------------------------------------------------------------
// Copy Constructor.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>::Field(const Field& field):
  FieldBase<Dimension>(field),
  mDataArray(field.mDataArray) {
}

//------------------------------------------------------------------------------
// The virtual clone method, allowing us to duplicate fields with just 
// FieldBase*.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
std::shared_ptr<FieldBase<Dimension> >
Field<Dimension, DataType>::clone() const {
  return std::shared_ptr<FieldBase<Dimension>>(new Field<Dimension, DataType>(*this));
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>::~Field() {
  mManagedData.free();
}

//------------------------------------------------------------------------------
// Assignment operator with FieldBase.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldBase<Dimension>&
Field<Dimension, DataType>::operator=(const FieldBase<Dimension>& rhs) {
  if (this != &rhs) {
    try {
      const auto* rhsPtr = dynamic_cast<const Field<Dimension, DataType>*>(&rhs);
      CHECK2(rhsPtr != 0, "Passed incorrect Field to operator=!");
      FieldBase<Dimension>::operator=(rhs);
      mDataArray = rhsPtr->mDataArray;
    } catch (const std::bad_cast &) {
      VERIFY2(false, "Attempt to assign a field to an incompatible field type.");
    }
  }
  return *this;
}

//------------------------------------------------------------------------------
// Assignment operator to another Field.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>&
Field<Dimension, DataType>::operator=(const Field<Dimension, DataType>& rhs) {
  if (this != &rhs) {
    FieldBase<Dimension>::operator=(rhs);
    mDataArray = rhs.mDataArray;
  }
  return *this;
}

//------------------------------------------------------------------------------
// Assigment operator with a vector.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>&
Field<Dimension, DataType>::operator=(const std::vector<DataType,DataAllocator<DataType>>& rhs) {
  REQUIRE(this->nodeList().numNodes() == rhs.size());
  mDataArray = rhs;
  return *this;
}

//------------------------------------------------------------------------------
// Assignment operator with a constant value of DataType
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>&
Field<Dimension, DataType>::operator=(const DataType& rhs) {
  std::fill(mDataArray.begin(), mDataArray.end(), rhs);
  return *this;
}

//------------------------------------------------------------------------------
// Test equivalence with a FieldBase.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
Field<Dimension, DataType>::operator==(const FieldBase<Dimension>& rhs) const {
  if (this->name() != rhs.name()) return false;
  if (this->nodeListPtr() != rhs.nodeListPtr()) return false;
  try {
    const auto* rhsPtr = dynamic_cast<const Field<Dimension, DataType>*>(&rhs);
    if (rhsPtr == nullptr) return false;
    return mDataArray == rhsPtr->mDataArray;
  } catch (const std::bad_cast &) {
    return false;
  }
}

//------------------------------------------------------------------------------
// Element access by integer index.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType&
Field<Dimension, DataType>::operator()(size_t index) {
  REQUIRE(index < numElements());
  return mDataArray[index];
}

template<typename Dimension, typename DataType>
inline
const DataType&
Field<Dimension, DataType>::operator()(size_t index) const {
  CHECK(index < numElements());
  return mDataArray[index];
}

template<typename Dimension, typename DataType>
inline
DataType&
Field<Dimension, DataType>::operator[](const size_t index) {
  return mDataArray[index];
}

template<typename Dimension, typename DataType>
inline
const DataType&
Field<Dimension, DataType>::operator[](const size_t index) const {
  return mDataArray[index];
}

//------------------------------------------------------------------------------
// at version, for consistency with STL interface.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType&
Field<Dimension, DataType>::at(size_t index) {
  return this->operator()(index);
}

template<typename Dimension, typename DataType>
inline
const DataType&
Field<Dimension, DataType>::at(size_t index) const {
  return this->operator()(index);
}

//------------------------------------------------------------------------------
// Element access by Node ID iterator.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType&
Field<Dimension, DataType>::
operator()(const NodeIteratorBase<Dimension>& itr) {
  CHECK(itr.nodeListPtr() == this->nodeListPtr());
  CHECK(itr.nodeID() >= 0 && itr.nodeID() < size());
  return mDataArray[itr.nodeID()];
}

template<typename Dimension, typename DataType>
inline
const DataType&
Field<Dimension, DataType>::
operator()(const NodeIteratorBase<Dimension>& itr) const {
  CHECK(itr.nodeListPtr() == this->nodeListPtr());
  CHECK(itr.nodeID() >= 0 && itr.nodeID() < size());
  return mDataArray[itr.nodeID()];
}

// //------------------------------------------------------------------------------
// // Element access by CoarseNodeIterator.
// //------------------------------------------------------------------------------
// template<typename Dimension, typename DataType>
// inline
// DataType&
// Field<Dimension, DataType>::
// operator()(const CoarseNodeIterator<Dimension>& itr) {
//   CHECK(itr.nodeListPtr() == nodeListPtr());
//   if (mNewCoarseNodes) {
//     cacheCoarseValues();
//     mNewCoarseNodes = false;
//   }
//   CHECK(itr.cacheID() >= 0 && itr.cacheID() < mCoarseCache.size());
//   return mCoarseCache[itr.cacheID()];
// }

// template<typename Dimension, typename DataType>
// inline
// const DataType&
// Field<Dimension, DataType>::
// operator()(const CoarseNodeIterator<Dimension>& itr) const {
//   CHECK(itr.nodeListPtr() == nodeListPtr());
//   if (mNewCoarseNodes) {
//     cacheCoarseValues();
//     mNewCoarseNodes = false;
//   }
//   CHECK(itr.cacheID() >= 0 && itr.cacheID() < mCoarseCache.size());
//   return mCoarseCache[itr.cacheID()];
// }

// //------------------------------------------------------------------------------
// // Element access by RefineNodeIterator.
// //------------------------------------------------------------------------------
// template<typename Dimension, typename DataType>
// inline
// DataType&
// Field<Dimension, DataType>::
// operator()(const RefineNodeIterator<Dimension>& itr) {
//   CHECK(itr.nodeListPtr() == nodeListPtr());
//   if (mNewRefineNodes) {
//     cacheRefineValues();
//     mNewRefineNodes = false;
//   }
//   CHECK(itr.cacheID() >= 0 && itr.cacheID() < mRefineCache.size());
//   return mRefineCache[itr.cacheID()];
// }

// template<typename Dimension, typename DataType>
// inline
// const DataType&
// Field<Dimension, DataType>::
// operator()(const RefineNodeIterator<Dimension>& itr) const {
//   CHECK(itr.nodeListPtr() == nodeListPtr());
//   if (mNewRefineNodes) {
//     cacheRefineValues();
//     mNewRefineNodes = false;
//   }
//   CHECK(itr.cacheID() >= 0 && itr.cacheID() < mRefineCache.size());
//   return mRefineCache[itr.cacheID()];
// }

//------------------------------------------------------------------------------
// Number of elements in the field.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
size_t
Field<Dimension, DataType>::numElements() const {
  return size();
}

template<typename Dimension, typename DataType>
inline
size_t 
Field<Dimension, DataType>::numInternalElements() const {
  return this->nodeList().numInternalNodes();
}

template<typename Dimension, typename DataType>
inline
size_t 
Field<Dimension, DataType>::numGhostElements() const {
  return this->nodeList().numGhostNodes();
}

//------------------------------------------------------------------------------
// Zero out the field elements.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
Field<Dimension, DataType>::Zero() {
  std::fill(mDataArray.begin(), mDataArray.end(), DataTypeTraits<DataType>::zero());
}

//------------------------------------------------------------------------------
// Apply a minimum value to the Field elements.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
Field<Dimension, DataType>::applyMin(const DataType& dataMin) {
  const size_t n = this->size();
  for (auto i = 0u; i < n; ++i) mDataArray[i] = std::max(mDataArray[i], dataMin);
}

//------------------------------------------------------------------------------
// Apply a maximum value to the Field elements.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
Field<Dimension, DataType>::applyMax(const DataType& dataMax) {
  const size_t n = this->size();
  for (auto i = 0u; i < n; ++i) mDataArray[i] = std::min(mDataArray[i], dataMax);
}

//------------------------------------------------------------------------------
// Apply a (double) minimum value  to the Field elements.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
Field<Dimension, DataType>::applyScalarMin(const double dataMin) {
  const size_t n = this->size();
  for (auto i = 0u; i < n; ++i) mDataArray[i] = std::max(mDataArray[i], dataMin);
}

//------------------------------------------------------------------------------
// Apply a (double) maximum value to the Field elements.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
Field<Dimension, DataType>::applyScalarMax(const double dataMax) {
  const size_t n = this->size();
  for (auto i = 0u; i < n; ++i) mDataArray[i] = std::min(mDataArray[i], dataMax);
}

//------------------------------------------------------------------------------
// Addition with another field.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>
Field<Dimension, DataType>::operator+(const Field<Dimension, DataType>& rhs) const {
  REQUIRE(this->nodeListPtr() == rhs.nodeListPtr());
  Field<Dimension, DataType> result(*this);
  const size_t n = this->size();
  for (size_t i = 0u; i < n; ++i) result(i) += rhs(i);
  return result;
}

//------------------------------------------------------------------------------
// Subtract another field.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>
Field<Dimension, DataType>::operator-(const Field<Dimension, DataType>& rhs) const {
  REQUIRE(this->nodeListPtr() == rhs.nodeListPtr());
  Field<Dimension, DataType> result(*this);
  const size_t n = this->size();
  for (auto i = 0u; i < n; ++i) result(i) -= rhs(i);
  return result;
}

//------------------------------------------------------------------------------
// Addition with another field in place
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>&
Field<Dimension, DataType>::operator+=(const Field<Dimension, DataType>& rhs) {
  REQUIRE(this->nodeListPtr() == rhs.nodeListPtr());
  const size_t n = this->size();
  for (auto i = 0u; i < n; ++i) (*this)(i) += rhs(i);
  return *this;
}

//------------------------------------------------------------------------------
// Subtract another field from this one in place.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>&
Field<Dimension, DataType>::operator-=(const Field<Dimension, DataType>& rhs) {
  REQUIRE(this->nodeListPtr() == rhs.nodeListPtr());
  const size_t n = this->size();
  for (auto i = 0u; i < n; ++i) (*this)(i) -= rhs(i);
  return *this;
}

//------------------------------------------------------------------------------
// Add a single value to every element of a field.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>
Field<Dimension, DataType>::operator+(const DataType& rhs) const {
  Field<Dimension, DataType> result(*this);
  const auto n = this->size();
  for (auto i = 0u; i < n; ++i) result(i) += rhs;
  return result;
}

//------------------------------------------------------------------------------
// Subtract a single value from every element of a field.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>
Field<Dimension, DataType>::operator-(const DataType& rhs) const {
  Field<Dimension, DataType> result(*this);
  const auto n = this->size();
  for (auto i = 0u; i < n; ++i) result(i) -= rhs;
  return result;
}

//------------------------------------------------------------------------------
// Addition with a single value in place
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>&
Field<Dimension, DataType>::operator+=(const DataType& rhs) {
  const auto n = this->numElements();
  for (auto i = 0u; i < n; ++i) (*this)(i) += rhs;
  return *this;
}

//------------------------------------------------------------------------------
// Subtract a single value in place
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>&
Field<Dimension, DataType>::operator-=(const DataType& rhs) {
  const auto n = this->numElements();
  for (auto i = 0u; i < n; ++i) (*this)(i) -= rhs;
  return *this;
}

//------------------------------------------------------------------------------
// Multiplication by a Scalar Field
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>
Field<Dimension, DataType>::
operator*(const Field<Dimension, Scalar>& rhs) const {
  REQUIRE(this->nodeListPtr() == rhs.nodeListPtr());
  Field<Dimension, DataType> result(*this);
  result *= rhs;
  return result;
}

//------------------------------------------------------------------------------
// Division by a Scalar Field
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>
Field<Dimension, DataType>::
operator/(const Field<Dimension, Scalar>& rhs) const {
  REQUIRE(this->nodeListPtr() == rhs.nodeListPtr());
  Field<Dimension, DataType> result(*this);
  result /= rhs;
  return result;
}

//------------------------------------------------------------------------------
// Multiplication by a Scalar Field in place.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>&
Field<Dimension, DataType>::
operator*=(const Field<Dimension, Scalar>& rhs) {
  REQUIRE(this->nodeListPtr() == rhs.nodeListPtr());
  const auto n = this->numElements();
  for (auto i = 0u; i < n; ++i) (*this)(i) *= rhs(i);
  return *this;
}

//------------------------------------------------------------------------------
// Division by a Scalar Field in place.
//-------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>&
Field<Dimension, DataType>::
operator/=(const Field<Dimension, typename Dimension::Scalar>& rhs) {
  REQUIRE(this->nodeListPtr() == rhs.nodeListPtr());
  const auto n = this->numElements();
  for (auto i = 0u; i < n; ++i) {
    (*this)(i) *= safeInvVar(rhs(i), 1.0e-60);
  }
  return *this;
}

//------------------------------------------------------------------------------
// Multiplication by a Scalar
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>
Field<Dimension, DataType>::
operator*(const Scalar& rhs) const {
  Field<Dimension, DataType> result(*this);
  result *= rhs;
  return result;
}

//------------------------------------------------------------------------------
// Division by a Scalar
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>
Field<Dimension, DataType>::
operator/(const Scalar& rhs) const {
  REQUIRE(rhs != 0.0);
  Field<Dimension, DataType> result(*this);
  result /= rhs;
  return result;
}

//------------------------------------------------------------------------------
// Multiplication by a Scalar in place
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>&
Field<Dimension, DataType>::
operator*=(const Scalar& rhs) {
  const auto n = this->numElements();
  for (auto i = 0u; i < n; ++i) (*this)(i) *= rhs;
  return *this;
}

//------------------------------------------------------------------------------
// Division by a Scalar value in place
//-------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>&
Field<Dimension, DataType>::
operator/=(const Scalar& rhs) {
  REQUIRE(rhs != 0.0);
  const auto n = this->numElements();
  for (auto i = 0u; i < n; ++i) (*this)(i) /= rhs;
  return *this;
}

//------------------------------------------------------------------------------
// Sum the elements of the field (assumes the DataType::operator+= is 
// available).
//-------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType
Field<Dimension, DataType>::
sumElements() const {
  return allReduce(this->localSumElements(), SPHERAL_OP_SUM);
}

//------------------------------------------------------------------------------
// Minimum.
//-------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType
Field<Dimension, DataType>::
min() const {
  return allReduce(this->localMin(), SPHERAL_OP_MIN);
}

//------------------------------------------------------------------------------
// Maximum.
//-------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType
Field<Dimension, DataType>::
max() const {
  return allReduce(this->localMax(), SPHERAL_OP_MAX);
}

//------------------------------------------------------------------------------
// Sum the elements of the field (assumes the DataType::operator+= is 
// available).  
// LOCAL to processor!
//-------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType
Field<Dimension, DataType>::
localSumElements() const {
  return std::accumulate(begin(), begin() + numInternalElements(), DataTypeTraits<DataType>::zero());
}

//------------------------------------------------------------------------------
// Minimum.
// LOCAL to processor!
//-------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType
Field<Dimension, DataType>::
localMin() const {
  DataType result;
  if (size() == 0u) {
    result = std::numeric_limits<DataType>::max();
  } else {
    result = *std::min_element(begin(), begin() + numInternalElements());
  }
  return result;
}

//------------------------------------------------------------------------------
// Maximum.
// LOCAL to processor!
//-------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType
Field<Dimension, DataType>::
localMax() const {
  DataType result;
  if (size() == 0u) {
    result = std::numeric_limits<DataType>::lowest();
  } else {
    result = *std::max_element(begin(), begin() + numInternalElements());
  }
  return result;
}

//------------------------------------------------------------------------------
// operator==(Field)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
Field<Dimension, DataType>::
operator==(const Field<Dimension, DataType>& rhs) const {
  return mDataArray == rhs.mDataArray;
}

//------------------------------------------------------------------------------
// operator!=(Field)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
Field<Dimension, DataType>::
operator!=(const Field<Dimension, DataType>& rhs) const {
  return !((*this) == rhs);
}

//------------------------------------------------------------------------------
// operator>(Field)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
Field<Dimension, DataType>::
operator>(const Field<Dimension, DataType>& rhs) const {
  if (this->size() != rhs.size()) return false;
  bool result = true;
  const_iterator lhsItr = this->begin();
  const_iterator rhsItr = rhs.begin();
  while (lhsItr < this->end() && result) {
    result = *lhsItr > *rhsItr;
    ++lhsItr;
    ++rhsItr;
  }
  return result;
}

//------------------------------------------------------------------------------
// operator<(Field)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
Field<Dimension, DataType>::
operator<(const Field<Dimension, DataType>& rhs) const {
  if (this->size() != rhs.size()) return false;
  bool result = true;
  const_iterator lhsItr = this->begin();
  const_iterator rhsItr = rhs.begin();
  while (lhsItr < this->end() && result) {
    result = *lhsItr < *rhsItr;
    ++lhsItr;
    ++rhsItr;
  }
  return result;
}

//------------------------------------------------------------------------------
// operator>=(Field)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
Field<Dimension, DataType>::
operator>=(const Field<Dimension, DataType>& rhs) const {
  if (this->size() != rhs.size()) return false;
  bool result = true;
  const_iterator lhsItr = this->begin();
  const_iterator rhsItr = rhs.begin();
  while (lhsItr < this->end() && result) {
    result = *lhsItr >= *rhsItr;
    ++lhsItr;
    ++rhsItr;
  }
  return result;
}

//------------------------------------------------------------------------------
// operator<=(Field)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
Field<Dimension, DataType>::
operator<=(const Field<Dimension, DataType>& rhs) const {
  if (this->size() != rhs.size()) return false;
  bool result = true;
  const_iterator lhsItr = this->begin();
  const_iterator rhsItr = rhs.begin();
  while (lhsItr < this->end() && result) {
    result = *lhsItr <= *rhsItr;
    ++lhsItr;
    ++rhsItr;
  }
  return result;
}

//------------------------------------------------------------------------------
// operator==(value)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
Field<Dimension, DataType>::
operator==(const DataType& rhs) const {
  bool result = true;
  const_iterator lhsItr = this->begin();
  while (lhsItr < this->end() && result) {
    result = *lhsItr == rhs;
    ++lhsItr;
  }
  return result;
}

//------------------------------------------------------------------------------
// operator!=(value)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
Field<Dimension, DataType>::
operator!=(const DataType& rhs) const {
  return !((*this) == rhs);
}

//------------------------------------------------------------------------------
// operator>(value)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
Field<Dimension, DataType>::
operator>(const DataType& rhs) const {
  bool result = true;
  const_iterator lhsItr = this->begin();
  while (lhsItr < this->end() && result) {
    result = *lhsItr > rhs;
    ++lhsItr;
  }
  return result;
}

//------------------------------------------------------------------------------
// operator<(value)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
Field<Dimension, DataType>::
operator<(const DataType& rhs) const {
  bool result = true;
  const_iterator lhsItr = this->begin();
  while (lhsItr < this->end() && result) {
    result = *lhsItr < rhs;
    ++lhsItr;
  }
  return result;
}

//------------------------------------------------------------------------------
// operator>=(value)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
Field<Dimension, DataType>::
operator>=(const DataType& rhs) const {
  bool result = true;
  const_iterator lhsItr = this->begin();
  while (lhsItr < this->end() && result) {
    result = *lhsItr >= rhs;
    ++lhsItr;
  }
  return result;
}

//------------------------------------------------------------------------------
// operator<=(value)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
Field<Dimension, DataType>::
operator<=(const DataType& rhs) const {
  bool result = true;
  const_iterator lhsItr = this->begin();
  while (lhsItr < this->end() && result) {
    result = *lhsItr <= rhs;
    ++lhsItr;
  }
  return result;
}

//------------------------------------------------------------------------------
// Iterator pointing to the beginning of the field values.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename Field<Dimension, DataType>::iterator
Field<Dimension, DataType>::begin() {
  return mDataArray.begin();
}

//------------------------------------------------------------------------------
// Iterator pointing to the end of the field values.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename Field<Dimension, DataType>::iterator
Field<Dimension, DataType>::end() {
  return mDataArray.end();
}

//------------------------------------------------------------------------------
// Iterator pointing to the beginning of the internal field values.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename Field<Dimension, DataType>::iterator
Field<Dimension, DataType>::internalBegin() {
  return mDataArray.begin();
}

//------------------------------------------------------------------------------
// Iterator pointing to the end of the internal field values.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename Field<Dimension, DataType>::iterator
Field<Dimension, DataType>::internalEnd() {
  CHECK(this->nodeList().firstGhostNode() <= this->nodeList().numNodes());
  return mDataArray.begin() + this->nodeList().firstGhostNode();
}

//------------------------------------------------------------------------------
// Iterator pointing to the beginning of the ghost field values.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename Field<Dimension, DataType>::iterator
Field<Dimension, DataType>::ghostBegin() {
  return this->internalEnd();
}

//------------------------------------------------------------------------------
// Iterator pointing to the end of the ghost field values.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename Field<Dimension, DataType>::iterator
Field<Dimension, DataType>::ghostEnd() {
  return mDataArray.end();
}

//------------------------------------------------------------------------------
// Const_iterator pointing to the beginning of the field values.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename Field<Dimension, DataType>::const_iterator
Field<Dimension, DataType>::begin() const {
  return mDataArray.begin();
}

//------------------------------------------------------------------------------
// Const_iterator pointing to the end of the field values.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename Field<Dimension, DataType>::const_iterator
Field<Dimension, DataType>::end() const {
  return mDataArray.end();
}

//------------------------------------------------------------------------------
// Const_iterator pointing to the beginning of the internal field values.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename Field<Dimension, DataType>::const_iterator
Field<Dimension, DataType>::internalBegin() const {
  return mDataArray.begin();
}

//------------------------------------------------------------------------------
// Const_iterator pointing to the end of the internal field values.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename Field<Dimension, DataType>::const_iterator
Field<Dimension, DataType>::internalEnd() const {
  REQUIRE(this->nodeList().firstGhostNode() <= this->nodeList().numNodes());
  return mDataArray.begin() + this->nodeList().firstGhostNode();
}

//------------------------------------------------------------------------------
// Const_iterator pointing to the beginning of the ghost field values.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename Field<Dimension, DataType>::const_iterator
Field<Dimension, DataType>::ghostBegin() const {
  return this->internalEnd();
}

//------------------------------------------------------------------------------
// Const_iterator pointing to the end of the ghost field values.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename Field<Dimension, DataType>::const_iterator
Field<Dimension, DataType>::ghostEnd() const {
  return mDataArray.end();
}

//------------------------------------------------------------------------------
// Set the NodeList with which this field is defined.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
Field<Dimension, DataType>::setNodeList(const NodeList<Dimension>& nodeList) {
  auto oldSize = this->size();
  this->setFieldBaseNodeList(nodeList);
  mDataArray.resize(nodeList.numNodes());
  if (this->size() > oldSize) {
    std::fill(mDataArray.begin() + oldSize, mDataArray.end(), DataTypeTraits<DataType>::zero());
  }
}

//------------------------------------------------------------------------------
// Pack the given buffer into the requested field elements.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
std::vector<char>
Field<Dimension, DataType>::
packValues(const std::vector<size_t>& nodeIDs) const {
  return packFieldValues(*this, nodeIDs);
}

//------------------------------------------------------------------------------
// Unpack the given buffer into the requested field elements.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
Field<Dimension, DataType>::
unpackValues(const std::vector<size_t>& nodeIDs,
             const std::vector<char>& buffer) {
  unpackFieldValues(*this, nodeIDs, buffer);
}

//------------------------------------------------------------------------------
// Copy values between sets of indices.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
Field<Dimension, DataType>::
copyElements(const std::vector<size_t>& fromIndices,
             const std::vector<size_t>& toIndices) {
  REQUIRE(fromIndices.size() == toIndices.size());
  REQUIRE(std::all_of(fromIndices.begin(), fromIndices.end(),
                      [&](const size_t i) { return i < this->size(); }));
  REQUIRE(std::all_of(toIndices.begin(), toIndices.end(),
                      [&](const size_t i) { return i < this->size(); }));
  const auto ni = fromIndices.size();
  for (auto k = 0u; k < ni; ++k) (*this)(toIndices[k]) = (*this)(fromIndices[k]);
}

//------------------------------------------------------------------------------
// fixedSizeDataType
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
Field<Dimension, DataType>::
fixedSizeDataType() const {
  return DataTypeTraits<DataType>::fixedSize();
}

//------------------------------------------------------------------------------
// numValsInDataType
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
size_t
Field<Dimension, DataType>::
numValsInDataType() const {
  return DataTypeTraits<DataType>::numElements(DataType());
}

//------------------------------------------------------------------------------
// sizeofDataType
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
size_t
Field<Dimension, DataType>::
sizeofDataType() const {
  return sizeof(DataTypeTraits<DataType>::zero());
}

//------------------------------------------------------------------------------
// computeCommBufferSize
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
size_t
Field<Dimension, DataType>::
computeCommBufferSize(const std::vector<size_t>& packIndices,
                      const int sendProc,
                      const int recvProc) const {
  return computeBufferSize(*this, packIndices, sendProc, recvProc);
}

//------------------------------------------------------------------------------
// Serialize the Field into a vector<char>.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
std::vector<char>
Field<Dimension, DataType>::
serialize() const {
  const size_t n = this->numInternalElements();
  vector<char> buf;
  packElement(this->name(), buf);
  packElement(n, buf);
  for (auto i = 0u; i < n; ++i) packElement((*this)[i], buf);
  return buf;
}

//------------------------------------------------------------------------------
// Deserialize the values from a vector<char> into this Field.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
Field<Dimension, DataType>::
deserialize(const std::vector<char>& buf) {
  auto itr = buf.begin();
  std::string nm;
  unpackElement(nm, itr, buf.end());
  this->name(nm);
  size_t n;
  unpackElement(n, itr, buf.end());
  VERIFY2(n == this->numInternalElements(),
          "Field ERROR: attempt to deserialize wrong number of elements: " << n << " != " << this->numInternalElements());
  for (auto i = 0u; i < n; ++i) unpackElement((*this)[i], itr, buf.end());
}

//------------------------------------------------------------------------------
// Construct std::vectors of the values.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
std::vector<DataType>
Field<Dimension, DataType>::
internalValues() const {
  const auto n = this->nodeList().numInternalNodes();
  std::vector<DataType> result(n);
  std::copy(internalBegin(), internalEnd(), result.begin());
  ENSURE(result.size() == n);
  return result;
}

template<typename Dimension, typename DataType>
inline
std::vector<DataType>
Field<Dimension, DataType>::
ghostValues() const {
  const auto n = this->nodeList().numGhostNodes();
  std::vector<DataType> result(n);
  std::copy(ghostBegin(), ghostEnd(), result.begin());
  ENSURE(result.size() == n);
  return result;
}

template<typename Dimension, typename DataType>
inline
std::vector<DataType>
Field<Dimension, DataType>::
allValues() const {
  const auto n = this->nodeList().numNodes();
  std::vector<DataType> result(n);
  std::copy(begin(), end(), result.begin());
  ENSURE(result.size() == n);
  return result;
}

//------------------------------------------------------------------------------
// Resize the field to the given number of nodes.  This operation ignores
// the distinction between internal and ghost nodes.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
Field<Dimension, DataType>::resizeField(size_t size) {
  REQUIRE(size == this->nodeList().numNodes());
  auto oldSize = this->size();
  mDataArray.resize(size);
  if (oldSize < size) {
    std::fill(mDataArray.begin() + oldSize,
              mDataArray.end(),
              DataTypeTraits<DataType>::zero());
  }
}

//------------------------------------------------------------------------------
// Resize the field to the given number of internal nodes, preserving any ghost
// values at the end.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
Field<Dimension, DataType>::resizeFieldInternal(const size_t newInternalSize,
                                                const size_t oldFirstGhostNode) {
  const auto currentSize = this->size();
  const auto currentInternalSize = oldFirstGhostNode;
  const auto numGhostNodes = this->nodeList().numGhostNodes();
  const auto newSize = newInternalSize + numGhostNodes;
  REQUIRE(numGhostNodes == currentSize - oldFirstGhostNode);
  REQUIRE(newSize == this->nodeList().numNodes());

  // If there is ghost data, we must preserve it.
  std::vector<DataType,DataAllocator<DataType>> oldGhostValues(numGhostNodes);
  if (numGhostNodes > 0u) {
    std::copy(mDataArray.begin() + oldFirstGhostNode, mDataArray.end(), oldGhostValues.begin());
  }

  // Resize the field data.
  mDataArray.resize(newSize);

  // Fill in any new internal values.
  if (newSize > currentSize) {
    CHECK(currentInternalSize < this->nodeList().firstGhostNode());
    std::fill(mDataArray.begin() + currentInternalSize,
              mDataArray.begin() + newInternalSize,
              DataTypeTraits<DataType>::zero());
  }

  // Fill the ghost data back in.
  if (numGhostNodes > 0u) {
    std::copy(oldGhostValues.begin(), oldGhostValues.end(), mDataArray.begin() + newInternalSize);
  }
}

//------------------------------------------------------------------------------
// Resize the field to the given number of ghost nodes.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
Field<Dimension, DataType>::resizeFieldGhost(const size_t size) {
  const auto currentSize = this->size();
  const auto numInternalNodes = this->nodeList().numInternalNodes();
  CHECK(currentSize >= numInternalNodes);
  const auto currentNumGhostNodes = currentSize - numInternalNodes;
  const auto newSize = numInternalNodes + size;
  REQUIRE(newSize == this->nodeList().numNodes());

  // Resize the field data.
  mDataArray.resize(newSize);
  CHECK(this->size() == newSize);

  // Fill in any new ghost values.
  if (newSize > currentSize) {
    std::fill(mDataArray.begin() + numInternalNodes + currentNumGhostNodes,
              mDataArray.end(),
              DataTypeTraits<DataType>::zero());
  }
}

//------------------------------------------------------------------------------
// Delete the given element id.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
Field<Dimension, DataType>::deleteElement(size_t nodeID) {
  const auto oldSize = this->size();
  CONTRACT_VAR(oldSize);
  REQUIRE(nodeID < oldSize);
  mDataArray.erase(mDataArray.begin() + nodeID);
  ENSURE(mDataArray.size() == oldSize - 1u);
}

//------------------------------------------------------------------------------
// Delete the given set of elements.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
Field<Dimension, DataType>::deleteElements(const std::vector<size_t>& nodeIDs) {
  // The standalone method does the actual work.
  removeElements(mDataArray, nodeIDs);
}

//****************************** Global Functions ******************************
//------------------------------------------------------------------------------
// Multiplication by another Field.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType, typename OtherDataType>
Field<Dimension, typename CombineTypes<DataType, OtherDataType>::ProductType>
operator*(const Field<Dimension, DataType>& lhs,
          const Field<Dimension, OtherDataType>& rhs) {
  CHECK(lhs.nodeList().numNodes() == rhs.nodeList().numNodes());
  Field<Dimension, typename CombineTypes<DataType, OtherDataType>::ProductType>
    result("product", const_cast<Field<Dimension, DataType>&>(lhs).nodeList());
  for (auto i = 0u; i < result.size(); ++i) {
    result(i) = lhs(i) * rhs(i);
  }
  return result;
}

//------------------------------------------------------------------------------
// Multiplication by a single value.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType, typename OtherDataType>
Field<Dimension, typename CombineTypes<DataType, OtherDataType>::ProductType>
operator*(const Field<Dimension, DataType>& lhs,
          const OtherDataType& rhs) {
  Field<Dimension, typename CombineTypes<DataType, OtherDataType>::ProductType>
    result("product", const_cast<Field<Dimension, DataType>&>(lhs).nodeList());
  for (auto i = 0u; i < result.size(); ++i) {
    result(i) = lhs(i) * rhs;
  }
  return result;
}

template<typename Dimension, typename DataType, typename OtherDataType>
Field<Dimension, typename CombineTypes<DataType, OtherDataType>::ProductType>
operator*(const DataType& lhs,
          const Field<Dimension, OtherDataType>& rhs) {
  Field<Dimension, typename CombineTypes<DataType, OtherDataType>::ProductType>
    result("product", const_cast<Field<Dimension, OtherDataType>&>(rhs).nodeList());
  for (auto i = 0u; i < result.size(); ++i) {
    result(i) = lhs * rhs(i);
  }
  return result;
}

// //------------------------------------------------------------------------------
// // Absolute value.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// Field<Dimension, typename Dimension::Scalar>
// abs(const Field<Dimension, typename Dimension::Scalar>& field) {
//   typedef typename Dimension::Scalar Scalar;
//   CHECK(field.valid());
//   Field<Dimension, Scalar> result(field);
//   for (int i = 0; i < result.size(); ++i) {
//     result(i) = std::abs(result(i));
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // Inverse cosine
// //------------------------------------------------------------------------------
// template<typename Dimension>
// Field<Dimension, typename Dimension::Scalar>
// acos(const Field<Dimension, typename Dimension::Scalar>& field) {
//   typedef typename Dimension::Scalar Scalar;
//   CHECK(field.valid());
//   Field<Dimension, Scalar> result(field);
//   for (int i = 0; i < result.size(); ++i) {
//     result(i) = std::acos(result(i));
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // Inverse sine
// //------------------------------------------------------------------------------
// template<typename Dimension>
// Field<Dimension, typename Dimension::Scalar>
// asin(const Field<Dimension, typename Dimension::Scalar>& field) {
//   typedef typename Dimension::Scalar Scalar;
//   CHECK(field.valid());
//   Field<Dimension, Scalar> result(field);
//   for (int i = 0; i < result.size(); ++i) {
//     result(i) = std::asin(result(i));
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // Inverse tangent
// //------------------------------------------------------------------------------
// template<typename Dimension>
// Field<Dimension, typename Dimension::Scalar>
// atan(const Field<Dimension, typename Dimension::Scalar>& field) {
//   typedef typename Dimension::Scalar Scalar;
//   CHECK(field.valid());
//   Field<Dimension, Scalar> result(field);
//   for (int i = 0; i < result.size(); ++i) {
//     result(i) = std::atan(result(i));
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // Inverse tangent2.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// Field<Dimension, typename Dimension::Scalar>
// atan2(const Field<Dimension, typename Dimension::Scalar>& field1,
//       const Field<Dimension, typename Dimension::Scalar>& field2) {
//   typedef typename Dimension::Scalar Scalar;
//   CHECK(field1.valid() && field2.valid());
//   CHECK(field1.nodeListPtr() == field2.nodeListPtr());
//   Field<Dimension, Scalar> result(field1);
//   for (int i = 0; i < result.size(); ++i) {
//     result(i) = std::atan2(field1(i), field2(i));
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // Ceiling -- smallest floating-point integer value not less than the argument.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// Field<Dimension, typename Dimension::Scalar>
// ceil(const Field<Dimension, typename Dimension::Scalar>& field) {
//   typedef typename Dimension::Scalar Scalar;
//   CHECK(field.valid());
//   Field<Dimension, Scalar> result(field);
//   for (int i = 0; i < result.size(); ++i) {
//     result(i) = std::ceil(result(i));
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // Cosine
// //------------------------------------------------------------------------------
// template<typename Dimension>
// Field<Dimension, typename Dimension::Scalar>
// cos(const Field<Dimension, typename Dimension::Scalar>& field) {
//   typedef typename Dimension::Scalar Scalar;
//   CHECK(field.valid());
//   Field<Dimension, Scalar> result(field);
//   for (int i = 0; i < result.size(); ++i) {
//     result(i) = std::cos(result(i));
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // Hyperbolic cosine
// //------------------------------------------------------------------------------
// template<typename Dimension>
// Field<Dimension, typename Dimension::Scalar>
// cosh(const Field<Dimension, typename Dimension::Scalar>& field) {
//   typedef typename Dimension::Scalar Scalar;
//   CHECK(field.valid());
//   Field<Dimension, Scalar> result(field);
//   for (int i = 0; i < result.size(); ++i) {
//     result(i) = std::cosh(result(i));
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // Exponential e^x
// //------------------------------------------------------------------------------
// template<typename Dimension>
// Field<Dimension, typename Dimension::Scalar>
// exp(const Field<Dimension, typename Dimension::Scalar>& field) {
//   typedef typename Dimension::Scalar Scalar;
//   CHECK(field.valid());
//   Field<Dimension, Scalar> result(field);
//   for (int i = 0; i < result.size(); ++i) {
//     result(i) = std::exp(result(i));
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // fabs -- same as abs, absolute value
// //------------------------------------------------------------------------------
// template<typename Dimension>
// Field<Dimension, typename Dimension::Scalar>
// fabs(const Field<Dimension, typename Dimension::Scalar>& field) {
//   typedef typename Dimension::Scalar Scalar;
//   CHECK(field.valid());
//   Field<Dimension, Scalar> result(field);
//   for (int i = 0; i < result.size(); ++i) {
//     result(i) = std::abs(result(i));
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // Floor -- largest floating-point integer value not greater than the argument
// //------------------------------------------------------------------------------
// template<typename Dimension>
// Field<Dimension, typename Dimension::Scalar>
// floor(const Field<Dimension, typename Dimension::Scalar>& field) {
//   typedef typename Dimension::Scalar Scalar;
//   CHECK(field.valid());
//   Field<Dimension, Scalar> result(field);
//   for (int i = 0; i < result.size(); ++i) {
//     result(i) = std::floor(result(i));
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // Natural logarithm
// //------------------------------------------------------------------------------
// template<typename Dimension>
// Field<Dimension, typename Dimension::Scalar>
// log(const Field<Dimension, typename Dimension::Scalar>& field) {
//   typedef typename Dimension::Scalar Scalar;
//   CHECK(field.valid());
//   Field<Dimension, Scalar> result(field);
//   for (int i = 0; i < result.size(); ++i) {
//     result(i) = std::log(result(i));
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // Log base 10
// //------------------------------------------------------------------------------
// template<typename Dimension>
// Field<Dimension, typename Dimension::Scalar>
// log10(const Field<Dimension, typename Dimension::Scalar>& field) {
//   typedef typename Dimension::Scalar Scalar;
//   CHECK(field.valid());
//   Field<Dimension, Scalar> result(field);
//   for (int i = 0; i < result.size(); ++i) {
//     result(i) = std::log10(result(i));
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // pow -- raise each element to an arbitrary power
// //------------------------------------------------------------------------------
// template<typename Dimension>
// Field<Dimension, typename Dimension::Scalar>
// pow(const Field<Dimension, typename Dimension::Scalar>& field,
//     const double exponent) {
//   typedef typename Dimension::Scalar Scalar;
//   CHECK(field.valid());
//   Field<Dimension, Scalar> result(field);
//   for (int i = 0; i < result.size(); ++i) {
//     result(i) = std::pow(result(i), exponent);
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // powN -- raise each element to the power N
// //------------------------------------------------------------------------------
// template<typename Dimension>
// Field<Dimension, typename Dimension::Scalar>
// pow2(const Field<Dimension, typename Dimension::Scalar>& field) {
//   typedef typename Dimension::Scalar Scalar;
//   CHECK(field.valid());
//   Field<Dimension, Scalar> result(field);
//   for (int i = 0; i < result.size(); ++i) {
//     result(i) = FastMath::pow2(result(i));
//   }
//   return result;
// }

// template<typename Dimension>
// Field<Dimension, typename Dimension::Scalar>
// pow3(const Field<Dimension, typename Dimension::Scalar>& field) {
//   typedef typename Dimension::Scalar Scalar;
//   CHECK(field.valid());
//   Field<Dimension, Scalar> result(field);
//   for (int i = 0; i < result.size(); ++i) {
//     result(i) = FastMath::pow3(result(i));
//   }
//   return result;
// }

// template<typename Dimension>
// Field<Dimension, typename Dimension::Scalar>
// pow4(const Field<Dimension, typename Dimension::Scalar>& field) {
//   typedef typename Dimension::Scalar Scalar;
//   CHECK(field.valid());
//   Field<Dimension, Scalar> result(field);
//   for (int i = 0; i < result.size(); ++i) {
//     result(i) = FastMath::pow4(result(i));
//   }
//   return result;
// }

// template<typename Dimension>
// Field<Dimension, typename Dimension::Scalar>
// pow5(const Field<Dimension, typename Dimension::Scalar>& field) {
//   typedef typename Dimension::Scalar Scalar;
//   CHECK(field.valid());
//   Field<Dimension, Scalar> result(field);
//   for (int i = 0; i < result.size(); ++i) {
//     result(i) = FastMath::pow5(result(i));
//   }
//   return result;
// }

// template<typename Dimension>
// Field<Dimension, typename Dimension::Scalar>
// pow6(const Field<Dimension, typename Dimension::Scalar>& field) {
//   typedef typename Dimension::Scalar Scalar;
//   CHECK(field.valid());
//   Field<Dimension, Scalar> result(field);
//   for (int i = 0; i < result.size(); ++i) {
//     result(i) = FastMath::pow6(result(i));
//   }
//   return result;
// }

// template<typename Dimension>
// Field<Dimension, typename Dimension::Scalar>
// pow7(const Field<Dimension, typename Dimension::Scalar>& field) {
//   typedef typename Dimension::Scalar Scalar;
//   CHECK(field.valid());
//   Field<Dimension, Scalar> result(field);
//   for (int i = 0; i < result.size(); ++i) {
//     result(i) = FastMath::pow7(result(i));
//   }
//   return result;
// }

// template<typename Dimension>
// Field<Dimension, typename Dimension::Scalar>
// pow8(const Field<Dimension, typename Dimension::Scalar>& field) {
//   typedef typename Dimension::Scalar Scalar;
//   Field<Dimension, Scalar> result(field);
//   for (int i = 0; i < result.size(); ++i) {
//     result(i) = FastMath::pow8(result(i));
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // Sine
// //------------------------------------------------------------------------------
// template<typename Dimension>
// Field<Dimension, typename Dimension::Scalar>
// sin(const Field<Dimension, typename Dimension::Scalar>& field) {
//   typedef typename Dimension::Scalar Scalar;
//   CHECK(field.valid());
//   Field<Dimension, Scalar> result(field);
//   for (int i = 0; i < result.size(); ++i) {
//     result(i) = std::sin(result(i));
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // Hyperbolic sine
// //------------------------------------------------------------------------------
// template<typename Dimension>
// Field<Dimension, typename Dimension::Scalar>
// sinh(const Field<Dimension, typename Dimension::Scalar>& field) {
//   typedef typename Dimension::Scalar Scalar;
//   CHECK(field.valid());
//   Field<Dimension, Scalar> result(field);
//   for (int i = 0; i < result.size(); ++i) {
//     result(i) = std::sinh(result(i));
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // Sqr -- square, same as pow2()
// //------------------------------------------------------------------------------
// template<typename Dimension>
// Field<Dimension, typename Dimension::Scalar>
// sqr(const Field<Dimension, typename Dimension::Scalar>& field) {
//   return field.pow2();
// }

// //------------------------------------------------------------------------------
// // Sqrt
// //------------------------------------------------------------------------------
// template<typename Dimension>
// Field<Dimension, typename Dimension::Scalar>
// sqrt(const Field<Dimension, typename Dimension::Scalar>& field) {
//   typedef typename Dimension::Scalar Scalar;
//   CHECK(field.valid());
//   Field<Dimension, Scalar> result(field);
//   for (int i = 0; i < result.size(); ++i) {
//     CHECK(result(i) >= 0.0);
//     result(i) = std::sqrt(result(i));
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // Minimum
// //------------------------------------------------------------------------------
// template<typename Dimension, typename DataType>
// Field<Dimension, DataType>
// min(const Field<Dimension, DataType>& field1,
//     const Field<Dimension, DataType>& field2) {
//   CHECK(field1.valid() && field2.valid());
//   CHECK(field1.size() == field2.size());
//   CHECK(field1.nodeListPtr() == field2.nodeListPtr());
//   Field<Dimension, DataType> result("min", const_cast<NodeList<Dimension>&>(field1.nodeList()));
//   for (int i = 0; i < field1.size(); ++i) {
//     result(i) = std::min(field1(i), field2(i));
//   }
//   return result;
// }

// template<typename Dimension, typename DataType>
// Field<Dimension, DataType>
// min(const DataType& value,
//     const Field<Dimension, DataType>& field) {
//   CHECK(field.valid());
//   Field<Dimension, DataType> result("min", const_cast<NodeList<Dimension>&>(field.nodeList()));
//   for (int i = 0; i < field.size(); ++i) {
//     result(i) = std::min(value, field(i));
//   }
//   return result;
// }

// template<typename Dimension, typename DataType>
// Field<Dimension, DataType>
// min(const Field<Dimension, DataType>& field, 
//     const DataType& value) {
//   CHECK(field.valid());
//   Field<Dimension, DataType> result("min", const_cast<NodeList<Dimension>&>(field.nodeList()));
//   for (int i = 0; i < field.size(); ++i) {
//     result(i) = std::min(field(i), value);
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // Maximum
// //------------------------------------------------------------------------------
// template<typename Dimension, typename DataType>
// Field<Dimension, DataType>
// max(const Field<Dimension, DataType>& field1,
//     const Field<Dimension, DataType>& field2) {
//   CHECK(field1.valid() && field2.valid());
//   CHECK(field1.size() == field2.size());
//   CHECK(field1.nodeListPtr() == field2.nodeListPtr());
//   Field<Dimension, DataType> result("max", const_cast<NodeList<Dimension>&>(field1.nodeList()));
//   for (int i = 0; i < field1.size(); ++i) {
//     result(i) = std::max(field1(i), field2(i));
//   }
//   return result;
// }

// template<typename Dimension, typename DataType>
// Field<Dimension, DataType>
// max(const DataType& value,
//     const Field<Dimension, DataType>& field) {
//   CHECK(field.valid());
//   Field<Dimension, DataType> result("max", const_cast<NodeList<Dimension>&>(field.nodeList()));
//   for (int i = 0; i < field.size(); ++i) {
//     result(i) = std::max(value, field(i));
//   }
//   return result;
// }

// template<typename Dimension, typename DataType>
// Field<Dimension, DataType>
// max(const Field<Dimension, DataType>& field, 
//     const DataType& value) {
//   CHECK(field.valid());
//   Field<Dimension, DataType> result("max", const_cast<NodeList<Dimension>&>(field.nodeList()));
//   for (int i = 0; i < field.size(); ++i) {
//     result(i) = std::max(field(i), value);
//   }
//   return result;
// }

//------------------------------------------------------------------------------
// Input (istream) operator.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
std::istream&
operator>>(std::istream& is, Field<Dimension, DataType>& field) {

  // Start by reading the number of elements.
  int numElementsInStream;
  is >> numElementsInStream;
  CHECK(numElementsInStream == (int)field.nodeList().numInternalNodes());

  // Read in the elements.
  for (auto itr = field.internalBegin();
       itr < field.internalEnd();
       ++itr) {
    is >> *itr;
  }
  return is;
}

//------------------------------------------------------------------------------
// Output (ostream) operator.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
std::ostream&
operator<<(std::ostream& os, const Field<Dimension, DataType>& field) {

  // Write the number of internal elements.
  os << field.nodeList().numInternalNodes() << " ";

  // Write the internal elements.
  for (auto itr = field.internalBegin();
       itr < field.internalEnd();
       ++itr) {
    os << *itr << " ";
  }
//   os << endl;
  return os;
}

//------------------------------------------------------------------------------
// getAxomType
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
axom::sidre::DataTypeId
Field<Dimension, DataType>::
getAxomTypeID() const {
  return DataTypeTraits<DataType>::axomTypeID();
}

// Default toView function call without an additionally defined callback.
template<typename Dimension, typename DataType>
inline
typename Field<Dimension, DataType>::ViewType
Field<Dimension, DataType>::
toView()
{
  auto func = [](
      const chai::PointerRecord *,
      chai::Action,
      chai::ExecutionSpace) {};

  return this->toView(func);
}

// The Primary toView() implementation. DataType MUST be implicitly copyable
// to call toView on a Field. Field::toView() passes the location of the
// std::vector allocation to a chai::ManagedArray. the MA does NOT own the
// host data. Field is still resposible for deallocation on destruction. On
// subsequent calls of toView() a check is made to see if the std::vector has
// been resized OR reallocated. If so the current MA object calls free -
// releasing any possible device memory that has been allocated. A new MA is
// then created with the std::vector pointer.
template<typename Dimension, typename DataType>
template<typename T, typename F>
std::enable_if_t<std::is_trivially_copyable<T>::value, typename Field<Dimension, DataType>::ViewType>
Field<Dimension, DataType>::
toView(F&& extension)
{
  if (mManagedData.size() != mDataArray.size() ||
      mManagedData.data(chai::CPU, false) != mDataArray.data()) {

    mManagedData.free();

    mManagedData = chai::makeManagedArray(
        mDataArray.data(), mDataArray.size(), chai::CPU, false);

    mManagedData.setUserCallback(
      getFieldCallback(std::forward<F>(extension)));
  }
  return ViewType(mManagedData);
}

// The inverse SFINAE of the above implementation. This should throw an error
// if it is ever called with a type that is not implicitly copyable.
template<typename Dimension, typename DataType>
template<typename T, typename F>
std::enable_if_t<!std::is_trivially_copyable<T>::value, typename Field<Dimension, DataType>::ViewType>
Field<Dimension, DataType>::
toView(F&&)
{
  ASSERT2(false, "Spheral::Field::toView() Is invalid when Field::DataType is not trivially copyable.");
  return ViewType(mManagedData);
}

// Default callback action to be used with chai Managed containers. An
// additional calback can be passed to extend this functionality. Useful for
// debuggin, testing and probing for performance counters / metrics.
template<typename Dimension, typename DataType>
template<typename F>
auto
Field<Dimension, DataType>::
getFieldCallback(F callback)
{
  return [n = this->name(), callback](
    const chai::PointerRecord * record,
    chai::Action action,
    chai::ExecutionSpace space) {
      if (action == chai::ACTION_MOVE) {
        if (space == chai::CPU)
          DEBUG_LOG << "Field :" << n << ": MOVED to the CPU";
        if (space == chai::GPU)
          DEBUG_LOG << "Field :" << n << ": MOVED to the GPU";
      }
      else if (action == chai::ACTION_ALLOC) {
        if (space == chai::CPU)
          DEBUG_LOG << "Field :" << n << ": ALLOC on the CPU";
        if (space == chai::GPU)
          DEBUG_LOG << "Field :" << n << ": ALLOC on the GPU";
      }
      else if (action == chai::ACTION_FREE) {
        if (space == chai::CPU)
          DEBUG_LOG << "Field :" << n << ": FREE on the CPU";
        if (space == chai::GPU)
          DEBUG_LOG << "Field :" << n << ": FREE on the GPU";
      }
      callback(record, action, space);
    };
}

} // namespace Spheral
