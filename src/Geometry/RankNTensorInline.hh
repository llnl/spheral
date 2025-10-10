#include <algorithm>
#include <limits.h>
#include <string>

#include "Utilities/SpheralFunctions.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Construct with the given value filling the tensor.
//------------------------------------------------------------------------------
template<int nDim, int rank, typename Descendant>
SPHERAL_HOST_DEVICE
inline
RankNTensor<nDim,rank, Descendant>::
RankNTensor(const double val) {
  for (size_t i = 0; i < numElements; ++i) mElements[i] = val;
}

//------------------------------------------------------------------------------
// Assignment (double).
//------------------------------------------------------------------------------
template<int nDim, int rank, typename Descendant>
SPHERAL_HOST_DEVICE
inline
RankNTensor<nDim,rank, Descendant>&
RankNTensor<nDim,rank, Descendant>::
operator=(const double rhs) {
  for (size_t i = 0; i < numElements; ++i) mElements[i] = rhs;
  return *this;
}

//------------------------------------------------------------------------------
// Return the (index) element using the bracket operator.
//------------------------------------------------------------------------------
template<int nDim, int rank, typename Descendant>
SPHERAL_HOST_DEVICE
inline
double
RankNTensor<nDim, rank, Descendant>::operator[](typename RankNTensor<nDim, rank, Descendant>::size_type index) const {
  REQUIRE(index < numElements);
  return *(begin() + index);
}

template<int nDim, int rank, typename Descendant>
SPHERAL_HOST_DEVICE
inline
double&
RankNTensor<nDim, rank, Descendant>::operator[](typename RankNTensor<nDim, rank, Descendant>::size_type index) {
  REQUIRE(index < numElements);
  return *(begin() + index);
}

//------------------------------------------------------------------------------
// Iterators.
//------------------------------------------------------------------------------
template<int nDim, int rank, typename Descendant>
SPHERAL_HOST_DEVICE
inline
typename RankNTensor<nDim,rank, Descendant>::iterator
RankNTensor<nDim,rank, Descendant>::
begin() {
  return mElements;
}

template<int nDim, int rank, typename Descendant>
SPHERAL_HOST_DEVICE
inline
typename RankNTensor<nDim,rank, Descendant>::iterator
RankNTensor<nDim,rank, Descendant>::
end() {
  return mElements + numElements;
}

template<int nDim, int rank, typename Descendant>
SPHERAL_HOST_DEVICE
inline
typename RankNTensor<nDim,rank, Descendant>::const_iterator
RankNTensor<nDim,rank, Descendant>::
begin() const {
  return mElements;
}

template<int nDim, int rank, typename Descendant>
SPHERAL_HOST_DEVICE
inline
typename RankNTensor<nDim,rank, Descendant>::const_iterator
RankNTensor<nDim,rank, Descendant>::
end() const {
  return mElements + numElements;
}

//------------------------------------------------------------------------------
// Zero out the tensor.
//------------------------------------------------------------------------------
template<int nDim, int rank, typename Descendant>
SPHERAL_HOST_DEVICE
inline
void
RankNTensor<nDim,rank, Descendant>::
Zero() {
  for (size_t i = 0; i < numElements; ++i) mElements[i] = 0.0;
}

//------------------------------------------------------------------------------
// Negative.
//------------------------------------------------------------------------------
template<int nDim, int rank, typename Descendant>
SPHERAL_HOST_DEVICE
inline
Descendant
RankNTensor<nDim, rank, Descendant>::
operator-() const {
  Descendant result(static_cast<const Descendant&>(*this));
  result *= -1.0;
  return result;
}

//------------------------------------------------------------------------------
// In place addition.
//------------------------------------------------------------------------------
template<int nDim, int rank, typename Descendant>
SPHERAL_HOST_DEVICE
inline
Descendant&
RankNTensor<nDim,rank, Descendant>::
operator+=(const RankNTensor& rhs) {
  for (size_type i = 0; i != numElements; ++i) mElements[i] += rhs.mElements[i];
  return static_cast<Descendant&>(*this);
}

//------------------------------------------------------------------------------
// In place subtraction.
//------------------------------------------------------------------------------
template<int nDim, int rank, typename Descendant>
SPHERAL_HOST_DEVICE
inline
Descendant&
RankNTensor<nDim,rank, Descendant>::
operator-=(const RankNTensor& rhs) {
  for (size_type i = 0; i != numElements; ++i) mElements[i] -= rhs.mElements[i];
  return static_cast<Descendant&>(*this);
}

//------------------------------------------------------------------------------
// Addition.
//------------------------------------------------------------------------------
template<int nDim, int rank, typename Descendant>
SPHERAL_HOST_DEVICE
inline
Descendant
RankNTensor<nDim,rank, Descendant>::
operator+(const RankNTensor& rhs) const {
  Descendant result(static_cast<const Descendant&>(*this));
  result += rhs;
  return result;
}

//------------------------------------------------------------------------------
// Subtraction.
//------------------------------------------------------------------------------
template<int nDim, int rank, typename Descendant>
SPHERAL_HOST_DEVICE
inline
Descendant
RankNTensor<nDim,rank, Descendant>::
operator-(const RankNTensor& rhs) const {
  Descendant result(static_cast<const Descendant&>(*this));
  result -= rhs;
  return result;
}

//------------------------------------------------------------------------------
// In place multiplication (double).
//------------------------------------------------------------------------------
template<int nDim, int rank, typename Descendant>
SPHERAL_HOST_DEVICE
inline
Descendant&
RankNTensor<nDim,rank, Descendant>::
operator*=(const double rhs) {
  for (size_type i = 0; i != numElements; ++i) mElements[i] *= rhs;
  return static_cast<Descendant&>(*this);
}

//------------------------------------------------------------------------------
// In place division (double).
//------------------------------------------------------------------------------
template<int nDim, int rank, typename Descendant>
SPHERAL_HOST_DEVICE
inline
Descendant&
RankNTensor<nDim,rank, Descendant>::
operator/=(const double rhs) {
  REQUIRE(rhs != 0.0);
  for (size_type i = 0; i != numElements; ++i) mElements[i] /= rhs;
  return static_cast<Descendant&>(*this);
}

//------------------------------------------------------------------------------
// Multiplication (double).
//------------------------------------------------------------------------------
template<int nDim, int rank, typename Descendant>
SPHERAL_HOST_DEVICE
inline
Descendant
RankNTensor<nDim,rank, Descendant>::
operator*(const double rhs) const {
  Descendant result(static_cast<const Descendant&>(*this));
  result *= rhs;
  return result;
}

//------------------------------------------------------------------------------
// Division (double).
//------------------------------------------------------------------------------
template<int nDim, int rank, typename Descendant>
SPHERAL_HOST_DEVICE
inline
Descendant
RankNTensor<nDim,rank, Descendant>::
operator/(const double rhs) const {
  REQUIRE(rhs != 0.0);
  Descendant result(static_cast<const Descendant&>(*this));
  result /= rhs;
  return result;
}

//------------------------------------------------------------------------------
// ==
//------------------------------------------------------------------------------
template<int nDim, int rank, typename Descendant>
SPHERAL_HOST_DEVICE
inline
bool
RankNTensor<nDim,rank, Descendant>::
operator==(const RankNTensor& rhs) const {
  bool result = mElements[0] == rhs.mElements[0];
  size_type i = 1;
  while (i != numElements and result) {
    result = result and (mElements[i] == rhs.mElements[i]);
    ++i;
  }
  return result;
}

//------------------------------------------------------------------------------
// !=
//------------------------------------------------------------------------------
template<int nDim, int rank, typename Descendant>
SPHERAL_HOST_DEVICE
inline
bool
RankNTensor<nDim,rank, Descendant>::
operator!=(const RankNTensor& rhs) const {
  return not this->operator==(rhs);
}

//------------------------------------------------------------------------------
// <
//------------------------------------------------------------------------------
template<int nDim, int rank, typename Descendant>
SPHERAL_HOST_DEVICE
inline
bool
RankNTensor<nDim,rank, Descendant>::
operator<(const RankNTensor& rhs) const {
  return this->doubledot(*this) < rhs.doubledot(rhs);
}

//------------------------------------------------------------------------------
// >
//------------------------------------------------------------------------------
template<int nDim, int rank, typename Descendant>
SPHERAL_HOST_DEVICE
inline
bool
RankNTensor<nDim,rank, Descendant>::
operator>(const RankNTensor& rhs) const {
  return this->doubledot(*this) > rhs.doubledot(rhs);
}

//------------------------------------------------------------------------------
// <=
//------------------------------------------------------------------------------
template<int nDim, int rank, typename Descendant>
SPHERAL_HOST_DEVICE
inline
bool
RankNTensor<nDim,rank, Descendant>::
operator<=(const RankNTensor& rhs) const {
  return this->operator<(rhs) or this->operator==(rhs);
}

//------------------------------------------------------------------------------
// >=
//------------------------------------------------------------------------------
template<int nDim, int rank, typename Descendant>
SPHERAL_HOST_DEVICE
inline
bool
RankNTensor<nDim,rank, Descendant>::
operator>=(const RankNTensor& rhs) const {
  return this->operator>(rhs) or this->operator==(rhs);
}

//------------------------------------------------------------------------------
// == (double)
//------------------------------------------------------------------------------
template<int nDim, int rank, typename Descendant>
SPHERAL_HOST_DEVICE
inline
bool
RankNTensor<nDim,rank, Descendant>::
operator==(const double rhs) const {
  bool result = mElements[0] == rhs;
  size_type i = 1;
  while (i != numElements and result) {
    result = result and (mElements[i] == rhs);
    ++i;
  }
  return result;
}

//------------------------------------------------------------------------------
// != (double)
//------------------------------------------------------------------------------
template<int nDim, int rank, typename Descendant>
SPHERAL_HOST_DEVICE
inline
bool
RankNTensor<nDim,rank, Descendant>::
operator!=(const double rhs) const {
  return not this->operator==(rhs);
}

//------------------------------------------------------------------------------
// doubledot
//------------------------------------------------------------------------------
template<int nDim, int rank, typename Descendant>
SPHERAL_HOST_DEVICE
inline
double
RankNTensor<nDim,rank, Descendant>::
doubledot(const RankNTensor& rhs) const {
  double result = mElements[0] * rhs.mElements[0];
  for (size_type i = 1; i != numElements; ++i) 
    result += mElements[i] * rhs.mElements[i];
  return result;
}

//------------------------------------------------------------------------------
// squareElements
//------------------------------------------------------------------------------
template<int nDim, int rank, typename Descendant>
SPHERAL_HOST_DEVICE
inline
Descendant
RankNTensor<nDim,rank, Descendant>::
squareElements() const {
  Descendant result(static_cast<const Descendant&>(*this));
  for (size_type i = 1; i != numElements; ++i) 
    *(result.begin() + i) *= mElements[i];
  return result;
}

//------------------------------------------------------------------------------
// maxAbsElement
//------------------------------------------------------------------------------
template<int nDim, int rank, typename Descendant>
SPHERAL_HOST_DEVICE
inline
double
RankNTensor<nDim,rank, Descendant>::
maxAbsElement() const {
  double result = std::abs(mElements[0]);
  for (size_type i = 1; i < numElements; ++i) result = std::max(result, std::abs(mElements[i]));
  return result;
}

//******************************************************************************
// Global methods.  (unbound to class)
//******************************************************************************

//------------------------------------------------------------------------------
// Multiplication with a double.
//------------------------------------------------------------------------------
template<int nDim, int rank, typename Descendant>
inline
Descendant
operator*(const double lhs,
          const RankNTensor<nDim, rank, Descendant>& rhs) {
  Descendant result(static_cast<const Descendant&>(rhs));
  result *= lhs;
  return result;
}

//------------------------------------------------------------------------------
// operator>> (input)
//------------------------------------------------------------------------------
template<int nDim, int rank, typename Descendant>
inline
std::istream&
operator>>(std::istream& is, RankNTensor<nDim, rank, Descendant>& ten) {
  std::string parenthesis;
  is >> parenthesis;
  for (typename RankNTensor<nDim,rank, Descendant>::iterator elementItr = ten.begin();
       elementItr != ten.end();
       ++elementItr) {
    is >> *elementItr;
  }
  is >> parenthesis;
  return is;
}

//------------------------------------------------------------------------------
// operator<< (output)
//------------------------------------------------------------------------------
template<int nDim, int rank, typename Descendant>
inline
std::ostream&
operator<<(std::ostream& os, const RankNTensor<nDim, rank,Descendant>& ten) {
  os << "( ";
  for (typename RankNTensor<nDim,rank,Descendant>::const_iterator itr = ten.begin();
       itr != ten.end(); ++itr) {
    os << *itr << " ";
  }
  os << ")";
  return os;
}

}
