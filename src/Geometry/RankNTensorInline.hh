#include <algorithm>
#include <limits.h>
#include <string>

#include "Utilities/SpheralFunctions.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

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
  for (size_type i = 0u; i < numElements(); ++i) mElements[i] += rhs.mElements[i];
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
  for (size_type i = 0u; i < numElements(); ++i) mElements[i] -= rhs.mElements[i];
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
  for (size_type i = 0u; i < numElements(); ++i) mElements[i] *= rhs;
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
  for (size_type i = 0u; i < numElements(); ++i) mElements[i] /= rhs;
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
operator==(const RankNTensor<nDim, rank, Descendant>& rhs) const {
  bool result = mElements[0u] == rhs.mElements[0u];
  size_type i = 1u;
  while (i != numElements() and result) {
    result = result and (mElements[i] == rhs.mElements[i]);
    ++i;
  }
  return result;
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
  bool result = mElements[0u] == rhs;
  size_type i = 1u;
  while (i != numElements() and result) {
    result = result and (mElements[i] == rhs);
    ++i;
  }
  return result;
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
  double result = 0.0;
  for (auto i = 0u; i < numElements(); ++i) result += mElements[i]*rhs.mElements[i];
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
  for (size_type i = 1u; i < numElements(); ++i) *(result.begin() + i) *= mElements[i];
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
  double result = std::abs(mElements[0u]);
  for (size_type i = 1u; i < numElements(); ++i) result = std::max(result, std::abs(mElements[i]));
  return result;
}

//******************************************************************************
// Global methods.  (unbound to class)
//******************************************************************************

//------------------------------------------------------------------------------
// Multiplication with a double.
//------------------------------------------------------------------------------
template<int nDim, int rank, typename Descendant>
SPHERAL_HOST_DEVICE
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
