#include <algorithm>
#include <limits.h>
#include <string>

#include "GeomThirdRankTensor.hh"
#include "Utilities/SpheralFunctions.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Construct with the given value filling the tensor.
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
GeomThirdRankTensor<nDim>::
GeomThirdRankTensor(const double val):
  RankNTensor<nDim, 3, GeomThirdRankTensor>(val) {
}

//------------------------------------------------------------------------------
// Assignment (double).
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
GeomThirdRankTensor<nDim>&
GeomThirdRankTensor<nDim>::
operator=(const double rhs) {
  RankNTensor<nDim, 3, GeomThirdRankTensor>::operator=(rhs);
  return *this;
}

//------------------------------------------------------------------------------
// Access the elements by indicies.
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
double
GeomThirdRankTensor<nDim>::
operator()(const GeomThirdRankTensor::size_type i,
           const GeomThirdRankTensor::size_type j,
           const GeomThirdRankTensor::size_type k) const {
  REQUIRE(i < nDim and j < nDim and k < nDim);
  return mElements[(i*nDim + j)*nDim + k];
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
double&
GeomThirdRankTensor<nDim>::
operator()(const GeomThirdRankTensor::size_type i,
           const GeomThirdRankTensor::size_type j,
           const GeomThirdRankTensor::size_type k) {
  REQUIRE(i < nDim and j < nDim and k < nDim);
  return mElements[(i*nDim + j)*nDim + k];
}

//------------------------------------------------------------------------------
// Have the global functions use the generic RankNTensor methods.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomThirdRankTensor<nDim>
operator*(const double lhs, const GeomThirdRankTensor<nDim>& rhs) {
  return lhs * static_cast<const RankNTensor<nDim, 3, GeomThirdRankTensor<nDim> >&>(rhs);
}

template<int nDim>
inline
::std::istream&
operator>>(std::istream& is, GeomThirdRankTensor<nDim>& rhs) {
  return operator>>(is, static_cast<RankNTensor<nDim, 3, GeomThirdRankTensor<nDim> >&>(rhs));
}

template<int nDim>
inline
::std::ostream&
operator<<(std::ostream& os, const GeomThirdRankTensor<nDim>& rhs) {
  return operator<<(os, static_cast<const RankNTensor<nDim, 3, GeomThirdRankTensor<nDim> >&>(rhs));
}

}
