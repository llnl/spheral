//---------------------------------Spheral++----------------------------------//
// A collection of useful helper methods to explicitly unroll Dimensional 
// loops for efficiency.
//----------------------------------------------------------------------------//
#ifndef __Spheral_GeometricUtilities__
#define __Spheral_GeometricUtilities__

#include "Geometry/Dimension.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Construct a Tensor with the given diagonal elements
//------------------------------------------------------------------------------
SPHERAL_HOST_DEVICE inline
Dim<1>::Tensor
constructTensorWithDiagonal(const Dim<1>::Vector& vec) {
  return Dim<1>::Tensor(vec[0]);
}

SPHERAL_HOST_DEVICE inline
Dim<2>::Tensor
constructTensorWithDiagonal(const Dim<2>::Vector& vec) {
  return Dim<2>::Tensor(vec[0], 0.0,
                        0.0,    vec[1]);
}

SPHERAL_HOST_DEVICE inline
Dim<3>::Tensor
constructTensorWithDiagonal(const Dim<3>::Vector& vec) {
  return Dim<3>::Tensor(vec[0], 0.0, 0.0,
                        0.0,    vec[1], 0.0,
                        0.0,    0.0,    vec[2]);
}

SPHERAL_HOST_DEVICE inline
Dim<1>::SymTensor
constructSymTensorWithDiagonal(const Dim<1>::Vector& vec) {
  return Dim<1>::SymTensor(vec[0]);
}

SPHERAL_HOST_DEVICE inline
Dim<2>::SymTensor
constructSymTensorWithDiagonal(const Dim<2>::Vector& vec) {
  return Dim<2>::SymTensor(vec[0], 0.0,
                           0.0,     vec[1]);
}

SPHERAL_HOST_DEVICE inline
Dim<3>::SymTensor
constructSymTensorWithDiagonal(const Dim<3>::Vector& vec) {
  return Dim<3>::SymTensor(vec[0], 0.0,    0.0,
                           0.0,    vec[1], 0.0,
                           0.0,    0.0,    vec[2]);
}

//------------------------------------------------------------------------------
// Construct a Tensor with the given diagonal elements, enforcing a minimum.
//------------------------------------------------------------------------------
SPHERAL_HOST_DEVICE inline
Dim<1>::Tensor
constructTensorWithMinDiagonal(const Dim<1>::Vector& vec,
                               const double maxval) {
  return Dim<1>::Tensor(std::min(maxval, vec[0]));
}

SPHERAL_HOST_DEVICE inline
Dim<2>::Tensor
constructTensorWithMinDiagonal(const Dim<2>::Vector& vec,
                               const double maxval) {
  return Dim<2>::Tensor(std::min(maxval, vec[0]), 0.0,
                        0.0, std::min(maxval, vec[1]));
}

SPHERAL_HOST_DEVICE inline
Dim<3>::Tensor
constructTensorWithMinDiagonal(const Dim<3>::Vector& vec,
                               const double maxval) {
  return Dim<3>::Tensor(std::min(maxval, vec[0]), 0.0, 0.0,
                        0.0, std::min(maxval, vec[1]), 0.0,
                        0.0, 0.0, std::min(maxval, vec[2]));
}

SPHERAL_HOST_DEVICE inline
Dim<1>::SymTensor
constructSymTensorWithMinDiagonal(const Dim<1>::Vector& vec,
                                  const double maxval) {
  return Dim<1>::SymTensor(std::min(maxval, vec[0]));
}

SPHERAL_HOST_DEVICE inline
Dim<2>::SymTensor
constructSymTensorWithMinDiagonal(const Dim<2>::Vector& vec,
                                  const double maxval) {
  return Dim<2>::SymTensor(std::min(maxval, vec[0]), 0.0,
                           0.0, std::min(maxval, vec[1]));
}

SPHERAL_HOST_DEVICE inline
Dim<3>::SymTensor
constructSymTensorWithMinDiagonal(const Dim<3>::Vector& vec,
                                  const double maxval) {
  return Dim<3>::SymTensor(std::min(maxval, vec[0]), 0.0, 0.0,
                           0.0, std::min(maxval, vec[1]), 0.0,
                           0.0, 0.0, std::min(maxval, vec[2]));
}

//------------------------------------------------------------------------------
// Construct a Tensor with the given diagonal elements, enforcing a maximum.
//------------------------------------------------------------------------------
SPHERAL_HOST_DEVICE inline
Dim<1>::Tensor
constructTensorWithMaxDiagonal(const Dim<1>::Vector& vec,
                               const double minval) {
  return Dim<1>::Tensor(std::max(minval, vec[0]));
}

SPHERAL_HOST_DEVICE inline
Dim<2>::Tensor
constructTensorWithMaxDiagonal(const Dim<2>::Vector& vec,
                               const double minval) {
  return Dim<2>::Tensor(std::max(minval, vec[0]), 0.0,
                        0.0, std::max(minval, vec[1]));
}

SPHERAL_HOST_DEVICE inline
Dim<3>::Tensor
constructTensorWithMaxDiagonal(const Dim<3>::Vector& vec,
                               const double minval) {
  return Dim<3>::Tensor(std::max(minval, vec[0]), 0.0, 0.0,
                        0.0, std::max(minval, vec[1]), 0.0,
                        0.0, 0.0, std::max(minval, vec[2]));
}

SPHERAL_HOST_DEVICE inline
Dim<1>::SymTensor
constructSymTensorWithMaxDiagonal(const Dim<1>::Vector& vec,
                                  const double minval) {
  return Dim<1>::SymTensor(std::max(minval, vec[0]));
}

SPHERAL_HOST_DEVICE inline
Dim<2>::SymTensor
constructSymTensorWithMaxDiagonal(const Dim<2>::Vector& vec,
                                  const double minval) {
  return Dim<2>::SymTensor(std::max(minval, vec[0]), 0.0,
                           0.0, std::max(minval, vec[1]));
}

SPHERAL_HOST_DEVICE inline
Dim<3>::SymTensor
constructSymTensorWithMaxDiagonal(const Dim<3>::Vector& vec,
                                  const double minval) {
  return Dim<3>::SymTensor(std::max(minval, vec[0]), 0.0, 0.0,
                           0.0, std::max(minval, vec[1]), 0.0,
                           0.0, 0.0, std::max(minval, vec[2]));
}

//------------------------------------------------------------------------------
// Construct a Tensor with the given diagonal elements, constrained to the 
// given range.
//------------------------------------------------------------------------------
SPHERAL_HOST_DEVICE inline
Dim<1>::Tensor
constructTensorWithBoundedDiagonal(const Dim<1>::Vector& vec,
                                   const double minval,
                                   const double maxval) {
  return Dim<1>::Tensor(std::max(minval, std::min(maxval, vec[0])));
}

SPHERAL_HOST_DEVICE inline
Dim<2>::Tensor
constructTensorWithBoundedDiagonal(const Dim<2>::Vector& vec,
                                   const double minval,
                                   const double maxval) {
  return Dim<2>::Tensor(std::max(minval, std::min(maxval, vec[0])), 0.0,
                        0.0, std::max(minval, std::min(maxval, vec[1])));
}

SPHERAL_HOST_DEVICE inline
Dim<3>::Tensor
constructTensorWithBoundedDiagonal(const Dim<3>::Vector& vec,
                                   const double minval,
                                   const double maxval) {
  return Dim<3>::Tensor(std::max(minval, std::min(maxval, vec[0])), 0.0, 0.0,
                        0.0, std::max(minval, std::min(maxval, vec[1])), 0.0,
                        0.0, 0.0, std::max(minval, std::min(maxval, vec[2])));
}

SPHERAL_HOST_DEVICE inline
Dim<1>::SymTensor
constructSymTensorWithBoundedDiagonal(const Dim<1>::Vector& vec,
                                      const double minval,
                                      const double maxval) {
  return Dim<1>::SymTensor(std::max(minval, std::min(maxval, vec[0])));
}

SPHERAL_HOST_DEVICE inline
Dim<2>::SymTensor
constructSymTensorWithBoundedDiagonal(const Dim<2>::Vector& vec,
                                      const double minval,
                                      const double maxval) {
  return Dim<2>::SymTensor(std::max(minval, std::min(maxval, vec[0])), 0.0,
                           0.0, std::max(minval, std::min(maxval, vec[1])));
}

SPHERAL_HOST_DEVICE inline
Dim<3>::SymTensor
constructSymTensorWithBoundedDiagonal(const Dim<3>::Vector& vec,
                                      const double minval,
                                      const double maxval) {
  return Dim<3>::SymTensor(std::max(minval, std::min(maxval, vec[0])), 0.0, 0.0,
                           0.0, std::max(minval, std::min(maxval, vec[1])), 0.0,
                           0.0, 0.0, std::max(minval, std::min(maxval, vec[2])));
}

//------------------------------------------------------------------------------
// Construct a SymTensor with the given diagonal elements, raised to a power.
//------------------------------------------------------------------------------
SPHERAL_HOST_DEVICE inline
Dim<1>::SymTensor
constructSymTensorWithPowDiagonal(const Dim<1>::Vector& vec,
                                  const double exponent) {
  return Dim<1>::SymTensor(pow(vec[0], exponent));
}

SPHERAL_HOST_DEVICE inline
Dim<2>::SymTensor
constructSymTensorWithPowDiagonal(const Dim<2>::Vector& vec,
                                  const double exponent) {
  return Dim<2>::SymTensor(pow(vec[0], exponent), 0.0,
                           0.0, pow(vec[1], exponent));
}

SPHERAL_HOST_DEVICE inline
Dim<3>::SymTensor
constructSymTensorWithPowDiagonal(const Dim<3>::Vector& vec,
                                  const double exponent) {
  return Dim<3>::SymTensor(pow(vec[0], exponent), 0.0, 0.0,
                           0.0, pow(vec[1], exponent), 0.0,
                           0.0, 0.0, pow(vec[2], exponent));
}

//------------------------------------------------------------------------------
// Construct a tensor with the given value in the given row.
//------------------------------------------------------------------------------
template<typename TensorType, size_t col>
SPHERAL_HOST_DEVICE inline
TensorType
constructTensorWithColumnValue(const double /*val*/) {
  VERIFY(false);
}

template<>
SPHERAL_HOST_DEVICE inline
Dim<1>::Tensor
constructTensorWithColumnValue<Dim<1>::Tensor, (size_t) 0>(const double val) {
  return Dim<1>::Tensor(val);
}

template<>
SPHERAL_HOST_DEVICE inline
Dim<2>::Tensor
constructTensorWithColumnValue<Dim<2>::Tensor, (size_t) 0>(const double val) {
  return Dim<2>::Tensor(val, 0.0,
                        val, 0.0);
}

template<>
SPHERAL_HOST_DEVICE inline
Dim<2>::Tensor
constructTensorWithColumnValue<Dim<2>::Tensor, (size_t) 1>(const double val) {
  return Dim<2>::Tensor(0.0, val,
                        0.0, val);
}

template<>
SPHERAL_HOST_DEVICE inline
Dim<3>::Tensor
constructTensorWithColumnValue<Dim<3>::Tensor, (size_t) 0>(const double val) {
  return Dim<3>::Tensor(val, 0.0, 0.0,
                        val, 0.0, 0.0,
                        val, 0.0, 0.0);
}

template<>
SPHERAL_HOST_DEVICE inline
Dim<3>::Tensor
constructTensorWithColumnValue<Dim<3>::Tensor, (size_t) 1>(const double val) {
  return Dim<3>::Tensor(0.0, val, 0.0,
                        0.0, val, 0.0,
                        0.0, val, 0.0);
}

template<>
SPHERAL_HOST_DEVICE inline
Dim<3>::Tensor
constructTensorWithColumnValue<Dim<3>::Tensor, (size_t) 2>(const double val) {
  return Dim<3>::Tensor(0.0, 0.0, val,
                        0.0, 0.0, val,
                        0.0, 0.0, val);
}

//------------------------------------------------------------------------------
// In place addition of the absolute value of the elements of a tensor.
//------------------------------------------------------------------------------
template<typename TensorType>
SPHERAL_HOST_DEVICE inline
void
inPlaceAbsAdd(TensorType& /*lhs*/, const TensorType& /*rhs*/) {
  VERIFY(false);
}

template<>
SPHERAL_HOST_DEVICE inline
void
inPlaceAbsAdd<Dim<1>::Tensor>(Dim<1>::Tensor& lhs, const Dim<1>::Tensor& rhs) {
  lhs.xx(lhs.xx() + std::abs(rhs.xx()));
}

template<>
SPHERAL_HOST_DEVICE inline
void
inPlaceAbsAdd<Dim<2>::Tensor>(Dim<2>::Tensor& lhs, const Dim<2>::Tensor& rhs) {
  lhs.xx(lhs.xx() + std::abs(rhs.xx()));
  lhs.xy(lhs.xy() + std::abs(rhs.xy()));
  lhs.yx(lhs.yx() + std::abs(rhs.yx()));
  lhs.yy(lhs.yy() + std::abs(rhs.yy()));
}

template<>
SPHERAL_HOST_DEVICE inline
void
inPlaceAbsAdd<Dim<3>::Tensor>(Dim<3>::Tensor& lhs, const Dim<3>::Tensor& rhs) {
  lhs.xx(lhs.xx() + std::abs(rhs.xx()));
  lhs.xy(lhs.xy() + std::abs(rhs.xy()));
  lhs.xz(lhs.xz() + std::abs(rhs.xz()));
  lhs.yx(lhs.yx() + std::abs(rhs.yx()));
  lhs.yy(lhs.yy() + std::abs(rhs.yy()));
  lhs.yz(lhs.yz() + std::abs(rhs.yz()));
  lhs.zx(lhs.zx() + std::abs(rhs.zx()));
  lhs.zy(lhs.zy() + std::abs(rhs.zy()));
  lhs.zz(lhs.zz() + std::abs(rhs.zz()));
}

//------------------------------------------------------------------------------
// In place element wise division of the elements of a tensor.
//------------------------------------------------------------------------------
template<typename TensorType>
SPHERAL_HOST_DEVICE inline
void
tensorElementWiseDivide(TensorType& /*lhs*/, const TensorType& /*rhs*/) {
  VERIFY(false);
}

template<>
SPHERAL_HOST_DEVICE inline
void
tensorElementWiseDivide<Dim<1>::Tensor>(Dim<1>::Tensor& lhs, const Dim<1>::Tensor& rhs) {
  REQUIRE(rhs.xx() != 0.0);
  lhs.xx(lhs.xx() / rhs.xx());
}

template<>
SPHERAL_HOST_DEVICE inline
void
tensorElementWiseDivide<Dim<2>::Tensor>(Dim<2>::Tensor& lhs, const Dim<2>::Tensor& rhs) {
  REQUIRE(std::find_if(rhs.begin(), rhs.end(),
                       [] (double val) { return val == 0.0; }) == rhs.end());
  lhs.xx(lhs.xx() / rhs.xx());
  lhs.xy(lhs.xy() / rhs.xy());
  lhs.yx(lhs.yx() / rhs.yx());
  lhs.yy(lhs.yy() / rhs.yy());
}

template<>
SPHERAL_HOST_DEVICE inline
void
tensorElementWiseDivide<Dim<3>::Tensor>(Dim<3>::Tensor& lhs, const Dim<3>::Tensor& rhs) {
  REQUIRE(std::find_if(rhs.begin(), rhs.end(),
                       [] (double val) { return val == 0.0; }) == rhs.end());
  lhs.xx(lhs.xx() / rhs.xx());
  lhs.xy(lhs.xy() / rhs.xy());
  lhs.xz(lhs.xz() / rhs.xz());
  lhs.yx(lhs.yx() / rhs.yx());
  lhs.yy(lhs.yy() / rhs.yy());
  lhs.yz(lhs.yz() / rhs.yz());
  lhs.zx(lhs.zx() / rhs.zx());
  lhs.zy(lhs.zy() / rhs.zy());
  lhs.zz(lhs.zz() / rhs.zz());
}

}

#endif
