#include "GeomVector.hh"
#include "GeomTensor.hh"
#include "EigenStruct.hh"
#include "Utilities/SpheralFunctions.hh"
#include "Utilities/FastMath.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/DBC.hh"

#include <algorithm>
#include <string>
#include <limits>

namespace Spheral {

//------------------------------------------------------------------------------
// Return the element index corresponding to the given (row,column)
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<1>::size_type
GeomSymmetricTensor<1>::elementIndex(const GeomSymmetricTensor<1>::size_type row,
                                     const GeomSymmetricTensor<1>::size_type column) const {
  CONTRACT_VAR(row);
  CONTRACT_VAR(column);
  REQUIRE(row < 3u);
  REQUIRE(column < 3u);
  REQUIRE(row == column);
  return row;
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<2>::size_type
GeomSymmetricTensor<2>::elementIndex(const GeomSymmetricTensor<2>::size_type row,
                                     const GeomSymmetricTensor<2>::size_type column) const {
  REQUIRE((row < 2u and column < 2u) or
          (row == 2u and column == 2u));
  return (row == 2u ? 3u :
          row == 1u ? column + 1u :
          column);
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<3>::size_type
GeomSymmetricTensor<3>::elementIndex(const GeomSymmetricTensor<3>::size_type row,
                                     const GeomSymmetricTensor<3>::size_type column) const {
  REQUIRE(row < 3u);
  REQUIRE(column < 3u);
  const auto i = std::min(row, column);
  const auto j = std::max(row, column);
  const auto result = (7u - i)*i/2u + j - i;
  ENSURE(result < numElements());
  return result;
}

//------------------------------------------------------------------------------
// Construct from a single scalar
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<nDim>::
GeomSymmetricTensor(const double a):
  GeomSymmetricTensorBase<nDim>() {
  (*this) = one() * a;
}

//------------------------------------------------------------------------------
// Construct with the given values for the elements.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<1>::
GeomSymmetricTensor(const double a11,
                    const double a22,
                    const double a33):
  GeomSymmetricTensorBase<1>(a11, a22, a33) {
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<2>::
GeomSymmetricTensor(const double a11, const double a12, 
                    const double a21, const double a22,
                    const double a33):
  GeomSymmetricTensorBase<2>(a11, a12,
                                  a22,
                                       a33) {
  CONTRACT_VAR(a21);
  REQUIRE(a12 == a21);
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<3>::
GeomSymmetricTensor(const double a11, const double a12, const double a13,
                    const double a21, const double a22, const double a23,
                    const double a31, const double a32, const double a33):
  GeomSymmetricTensorBase<3>(a11, a12, a13,
                                  a22, a23,
                                       a33) {
  CONTRACT_VAR(a21);
  CONTRACT_VAR(a31);
  CONTRACT_VAR(a32);
  REQUIRE(a12 == a21);
  REQUIRE(a13 == a31);
  REQUIRE(a23 == a32);
}

//------------------------------------------------------------------------------
// Copy constructor (tensor)
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<nDim>::
GeomSymmetricTensor(const GeomTensor<nDim>& rhs) {
  this->operator=(rhs.Symmetric());
}

//------------------------------------------------------------------------------
// Copy constructor (3D tensors)
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<1>::
GeomSymmetricTensor(const GeomSymmetricTensor<3>& rhs):
  GeomSymmetricTensorBase<1>(rhs.xx(), rhs.yy(), rhs.zz()) {
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<2>::
GeomSymmetricTensor(const GeomSymmetricTensor<3>& rhs):
  GeomSymmetricTensorBase<2>(rhs.xx(), rhs.xy(),
                                       rhs.yy(), rhs.zz()) {
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<1>::
GeomSymmetricTensor(const GeomTensor<3>& rhs):
  GeomSymmetricTensorBase<1>(rhs.xx(), rhs.yy(), rhs.zz()) {
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<2>::
GeomSymmetricTensor(const GeomTensor<3>& rhs):
  GeomSymmetricTensorBase<2>(rhs.xx(), 0.5*(rhs.xy() + rhs.yx()),
                                       rhs.yy(), rhs.zz()) {
}

//------------------------------------------------------------------------------
// Construct from an Eigen Tensor
//------------------------------------------------------------------------------
template<>
template<typename Derived>
inline
GeomSymmetricTensor<1>::GeomSymmetricTensor(const Eigen::MatrixBase<Derived>& ten):
  GeomSymmetricTensorBase<1>(ten(0,0)) {
}

template<>
template<typename Derived>
inline
GeomSymmetricTensor<2>::GeomSymmetricTensor(const Eigen::MatrixBase<Derived>& ten):
  GeomSymmetricTensorBase<2>(ten(0,0), ten(0,1),
                                       ten(1,1)) {
}

template<>
template<typename Derived>
inline
GeomSymmetricTensor<3>::GeomSymmetricTensor(const Eigen::MatrixBase<Derived>& ten):
  GeomSymmetricTensorBase<3>(ten(0,0), ten(0,1), ten(0,2),
                                       ten(1,1), ten(1,2),
                                                 ten(2,2)) {
}

//------------------------------------------------------------------------------
// Construct from diagonal elements
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<nDim>::GeomSymmetricTensor(const VectorType& diagonal):
  GeomSymmetricTensorBase<nDim>() {
  this->mxx = diagonal.x();
  this->myy = diagonal.y();
  this->mzz = diagonal.z();
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<nDim>::GeomSymmetricTensor(const GeomVector<3>& diagonal) requires (nDim < 3):
  GeomSymmetricTensorBase<nDim>() {
  this->mxx = diagonal.x();
  this->myy = diagonal.y();
  this->mzz = diagonal.z();
}

//------------------------------------------------------------------------------
// Construct from diagonal elements with bounding values
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<nDim>::GeomSymmetricTensor(const VectorType& diagonal,
                                               const double minvalue,
                                               const double maxvalue):
  GeomSymmetricTensorBase<nDim>() {
  this->mxx = std::clamp(diagonal.x(), minvalue, maxvalue);
  this->myy = std::clamp(diagonal.y(), minvalue, maxvalue);
  this->mzz = std::clamp(diagonal.z(), minvalue, maxvalue);
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<nDim>::GeomSymmetricTensor(const GeomVector<3>& diagonal,
                                               const double minvalue,
                                               const double maxvalue) requires (nDim < 3):
  GeomSymmetricTensorBase<nDim>() {
  this->mxx = std::clamp(diagonal.x(), minvalue, maxvalue);
  this->myy = std::clamp(diagonal.y(), minvalue, maxvalue);
  this->mzz = std::clamp(diagonal.z(), minvalue, maxvalue);
}

//------------------------------------------------------------------------------
// Assignment operators.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<1>&
GeomSymmetricTensor<1>::
operator=(const GeomTensor<1>& rhs) {
  this->mxx = rhs.xx();
  this->myy = rhs.yy();
  this->mzz = rhs.zz();
  return *this;
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<2>&
GeomSymmetricTensor<2>::
operator=(const GeomTensor<2>& rhs) {
  REQUIRE(rhs.xy() == rhs.yx());
  this->mxx = rhs.xx();
  this->mxy = 0.5*(rhs.xy() + rhs.yx());
  this->myy = rhs.yy();
  this->mzz = rhs.zz();
  return *this;
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<3>&
GeomSymmetricTensor<3>::
operator=(const GeomTensor<3>& rhs) {
  this->mxx = rhs.xx();
  this->mxy = 0.5*(rhs.xy() + rhs.yx());
  this->mxz = 0.5*(rhs.xz() + rhs.zx());
  this->myy = rhs.yy();
  this->myz = 0.5*(rhs.yz() + rhs.zy());
  this->mzz = rhs.zz();
  return *this;
}

//------------------------------------------------------------------------------
// The assignment operator (Eigen Tensor).
//------------------------------------------------------------------------------
template<>
template<typename Derived>
inline
GeomSymmetricTensor<1>&
GeomSymmetricTensor<1>::operator=(const Eigen::MatrixBase<Derived>& ten) {
  this->mxx = ten(0,0);
  return *this;
}

template<>
template<typename Derived>
inline
GeomSymmetricTensor<2>&
GeomSymmetricTensor<2>::operator=(const Eigen::MatrixBase<Derived>& ten) {
  this->mxx = ten(0,0);
  this->mxy = ten(0,1);
  this->myy = ten(1,1);
  return *this;
}

template<>
template<typename Derived>
inline
GeomSymmetricTensor<3>&
GeomSymmetricTensor<3>::operator=(const Eigen::MatrixBase<Derived>& ten) {
  this->mxx = ten(0,0);
  this->mxy = ten(0,1);
  this->mxz = ten(0,2);
  this->myy = ten(1,1);
  this->myz = ten(1,2);
  this->mzz = ten(2,2);
  return *this;
}

//------------------------------------------------------------------------------
// Access the elements by indicies.
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
double
GeomSymmetricTensor<nDim>::operator()(const typename GeomSymmetricTensor<nDim>::size_type row,
                                      const typename GeomSymmetricTensor<nDim>::size_type column) const {
  REQUIRE(row < nDim or row == nDim);
  REQUIRE(column < nDim or row == nDim);
  return *(begin() + elementIndex(row, column));
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
double&
GeomSymmetricTensor<nDim>::operator()(const typename GeomSymmetricTensor<nDim>::size_type row,
                                      const typename GeomSymmetricTensor<nDim>::size_type column) {
  REQUIRE(row < nDim or row == nDim);
  REQUIRE(column < nDim or row == nDim);
  return *(begin() + elementIndex(row, column));
}

//------------------------------------------------------------------------------
// Return the (index) element using the bracket operator.
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
double
GeomSymmetricTensor<nDim>::operator[](typename GeomSymmetricTensor<nDim>::size_type index) const {
  REQUIRE(index < numElements());
  return *(begin() + index);
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
double&
GeomSymmetricTensor<nDim>::operator[](typename GeomSymmetricTensor<nDim>::size_type index) {
  REQUIRE(index < numElements());
  return *(begin() + index);
}

//------------------------------------------------------------------------------
// Return the individual elements, mapped as:
//    xx, xy, xz     11, 12, 13
//    yx, yy, yz  =  21, 22, 23
//    zx, zy, zz     31, 32, 33
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
double
GeomSymmetricTensor<nDim>::xx() const {
  return this->mxx;
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
double
GeomSymmetricTensor<nDim>::xy() const {
  return this->mxy;
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
double
GeomSymmetricTensor<nDim>::xz() const {
  return this->mxz;
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
double
GeomSymmetricTensor<nDim>::yx() const {
  return this->mxy;
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
double
GeomSymmetricTensor<nDim>::yy() const {
  return this->myy;
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
double
GeomSymmetricTensor<nDim>::yz() const {
  return this->myz;
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
double
GeomSymmetricTensor<nDim>::zx() const {
  return this->mxz;
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
double
GeomSymmetricTensor<nDim>::zy() const {
  return this->myz;
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
double
GeomSymmetricTensor<nDim>::zz() const {
  return this->mzz;
}

//------------------------------------------------------------------------------
// 1D dummy elements
template<> SPHERAL_HOST_DEVICE inline double GeomSymmetricTensor<1>::xy() const { return 0.0; }
template<> SPHERAL_HOST_DEVICE inline double GeomSymmetricTensor<1>::xz() const { return 0.0; }
template<> SPHERAL_HOST_DEVICE inline double GeomSymmetricTensor<1>::yx() const { return 0.0; }
template<> SPHERAL_HOST_DEVICE inline double GeomSymmetricTensor<1>::yz() const { return 0.0; }
template<> SPHERAL_HOST_DEVICE inline double GeomSymmetricTensor<1>::zx() const { return 0.0; }
template<> SPHERAL_HOST_DEVICE inline double GeomSymmetricTensor<1>::zy() const { return 0.0; }

//------------------------------------------------------------------------------
// 2D dummy elements
template<> SPHERAL_HOST_DEVICE inline double GeomSymmetricTensor<2>::xz() const { return 0.0; }
template<> SPHERAL_HOST_DEVICE inline double GeomSymmetricTensor<2>::yz() const { return 0.0; }
template<> SPHERAL_HOST_DEVICE inline double GeomSymmetricTensor<2>::zx() const { return 0.0; }
template<> SPHERAL_HOST_DEVICE inline double GeomSymmetricTensor<2>::zy() const { return 0.0; }

//------------------------------------------------------------------------------
// Set the individual elements, as above.
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
void
GeomSymmetricTensor<nDim>::xx(const double val) {
  this->mxx = val;
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
void
GeomSymmetricTensor<nDim>::xy(const double val) {
  this->mxy = val;
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
void
GeomSymmetricTensor<nDim>::xz(const double val) {
  this->mxz = val;
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
void
GeomSymmetricTensor<nDim>::yx(const double val) {
  this->mxy = val;
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
void
GeomSymmetricTensor<nDim>::yy(const double val) {
  this->myy = val;
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
void
GeomSymmetricTensor<nDim>::yz(const double val) {
  this->myz = val;
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
void
GeomSymmetricTensor<nDim>::zx(const double val) {
  this->mxz = val;
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
void
GeomSymmetricTensor<nDim>::zy(const double val) {
  this->myz = val;
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
void
GeomSymmetricTensor<nDim>::zz(double val) {
  this->mzz = val;
}

//------------------------------------------------------------------------------
// 1D dummy elements
template<> SPHERAL_HOST_DEVICE inline void GeomSymmetricTensor<1>::xy(const double /*val*/) {}
template<> SPHERAL_HOST_DEVICE inline void GeomSymmetricTensor<1>::xz(const double /*val*/) {}
template<> SPHERAL_HOST_DEVICE inline void GeomSymmetricTensor<1>::yx(const double /*val*/) {}
template<> SPHERAL_HOST_DEVICE inline void GeomSymmetricTensor<1>::yz(const double /*val*/) {}
template<> SPHERAL_HOST_DEVICE inline void GeomSymmetricTensor<1>::zx(const double /*val*/) {}
template<> SPHERAL_HOST_DEVICE inline void GeomSymmetricTensor<1>::zy(const double /*val*/) {}

//------------------------------------------------------------------------------
// 2D dummy elements
template<> SPHERAL_HOST_DEVICE inline void GeomSymmetricTensor<2>::xz(const double /*val*/) {}
template<> SPHERAL_HOST_DEVICE inline void GeomSymmetricTensor<2>::yz(const double /*val*/) {}
template<> SPHERAL_HOST_DEVICE inline void GeomSymmetricTensor<2>::zx(const double /*val*/) {}
template<> SPHERAL_HOST_DEVICE inline void GeomSymmetricTensor<2>::zy(const double /*val*/) {}

//------------------------------------------------------------------------------
// Access the individual rows of the GeomSymmetricTensor.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<1>
GeomSymmetricTensor<1>::
getRow(const GeomSymmetricTensor<2>::size_type index) const {
  REQUIRE(index < 1);
  return GeomVector<1>((*this)(index, 0));
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<2>
GeomSymmetricTensor<2>::
getRow(const GeomSymmetricTensor<2>::size_type index) const {
  REQUIRE(index < 2);
  return GeomVector<2>((*this)(index, 0), (*this)(index, 1));
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<3>
GeomSymmetricTensor<3>::
getRow(const GeomSymmetricTensor<3>::size_type index) const {
  REQUIRE(index < 3);
  return GeomVector<3>((*this)(index, 0), (*this)(index, 1), (*this)(index, 2));
}

//------------------------------------------------------------------------------
// Access the individual columns of the GeomSymmetricTensor.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<1>
GeomSymmetricTensor<1>::
getColumn(const GeomSymmetricTensor<2>::size_type index) const {
  REQUIRE(index < 1);
  return GeomVector<1>((*this)(0, index));
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<2>
GeomSymmetricTensor<2>::
getColumn(const GeomSymmetricTensor<2>::size_type index) const {
  REQUIRE(index < 2);
  return GeomVector<2>((*this)(0, index), (*this)(1, index));
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<3>
GeomSymmetricTensor<3>::
getColumn(const GeomSymmetricTensor<3>::size_type index) const {
  REQUIRE(index < 3);
  return GeomVector<3>((*this)(0, index), (*this)(1, index), (*this)(2, index));
}

//------------------------------------------------------------------------------
// Set a row of the GeomSymmetricTensor to a GeomVector.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
void
GeomSymmetricTensor<1>::
setRow(const GeomSymmetricTensor<1>::size_type index,
       const GeomVector<1>& vec) {
  REQUIRE(index < 1);
  (*this)(index, 0) = vec(0);
}

template<>
SPHERAL_HOST_DEVICE
inline
void
GeomSymmetricTensor<2>::
setRow(const GeomSymmetricTensor<2>::size_type index,
       const GeomVector<2>& vec) {
  REQUIRE(index < 2);
  (*this)(index, 0) = vec(0);
  (*this)(index, 1) = vec(1);
}

template<>
SPHERAL_HOST_DEVICE
inline
void
GeomSymmetricTensor<3>::
setRow(const GeomSymmetricTensor<3>::size_type index,
       const GeomVector<3>& vec) {
  REQUIRE(index < 3);
  (*this)(index, 0) = vec(0);
  (*this)(index, 1) = vec(1);
  (*this)(index, 2) = vec(2);
}

//------------------------------------------------------------------------------
// Set a column of the GeomSymmetricTensor to a GeomVector.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
void
GeomSymmetricTensor<1>::
setColumn(const GeomSymmetricTensor<1>::size_type index,
          const GeomVector<1>& vec) {
  REQUIRE(index < 1);
  (*this)(0, index) = vec(0);
}

template<>
SPHERAL_HOST_DEVICE
inline
void
GeomSymmetricTensor<2>::
setColumn(const GeomSymmetricTensor<2>::size_type index,
          const GeomVector<2>& vec) {
  REQUIRE(index < 2);
  (*this)(0, index) = vec(0);
  (*this)(1, index) = vec(1);
}

template<>
SPHERAL_HOST_DEVICE
inline
void
GeomSymmetricTensor<3>::
setColumn(const GeomSymmetricTensor<3>::size_type index,
          const GeomVector<3>& vec) {
  REQUIRE(index < 3);
  (*this)(0, index) = vec(0);
  (*this)(1, index) = vec(1);
  (*this)(2, index) = vec(2);
}

//------------------------------------------------------------------------------
// Iterators to the raw data.
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
typename GeomSymmetricTensor<nDim>::iterator
GeomSymmetricTensor<nDim>::begin() {
  return &(this->mxx);
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
typename GeomSymmetricTensor<nDim>::iterator
GeomSymmetricTensor<nDim>::end() {
  return &(this->mxx) + numElements();
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
typename GeomSymmetricTensor<nDim>::const_iterator
GeomSymmetricTensor<nDim>::begin() const{
  return &(this->mxx);
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
typename GeomSymmetricTensor<nDim>::const_iterator
GeomSymmetricTensor<nDim>::end() const {
  return &(this->mxx) + numElements();
}

//------------------------------------------------------------------------------
// Zero out the GeomSymmetricTensor.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
void
GeomSymmetricTensor<1>::Zero() {
  this->mxx = 0.0;
  this->myy = 0.0;
  this->mzz = 0.0;
}

template<>
SPHERAL_HOST_DEVICE
inline
void
GeomSymmetricTensor<2>::Zero() {
  this->mxx = 0.0;
  this->mxy = 0.0;
  this->myy = 0.0;
  this->mzz = 0.0;
}

template<>
SPHERAL_HOST_DEVICE
inline
void
GeomSymmetricTensor<3>::Zero() {
  this->mxx = 0.0;
  this->mxy = 0.0;
  this->mxz = 0.0;
  this->myy = 0.0;
  this->myz = 0.0;
  this->mzz = 0.0;
}

//------------------------------------------------------------------------------
// Force the tensor to be the identity tensor.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
void
GeomSymmetricTensor<1>::Identity() {
  this->mxx = 1.0;
  this->myy = 1.0;
  this->mzz = 1.0;
}

template<>
SPHERAL_HOST_DEVICE
inline
void
GeomSymmetricTensor<2>::Identity() {
  this->mxx = 1.0;
  this->mxy = 0.0;
  this->myy = 1.0;
  this->mzz = 1.0;
}

template<>
SPHERAL_HOST_DEVICE
inline
void
GeomSymmetricTensor<3>::Identity() {
  this->mxx = 1.0;
  this->mxy = 0.0;
  this->mxz = 0.0;
  this->myy = 1.0;
  this->myz = 0.0;
  this->mzz = 1.0;
}

//------------------------------------------------------------------------------
// Return the negative of a tensor.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<1>
GeomSymmetricTensor<1>::operator-() const {
  return GeomSymmetricTensor<1>(-(this->mxx),
                                -(this->myy),
                                -(this->mzz));
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<2>
GeomSymmetricTensor<2>::operator-() const {
  return GeomSymmetricTensor<2>(-(this->mxx), -(this->mxy),
                                -(this->mxy), -(this->myy),
                                -(this->mzz));
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<3>
GeomSymmetricTensor<3>::operator-() const {
  return GeomSymmetricTensor<3>(-(this->mxx), -(this->mxy), -(this->mxz),
                                -(this->mxy), -(this->myy), -(this->myz),
                                -(this->mxz), -(this->myz), -(this->mzz));
}

//------------------------------------------------------------------------------
// Add two tensors.
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
GeomTensor<nDim>
GeomSymmetricTensor<nDim>::
operator+(const GeomTensor<nDim>& rhs) const {
  GeomTensor<nDim> result(*this);
  result += rhs;
  return result;
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<nDim>
GeomSymmetricTensor<nDim>::
operator+(const GeomSymmetricTensor<nDim>& rhs) const {
  GeomSymmetricTensor<nDim> result(*this);
  result += rhs;
  return result;
}

//------------------------------------------------------------------------------
// Subtract a tensor from another.
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
GeomTensor<nDim>
GeomSymmetricTensor<nDim>::
operator-(const GeomTensor<nDim>& rhs) const {
  GeomTensor<nDim> result(*this);
  result -= rhs;
  return result;
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<nDim>
GeomSymmetricTensor<nDim>::
operator-(const GeomSymmetricTensor<nDim>& rhs) const {
  GeomSymmetricTensor<nDim> result(*this);
  result -= rhs;
  return result;
}

//------------------------------------------------------------------------------
// Multiply two tensors.
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
GeomTensor<nDim>
GeomSymmetricTensor<nDim>::
operator*(const GeomTensor<nDim>& rhs) const {
  return this->dot(rhs);
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
GeomTensor<nDim>
GeomSymmetricTensor<nDim>::
operator*(const GeomSymmetricTensor<nDim>& rhs) const {
  return this->dot(rhs);
}

//------------------------------------------------------------------------------
// Multiply a tensor with a vector.
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
GeomVector<nDim>
GeomSymmetricTensor<nDim>::operator*(const GeomVector<nDim>& rhs) const {
  return this->dot(rhs);
}

//------------------------------------------------------------------------------
// Multiply a tensor by a scalar
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<1>
GeomSymmetricTensor<1>::operator*(const double rhs) const {
  return GeomSymmetricTensor<1>(this->mxx * rhs,
                                this->myy * rhs,
                                this->mzz * rhs);
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<2>
GeomSymmetricTensor<2>::operator*(const double rhs) const {
  GeomSymmetricTensor<2> result(*this);
  result.mxx *= rhs;
  result.mxy *= rhs;
  result.myy *= rhs;
  result.mzz *= rhs;
  return result;
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<3>
GeomSymmetricTensor<3>::operator*(const double rhs) const {
  GeomSymmetricTensor<3> result(*this);
  result.mxx *= rhs;
  result.mxy *= rhs;
  result.mxz *= rhs;
  result.myy *= rhs;
  result.myz *= rhs;
  result.mzz *= rhs;
  return result;
}

//------------------------------------------------------------------------------
// Divide a tensor by a scalar
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<1>
GeomSymmetricTensor<1>::operator/(const double rhs) const {
  REQUIRE(rhs != 0.0);
  const double rhsInv = safeInvVar(rhs);
  return GeomSymmetricTensor<1>(this->mxx * rhsInv,
                                this->myy * rhsInv,
                                this->mzz * rhsInv);
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<2>
GeomSymmetricTensor<2>::operator/(const double rhs) const {
  REQUIRE(rhs != 0.0);
  const double rhsInv = safeInvVar(rhs);
  GeomSymmetricTensor<2> result(*this);
  result.mxx *= rhsInv;
  result.mxy *= rhsInv;
  result.myy *= rhsInv;
  result.mzz *= rhsInv;
  return result;
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<3>
GeomSymmetricTensor<3>::operator/(const double rhs) const {
  REQUIRE(rhs != 0.0);
  const double rhsInv = safeInvVar(rhs);
  GeomSymmetricTensor<3> result(*this);
  result.mxx *= rhsInv;
  result.mxy *= rhsInv;
  result.mxz *= rhsInv;
  result.myy *= rhsInv;
  result.myz *= rhsInv;
  result.mzz *= rhsInv;
  return result;
}

//------------------------------------------------------------------------------
// += symmetric tensor
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<1>&
GeomSymmetricTensor<1>::operator+=(const GeomSymmetricTensor<1>& rhs) {
  this->mxx += rhs.mxx;
  this->myy += rhs.myy;
  this->mzz += rhs.mzz;
  return *this;
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<2>&
GeomSymmetricTensor<2>::operator+=(const GeomSymmetricTensor<2>& rhs) {
  this->mxx += rhs.mxx;
  this->mxy += rhs.mxy;
  this->myy += rhs.myy;
  this->mzz += rhs.mzz;
  return *this;
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<3>&
GeomSymmetricTensor<3>::operator+=(const GeomSymmetricTensor<3>& rhs) {
  this->mxx += rhs.mxx;
  this->mxy += rhs.mxy;
  this->mxz += rhs.mxz;
  this->myy += rhs.myy;
  this->myz += rhs.myz;
  this->mzz += rhs.mzz;
  return *this;
}

//------------------------------------------------------------------------------
// -= symmetric tensor
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<1>&
GeomSymmetricTensor<1>::operator-=(const GeomSymmetricTensor<1>& rhs) {
  this->mxx -= rhs.mxx;
  this->myy -= rhs.myy;
  this->mzz -= rhs.mzz;
  return *this;
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<2>&
GeomSymmetricTensor<2>::operator-=(const GeomSymmetricTensor<2>& rhs) {
  this->mxx -= rhs.mxx;
  this->mxy -= rhs.mxy;
  this->myy -= rhs.myy;
  this->mzz -= rhs.mzz;
  return *this;
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<3>&
GeomSymmetricTensor<3>::operator-=(const GeomSymmetricTensor<3>& rhs) {
  this->mxx -= rhs.mxx;
  this->mxy -= rhs.mxy;
  this->mxz -= rhs.mxz;
  this->myy -= rhs.myy;
  this->myz -= rhs.myz;
  this->mzz -= rhs.mzz;
  return *this;
}

//------------------------------------------------------------------------------
// += eigen tensor
//------------------------------------------------------------------------------
template<>
template<typename Derived>
inline
GeomSymmetricTensor<1>&
GeomSymmetricTensor<1>::operator+=(const Eigen::MatrixBase<Derived>& rhs) {
  this->mxx += rhs(0,0);
  return *this;
}

template<>
template<typename Derived>
inline
GeomSymmetricTensor<2>&
GeomSymmetricTensor<2>::operator+=(const Eigen::MatrixBase<Derived>& rhs) {
  REQUIRE(fuzzyEqual(rhs(0,1), rhs(1,0), 1.e-10));
  this->mxx += rhs(0,0);
  this->mxy += rhs(0,1);
  this->myy += rhs(1,1);
  return *this;
}

template<>
template<typename Derived>
inline
GeomSymmetricTensor<3>&
GeomSymmetricTensor<3>::operator+=(const Eigen::MatrixBase<Derived>& rhs) {
  REQUIRE(fuzzyEqual(rhs(0,1), rhs(1,0), 1.e-10));
  REQUIRE(fuzzyEqual(rhs(0,2), rhs(2,0), 1.e-10));
  REQUIRE(fuzzyEqual(rhs(1,2), rhs(2,1), 1.e-10));
  this->mxx += rhs(0,0);
  this->mxy += rhs(0,1);
  this->mxz += rhs(0,2);
  this->myy += rhs(1,1);
  this->myz += rhs(1,2);
  this->mzz += rhs(2,2);
  return *this;
}

//------------------------------------------------------------------------------
// -= eigen tensor
//------------------------------------------------------------------------------
template<>
template<typename Derived>
inline
GeomSymmetricTensor<1>&
GeomSymmetricTensor<1>::operator-=(const Eigen::MatrixBase<Derived>& rhs) {
  this->mxx -= rhs(0,0);
  return *this;
}

template<>
template<typename Derived>
inline
GeomSymmetricTensor<2>&
GeomSymmetricTensor<2>::operator-=(const Eigen::MatrixBase<Derived>& rhs) {
  REQUIRE(fuzzyEqual(rhs(0,1), rhs(1,0), 1.e-10));
  this->mxx -= rhs(0,0);
  this->mxy -= rhs(0,1);
  this->myy -= rhs(1,1);
  return *this;
}

template<>
template<typename Derived>
inline
GeomSymmetricTensor<3>&
GeomSymmetricTensor<3>::operator-=(const Eigen::MatrixBase<Derived>& rhs) {
  REQUIRE(fuzzyEqual(rhs(0,1), rhs(1,0), 1.e-10));
  REQUIRE(fuzzyEqual(rhs(0,2), rhs(2,0), 1.e-10));
  REQUIRE(fuzzyEqual(rhs(1,2), rhs(2,1), 1.e-10));
  this->mxx -= rhs(0,0);
  this->mxy -= rhs(0,1);
  this->mxz -= rhs(0,2);
  this->myy -= rhs(1,1);
  this->myz -= rhs(1,2);
  this->mzz -= rhs(2,2);
  return *this;
}

//------------------------------------------------------------------------------
// Multiply this tensor by a scalar in place.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<1>&
GeomSymmetricTensor<1>::operator*=(const double rhs) {
  this->mxx *= rhs;
  this->myy *= rhs;
  this->mzz *= rhs;
  return *this;
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<2>&
GeomSymmetricTensor<2>::operator*=(const double rhs) {
  this->mxx *= rhs;
  this->mxy *= rhs;
  this->myy *= rhs;
  this->mzz *= rhs;
  return *this;
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<3>&
GeomSymmetricTensor<3>::operator*=(const double rhs) {
  this->mxx *= rhs;
  this->mxy *= rhs;
  this->mxz *= rhs;
  this->myy *= rhs;
  this->myz *= rhs;
  this->mzz *= rhs;
  return *this;
}

//------------------------------------------------------------------------------
// Divide this tensor by a scalar in place
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<1>&
GeomSymmetricTensor<1>::operator/=(const double rhs) {
  REQUIRE(rhs != 0.0);
  const auto rhsInv = safeInvVar(rhs);
  this->mxx *= rhsInv;
  this->myy *= rhsInv;
  this->mzz *= rhsInv;
  return *this;
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<2>&
GeomSymmetricTensor<2>::operator/=(const double rhs) {
  REQUIRE(rhs != 0.0);
  const auto rhsInv = safeInvVar(rhs);
  this->mxx *= rhsInv;
  this->mxy *= rhsInv;
  this->myy *= rhsInv;
  this->mzz *= rhsInv;
  return *this;
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<3>&
GeomSymmetricTensor<3>::operator/=(const double rhs) {
  REQUIRE(rhs != 0.0);
  const double rhsInv = safeInvVar(rhs);
  this->mxx *= rhsInv;
  this->mxy *= rhsInv;
  this->mxz *= rhsInv;
  this->myy *= rhsInv;
  this->myz *= rhsInv;
  this->mzz *= rhsInv;
  return *this;
}

//------------------------------------------------------------------------------
// Define the equivalence operator.
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
bool
GeomSymmetricTensor<nDim>::
operator==(const GeomTensor<nDim>& rhs) const {
  return *this == rhs.Symmetric();
}

template<>
SPHERAL_HOST_DEVICE
inline
bool
GeomSymmetricTensor<1>::
operator==(const GeomSymmetricTensor<1>& rhs) const {
  return (this->mxx == rhs.mxx and
          this->myy == rhs.myy and
          this->mzz == rhs.mzz);
}

template<>
SPHERAL_HOST_DEVICE
inline
bool
GeomSymmetricTensor<2>::
operator==(const GeomSymmetricTensor<2>& rhs) const {
  return (this->mxx == rhs.mxx and
          this->mxy == rhs.mxy and
          this->myy == rhs.myy and
          this->mzz == rhs.mzz);
}

template<>
SPHERAL_HOST_DEVICE
inline
bool
GeomSymmetricTensor<3>::
operator==(const GeomSymmetricTensor<3>& rhs) const {
  return (this->mxx == rhs.mxx and
          this->mxy == rhs.mxy and
          this->mxz == rhs.mxz and
          this->myy == rhs.myy and
          this->myz == rhs.myz and
          this->mzz == rhs.mzz);
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
bool
GeomSymmetricTensor<nDim>::
operator==(const double rhs) const {
  return *this == (one() * rhs);
}

//------------------------------------------------------------------------------
// Define the not equivalence than comparitor.
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
bool
GeomSymmetricTensor<nDim>::
operator!=(const GeomTensor<nDim>& rhs) const {
  return !(*this == rhs);
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
bool
GeomSymmetricTensor<nDim>::
operator!=(const GeomSymmetricTensor<nDim>& rhs) const {
  return !(*this == rhs);
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
bool
GeomSymmetricTensor<nDim>::
operator!=(const double rhs) const {
  return !(*this == rhs);
}

//------------------------------------------------------------------------------
// Define the less than operator.
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
bool
GeomSymmetricTensor<nDim>::
operator<(const GeomTensor<nDim>& rhs) const {
  return this->Determinant() < rhs.Determinant();
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
bool
GeomSymmetricTensor<nDim>::
operator<(const GeomSymmetricTensor<nDim>& rhs) const {
  return this->Determinant() < rhs.Determinant();
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
bool
GeomSymmetricTensor<nDim>::
operator<(const double rhs) const {
  const auto ev = this->eigenValues();
  return ev.maxElement() < rhs;
}

//------------------------------------------------------------------------------
// Define the greater than operator.
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
bool
GeomSymmetricTensor<nDim>::
operator>(const GeomTensor<nDim>& rhs) const {
  return this->Determinant() > rhs.Determinant();
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
bool
GeomSymmetricTensor<nDim>::
operator>(const GeomSymmetricTensor<nDim>& rhs) const {
  return this->Determinant() > rhs.Determinant();
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
bool
GeomSymmetricTensor<nDim>::
operator>(const double rhs) const {
  const auto ev = this->eigenValues();
  return ev.minElement() > rhs;
}

//------------------------------------------------------------------------------
// Define the less than or equal operator.
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
bool
GeomSymmetricTensor<nDim>::
operator<=(const GeomTensor<nDim>& rhs) const {
  return (*this < rhs) or (*this == rhs);
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
bool
GeomSymmetricTensor<nDim>::
operator<=(const GeomSymmetricTensor<nDim>& rhs) const {
  return (*this < rhs) or (*this == rhs);
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
bool
GeomSymmetricTensor<nDim>::
operator<=(const double rhs) const {
  return (*this < rhs) or (*this == rhs);
}

//------------------------------------------------------------------------------
// Define the greater than or equal operator.
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
bool
GeomSymmetricTensor<nDim>::
operator>=(const GeomTensor<nDim>& rhs) const {
  return (*this > rhs) or (*this == rhs);
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
bool
GeomSymmetricTensor<nDim>::
operator>=(const GeomSymmetricTensor<nDim>& rhs) const {
  return (*this > rhs) or (*this == rhs);
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
bool
GeomSymmetricTensor<nDim>::
operator>=(const double rhs) const {
  return (*this > rhs) or (*this == rhs);
}

//------------------------------------------------------------------------------
// Return the symmetric part.  A no-op for this class.
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<nDim>
GeomSymmetricTensor<nDim>::Symmetric() const {
  return *this;
}

//------------------------------------------------------------------------------
// Return the skew-symmetric part of a GeomSymmetricTensor.
//   Bij = 0.5*(Aij - Aji)
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
GeomTensor<nDim>
GeomSymmetricTensor<nDim>::SkewSymmetric() const {
  return GeomTensor<nDim>::zero();
}

//------------------------------------------------------------------------------
// Return the transpose of the GeomSymmetricTensor.
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<nDim>
GeomSymmetricTensor<nDim>::
Transpose() const {
  return *this;
}

//------------------------------------------------------------------------------
// Return the inverse of the tensor.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<1>
GeomSymmetricTensor<1>::Inverse() const {
  CHECK(this->mxx != 0.0);
  return GeomSymmetricTensor<1>(safeInvVar(this->mxx),
                                safeInvVar(this->myy),
                                safeInvVar(this->mzz));
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<2>
GeomSymmetricTensor<2>::Inverse() const {
  REQUIRE(Determinant() != 0.0);
  return GeomSymmetricTensor<2>( (this->myy), -(this->mxy),
                                -(this->mxy),  (this->mxx))*safeInvVar(this->Determinant());
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<3>
GeomSymmetricTensor<3>::Inverse() const {
  REQUIRE(Determinant() != 0.0);
  return GeomSymmetricTensor<3>((this->myy)*(this->mzz) - (this->myz)*(this->myz),
                                (this->myz)*(this->mxz) - (this->mxy)*(this->mzz),
                                (this->mxy)*(this->myz) - (this->myy)*(this->mxz),
                                (this->mxz)*(this->myz) - (this->mxy)*(this->mzz),
                                (this->mxx)*(this->mzz) - (this->mxz)*(this->mxz),
                                (this->mxy)*(this->mxz) - (this->mxx)*(this->myz),
                                (this->mxy)*(this->myz) - (this->mxz)*(this->myy),
                                (this->mxz)*(this->mxy) - (this->mxx)*(this->myz),
                                (this->mxx)*(this->myy) - (this->mxy)*(this->mxy))*safeInvVar(this->Determinant());
}

//------------------------------------------------------------------------------
// Return the diagonal elements of the GeomSymmetricTensor.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<1>
GeomSymmetricTensor<1>::diagonalElements() const {
  return GeomVector<1>(this->mxx);
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<2>
GeomSymmetricTensor<2>::diagonalElements() const {
  return GeomVector<2>(this->mxx, this->myy);
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<3>
GeomSymmetricTensor<3>::diagonalElements() const {
  return GeomVector<3>(this->mxx, this->myy, this->mzz);
}

//------------------------------------------------------------------------------
// Return the trace of the tensor.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
double
GeomSymmetricTensor<1>::Trace() const {
  return this->mxx;
}

template<>
SPHERAL_HOST_DEVICE
inline
double
GeomSymmetricTensor<2>::Trace() const {
  return this->mxx + this->myy;
}

template<>
SPHERAL_HOST_DEVICE
inline
double
GeomSymmetricTensor<3>::Trace() const {
  return this->mxx + this->myy + this->mzz;
}

//------------------------------------------------------------------------------
// Return the determinant of the tensor.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
double
GeomSymmetricTensor<1>::Determinant() const {
  return this->mxx;
}

template<>
SPHERAL_HOST_DEVICE
inline
double
GeomSymmetricTensor<2>::Determinant() const {
  return (this->mxx)*(this->myy) - (this->mxy)*(this->mxy);
}

template<>
SPHERAL_HOST_DEVICE
inline
double
GeomSymmetricTensor<3>::Determinant() const {
  return ((this->mxx)*(this->myy)*(this->mzz) + 
	  (this->mxy)*(this->myz)*(this->mxz) + 
	  (this->mxz)*(this->mxy)*(this->myz) - 
	  (this->mxx)*(this->myz)*(this->myz) - 
	  (this->mxy)*(this->mxy)*(this->mzz) - 
	  (this->mxz)*(this->myy)*(this->mxz));
}

//------------------------------------------------------------------------------
// Multiply a tensor with a vector.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<1>
GeomSymmetricTensor<1>::dot(const GeomVector<1>& rhs) const {
  return GeomVector<1>((this->mxx)*(rhs.x()));
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<2>
GeomSymmetricTensor<2>::dot(const GeomVector<2>& rhs) const {
  return GeomVector<2>((this->mxx)*(rhs.x()) + (this->mxy)*(rhs.y()),
                       (this->mxy)*(rhs.x()) + (this->myy)*(rhs.y()));
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<3>
GeomSymmetricTensor<3>::dot(const GeomVector<3>& rhs) const {
  return GeomVector<3>((this->mxx)*(rhs.x()) + (this->mxy)*(rhs.y()) + (this->mxz)*(rhs.z()),
                       (this->mxy)*(rhs.x()) + (this->myy)*(rhs.y()) + (this->myz)*(rhs.z()),
                       (this->mxz)*(rhs.x()) + (this->myz)*(rhs.y()) + (this->mzz)*(rhs.z()));
}

//------------------------------------------------------------------------------
// Multiply two tensors.  This is just the linear algebra definition for matrix
// multiplication.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<1>
GeomSymmetricTensor<1>::dot(const GeomTensor<1>& rhs) const {
  return GeomTensor<1>(this->mxx * rhs.xx(),
                       this->myy * rhs.yy(),
                       this->mzz + rhs.zz());
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<2>
GeomSymmetricTensor<2>::dot(const GeomTensor<2>& rhs) const {
  return GeomTensor<2>(this->mxx * rhs.xx() + this->mxy * rhs.yx(),
                       this->mxx * rhs.xy() + this->mxy * rhs.yy(),
                       this->mxy * rhs.xx() + this->myy * rhs.yx(),
                       this->mxy * rhs.xy() + this->myy * rhs.yy(),
                       this->mzz * rhs.zz());
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<3>
GeomSymmetricTensor<3>::dot(const GeomTensor<3>& rhs) const {
  return GeomTensor<3>(this->mxx * rhs.xx() + this->mxy * rhs.yx() + this->mxz * rhs.zx(),
                       this->mxx * rhs.xy() + this->mxy * rhs.yy() + this->mxz * rhs.zy(),
                       this->mxx * rhs.xz() + this->mxy * rhs.yz() + this->mxz * rhs.zz(),
                       this->mxy * rhs.xx() + this->myy * rhs.yx() + this->myz * rhs.zx(),
                       this->mxy * rhs.xy() + this->myy * rhs.yy() + this->myz * rhs.zy(),
                       this->mxy * rhs.xz() + this->myy * rhs.yz() + this->myz * rhs.zz(),
                       this->mxz * rhs.xx() + this->myz * rhs.yx() + this->mzz * rhs.zx(),
                       this->mxz * rhs.xy() + this->myz * rhs.yy() + this->mzz * rhs.zy(),
                       this->mxz * rhs.xz() + this->myz * rhs.yz() + this->mzz * rhs.zz());
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<1>
GeomSymmetricTensor<1>::dot(const GeomSymmetricTensor<1>& rhs) const {
  return GeomTensor<1>(this->mxx * rhs.xx(),
                       this->myy * rhs.yy(),
                       this->mzz * rhs.zz());
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<2>
GeomSymmetricTensor<2>::dot(const GeomSymmetricTensor<2>& rhs) const {
  return GeomTensor<2>(this->mxx * rhs.xx() + this->mxy * rhs.yx(),
                       this->mxx * rhs.xy() + this->mxy * rhs.yy(),
                       this->mxy * rhs.xx() + this->myy * rhs.yx(),
                       this->mxy * rhs.xy() + this->myy * rhs.yy(),
                       this->mzz * rhs.zz());
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<3>
GeomSymmetricTensor<3>::dot(const GeomSymmetricTensor<3>& rhs) const {
  return GeomTensor<3>(this->mxx * rhs.xx() + this->mxy * rhs.yx() + this->mxz * rhs.zx(),
                       this->mxx * rhs.xy() + this->mxy * rhs.yy() + this->mxz * rhs.zy(),
                       this->mxx * rhs.xz() + this->mxy * rhs.yz() + this->mxz * rhs.zz(),
                       this->mxy * rhs.xx() + this->myy * rhs.yx() + this->myz * rhs.zx(),
                       this->mxy * rhs.xy() + this->myy * rhs.yy() + this->myz * rhs.zy(),
                       this->mxy * rhs.xz() + this->myy * rhs.yz() + this->myz * rhs.zz(),
                       this->mxz * rhs.xx() + this->myz * rhs.yx() + this->mzz * rhs.zx(),
                       this->mxz * rhs.xy() + this->myz * rhs.yy() + this->mzz * rhs.zy(),
                       this->mxz * rhs.xz() + this->myz * rhs.yz() + this->mzz * rhs.zz());
}


//------------------------------------------------------------------------------
// Return the doubledot product.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
double
GeomSymmetricTensor<1>::
doubledot(const GeomTensor<1>& rhs) const {
  return this->mxx * rhs.xx();
}

template<>
SPHERAL_HOST_DEVICE
inline
double
GeomSymmetricTensor<2>::
doubledot(const GeomTensor<2>& rhs) const {
  return (this->mxx * rhs.xx() +this->mxy * rhs.yx() + 
          this->mxy * rhs.xy() +this->myy * rhs.yy());
}

template<>
SPHERAL_HOST_DEVICE
inline
double
GeomSymmetricTensor<3>::
doubledot(const GeomTensor<3>& rhs) const {
  return (this->mxx * rhs.xx() + this->mxy * rhs.yx() + this->mxz * rhs.zx() +
          this->mxy * rhs.xy() + this->myy * rhs.yy() + this->myz * rhs.zy() +
          this->mxz * rhs.xz() + this->myz * rhs.yz() + this->mzz * rhs.zz());
}

//------------------------------------------------------------------------------
// Return the doubledot product with a symmetric tensor.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
double
GeomSymmetricTensor<1>::
doubledot(const GeomSymmetricTensor<1>& rhs) const {
  return this->mxx * rhs.xx();
}

template<>
SPHERAL_HOST_DEVICE
inline
double
GeomSymmetricTensor<2>::
doubledot(const GeomSymmetricTensor<2>& rhs) const {
  return (this->mxx * rhs.xx() + this->mxy * rhs.yx() + 
          this->mxy * rhs.xy() + this->myy * rhs.yy());
}

template<>
SPHERAL_HOST_DEVICE
inline
double
GeomSymmetricTensor<3>::
doubledot(const GeomSymmetricTensor<3>& rhs) const {
  return (this->mxx * rhs.xx() + this->mxy * rhs.yx() + this->mxz * rhs.zx() +
          this->mxy * rhs.xy() + this->myy * rhs.yy() + this->myz * rhs.zy() +
          this->mxz * rhs.xz() + this->myz * rhs.yz() + this->mzz * rhs.zz());
}

//------------------------------------------------------------------------------
// Return the doubledot product with ourself.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
double
GeomSymmetricTensor<1>::
selfDoubledot() const {
  return FastMath::square(this->mxx);
}

template<>
SPHERAL_HOST_DEVICE
inline
double
GeomSymmetricTensor<2>::
selfDoubledot() const {
  return (this->mxx * this->mxx + this->mxy * this->mxy + 
          this->mxy * this->mxy + this->myy * this->myy);
}

template<>
SPHERAL_HOST_DEVICE
inline
double
GeomSymmetricTensor<3>::
selfDoubledot() const {
  return (this->mxx * this->mxx + this->mxy * this->mxy + this->mxz * this->mxz +
          this->mxy * this->mxy + this->myy * this->myy + this->myz * this->myz +
          this->mxz * this->mxz + this->myz * this->myz + this->mzz * this->mzz);
}

//------------------------------------------------------------------------------
// Return the square of this tensor.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<1>
GeomSymmetricTensor<1>::
square() const {
  return GeomSymmetricTensor<1>(this->mxx * this->mxx,
                                this->myy * this->myy,
                                this->mzz * this->mzz);
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<2>
GeomSymmetricTensor<2>::
square() const {
  GeomSymmetricTensor<2> result;
  result.mxx = this->mxx * this->mxx +  this->mxy * this->mxy;
  result.mxy = this->mxy * (this->mxx + this->myy);
  result.myy = this->mxy * this->mxy  +  this->myy * this->myy;
  result.mzz = this->mzz * this->mzz;
  return result;
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<3>
GeomSymmetricTensor<3>::
square() const {
  GeomSymmetricTensor<3> result;
  result.mxx = this->mxx * this->mxx +  this->mxy * this->mxy + this->mxz * this->mxz;
  result.mxy = this->mxy * (this->mxx + this->myy) + this->mxz * this->myz;
  result.mxz = this->mxz * (this->mxx + this->mzz) + this->mxy * this->myz;
  result.myy = this->mxy * this->mxy + this->myy * this->myy + this->myz * this->myz;
  result.myz = this->myz * (this->myy + this->mzz) + this->mxy * this->mxz;
  result.mzz = this->mxz * this->mxz + this->myz * this->myz + this->mzz * this->mzz;
  return result;
}

//------------------------------------------------------------------------------
// Compute the square root of the tensor.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<1>
GeomSymmetricTensor<1>::
sqrt() const {
  return GeomSymmetricTensor<1>(std::sqrt(std::max(0.0, this->mxx)),
                                std::sqrt(std::max(0.0, this->myy)),
                                std::sqrt(std::max(0.0, this->mzz)));
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<2>
GeomSymmetricTensor<2>::
sqrt() const {
  const auto eigen = this->eigenVectors();
  GeomSymmetricTensor<2> result(std::sqrt(std::max(0.0, eigen.eigenValues(0))), 0.0,
                                0.0, std::sqrt(std::max(0.0, eigen.eigenValues(1))),
                                std::sqrt(std::max(0.0, this->mzz)));
  result.rotationalTransform(eigen.eigenVectors);
  return result;
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<3>
GeomSymmetricTensor<3>::
sqrt() const {
  const auto eigen = this->eigenVectors();
  GeomSymmetricTensor<3> result(std::sqrt(std::max(0.0, eigen.eigenValues(0))), 0.0, 0.0,
                                0.0, std::sqrt(std::max(0.0, eigen.eigenValues(1))), 0.0,
                                0.0, 0.0, std::sqrt(std::max(0.0, eigen.eigenValues(2))));
  result.rotationalTransform(eigen.eigenVectors);
  return result;
}

//------------------------------------------------------------------------------
// The general version, raise to an arbitrary power.
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<nDim>
GeomSymmetricTensor<nDim>::
pow(const double p) const {
  const auto eigen = this->eigenVectors3D();
  GeomSymmetricTensor<nDim> result;
  for (auto i = 0u; i < 3u; ++i) {
    result(i,i) = std::pow(std::abs(eigen.eigenValues(i)), p) * sgn(eigen.eigenValues(i));
  }
  result.rotationalTransform(eigen.eigenVectors);
  return result;
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<1>
GeomSymmetricTensor<1>::
pow(const double p) const {
  return GeomSymmetricTensor<1>(std::pow(this->mxx, p),
                                std::pow(this->myy, p),
                                std::pow(this->mzz, p));
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<2>
GeomSymmetricTensor<2>::
pow(const double p) const {
  const auto eigen = this->eigenVectors();
  GeomSymmetricTensor<2> result(std::pow(eigen.eigenValues(0), p), 0.0,
                                0.0, std::pow(eigen.eigenValues(1), p),
                                std::pow(this->mzz, p));
  result.rotationalTransform(eigen.eigenVectors);
  return result;
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<3>
GeomSymmetricTensor<3>::
pow(const double p) const {
  const auto eigen = this->eigenVectors();
  GeomSymmetricTensor<3> result(std::pow(eigen.eigenValues(0), p), 0.0, 0.0,
                                0.0, std::pow(eigen.eigenValues(1), p), 0.0,
                                0.0, 0.0, std::pow(eigen.eigenValues(2), p));
  result.rotationalTransform(eigen.eigenVectors);
  return result;
}

//------------------------------------------------------------------------------
// Element-wise square
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<1>
GeomSymmetricTensor<1>::
squareElements() const {
  return GeomSymmetricTensor<1>(this->mxx * this->mxx,
                                this->myy * this->myy,
                                this->mzz * this->mzz);
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<2>
GeomSymmetricTensor<2>::
squareElements() const {
  GeomSymmetricTensor<2> result;
  result.mxx = this->mxx * this->mxx;
  result.mxy = this->mxy * this->mxy;
  result.myy = this->myy * this->myy;
  result.mzz = this->mzz * this->mzz;
  return result;
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<3>
GeomSymmetricTensor<3>::
squareElements() const {
  GeomSymmetricTensor<3> result;
  result.mxx = this->mxx * this->mxx;
  result.mxy = this->mxy * this->mxy;
  result.mxz = this->mxz * this->mxz;
  result.myy = this->myy * this->myy;
  result.myz = this->myz * this->myz;
  result.mzz = this->mzz * this->mzz;
  return result;
}

//------------------------------------------------------------------------------
// Find the eigenvalues of a tensor.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<1>
GeomSymmetricTensor<1>::eigenValues() const {
  return GeomVector<1>(this->mxx);
}

//----------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<2>
GeomSymmetricTensor<2>::eigenValues() const {
  if (std::abs(xy()) < 1.0e-50) {
    return diagonalElements();
  } else {
    const auto b = Trace();
    const auto c = Determinant();
    const auto q = 0.5*(b + sgn(b)*std::sqrt(std::max(0.0, b*b - 4.0*c)));
    CHECK(q != 0.0);
    return GeomVector<2>(q, c/q); // (q + 1.0e-50*sgn(q)));
  }
}

//----------------------------------------------------------------------
// This 3-D version is based on the ideas from David Eberly at
// www.geometrictools.com
//----------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<3>
GeomSymmetricTensor<3>::eigenValues() const {
  const auto fscale = std::max(10.0*std::numeric_limits<double>::epsilon(), this->maxAbsElement()); 
  CHECK(fscale > 0.0);
  const auto fscalei = safeInvVar(fscale);
  const auto a00 = this->mxx*fscalei;
  const auto a01 = this->mxy*fscalei;
  const auto a02 = this->mxz*fscalei;
  const auto a11 = this->myy*fscalei;
  const auto a12 = this->myz*fscalei;
  const auto a22 = this->mzz*fscalei;
  const auto c0 = a00*a11*a22 + 2.0*a01*a02*a12 - a00*a12*a12 - a11*a02*a02 - a22*a01*a01;
  const auto c1 = a00*a11 - a01*a01 + a00*a22 - a02*a02 + a11*a22 - a12*a12;
  const auto c2 = a00 + a11 + a22;
  const auto c2Div3 = c2*onethird();
  const auto aDiv3 = std::min(0.0, onethird()*(c1 - c2*c2Div3));
  const auto mbDiv2 = 0.5*(c0 + c2Div3*(2.0*c2Div3*c2Div3 - c1));
  const auto q = std::min(0.0, mbDiv2*mbDiv2 + aDiv3*aDiv3*aDiv3);
  CHECK(-aDiv3 >= 0.0);
  CHECK(-q >= 0.0);
  const auto mag = std::sqrt(-aDiv3);
  const auto angle = atan2(std::sqrt(-q), mbDiv2)*onethird();
  const auto cs = cos(angle);
  const auto sn = sin(angle);
  return GeomVector<3>(fscale*(c2Div3 + 2.0*mag*cs),
                       fscale*(c2Div3 - mag*(cs + sqrt3()*sn)),
                       fscale*(c2Div3 - mag*(cs - sqrt3()*sn)));
}

//------------------------------------------------------------------------------
// Return the eigen values and eigen vectors of a symmetric tensor
//------------------------------------------------------------------------------
// 1-D.
template<>
SPHERAL_HOST_DEVICE
inline
EigenStruct<1>
GeomSymmetricTensor<1>::eigenVectors() const {
  EigenStruct<1> result;
  result.eigenValues.x(mxx);
  result.eigenVectors = GeomTensor<1>(1.0, 1.0, 1.0);
  return result;
}

//------------------------------------------------------------------------------
// 2-D.
template<>
SPHERAL_HOST_DEVICE
inline
EigenStruct<2>
GeomSymmetricTensor<2>::eigenVectors() const {
  const auto fscale = std::max(10.0*std::numeric_limits<double>::epsilon(), this->maxAbsElement()); 
  CHECK(fscale > 0.0);
  const auto fscalei = safeInvVar(fscale);
  const auto axx = this->mxx*fscalei;
  const auto axy = this->mxy*fscalei;
  const auto ayy = this->myy*fscalei;
  EigenStruct<2> result;
  if (std::abs(axy) < 1.0e-50) {
    result.eigenValues = diagonalElements();
    result.eigenVectors = one();
  } else {
    const auto theta = 0.5*atan2(2.0*axy, ayy - axx);
    const auto xhat = cos(theta);
    const auto yhat = sin(theta);
    result.eigenValues.x(xhat*(axx*xhat - axy*yhat) - yhat*(axy*xhat - ayy*yhat));
    result.eigenValues.y(yhat*(axx*yhat + axy*xhat) + xhat*(axy*yhat + ayy*xhat));
    result.eigenValues *= fscale;
    result.eigenVectors = GeomTensor<2>( xhat, yhat,
                                        -yhat, xhat,
                                         1.0);
  }
  return result;
}

//------------------------------------------------------------------------------
// 3-D.
template<>
SPHERAL_HOST_DEVICE
inline
EigenStruct<3>
GeomSymmetricTensor<3>::eigenVectors() const {

  // Some useful typedefs.
  using Vector =GeomVector<3>;
  using Tensor = GeomTensor<3>;
  using SymTensor = GeomSymmetricTensor<3>;

  // Tolerances for fuzzy math.
  const auto degenerate = 1.0e-20;
  const auto tolerance = 5.0e-5;

  // Prepare the result.
  EigenStruct<3> result;

  // Create a scaled version of this tensor, with all elements in the range [-1,1].
  const auto fscale = std::max(10.0*std::numeric_limits<double>::epsilon(), this->maxAbsElement());
  CHECK(fscale > 0.0);
  const auto fscalei = safeInvVar(fscale);
  SymTensor A = (*this)*fscalei;

  // Check for any degenerate elements, and just zero 'em out.
  A.xx(std::abs(A.xx()) < degenerate ? 0.0 : A.xx());
  A.xy(std::abs(A.xy()) < degenerate ? 0.0 : A.xy());
  A.xz(std::abs(A.xz()) < degenerate ? 0.0 : A.xz());
  A.yy(std::abs(A.yy()) < degenerate ? 0.0 : A.yy());
  A.yz(std::abs(A.yz()) < degenerate ? 0.0 : A.yz());
  A.zz(std::abs(A.zz()) < degenerate ? 0.0 : A.zz());

// #ifdef USEJACOBI

//   // Use the Jacobi iterative diagonalization method to determine
//   // the eigen values/vectors.
//   const int nrot = jacobiDiagonalize<Dim<3> >(A,
//                                               result.eigenVectors,
//                                               result.eigenValues);
//   result.eigenValues *= fscale;

// #elif USEEIGEN

  // Use the Eigen library to determine the eigen values/vectors.
  {
    Eigen::Matrix3d B;
    B << 
      A.xx(), A.xy(), A.xz(),
      A.yx(), A.yy(), A.yz(),
      A.zx(), A.zy(), A.zz();
    const Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(B);
    const Eigen::Vector3d& Bvals = eigensolver.eigenvalues();
    const Eigen::Matrix3d& Bvecs = eigensolver.eigenvectors();
    result.eigenValues = Vector(Bvals(0), Bvals(1), Bvals(2)) * fscale;
    const auto x1 = 1.0/std::sqrt(Bvecs(0,0)*Bvecs(0,0) + Bvecs(1,0)*Bvecs(1,0) + Bvecs(2,0)*Bvecs(2,0));
    const auto x2 = 1.0/std::sqrt(Bvecs(0,1)*Bvecs(0,1) + Bvecs(1,1)*Bvecs(1,1) + Bvecs(2,1)*Bvecs(2,1));
    const auto x3 = 1.0/std::sqrt(Bvecs(0,2)*Bvecs(0,2) + Bvecs(1,2)*Bvecs(1,2) + Bvecs(2,2)*Bvecs(2,2));
    result.eigenVectors = Tensor(Bvecs(0,0)*x1, Bvecs(0,1)*x2, Bvecs(0,2)*x3,
                                 Bvecs(1,0)*x1, Bvecs(1,1)*x2, Bvecs(1,2)*x3,
                                 Bvecs(2,0)*x1, Bvecs(2,1)*x2, Bvecs(2,2)*x3);
  }

// #else

//   // Compute the scaled eigen-values, and sort them.
//   Vector lambdaVec = A.eigenValues();
//   sort(lambdaVec.begin(), lambdaVec.end());
//   CHECK(lambdaVec.x() <= lambdaVec.y() and
//         lambdaVec.y() <= lambdaVec.z());

//   // Assign the true eigen-values in the result.
//   result.eigenValues = fscale*lambdaVec;
//   result.eigenVectors = SymTensor::one();

//   // If any of the eigen-values result in a tensor that is not positive-rank 
//   // (all zero elements), we assume the eigen-values are equal and punt
//   // with the identity tensor for the eigen-vectors.
//   // We simultaneously compute the row containing the maximum absolute value 
//   // element for each eigen-value.
//   bool punt = false;
//   double maxEVelement = -1.0;
//   Vector maxEVrow;
//   int iFirst = -1;
//   for (int ivalue = 0; ivalue != 3; ++ivalue) {
//     const SymTensor M = A - lambdaVec(ivalue)*SymTensor::one();
//     if (M.maxAbsElement() < degenerate) punt = true;
//     for (int irow = 0; irow != 3; ++irow) {
//       const Vector Mvec = M.getRow(irow);
//       const double thpt = Mvec.maxAbsElement();
//       if (thpt > maxEVelement) {
//         maxEVelement = thpt;
//         maxEVrow = Mvec;
//         iFirst = ivalue;
//       }
//     }
//   }

//   // If we found an all zero M (= A - lambda*I) matrix, we punt and accept the identity
//   // tensor as our eigen-vectors.  Otherwise, continue the compuation.
//   if (!punt) {
//     CHECK(iFirst >= 0 and iFirst < 3);

//     // Select the ordering we'll go through the eigen-values in, starting
//     // with the row with the largest absolute value element.
//     const int iSecond = (iFirst + 1) % 3;
//     const int iThird = (iSecond + 1) % 3;
//     CHECK(iFirst + iSecond + iThird == 3);

//     // We need two orthogonal unit vectors in the plane perpendicular to
//     // the maximum row selected previously.  We can do this by finding the
//     // rotational transformation wherein x' axis is aligned with this row, and 
//     // taking our two vectors as the other two rows of this transform.
//     const Vector R = maxEVrow.unitVector();
//     const Tensor Tr = rotationMatrix(R);
//     const Vector U0 = Tr.getRow(1);
//     const Vector U1 = Tr.getRow(2);
    
//     // Now we can compute the eigen-vector corresponding the first eigen-value
//     // selected previously.
//     const Vector V0 = buildUniqueEigenVector(A, 
//                                              lambdaVec(iFirst),
//                                              U0,
//                                              U1);
//     result.eigenVectors.setColumn(iFirst, V0);

//     // Now we know the remaining eigen-vectors are in the plane perpendicular to
//     // V0.  We know R is in that plane, and so is R x V0.  With that knowledge
//     // we can basically repeat the same procedure for the next eigen-vector.
//     Vector S = R.cross(V0);
//     CHECK(fuzzyEqual(S.magnitude2(), 1.0, tolerance));
//     const Vector V1 = buildUniqueEigenVector(A,
//                                              lambdaVec(iSecond),
//                                              R,
//                                              S);
//     result.eigenVectors.setColumn(iSecond, V1);
    
//     // The last eigen-vector is orthogonal to the first two, so we can find it
//     // simply by taking the cross-product of the previous eigen-vectors.
//     const Vector V2 = V0.cross(V1);
//     CHECK(fuzzyEqual(V2.magnitude2(), 1.0, tolerance));
//     CHECK(fuzzyEqual(((A - lambdaVec(iThird)*SymTensor::one())*V2).maxAbsElement(), 0.0, tolerance));
//     result.eigenVectors.setColumn(iThird, V2);
//   }

// #endif

  BEGIN_CONTRACT_SCOPE
  // Check the result.
  const auto lambda1 = result.eigenValues.x();
  const auto lambda2 = result.eigenValues.y();
  const auto lambda3 = result.eigenValues.z();
  CONTRACT_VAR(lambda1);
  CONTRACT_VAR(lambda2);
  CONTRACT_VAR(lambda3);
  const Vector v1 = result.eigenVectors.getColumn(0);
  const Vector v2 = result.eigenVectors.getColumn(1);
  const Vector v3 = result.eigenVectors.getColumn(2);
  CONTRACT_VAR(v1);
  CONTRACT_VAR(v2);
  CONTRACT_VAR(v3);
  ENSURE2(fuzzyEqual(v1.dot(v2), 0.0, tolerance) and 
          fuzzyEqual(v1.dot(v3), 0.0, tolerance) and 
          fuzzyEqual(v2.dot(v3), 0.0, tolerance),
          v1 << " " << v2 << " " << v3 << " : " << *this);
  ENSURE2(fuzzyEqual(v1.magnitude2(), 1.0, tolerance) and
          fuzzyEqual(v2.magnitude2(), 1.0, tolerance) and
          fuzzyEqual(v3.magnitude2(), 1.0, tolerance),
          v1 << " " << v2 << " " << v3);
  const auto tol = tolerance*std::max(1.0, this->maxAbsElement());
  CONTRACT_VAR(tol);
  ENSURE2(fuzzyEqual((SymTensor(xx() - lambda1, xy(), xz(),
                                yx(), yy() - lambda1, yz(),
                                zx(), zy(), zz() - lambda1)*v1).maxAbsElement(), 0.0, tol),
          *this << " " << A << " " << lambda1 << " " << v1 << " " << tol << " "
          << SymTensor(xx() - lambda1, xy(), xz(),
                       yx(), yy() - lambda1, yz(),
                       zx(), zy(), zz() - lambda1)*v1);
  ENSURE(fuzzyEqual((SymTensor(xx() - lambda2, xy(), xz(),
                               yx(), yy() - lambda2, yz(),
                               zx(), zy(), zz() - lambda2)*v2).maxAbsElement(), 0.0, tol));
  ENSURE(fuzzyEqual((SymTensor(xx() - lambda3, xy(), xz(),
                               yx(), yy() - lambda3, yz(),
                               zx(), zy(), zz() - lambda3)*v3).maxAbsElement(), 0.0, tol));
  ENSURE(fuzzyEqual(abs(result.eigenVectors.Determinant()), 1.0, tolerance));
  END_CONTRACT_SCOPE

  return result;
}

//------------------------------------------------------------------------------
// Apply a rotational transform to this tensor.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
void
GeomSymmetricTensor<1>::
rotationalTransform(const GeomTensor<1>& R) {
  CONTRACT_VAR(R);
  REQUIRE2(fuzzyEqual(std::abs(R.Determinant()), 1.0, 1.0e-5), R);
}

template<>
SPHERAL_HOST_DEVICE
inline
void
GeomSymmetricTensor<2>::
rotationalTransform(const GeomTensor<2>& R) {
  REQUIRE2(fuzzyEqual(std::abs(R.Determinant()), 1.0, 1.0e-5), R);

  const auto A0 = this->mxx;
  const auto A1 = this->mxy;
  const auto A2 = this->myy;

  const auto R0 = R.xx();
  const auto R1 = R.xy();
  const auto R2 = R.yx();
  const auto R3 = R.yy();

  const auto T0 = A0*R0 + A1*R1;
  const auto T1 = A1*R0 + A2*R1;

  this->mxx = R0*T0 + R1*T1;
  this->mxy = R2*T0 + R3*T1;
  this->myy = R2*(A0*R2 + A1*R3) + R3*(A1*R2 + A2*R3);
}

template<>
SPHERAL_HOST_DEVICE
inline
void
GeomSymmetricTensor<3>::
rotationalTransform(const GeomTensor<3>& R) {
  REQUIRE2(fuzzyEqual(std::abs(R.Determinant()), 1.0, 1.0e-5), R);

  const auto A0 = this->mxx;
  const auto A1 = this->mxy;
  const auto A2 = this->mxz;
  const auto A3 = this->myy;
  const auto A4 = this->myz;
  const auto A5 = this->mzz;

  const auto R0 = R.xx();
  const auto R1 = R.xy();
  const auto R2 = R.xz();
  const auto R3 = R.yx();
  const auto R4 = R.yy();
  const auto R5 = R.yz();
  const auto R6 = R.zx();
  const auto R7 = R.zy();
  const auto R8 = R.zz();

  const auto T0 = A0*R0 + A1*R1 + A2*R2;
  const auto T1 = A0*R3 + A1*R4 + A2*R5;
  const auto T2 = A1*R0 + A3*R1 + A4*R2;
  const auto T3 = A1*R3 + A3*R4 + A4*R5;
  const auto T4 = A2*R0 + A4*R1 + A5*R2;
  const auto T5 = A2*R3 + A4*R4 + A5*R5;

  this->mxx = R0*T0 + R1*T2 + R2*T4;
  this->mxy = R3*T0 + R4*T2 + R5*T4;
  this->mxz = R6*T0 + R7*T2 + R8*T4;
  this->myy = R3*T1 + R4*T3 + R5*T5;
  this->myz = R6*T1 + R7*T3 + R8*T5;
  this->mzz = R6*(A0*R6 + A1*R7 + A2*R8) + R7*(A1*R6 + A3*R7 + A4*R8) + R8*(A2*R6 + A4*R7 + A5*R8);
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
void
GeomSymmetricTensor<nDim>::
rotationalTransform(const GeomTensor<3>& R) requires (nDim != 3) {
  this->rotationalTransform(TensorType(R));
}

//------------------------------------------------------------------------------
// Return the maximum absolute value of the elements.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
double
GeomSymmetricTensor<1>::
maxAbsElement() const {
  return std::abs(this->mxx);
}

template<>
SPHERAL_HOST_DEVICE
inline
double
GeomSymmetricTensor<2>::
maxAbsElement() const {
  return std::max({std::abs(this->mxx), std::abs(this->mxy), std::abs(this->myy)});
}

template<>
SPHERAL_HOST_DEVICE
inline
double
GeomSymmetricTensor<3>::
maxAbsElement() const {
  return std::max({std::abs(this->mxx), std::abs(this->mxy), std::abs(this->mxz), 
                   std::abs(this->myy), std::abs(this->myz), std::abs(this->mzz)});
}

//------------------------------------------------------------------------------
// Enforce a minimum eigen value
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<1>
GeomSymmetricTensor<1>::
enforceMinEigenValue(const double& x) const {
  return GeomSymmetricTensor(std::max(this->mxx, x),
                             std::max(this->myy, x),
                             std::max(this->mzz, x));
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<2>
GeomSymmetricTensor<2>::
enforceMinEigenValue(const double& x) const {
  auto ev = this->eigenVectors();
  if (ev.eigenValues.minElement() < x) {
    GeomSymmetricTensor result;
    result.xx(std::max(ev.eigenValues(0), x));
    result.yy(std::max(ev.eigenValues(1), x));
    result.zz(std::max(this->mzz, x));
    result.rotationalTransform(ev.eigenVectors);
    return result;
  } else if (this->mzz < x) {
    GeomSymmetricTensor result(*this);
    result.zz(x);
    return result;
  } else {
    return *this;
  }
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<3>
GeomSymmetricTensor<3>::
enforceMinEigenValue(const double& x) const {
  auto ev = this->eigenVectors();
  if (ev.eigenValues.minElement() < x) {
    GeomSymmetricTensor result;
    result.xx(std::max(ev.eigenValues(0), x));
    result.yy(std::max(ev.eigenValues(1), x));
    result.zz(std::max(ev.eigenValues(2), x));
    result.rotationalTransform(ev.eigenVectors);
    return result;
  } else {
    return *this;
  }
}

//------------------------------------------------------------------------------
// Enforce a maximum eigen value
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<1>
GeomSymmetricTensor<1>::
enforceMaxEigenValue(const double& x) const {
  return GeomSymmetricTensor(std::min(this->mxx, x),
                             std::min(this->myy, x),
                             std::min(this->mzz, x));
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<2>
GeomSymmetricTensor<2>::
enforceMaxEigenValue(const double& x) const {
  auto ev = this->eigenVectors();
  if (ev.eigenValues.maxElement() > x) {
    GeomSymmetricTensor result;
    result.xx(std::min(ev.eigenValues(0), x));
    result.yy(std::min(ev.eigenValues(1), x));
    result.zz(std::min(this->mzz, x));
    result.rotationalTransform(ev.eigenVectors);
    return result;
  } else if (this->mzz > x) {
    GeomSymmetricTensor result(*this);
    result.zz(x);
    return result;
  } else {
    return *this;
  }
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<3>
GeomSymmetricTensor<3>::
enforceMaxEigenValue(const double& x) const {
  auto ev = this->eigenVectors();
  if (ev.eigenValues.maxElement() > x) {
    GeomSymmetricTensor result;
    result.xx(std::min(ev.eigenValues(0), x));
    result.yy(std::min(ev.eigenValues(1), x));
    result.zz(std::min(ev.eigenValues(2), x));
    result.rotationalTransform(ev.eigenVectors);
    return result;
  } else {
    return *this;
  }
}

//------------------------------------------------------------------------------
// Clamp the min/max eigen values
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<1>
GeomSymmetricTensor<1>::
clampEigenValues(const double& minvalue,
                 const double& maxvalue) const {
  return GeomSymmetricTensor(std::clamp(this->mxx, minvalue, maxvalue),
                             std::clamp(this->myy, minvalue, maxvalue),
                             std::clamp(this->mzz, minvalue, maxvalue));
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<2>
GeomSymmetricTensor<2>::
clampEigenValues(const double& minvalue,
                 const double& maxvalue) const {
  auto ev = this->eigenVectors();
  if (ev.eigenValues.minElement() < minvalue or
      ev.eigenValues.maxElement() > maxvalue) {
    GeomSymmetricTensor result;
    result.xx(std::clamp(ev.eigenValues(0), minvalue, maxvalue));
    result.yy(std::clamp(ev.eigenValues(1), minvalue, maxvalue));
    result.zz(std::clamp(this->mzz, minvalue, maxvalue));
    result.rotationalTransform(ev.eigenVectors);
    return result;
  } else if (this->mzz < minvalue or
             this->mzz > maxvalue) {
    GeomSymmetricTensor result(*this);
    result.zz(std::clamp(this->mzz, minvalue, maxvalue));
    return result;
  } else {
    return *this;
  }
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<3>
GeomSymmetricTensor<3>::
clampEigenValues(const double& minvalue,
                 const double& maxvalue) const {
  auto ev = this->eigenVectors();
  if (ev.eigenValues.minElement() < minvalue or
      ev.eigenValues.maxElement() > maxvalue) {
    GeomSymmetricTensor result;
    result.xx(std::clamp(ev.eigenValues(0), minvalue, maxvalue));
    result.yy(std::clamp(ev.eigenValues(1), minvalue, maxvalue));
    result.zz(std::clamp(ev.eigenValues(2), minvalue, maxvalue));
    result.rotationalTransform(ev.eigenVectors);
    return result;
  } else {
    return *this;
  }
}

//------------------------------------------------------------------------------
// diagonal (strictly 3D)
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
GeomVector<3>
GeomSymmetricTensor<nDim>::diagonalElements3D() const {
  return GeomVector<3>(this->mxx, this->myy, this->mzz);
}

//------------------------------------------------------------------------------
// trace (strictly 3D)
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
double
GeomSymmetricTensor<nDim>::Trace3D() const {
  return this->mxx + this->myy + this->mzz;
}

//------------------------------------------------------------------------------
// determinant (strictly 3D)
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
double
GeomSymmetricTensor<nDim>::Determinant3D() const {
  return (xx()*yy()*zz() + xy()*yz()*zx() + xz()*yx()*zy() -
          xx()*yz()*zy() - xy()*yx()*zz() - xz()*yy()*zx());
}

//------------------------------------------------------------------------------
// doubledot product (a_ij * b_ji -- strictly 3D)
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
double
GeomSymmetricTensor<nDim>::
doubledot3D(const GeomTensor<nDim>& rhs) const {
  return (xx()*rhs.xx() + xy()*rhs.yx() + xz()*rhs.zx() +
          yx()*rhs.xy() + yy()*rhs.yy() + yz()*rhs.zy() +
          zx()*rhs.xz() + zy()*rhs.yz() + zz()*rhs.zz());
}

//------------------------------------------------------------------------------
// doubledot product with a symmetric tensor (a_ij * b_ji -- strictly 3D)
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
double
GeomSymmetricTensor<nDim>::
doubledot3D(const GeomSymmetricTensor<nDim>& rhs) const {
  return (xx()*rhs.xx() + xy()*rhs.yx() + xz()*rhs.zx() +
          yx()*rhs.xy() + yy()*rhs.yy() + yz()*rhs.zy() +
          zx()*rhs.xz() + zy()*rhs.yz() + zz()*rhs.zz());
}

//------------------------------------------------------------------------------
// Return the doubledot product with ourself -- strictly 3D
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
double
GeomSymmetricTensor<nDim>::
selfDoubledot3D() const {
  return (xx()*xx() + 2.0*xy()*yx() + 2.0*xz()*zx() +
          yy()*yy() + 2.0*yz()*zy() + zz()*zz());
}

//------------------------------------------------------------------------------
// Return the maximum absolute value of the elements (3D)
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
double
GeomSymmetricTensor<1>::
maxAbsElement3D() const {
  return std::max({std::abs(this->mxx), std::abs(this->myy), std::abs(this->mzz)});
}

template<>
SPHERAL_HOST_DEVICE
inline
double
GeomSymmetricTensor<2>::
maxAbsElement3D() const {
  return std::max({std::abs(this->mxx), std::abs(this->mxy), std::abs(this->myy),
                   std::abs(this->mzz)});
}

template<>
SPHERAL_HOST_DEVICE
inline
double
GeomSymmetricTensor<3>::
maxAbsElement3D() const {
  return maxAbsElement();
}

//------------------------------------------------------------------------------
// Find the eigenvalues (3D)
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<3>
GeomSymmetricTensor<1>::eigenValues3D() const {
  return diagonalElements3D();
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<3>
GeomSymmetricTensor<2>::eigenValues3D() const {
  const auto ev2D = eigenValues();
  return GeomVector<3>(ev2D.x(), ev2D.y(), this->mzz);
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<3>
GeomSymmetricTensor<3>::eigenValues3D() const {
  return eigenValues();
}

//------------------------------------------------------------------------------
// Find the eigenvectors (3D)
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
EigenStruct<3>
GeomSymmetricTensor<1>::eigenVectors3D() const {
  EigenStruct<3> result;
  result.eigenValues = eigenValues3D();
  result.eigenVectors = GeomTensor<3>::one();
  return result;
}

template<>
SPHERAL_HOST_DEVICE
inline
EigenStruct<3>
GeomSymmetricTensor<2>::eigenVectors3D() const {
  EigenStruct<3> result;
  const auto ev2D = eigenVectors();
  result.eigenValues.x(ev2D.eigenValues.x());
  result.eigenValues.y(ev2D.eigenValues.y());
  result.eigenValues.z(this->mzz);
  result.eigenVectors(0,0) = ev2D.eigenVectors(0,0);
  result.eigenVectors(1,0) = ev2D.eigenVectors(1,0);
  result.eigenVectors(0,1) = ev2D.eigenVectors(0,1);
  result.eigenVectors(1,1) = ev2D.eigenVectors(1,1);
  result.eigenVectors(2,2) = 1.0;
  return result;
}

template<>
SPHERAL_HOST_DEVICE
inline
EigenStruct<3>
GeomSymmetricTensor<3>::eigenVectors3D() const {
  return eigenVectors();
}

//------------------------------------------------------------------------------
// Multiply a tensor with a vector (strictly 3D)
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
GeomVector<3>
GeomSymmetricTensor<nDim>::dot(const GeomVector<3>& rhs) const requires (nDim != 3) {
  return GeomVector<3>(xx()*rhs.x() + xy()*rhs.y() + xz()*rhs.z(),
                       yx()*rhs.x() + yy()*rhs.y() + yz()*rhs.z(),
                       zx()*rhs.x() + zy()*rhs.y() + zz()*rhs.z());
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
GeomVector<3>
GeomSymmetricTensor<nDim>::operator*(const GeomVector<3>& rhs) const requires (nDim != 3) {
  return dot(rhs);
}

//------------------------------------------------------------------------------
// Addition with a 3D tensor
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<1>&
GeomSymmetricTensor<1>::operator+=(const GeomSymmetricTensor<3>& rhs) {
  REQUIRE(rhs.xy() == 0.0 and
          rhs.xz() == 0.0 and
          rhs.yz() == 0.0);
  this->mxx += rhs.xx();
  this->myy += rhs.yy();
  this->mzz += rhs.zz();
  return *this;
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<2>&
GeomSymmetricTensor<2>::operator+=(const GeomSymmetricTensor<3>& rhs) {
  REQUIRE(rhs.xz() == 0.0 and
          rhs.yz() == 0.0);
  this->mxx += rhs.xx();
  this->mxy += rhs.xy();
  this->myy += rhs.yy();
  this->mzz += rhs.zz();
  return *this;
}

//------------------------------------------------------------------------------
// Generate an Eigen Tensor.
//------------------------------------------------------------------------------
template<>
inline
GeomSymmetricTensor<1>::EigenType
GeomSymmetricTensor<1>::eigen() const {
  return EigenType(this->mxx);
}

template<>
inline
GeomSymmetricTensor<2>::EigenType
GeomSymmetricTensor<2>::eigen() const {
  EigenType result;
  result << this->mxx, this->mxy,
            this->mxy, this->myy;
  return result;
}

template<>
inline
GeomSymmetricTensor<3>::EigenType
GeomSymmetricTensor<3>::eigen() const {
  EigenType result;
  result << this->mxx, this->mxy, this->mxz,
            this->mxy, this->myy, this->myz,
            this->mxz, this->myz, this->mzz;
  return result;
}

//********************************************************************************
// Global Functions.
//********************************************************************************

//------------------------------------------------------------------------------
// Multiplication by a scalar
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
Spheral::GeomSymmetricTensor<nDim>
operator*(double lhs, const Spheral::GeomSymmetricTensor<nDim>& rhs) {
  return rhs*lhs;
}

//------------------------------------------------------------------------------
// Input (istream) operator.
//------------------------------------------------------------------------------
template<int nDim>
inline
std::istream&
operator>>(std::istream& is, Spheral::GeomSymmetricTensor<nDim>& ten) {
  std::string parenthesis;
  is >> parenthesis;
  for (auto itr = ten.begin(); itr < ten.end(); ++itr) {
    is >> *itr;
  }
  is >> parenthesis;
  return is;
}

//------------------------------------------------------------------------------
// Output (ostream) operator.
//------------------------------------------------------------------------------
template<int nDim>
inline
std::ostream&
operator<<(std::ostream& os, const Spheral::GeomSymmetricTensor<nDim>& ten) {
  os << "( ";
  for (auto itr = ten.begin(); itr < ten.end(); ++itr) {
    os << *itr << " ";
  }
  os << ")";
  return os;
}

//------------------------------------------------------------------------------
// Comparison with doubles as first argument
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
bool
operator<(const double& lhs, const GeomSymmetricTensor<nDim>& rhs) {
  return rhs > lhs;
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
bool
operator>(const double& lhs, const GeomSymmetricTensor<nDim>& rhs) {
  return rhs < lhs;
}

//------------------------------------------------------------------------------
// Symmetric tensor specializations for min/max
//------------------------------------------------------------------------------
template<int ndim>
SPHERAL_HOST_DEVICE
GeomSymmetricTensor<ndim>
min(const double& lhs, const GeomSymmetricTensor<ndim>& rhs) {
  return rhs.enforceMaxEigenValue(lhs);
}

template<int ndim>
SPHERAL_HOST_DEVICE
GeomSymmetricTensor<ndim>
min(const GeomSymmetricTensor<ndim>& lhs, const double& rhs) {
  return lhs.enforceMaxEigenValue(rhs);
}

template<int ndim>
SPHERAL_HOST_DEVICE
GeomSymmetricTensor<ndim>
max(const double& lhs, const GeomSymmetricTensor<ndim>& rhs) {
  return rhs.enforceMinEigenValue(lhs);
}

template<int ndim>
SPHERAL_HOST_DEVICE
GeomSymmetricTensor<ndim>
max(const GeomSymmetricTensor<ndim>& lhs, const double& rhs) {
  return lhs.enforceMinEigenValue(rhs);
}

}
