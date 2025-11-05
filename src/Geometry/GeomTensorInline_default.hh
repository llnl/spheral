#include <algorithm>
#include <limits.h>
#include <string>

#include "GeomVector.hh"
#include "GeomSymmetricTensor.hh"
#include "Utilities/SpheralFunctions.hh"
#include "Utilities/FastMath.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Return the element index corresponding to the given (row,column)
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
typename GeomTensor<nDim>::size_type
GeomTensor<nDim>::elementIndex(const typename GeomTensor<nDim>::size_type row,
                               const typename GeomTensor<nDim>::size_type column) const {
  REQUIRE(row < nDim);
  REQUIRE(column < nDim);
  return row*nDim + column;
}

//------------------------------------------------------------------------------
// Construct with the given values for the elements.
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
GeomTensor<nDim>::
GeomTensor(const double a11):
  GeomTensorBase<nDim>(a11) {
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<2>::
GeomTensor(const double a11, const double a12, 
	   const double a21, const double a22):
  GeomTensorBase<2>(a11, a12,
                    a21, a22) {
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<3>::
GeomTensor(const double a11, const double a12, const double a13,
	   const double a21, const double a22, const double a23,
	   const double a31, const double a32, const double a33):
  GeomTensorBase<3>(a11, a12, a13,
                    a21, a22, a23,
                    a31, a32, a33) {
}

//------------------------------------------------------------------------------
// Override the generic constructors to throw if they're called in the wrong 
// dimensions.
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
GeomTensor<nDim>::
GeomTensor(const double /*a11*/, const double /*a12*/,
           const double /*a21*/, const double /*a22*/):
  GeomTensorBase<nDim>(0.0) {
  VERIFY2(false, "GeomTensor(a11, a12, a21, a22): wrong number of dimensions.");
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
GeomTensor<nDim>::
GeomTensor(const double /*a11*/, const double /*a12*/, const double /*a13*/,
           const double /*a21*/, const double /*a22*/, const double /*a23*/,
           const double /*a31*/, const double /*a32*/, const double /*a33*/):
  GeomTensorBase<nDim>(0.0) {
  VERIFY2(false, "GeomTensor(a11, a12, a13, a21, a22, a23, a31, a32, a33): wrong number of dimensions.");
}

//------------------------------------------------------------------------------
// Copy constructors.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<1>::
GeomTensor(const GeomSymmetricTensor<1>& ten):
  GeomTensorBase<1>(ten.xx()) {
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<2>::
GeomTensor(const GeomSymmetricTensor<2>& ten):
  GeomTensorBase<2>(ten.xx(), ten.xy(),
                    ten.yx(), ten.yy()) {
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<3>::
GeomTensor(const GeomSymmetricTensor<3>& ten):
  GeomTensorBase<3>(ten.xx(), ten.xy(), ten.xz(),
                    ten.yx(), ten.yy(), ten.yz(),
                    ten.zx(), ten.zy(), ten.zz()) {
}

//------------------------------------------------------------------------------
// Construct from an Eigen Tensor.
//------------------------------------------------------------------------------
template<>
template<typename Derived>
inline
GeomTensor<1>::GeomTensor(const Eigen::MatrixBase<Derived>& ten):
  GeomTensorBase<1>(ten(0,0)) {
}

template<>
template<typename Derived>
inline
GeomTensor<2>::GeomTensor(const Eigen::MatrixBase<Derived>& ten):
  GeomTensorBase<2>(ten(0,0), ten(0,1),
                    ten(1,0), ten(1,1)) {
}

template<>
template<typename Derived>
inline
GeomTensor<3>::GeomTensor(const Eigen::MatrixBase<Derived>& ten):
  GeomTensorBase<3>(ten(0,0), ten(0,1), ten(0,2),
                    ten(1,0), ten(1,1), ten(1,2),
                    ten(2,0), ten(2,1), ten(2,2)) {
}

//------------------------------------------------------------------------------
// Assignment operators.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<1>&
GeomTensor<1>::
operator=(const GeomSymmetricTensor<1>& ten) {
  this->mxx = ten.xx();
  return *this;
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<2>&
GeomTensor<2>::
operator=(const GeomSymmetricTensor<2>& ten) {
  this->mxx = ten.xx();
  this->mxy = ten.xy();
  this->myx = ten.yx();
  this->myy = ten.yy();
  return *this;
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<3>&
GeomTensor<3>::
operator=(const GeomSymmetricTensor<3>& ten) {
  this->mxx = ten.xx();
  this->mxy = ten.xy();
  this->mxz = ten.xz();
  this->myx = ten.yx();
  this->myy = ten.yy();
  this->myz = ten.yz();
  this->mzx = ten.zx();
  this->mzy = ten.zy();
  this->mzz = ten.zz();
  return *this;
}

//------------------------------------------------------------------------------
// The assignment operator (Eigen Tensor).
//------------------------------------------------------------------------------
template<>
template<typename Derived>
inline
GeomTensor<1>&
GeomTensor<1>::operator=(const Eigen::MatrixBase<Derived>& ten) {
  this->mxx = ten(0,0);
  return *this;
}

template<>
template<typename Derived>
inline
GeomTensor<2>&
GeomTensor<2>::operator=(const Eigen::MatrixBase<Derived>& ten) {
  this->mxx = ten(0,0);
  this->mxy = ten(0,1);
  this->myx = ten(1,0);
  this->myy = ten(1,1);
  return *this;
}

template<>
template<typename Derived>
inline
GeomTensor<3>&
GeomTensor<3>::operator=(const Eigen::MatrixBase<Derived>& ten) {
  this->mxx = ten(0,0);
  this->mxy = ten(0,1);
  this->mxz = ten(0,2);
  this->myx = ten(1,0);
  this->myy = ten(1,1);
  this->myz = ten(1,2);
  this->mzx = ten(2,0);
  this->mzy = ten(2,1);
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
GeomTensor<nDim>::operator()(const typename GeomTensor<nDim>::size_type row,
                             const typename GeomTensor<nDim>::size_type column) const {
  REQUIRE(row < nDim);
  REQUIRE(column < nDim);
  return *(begin() + elementIndex(row, column));
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
double&
GeomTensor<nDim>::operator()(const typename GeomTensor<nDim>::size_type row,
                             const typename GeomTensor<nDim>::size_type column) {
  REQUIRE(row < nDim);
  REQUIRE(column < nDim);
  return *(begin() + elementIndex(row, column));
}

//------------------------------------------------------------------------------
// Return the (index) element using the bracket operator.
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
double
GeomTensor<nDim>::operator[](typename GeomTensor<nDim>::size_type index) const {
  REQUIRE(index < numElements);
  return *(begin() + index);
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
double&
GeomTensor<nDim>::operator[](typename GeomTensor<nDim>::size_type index) {
  REQUIRE(index < numElements);
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
GeomTensor<nDim>::xx() const {
  return this->mxx;
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
double
GeomTensor<nDim>::xy() const {
  return this->mxy;
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
double
GeomTensor<nDim>::xz() const {
  return this->mxz;
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
double
GeomTensor<nDim>::yx() const {
  return this->myx;
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
double
GeomTensor<nDim>::yy() const {
  return this->myy;
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
double
GeomTensor<nDim>::yz() const {
  return this->myz;
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
double
GeomTensor<nDim>::zx() const {
  return this->mzx;
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
double
GeomTensor<nDim>::zy() const {
  return this->mzy;
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
double
GeomTensor<nDim>::zz() const {
  return this->mzz;
}

//------------------------------------------------------------------------------
// 1D dummy elements
template<> SPHERAL_HOST_DEVICE inline double GeomTensor<1>::xy() const { return 0.0; }
template<> SPHERAL_HOST_DEVICE inline double GeomTensor<1>::xz() const { return 0.0; }
template<> SPHERAL_HOST_DEVICE inline double GeomTensor<1>::yx() const { return 0.0; }
template<> SPHERAL_HOST_DEVICE inline double GeomTensor<1>::yy() const { return 0.0; }
template<> SPHERAL_HOST_DEVICE inline double GeomTensor<1>::yz() const { return 0.0; }
template<> SPHERAL_HOST_DEVICE inline double GeomTensor<1>::zx() const { return 0.0; }
template<> SPHERAL_HOST_DEVICE inline double GeomTensor<1>::zy() const { return 0.0; }
template<> SPHERAL_HOST_DEVICE inline double GeomTensor<1>::zz() const { return 0.0; }

//------------------------------------------------------------------------------
// 2D dummy elements
template<> SPHERAL_HOST_DEVICE inline double GeomTensor<2>::xz() const { return 0.0; }
template<> SPHERAL_HOST_DEVICE inline double GeomTensor<2>::yz() const { return 0.0; }
template<> SPHERAL_HOST_DEVICE inline double GeomTensor<2>::zx() const { return 0.0; }
template<> SPHERAL_HOST_DEVICE inline double GeomTensor<2>::zy() const { return 0.0; }
template<> SPHERAL_HOST_DEVICE inline double GeomTensor<2>::zz() const { return 0.0; }

//------------------------------------------------------------------------------
// Set the individual elements, as above.
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
void
GeomTensor<nDim>::xx(const double val) {
  this->mxx = val;
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
void
GeomTensor<nDim>::xy(const double val) {
  this->mxy = val;
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
void
GeomTensor<nDim>::xz(const double val) {
  this->mxz = val;
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
void
GeomTensor<nDim>::yx(const double val) {
  this->myx = val;
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
void
GeomTensor<nDim>::yy(const double val) {
  this->myy = val;
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
void
GeomTensor<nDim>::yz(const double val) {
  this->myz = val;
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
void
GeomTensor<nDim>::zx(const double val) {
  this->mzx = val;
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
void
GeomTensor<nDim>::zy(const double val) {
  this->mzy = val;
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
void
GeomTensor<nDim>::zz(double val) {
  this->mzz = val;
}

//------------------------------------------------------------------------------
// 1D dummy elements
template<> SPHERAL_HOST_DEVICE inline void GeomTensor<1>::xy(const double /*val*/) {}
template<> SPHERAL_HOST_DEVICE inline void GeomTensor<1>::xz(const double /*val*/) {}
template<> SPHERAL_HOST_DEVICE inline void GeomTensor<1>::yx(const double /*val*/) {}
template<> SPHERAL_HOST_DEVICE inline void GeomTensor<1>::yy(const double /*val*/) {}
template<> SPHERAL_HOST_DEVICE inline void GeomTensor<1>::yz(const double /*val*/) {}
template<> SPHERAL_HOST_DEVICE inline void GeomTensor<1>::zx(const double /*val*/) {}
template<> SPHERAL_HOST_DEVICE inline void GeomTensor<1>::zy(const double /*val*/) {}
template<> SPHERAL_HOST_DEVICE inline void GeomTensor<1>::zz(const double /*val*/) {}

//------------------------------------------------------------------------------
// 2D dummy elements
template<> SPHERAL_HOST_DEVICE inline void GeomTensor<2>::xz(const double /*val*/) {}
template<> SPHERAL_HOST_DEVICE inline void GeomTensor<2>::yz(const double /*val*/) {}
template<> SPHERAL_HOST_DEVICE inline void GeomTensor<2>::zx(const double /*val*/) {}
template<> SPHERAL_HOST_DEVICE inline void GeomTensor<2>::zy(const double /*val*/) {}
template<> SPHERAL_HOST_DEVICE inline void GeomTensor<2>::zz(const double /*val*/) {}

//------------------------------------------------------------------------------
// Access the individual rows of the GeomTensor.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<1>
GeomTensor<1>::getRow(const GeomTensor<1>::size_type index) const {
  REQUIRE(index < 1);
  return GeomVector<1>((*this)(index, 0));
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<2>
GeomTensor<2>::getRow(const GeomTensor<2>::size_type index) const {
  REQUIRE(index < 2);
  return GeomVector<2>((*this)(index, 0), (*this)(index, 1));
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<3>
GeomTensor<3>::getRow(const GeomTensor<3>::size_type index) const {
  REQUIRE(index < 3);
  return GeomVector<3>((*this)(index, 0), (*this)(index, 1), (*this)(index, 2));
}

//------------------------------------------------------------------------------
// Access the individual columns of the GeomTensor.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<1>
GeomTensor<1>::getColumn(const GeomTensor<2>::size_type index) const {
  REQUIRE(index < 1);
  return GeomVector<1>((*this)(0, index));
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<2>
GeomTensor<2>::getColumn(const GeomTensor<2>::size_type index) const {
  REQUIRE(index < 2);
  return GeomVector<2>((*this)(0, index), (*this)(1, index));
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<3>
GeomTensor<3>::getColumn(const GeomTensor<3>::size_type index) const {
  REQUIRE(index < 3);
  return GeomVector<3>((*this)(0, index), (*this)(1, index), (*this)(2, index));
}

//------------------------------------------------------------------------------
// Set a row of the GeomTensor to a GeomVector.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
void
GeomTensor<1>::setRow(const GeomTensor<1>::size_type index,
                      const GeomVector<1>& vec) {
  REQUIRE(index < 1);
  (*this)(index, 0) = vec(0);
}

template<>
SPHERAL_HOST_DEVICE
inline
void
GeomTensor<2>::setRow(const GeomTensor<2>::size_type index,
                      const GeomVector<2>& vec) {
  REQUIRE(index < 2);
  (*this)(index, 0) = vec(0);
  (*this)(index, 1) = vec(1);
}

template<>
SPHERAL_HOST_DEVICE
inline
void
GeomTensor<3>::setRow(const GeomTensor<3>::size_type index,
                      const GeomVector<3>& vec) {
  REQUIRE(index < 3);
  (*this)(index, 0) = vec(0);
  (*this)(index, 1) = vec(1);
  (*this)(index, 2) = vec(2);
}

//------------------------------------------------------------------------------
// Set a column of the GeomTensor to a GeomVector.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
void
GeomTensor<1>::setColumn(const GeomTensor<1>::size_type index,
                      const GeomVector<1>& vec) {
  REQUIRE(index < 1);
  (*this)(0, index) = vec.x();
}

template<>
SPHERAL_HOST_DEVICE
inline
void
GeomTensor<2>::setColumn(const GeomTensor<2>::size_type index,
                      const GeomVector<2>& vec) {
  REQUIRE(index < 2);
  (*this)(0, index) = vec.x();
  (*this)(1, index) = vec.y();
}

template<>
SPHERAL_HOST_DEVICE
inline
void
GeomTensor<3>::setColumn(const GeomTensor<3>::size_type index,
                      const GeomVector<3>& vec) {
  REQUIRE(index < 3);
  (*this)(0, index) = vec.x();
  (*this)(1, index) = vec.y();
  (*this)(2, index) = vec.z();
}

//------------------------------------------------------------------------------
// Iterators to the raw data.
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
typename GeomTensor<nDim>::iterator
GeomTensor<nDim>::begin() {
  return &(this->mxx);
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
typename GeomTensor<nDim>::iterator
GeomTensor<nDim>::end() {
  return &(this->mxx) + nDim*nDim;
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
typename GeomTensor<nDim>::const_iterator
GeomTensor<nDim>::begin() const{
  return &(this->mxx);
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
typename GeomTensor<nDim>::const_iterator
GeomTensor<nDim>::end() const {
  return &(this->mxx) + nDim*nDim;
}

//------------------------------------------------------------------------------
// Zero out the GeomTensor.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
void
GeomTensor<1>::Zero() {
  this->mxx = 0.0;
}

template<>
SPHERAL_HOST_DEVICE
inline
void
GeomTensor<2>::Zero() {
  this->mxx = 0.0;
  this->mxy = 0.0;
  this->myx = 0.0;
  this->myy = 0.0;
}

template<>
SPHERAL_HOST_DEVICE
inline
void
GeomTensor<3>::Zero() {
  this->mxx = 0.0;
  this->mxy = 0.0;
  this->mxz = 0.0;
  this->myx = 0.0;
  this->myy = 0.0;
  this->myz = 0.0;
  this->mzx = 0.0;
  this->mzy = 0.0;
  this->mzz = 0.0;
}

//------------------------------------------------------------------------------
// Force the tensor to be the identity tensor.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
void
GeomTensor<1>::Identity() {
  this->mxx = 1.0;
}

template<>
SPHERAL_HOST_DEVICE
inline
void
GeomTensor<2>::Identity() {
  this->mxx = 1.0;
  this->mxy = 0.0;
  this->myx = 0.0;
  this->myy = 1.0;
}

template<>
SPHERAL_HOST_DEVICE
inline
void
GeomTensor<3>::Identity() {
  this->mxx = 1.0;
  this->mxy = 0.0;
  this->mxz = 0.0;
  this->myx = 0.0;
  this->myy = 1.0;
  this->myz = 0.0;
  this->mzx = 0.0;
  this->mzy = 0.0;
  this->mzz = 1.0;
}

//------------------------------------------------------------------------------
// Return the negative of a tensor.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<1>
GeomTensor<1>::operator-() const {
  GeomTensor<1> result;
  result.mxx = -(this->mxx);
  return result;
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<2>
GeomTensor<2>::operator-() const {
  GeomTensor<2> result;
  result.mxx = -(this->mxx);
  result.mxy = -(this->mxy);
  result.myx = -(this->myx);
  result.myy = -(this->myy);
  return result;
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<3>
GeomTensor<3>::operator-() const {
  GeomTensor<3> result;
  result.mxx = -(this->mxx);
  result.mxy = -(this->mxy);
  result.mxz = -(this->mxz);
  result.myx = -(this->myx);
  result.myy = -(this->myy);
  result.myz = -(this->myz);
  result.mzx = -(this->mzx);
  result.mzy = -(this->mzy);
  result.mzz = -(this->mzz);
  return result;
}

//------------------------------------------------------------------------------
// Add two tensors.
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
GeomTensor<nDim>
GeomTensor<nDim>::
operator+(const GeomTensor<nDim>& rhs) const {
  GeomTensor<nDim> result(*this);
  result += rhs;
  return result;
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
GeomTensor<nDim>
GeomTensor<nDim>::
operator+(const GeomSymmetricTensor<nDim>& rhs) const {
  GeomTensor<nDim> result(*this);
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
GeomTensor<nDim>::
operator-(const GeomTensor<nDim>& rhs) const {
  GeomTensor<nDim> result(*this);
  result -= rhs;
  return result;
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
GeomTensor<nDim>
GeomTensor<nDim>::
operator-(const GeomSymmetricTensor<nDim>& rhs) const {
  GeomTensor<nDim> result(*this);
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
GeomTensor<nDim>::
operator*(const GeomTensor<nDim>& rhs) const {
  return this->dot(rhs);
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
GeomTensor<nDim>
GeomTensor<nDim>::
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
GeomTensor<nDim>::operator*(const GeomVector<nDim>& rhs) const {
  return this->dot(rhs);
}

// //------------------------------------------------------------------------------
// // Add a scalar to a tensor.
// //------------------------------------------------------------------------------
// template<>
// SPHERAL_HOST_DEVICE
// inline
// GeomTensor<1>
// GeomTensor<1>::operator+(const double rhs) const {
//   return GeomTensor<1>((this->mxx) + rhs);
// }

// template<>
// SPHERAL_HOST_DEVICE
// inline
// GeomTensor<2>
// GeomTensor<2>::operator+(const double rhs) const {
//   return GeomTensor<2>(this->mxx + rhs, this->mxy + rhs, 
//                        this->myx + rhs, this->myy+ rhs);
// }

// template<>
// SPHERAL_HOST_DEVICE
// inline
// GeomTensor<3>
// GeomTensor<3>::operator+(const double rhs) const {
//   return GeomTensor<3>(this->mxx + rhs, this->mxy + rhs, this->mxz + rhs, 
//                        this->myx + rhs, this->myy + rhs, this->myz + rhs,
//                        this->mzx + rhs, this->mzy + rhs, this->mzz + rhs);
// }

// //------------------------------------------------------------------------------
// // Subtract a scalar from a tensor.
// //------------------------------------------------------------------------------
// template<>
// SPHERAL_HOST_DEVICE
// inline
// GeomTensor<1>
// GeomTensor<1>::operator-(const double rhs) const {
//   return GeomTensor<1>(this->mxx - rhs);
// }

// template<>
// SPHERAL_HOST_DEVICE
// inline
// GeomTensor<2>
// GeomTensor<2>::operator-(const double rhs) const {
//   return GeomTensor<2>(this->mxx - rhs, this->mxy - rhs, 
//                        this->myx - rhs, this->myy- rhs);
// }

// template<>
// SPHERAL_HOST_DEVICE
// inline
// GeomTensor<3>
// GeomTensor<3>::operator-(const double rhs) const {
//   return GeomTensor<3>(this->mxx - rhs, this->mxy - rhs, this->mxz - rhs, 
//                        this->myx - rhs, this->myy - rhs, this->myz - rhs,
//                        this->mzx - rhs, this->mzy - rhs, this->mzz - rhs);
// }


//------------------------------------------------------------------------------
// Multiply a tensor by a scalar
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<1>
GeomTensor<1>::operator*(const double rhs) const {
  return GeomTensor<1>(this->mxx * rhs);
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<2>
GeomTensor<2>::operator*(const double rhs) const {
  return GeomTensor<2>(this->mxx * rhs, this->mxy * rhs, 
                       this->myx * rhs, this->myy* rhs);
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<3>
GeomTensor<3>::operator*(const double rhs) const {
  return GeomTensor<3>(this->mxx * rhs, this->mxy * rhs, this->mxz * rhs, 
                       this->myx * rhs, this->myy * rhs, this->myz * rhs,
                       this->mzx * rhs, this->mzy * rhs, this->mzz * rhs);
}


//------------------------------------------------------------------------------
// Divide a tensor by a scalar
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<1>
GeomTensor<1>::operator/(const double rhs) const {
  REQUIRE(rhs != 0.0);
  return GeomTensor<1>(this->mxx / rhs);
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<2>
GeomTensor<2>::operator/(const double rhs) const {
  REQUIRE(rhs != 0.0);
  const double rhsInv = 1.0/rhs;
  return GeomTensor<2>(this->mxx * rhsInv, this->mxy * rhsInv, 
                       this->myx * rhsInv, this->myy* rhsInv);
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<3>
GeomTensor<3>::operator/(const double rhs) const {
  REQUIRE(rhs != 0.0);
  const double rhsInv = 1.0/rhs;
  return GeomTensor<3>(this->mxx * rhsInv, this->mxy * rhsInv, this->mxz * rhsInv, 
                       this->myx * rhsInv, this->myy * rhsInv, this->myz * rhsInv,
                       this->mzx * rhsInv, this->mzy * rhsInv, this->mzz * rhsInv);
}


//------------------------------------------------------------------------------
// += tensor
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<1>&
GeomTensor<1>::operator+=(const GeomTensor<1>& rhs) {
  this->mxx += rhs.mxx;
  return *this;
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<2>&
GeomTensor<2>::operator+=(const GeomTensor<2>& rhs) {
  this->mxx += rhs.mxx;
  this->mxy += rhs.mxy;
  this->myx += rhs.myx;
  this->myy += rhs.myy;
  return *this;
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<3>&
GeomTensor<3>::operator+=(const GeomTensor<3>& rhs) {
  this->mxx += rhs.mxx;
  this->mxy += rhs.mxy;
  this->mxz += rhs.mxz;
  this->myx += rhs.myx;
  this->myy += rhs.myy;
  this->myz += rhs.myz;
  this->mzx += rhs.mzx;
  this->mzy += rhs.mzy;
  this->mzz += rhs.mzz;
  return *this;
}

//------------------------------------------------------------------------------
// += symmetric tensor
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<1>&
GeomTensor<1>::operator+=(const GeomSymmetricTensor<1>& rhs) {
  this->mxx += rhs.xx();
  return *this;
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<2>&
GeomTensor<2>::operator+=(const GeomSymmetricTensor<2>& rhs) {
  this->mxx += rhs.xx();
  this->mxy += rhs.xy();
  this->myx += rhs.yx();
  this->myy += rhs.yy();
  return *this;
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<3>&
GeomTensor<3>::operator+=(const GeomSymmetricTensor<3>& rhs) {
  this->mxx += rhs.xx();
  this->mxy += rhs.xy();
  this->mxz += rhs.xz();
  this->myx += rhs.yx();
  this->myy += rhs.yy();
  this->myz += rhs.yz();
  this->mzx += rhs.zx();
  this->mzy += rhs.zy();
  this->mzz += rhs.zz();
  return *this;
}

//------------------------------------------------------------------------------
// += eigen tensor
//------------------------------------------------------------------------------
template<>
template<typename Derived>
inline
GeomTensor<1>&
GeomTensor<1>::operator+=(const Eigen::MatrixBase<Derived>& rhs) {
  this->mxx += rhs(0,0);
  return *this;
}

template<>
template<typename Derived>
inline
GeomTensor<2>&
GeomTensor<2>::operator+=(const Eigen::MatrixBase<Derived>& rhs) {
  this->mxx += rhs(0,0);
  this->mxy += rhs(0,1);
  this->myx += rhs(1,0);
  this->myy += rhs(1,1);
  return *this;
}

template<>
template<typename Derived>
inline
GeomTensor<3>&
GeomTensor<3>::operator+=(const Eigen::MatrixBase<Derived>& rhs) {
  this->mxx += rhs(0,0);
  this->mxy += rhs(0,1);
  this->mxz += rhs(0,2);
  this->myx += rhs(1,0);
  this->myy += rhs(1,1);
  this->myz += rhs(1,2);
  this->mzx += rhs(2,0);
  this->mzy += rhs(2,1);
  this->mzz += rhs(2,2);
  return *this;
}

//------------------------------------------------------------------------------
// -= tensor
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<1>&
GeomTensor<1>::operator-=(const GeomTensor<1>& rhs) {
  this->mxx -= rhs.mxx;
  return *this;
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<2>&
GeomTensor<2>::operator-=(const GeomTensor<2>& rhs) {
  this->mxx -= rhs.mxx;
  this->mxy -= rhs.mxy;
  this->myx -= rhs.myx;
  this->myy -= rhs.myy;
  return *this;
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<3>&
GeomTensor<3>::operator-=(const GeomTensor<3>& rhs) {
  this->mxx -= rhs.mxx;
  this->mxy -= rhs.mxy;
  this->mxz -= rhs.mxz;
  this->myx -= rhs.myx;
  this->myy -= rhs.myy;
  this->myz -= rhs.myz;
  this->mzx -= rhs.mzx;
  this->mzy -= rhs.mzy;
  this->mzz -= rhs.mzz;
  return *this;
}

//------------------------------------------------------------------------------
// -= symmetric tensor
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<1>&
GeomTensor<1>::operator-=(const GeomSymmetricTensor<1>& rhs) {
  this->mxx -= rhs.xx();
  return *this;
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<2>&
GeomTensor<2>::operator-=(const GeomSymmetricTensor<2>& rhs) {
  this->mxx -= rhs.xx();
  this->mxy -= rhs.xy();
  this->myx -= rhs.yx();
  this->myy -= rhs.yy();
  return *this;
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<3>&
GeomTensor<3>::operator-=(const GeomSymmetricTensor<3>& rhs) {
  this->mxx -= rhs.xx();
  this->mxy -= rhs.xy();
  this->mxz -= rhs.xz();
  this->myx -= rhs.yx();
  this->myy -= rhs.yy();
  this->myz -= rhs.yz();
  this->mzx -= rhs.zx();
  this->mzy -= rhs.zy();
  this->mzz -= rhs.zz();
  return *this;
}

//------------------------------------------------------------------------------
// -= eigen tensor
//------------------------------------------------------------------------------
template<>
template<typename Derived>
inline
GeomTensor<1>&
GeomTensor<1>::operator-=(const Eigen::MatrixBase<Derived>& rhs) {
  this->mxx -= rhs(0,0);
  return *this;
}

template<>
template<typename Derived>
inline
GeomTensor<2>&
GeomTensor<2>::operator-=(const Eigen::MatrixBase<Derived>& rhs) {
  this->mxx -= rhs(0,0);
  this->mxy -= rhs(0,1);
  this->myx -= rhs(1,0);
  this->myy -= rhs(1,1);
  return *this;
}

template<>
template<typename Derived>
inline
GeomTensor<3>&
GeomTensor<3>::operator-=(const Eigen::MatrixBase<Derived>& rhs) {
  this->mxx -= rhs(0,0);
  this->mxy -= rhs(0,1);
  this->mxz -= rhs(0,2);
  this->myx -= rhs(1,0);
  this->myy -= rhs(1,1);
  this->myz -= rhs(1,2);
  this->mzx -= rhs(2,0);
  this->mzy -= rhs(2,1);
  this->mzz -= rhs(2,2);
  return *this;
}

//------------------------------------------------------------------------------
// Multiply by a tensor in place.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<1>&
GeomTensor<1>::operator*=(const GeomTensor<1>& rhs) {
  this->mxx *= rhs.mxx;
  return *this;
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<2>&
GeomTensor<2>::operator*=(const GeomTensor<2>& rhs) {
  const double A0 = this->mxx;
  const double A1 = this->mxy;
  const double A2 = this->myx;
  const double A3 = this->myy;
  this->mxx = A0*rhs.mxx + A1*rhs.myx;
  this->mxy = A0*rhs.mxy + A1*rhs.myy;
  this->myx = A2*rhs.mxx + A3*rhs.myx;
  this->myy = A2*rhs.mxy + A3*rhs.myy;
  return *this;
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<3>&
GeomTensor<3>::operator*=(const GeomTensor<3>& rhs) {
  const double A0 = this->mxx;
  const double A1 = this->mxy;
  const double A2 = this->mxz;
  const double A3 = this->myx;
  const double A4 = this->myy;
  const double A5 = this->myz;
  const double A6 = this->mzx;
  const double A7 = this->mzy;
  const double A8 = this->mzz;
  this->mxx = A0*rhs.mxx + A1*rhs.myx + A2*rhs.mzx;
  this->mxy = A0*rhs.mxy + A1*rhs.myy + A2*rhs.mzy;
  this->mxz = A0*rhs.mxz + A1*rhs.myz + A2*rhs.mzz;
  this->myx = A3*rhs.mxx + A4*rhs.myx + A5*rhs.mzx;
  this->myy = A3*rhs.mxy + A4*rhs.myy + A5*rhs.mzy;
  this->myz = A3*rhs.mxz + A4*rhs.myz + A5*rhs.mzz;
  this->mzx = A6*rhs.mxx + A7*rhs.myx + A8*rhs.mzx;
  this->mzy = A6*rhs.mxy + A7*rhs.myy + A8*rhs.mzy;
  this->mzz = A6*rhs.mxz + A7*rhs.myz + A8*rhs.mzz;
  return *this;
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<1>&
GeomTensor<1>::operator*=(const GeomSymmetricTensor<1>& rhs) {
  this->mxx *= rhs.xx();
  return *this;
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<2>&
GeomTensor<2>::operator*=(const GeomSymmetricTensor<2>& rhs) {
  const double A0 = this->mxx;
  const double A1 = this->mxy;
  const double A2 = this->myx;
  const double A3 = this->myy;
  this->mxx = A0*rhs.xx() + A1*rhs.yx();
  this->mxy = A0*rhs.xy() + A1*rhs.yy();
  this->myx = A2*rhs.xx() + A3*rhs.yx();
  this->myy = A2*rhs.xy() + A3*rhs.yy();
  return *this;
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<3>&
GeomTensor<3>::operator*=(const GeomSymmetricTensor<3>& rhs) {
  const double A0 = this->mxx;
  const double A1 = this->mxy;
  const double A2 = this->mxz;
  const double A3 = this->myx;
  const double A4 = this->myy;
  const double A5 = this->myz;
  const double A6 = this->mzx;
  const double A7 = this->mzy;
  const double A8 = this->mzz;
  this->mxx = A0*rhs.xx() + A1*rhs.yx() + A2*rhs.zx();
  this->mxy = A0*rhs.xy() + A1*rhs.yy() + A2*rhs.zy();
  this->mxz = A0*rhs.xz() + A1*rhs.yz() + A2*rhs.zz();
  this->myx = A3*rhs.xx() + A4*rhs.yx() + A5*rhs.zx();
  this->myy = A3*rhs.xy() + A4*rhs.yy() + A5*rhs.zy();
  this->myz = A3*rhs.xz() + A4*rhs.yz() + A5*rhs.zz();
  this->mzx = A6*rhs.xx() + A7*rhs.yx() + A8*rhs.zx();
  this->mzy = A6*rhs.xy() + A7*rhs.yy() + A8*rhs.zy();
  this->mzz = A6*rhs.xz() + A7*rhs.yz() + A8*rhs.zz();
  return *this;
}

// //------------------------------------------------------------------------------
// // Add a scalar to this tensor in place.
// //------------------------------------------------------------------------------
// template<>
// SPHERAL_HOST_DEVICE
// inline
// GeomTensor<1>&
// GeomTensor<1>::operator+=(const double rhs) {
//   this->mxx += rhs;
//   return *this;
// }

// template<>
// SPHERAL_HOST_DEVICE
// inline
// GeomTensor<2>&
// GeomTensor<2>::operator+=(const double rhs) {
//   this->mxx += rhs;
//   this->mxy += rhs;
//   this->myx += rhs;
//   this->myy += rhs;
//   return *this;
// }

// template<>
// SPHERAL_HOST_DEVICE
// inline
// GeomTensor<3>&
// GeomTensor<3>::operator+=(const double rhs) {
//   this->mxx += rhs;
//   this->mxy += rhs;
//   this->mxz += rhs;
//   this->myx += rhs;
//   this->myy += rhs;
//   this->myz += rhs;
//   this->mzx += rhs;
//   this->mzy += rhs;
//   this->mzz += rhs;
//   return *this;
// }

// //------------------------------------------------------------------------------
// // Subtract a scalar from this tensor in place.
// //------------------------------------------------------------------------------
// template<>
// SPHERAL_HOST_DEVICE
// inline
// GeomTensor<1>&
// GeomTensor<1>::operator-=(const double rhs) {
//   this->mxx -= rhs;
//   return *this;
// }

// template<>
// SPHERAL_HOST_DEVICE
// inline
// GeomTensor<2>&
// GeomTensor<2>::operator-=(const double rhs) {
//   this->mxx -= rhs;
//   this->mxy -= rhs;
//   this->myx -= rhs;
//   this->myy -= rhs;
//   return *this;
// }

// template<>
// SPHERAL_HOST_DEVICE
// inline
// GeomTensor<3>&
// GeomTensor<3>::operator-=(const double rhs) {
//   this->mxx -= rhs;
//   this->mxy -= rhs;
//   this->mxz -= rhs;
//   this->myx -= rhs;
//   this->myy -= rhs;
//   this->myz -= rhs;
//   this->mzx -= rhs;
//   this->mzy -= rhs;
//   this->mzz -= rhs;
//   return *this;
// }

//------------------------------------------------------------------------------
// Multiply this tensor by a scalar in place.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<1>&
GeomTensor<1>::operator*=(const double rhs) {
  this->mxx *= rhs;
  return *this;
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<2>&
GeomTensor<2>::operator*=(const double rhs) {
  this->mxx *= rhs;
  this->mxy *= rhs;
  this->myx *= rhs;
  this->myy *= rhs;
  return *this;
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<3>&
GeomTensor<3>::operator*=(const double rhs) {
  this->mxx *= rhs;
  this->mxy *= rhs;
  this->mxz *= rhs;
  this->myx *= rhs;
  this->myy *= rhs;
  this->myz *= rhs;
  this->mzx *= rhs;
  this->mzy *= rhs;
  this->mzz *= rhs;
  return *this;
}

//------------------------------------------------------------------------------
// Divide this tensor by a scalar in place
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<1>&
GeomTensor<1>::operator/=(const double rhs) {
  REQUIRE(rhs != 0.0);
  this->mxx /= rhs;
  return *this;
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<2>&
GeomTensor<2>::operator/=(const double rhs) {
  REQUIRE(rhs != 0.0);
  const double rhsInv = 1.0/rhs;
  this->mxx *= rhsInv;
  this->mxy *= rhsInv;
  this->myx *= rhsInv;
  this->myy *= rhsInv;
  return *this;
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<3>&
GeomTensor<3>::operator/=(const double rhs) {
  REQUIRE(rhs != 0.0);
  const double rhsInv = 1.0/rhs;
  this->mxx *= rhsInv;
  this->mxy *= rhsInv;
  this->mxz *= rhsInv;
  this->myx *= rhsInv;
  this->myy *= rhsInv;
  this->myz *= rhsInv;
  this->mzx *= rhsInv;
  this->mzy *= rhsInv;
  this->mzz *= rhsInv;
  return *this;
}

//------------------------------------------------------------------------------
// Define the equivalence operator.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
bool
GeomTensor<1>::
operator==(const GeomTensor<1>& rhs) const {
  return this->mxx == rhs.xx();
}

template<>
SPHERAL_HOST_DEVICE
inline
bool
GeomTensor<2>::
operator==(const GeomTensor<2>& rhs) const {
  return (this->mxx == rhs.xx() and
          this->mxy == rhs.xy() and
          this->myx == rhs.yx() and
          this->myy == rhs.yy());
}

template<>
SPHERAL_HOST_DEVICE
inline
bool
GeomTensor<3>::
operator==(const GeomTensor<3>& rhs) const {
  return (this->mxx == rhs.xx() and
          this->mxy == rhs.xy() and
          this->mxz == rhs.xz() and
          this->myx == rhs.yx() and
          this->myy == rhs.yy() and
          this->myz == rhs.yz() and
          this->mzx == rhs.zx() and
          this->mzy == rhs.zy() and
          this->mzz == rhs.zz());
}

template<>
SPHERAL_HOST_DEVICE
inline
bool
GeomTensor<1>::
operator==(const GeomSymmetricTensor<1>& rhs) const {
  return this->mxx == rhs.xx();
}

template<>
SPHERAL_HOST_DEVICE
inline
bool
GeomTensor<2>::
operator==(const GeomSymmetricTensor<2>& rhs) const {
  return (this->mxx == rhs.xx() and
          this->mxy == rhs.xy() and
          this->myx == rhs.yx() and
          this->myy == rhs.yy());
}

template<>
SPHERAL_HOST_DEVICE
inline
bool
GeomTensor<3>::
operator==(const GeomSymmetricTensor<3>& rhs) const {
  return (this->mxx == rhs.xx() and
          this->mxy == rhs.xy() and
          this->mxz == rhs.xz() and
          this->myx == rhs.yx() and
          this->myy == rhs.yy() and
          this->myz == rhs.yz() and
          this->mzx == rhs.zx() and
          this->mzy == rhs.zy() and
          this->mzz == rhs.zz());
}

// template<>
// SPHERAL_HOST_DEVICE
// inline
// bool
// GeomTensor<1>::
// operator==(const double rhs) const {
//   return this->mxx == rhs;
// }

// template<>
// SPHERAL_HOST_DEVICE
// inline
// bool
// GeomTensor<2>::
// operator==(const double rhs) const {
//   return (this->mxx == rhs and
//           this->mxy == rhs and
//           this->myx == rhs and
//           this->myy == rhs);
// }

// template<>
// SPHERAL_HOST_DEVICE
// inline
// bool
// GeomTensor<3>::
// operator==(const double rhs) const {
//   return (this->mxx == rhs and
//           this->mxy == rhs and
//           this->mxz == rhs and
//           this->myx == rhs and
//           this->myy == rhs and
//           this->myz == rhs and
//           this->mzx == rhs and
//           this->mzy == rhs and
//           this->mzz == rhs);
// }

//------------------------------------------------------------------------------
// Define the not equivalence than comparitor.
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
bool
GeomTensor<nDim>::
operator!=(const GeomTensor<nDim>& rhs) const {
  return !(*this == rhs);
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
bool
GeomTensor<nDim>::
operator!=(const GeomSymmetricTensor<nDim>& rhs) const {
  return !(*this == rhs);
}

// template<int nDim>
// SPHERAL_HOST_DEVICE
// inline
// bool
// GeomTensor<nDim>::
// operator!=(const double rhs) const {
//   return !(*this == rhs);
// }

//------------------------------------------------------------------------------
// Define the less than operator.
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
bool
GeomTensor<nDim>::
operator<(const GeomTensor<nDim>& rhs) const {
  return this->Determinant() < rhs.Determinant();
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
bool
GeomTensor<nDim>::
operator<(const GeomSymmetricTensor<nDim>& rhs) const {
  return this->Determinant() < rhs.Determinant();
}

// template<>
// SPHERAL_HOST_DEVICE
// inline
// bool
// GeomTensor<1>::
// operator<(const double rhs) const {
//   return this->mxx < rhs;
// }

// template<>
// SPHERAL_HOST_DEVICE
// inline
// bool
// GeomTensor<2>::
// operator<(const double rhs) const {
//   return (this->mxx < rhs and
//           this->mxy < rhs and
//           this->myx < rhs and
//           this->myy < rhs);
// }

// template<>
// SPHERAL_HOST_DEVICE
// inline
// bool
// GeomTensor<3>::
// operator<(const double rhs) const {
//   return (this->mxx < rhs and
//           this->mxy < rhs and
//           this->mxz < rhs and
//           this->myx < rhs and
//           this->myy < rhs and
//           this->myz < rhs and
//           this->mzx < rhs and
//           this->mzy < rhs and
//           this->mzz < rhs);
// }

//------------------------------------------------------------------------------
// Define the greater than operator.
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
bool
GeomTensor<nDim>::
operator>(const GeomTensor<nDim>& rhs) const {
  return this->Determinant() > rhs.Determinant();
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
bool
GeomTensor<nDim>::
operator>(const GeomSymmetricTensor<nDim>& rhs) const {
  return this->Determinant() > rhs.Determinant();
}

// template<int nDim>
// SPHERAL_HOST_DEVICE
// inline
// bool
// GeomTensor<nDim>::
// operator>(const double rhs) const {
//   return this->mxx > rhs;
// }

// template<>
// SPHERAL_HOST_DEVICE
// inline
// bool
// GeomTensor<2>::
// operator>(const double rhs) const {
//   return (this->mxx > rhs and
//           this->mxy > rhs and
//           this->myx > rhs and
//           this->myy > rhs);
// }

// template<>
// SPHERAL_HOST_DEVICE
// inline
// bool
// GeomTensor<3>::
// operator>(const double rhs) const {
//   return (this->mxx > rhs and
//           this->mxy > rhs and
//           this->mxz > rhs and
//           this->myx > rhs and
//           this->myy > rhs and
//           this->myz > rhs and
//           this->mzx > rhs and
//           this->mzy > rhs and
//           this->mzz > rhs);
// }

//------------------------------------------------------------------------------
// Define the less than or equal operator.
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
bool
GeomTensor<nDim>::
operator<=(const GeomTensor<nDim>& rhs) const {
  return (*this < rhs) or (*this == rhs);
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
bool
GeomTensor<nDim>::
operator<=(const GeomSymmetricTensor<nDim>& rhs) const {
  return (*this < rhs) or (*this == rhs);
}

// template<>
// SPHERAL_HOST_DEVICE
// inline
// bool
// GeomTensor<1>::
// operator<=(const double rhs) const {
//   return this->mxx <= rhs;
// }

// template<>
// SPHERAL_HOST_DEVICE
// inline
// bool
// GeomTensor<2>::
// operator<=(const double rhs) const {
//   return (this->mxx <= rhs and
//           this->mxy <= rhs and
//           this->myx <= rhs and
//           this->myy <= rhs);
// }

// template<>
// SPHERAL_HOST_DEVICE
// inline
// bool
// GeomTensor<3>::
// operator<=(const double rhs) const {
//   return (this->mxx <= rhs and
//           this->mxy <= rhs and
//           this->mxz <= rhs and
//           this->myx <= rhs and
//           this->myy <= rhs and
//           this->myz <= rhs and
//           this->mzx <= rhs and
//           this->mzy <= rhs and
//           this->mzz <= rhs);
// }

//------------------------------------------------------------------------------
// Define the greater than or equal operator.
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
bool
GeomTensor<nDim>::
operator>=(const GeomTensor<nDim>& rhs) const {
  return (*this > rhs) or (*this == rhs);
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
bool
GeomTensor<nDim>::
operator>=(const GeomSymmetricTensor<nDim>& rhs) const {
  return (*this > rhs) or (*this == rhs);
}

// template<int nDim>
// SPHERAL_HOST_DEVICE
// inline
// bool
// GeomTensor<nDim>::
// operator>=(const double rhs) const {
//   return this->mxx >= rhs;
// }

// template<>
// SPHERAL_HOST_DEVICE
// inline
// bool
// GeomTensor<2>::
// operator>=(const double rhs) const {
//   return (this->mxx >= rhs and
//           this->mxy >= rhs and
//           this->myx >= rhs and
//           this->myy >= rhs);
// }

// template<>
// SPHERAL_HOST_DEVICE
// inline
// bool
// GeomTensor<3>::
// operator>=(const double rhs) const {
//   return (this->mxx >= rhs and
//           this->mxy >= rhs and
//           this->mxz >= rhs and
//           this->myx >= rhs and
//           this->myy >= rhs and
//           this->myz >= rhs and
//           this->mzx >= rhs and
//           this->mzy >= rhs and
//           this->mzz >= rhs);
// }

//------------------------------------------------------------------------------
// Return the symmetric part of a GeomTensor.
//   Bij = 0.5*(Aij + Aji)
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<1>
GeomTensor<1>::Symmetric() const {
  return GeomSymmetricTensor<1>(this->mxx);
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<2>
GeomTensor<2>::Symmetric() const {
  GeomSymmetricTensor<2> result;
  result.xx(this->mxx);
  result.xy(0.5*(this->mxy + this->myx));
  result.yy(myy);
  return result;
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<3>
GeomTensor<3>::Symmetric() const {
  GeomSymmetricTensor<3> result;
  result.xx(this->mxx);
  result.xy(0.5*(this->mxy + this->myx));
  result.xz(0.5*(this->mxz + this->mzx));
  result.yy(this->myy);
  result.yz(0.5*(this->myz + this->mzy));
  result.zz(this->mzz);
  return result;
}

//------------------------------------------------------------------------------
// Return the skew-symmetric part of a GeomTensor.
//   Bij = 0.5*(Aij - Aji)
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<1>
GeomTensor<1>::SkewSymmetric() const {
  return GeomTensor<1>(0.0);
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<2>
GeomTensor<2>::SkewSymmetric() const {
  GeomTensor<2> result;
  result.xy(0.5*(this->mxy - this->myx));
  result.yx(-(result.xy()));
  return result;
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<3>
GeomTensor<3>::SkewSymmetric() const {
  GeomTensor<3> result;
  result.xy(0.5*(this->mxy - this->myx));
  result.xz(0.5*(this->mxz - this->mzx));
  result.yz(0.5*(this->myz - this->mzy));
  result.yx(-(result.xy()));
  result.zx(-(result.xz()));
  result.zy(-(result.yz()));
  return result;
}

//------------------------------------------------------------------------------
// Return the transpose of the GeomTensor.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<1>
GeomTensor<1>::
Transpose() const {
  return *this;
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<2>
GeomTensor<2>::
Transpose() const {
  return GeomTensor<2>(this->mxx, this->myx,
                       this->mxy, this->myy);
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<3>
GeomTensor<3>::
Transpose() const {
  return GeomTensor<3>(this->mxx, this->myx, this->mzx,
                       this->mxy, this->myy, this->mzy,
                       this->mxz, this->myz, this->mzz);
}

//------------------------------------------------------------------------------
// Return the inverse of the tensor.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<1>
GeomTensor<1>::Inverse() const {
  CHECK(this->mxx != 0.0);
  return GeomTensor<1>(1.0/(this->mxx));
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<2>
GeomTensor<2>::Inverse() const {
  REQUIRE(Determinant() != 0.0);
  return GeomTensor<2>((this->myy), -(this->mxy),
                       -(this->myx), (this->mxx))/this->Determinant();
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<3>
GeomTensor<3>::Inverse() const {
  REQUIRE(Determinant() != 0.0);
  return GeomTensor<3>((this->myy)*(this->mzz) - (this->myz)*(this->mzy), (this->mxz)*(this->mzy) - (this->mxy)*(this->mzz), (this->mxy)*(this->myz) - (this->mxz)*(this->myy),
                       (this->myz)*(this->mzx) - (this->myx)*(this->mzz), (this->mxx)*(this->mzz) - (this->mxz)*(this->mzx), (this->mxz)*(this->myx) - (this->mxx)*(this->myz),
                       (this->myx)*(this->mzy) - (this->myy)*(this->mzx), (this->mxy)*(this->mzx) - (this->mxx)*(this->mzy), (this->mxx)*(this->myy) - (this->mxy)*(this->myx))/this->Determinant();
}

//------------------------------------------------------------------------------
// Return the diagonal elements of the GeomTensor.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<1>
GeomTensor<1>::diagonalElements() const {
  return GeomVector<1>(this->mxx);
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<2>
GeomTensor<2>::diagonalElements() const {
  return GeomVector<2>(this->mxx, this->myy);
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<3>
GeomTensor<3>::diagonalElements() const {
  return GeomVector<3>(this->mxx, this->myy, this->mzz);
}

//------------------------------------------------------------------------------
// Return the trace of the tensor.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
double
GeomTensor<1>::Trace() const {
  return this->mxx;
}

template<>
SPHERAL_HOST_DEVICE
inline
double
GeomTensor<2>::Trace() const {
  return this->mxx + this->myy;
}

template<>
SPHERAL_HOST_DEVICE
inline
double
GeomTensor<3>::Trace() const {
  return this->mxx + this->myy + this->mzz;
}

//------------------------------------------------------------------------------
// Return the determinant of the tensor.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
double
GeomTensor<1>::Determinant() const {
  return this->mxx;
}

template<>
SPHERAL_HOST_DEVICE
inline
double
GeomTensor<2>::Determinant() const {
  return (this->mxx)*(this->myy) - (this->mxy)*(this->myx);
}

template<>
SPHERAL_HOST_DEVICE
inline
double
GeomTensor<3>::Determinant() const {
  return ((this->mxx)*(this->myy)*(this->mzz) + (this->mxy)*(this->myz)*(this->mzx) + (this->mxz)*(this->myx)*(this->mzy) -
	  (this->mxx)*(this->myz)*(this->mzy) - (this->mxy)*(this->myx)*(this->mzz) - (this->mxz)*(this->myy)*(this->mzx));
}

//------------------------------------------------------------------------------
// Multiply a tensor with a vector.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<1>
GeomTensor<1>::dot(const GeomVector<1>& rhs) const {
  return GeomVector<1>((this->mxx)*rhs.x());
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<2>
GeomTensor<2>::dot(const GeomVector<2>& rhs) const {
  return GeomVector<2>((this->mxx)*rhs.x() + (this->mxy)*rhs.y(),
                       (this->myx)*rhs.x() + (this->myy)*rhs.y());
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<3>
GeomTensor<3>::dot(const GeomVector<3>& rhs) const {
  return GeomVector<3>((this->mxx)*rhs.x() + (this->mxy)*rhs.y() + (this->mxz)*rhs.z(),
                       (this->myx)*rhs.x() + (this->myy)*rhs.y() + (this->myz)*rhs.z(),
                       (this->mzx)*rhs.x() + (this->mzy)*rhs.y() + (this->mzz)*rhs.z());
}

//------------------------------------------------------------------------------
// Multiply two tensors.  This is just the linear algebra definition for matrix
// multiplication.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<1>
GeomTensor<1>::dot(const GeomTensor<1>& rhs) const {
  return GeomTensor<1>((this->mxx) * rhs.mxx);
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<2>
GeomTensor<2>::dot(const GeomTensor<2>& rhs) const {
  return GeomTensor<2>((this->mxx)*(rhs.mxx) + (this->mxy)*(rhs.myx),
                       (this->mxx)*(rhs.mxy) + (this->mxy)*(rhs.myy),
                       (this->myx)*(rhs.mxx) + (this->myy)*(rhs.myx),
                       (this->myx)*(rhs.mxy) + (this->myy)*(rhs.myy));
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<3>
GeomTensor<3>::dot(const GeomTensor<3>& rhs) const {
  return GeomTensor<3>((this->mxx)*(rhs.mxx) + (this->mxy)*(rhs.myx) + (this->mxz)*(rhs.mzx),
                       (this->mxx)*(rhs.mxy) + (this->mxy)*(rhs.myy) + (this->mxz)*(rhs.mzy),
                       (this->mxx)*(rhs.mxz) + (this->mxy)*(rhs.myz) + (this->mxz)*(rhs.mzz),
                       (this->myx)*(rhs.mxx) + (this->myy)*(rhs.myx) + (this->myz)*(rhs.mzx),
                       (this->myx)*(rhs.mxy) + (this->myy)*(rhs.myy) + (this->myz)*(rhs.mzy),
                       (this->myx)*(rhs.mxz) + (this->myy)*(rhs.myz) + (this->myz)*(rhs.mzz),
                       (this->mzx)*(rhs.mxx) + (this->mzy)*(rhs.myx) + (this->mzz)*(rhs.mzx),
                       (this->mzx)*(rhs.mxy) + (this->mzy)*(rhs.myy) + (this->mzz)*(rhs.mzy),
                       (this->mzx)*(rhs.mxz) + (this->mzy)*(rhs.myz) + (this->mzz)*(rhs.mzz));
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<1>
GeomTensor<1>::dot(const GeomSymmetricTensor<1>& rhs) const {
  return GeomTensor<1>((this->mxx) * rhs.xx());
  return *this;
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<2>
GeomTensor<2>::dot(const GeomSymmetricTensor<2>& rhs) const {
  return GeomTensor<2>((this->mxx)*(rhs.xx()) + (this->mxy)*(rhs.yx()),
                       (this->mxx)*(rhs.xy()) + (this->mxy)*(rhs.yy()),
                       (this->myx)*(rhs.xx()) + (this->myy)*(rhs.yx()),
                       (this->myx)*(rhs.xy()) + (this->myy)*(rhs.yy()));
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<3>
GeomTensor<3>::dot(const GeomSymmetricTensor<3>& rhs) const {
  return GeomTensor<3>((this->mxx)*(rhs.xx()) + (this->mxy)*(rhs.yx()) + (this->mxz)*(rhs.zx()),
                       (this->mxx)*(rhs.xy()) + (this->mxy)*(rhs.yy()) + (this->mxz)*(rhs.zy()),
                       (this->mxx)*(rhs.xz()) + (this->mxy)*(rhs.yz()) + (this->mxz)*(rhs.zz()),
                       (this->myx)*(rhs.xx()) + (this->myy)*(rhs.yx()) + (this->myz)*(rhs.zx()),
                       (this->myx)*(rhs.xy()) + (this->myy)*(rhs.yy()) + (this->myz)*(rhs.zy()),
                       (this->myx)*(rhs.xz()) + (this->myy)*(rhs.yz()) + (this->myz)*(rhs.zz()),
                       (this->mzx)*(rhs.xx()) + (this->mzy)*(rhs.yx()) + (this->mzz)*(rhs.zx()),
                       (this->mzx)*(rhs.xy()) + (this->mzy)*(rhs.yy()) + (this->mzz)*(rhs.zy()),
                       (this->mzx)*(rhs.xz()) + (this->mzy)*(rhs.yz()) + (this->mzz)*(rhs.zz()));
}

//------------------------------------------------------------------------------
// Return the doubledot product.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
double
GeomTensor<1>::
doubledot(const GeomTensor<1>& rhs) const {
  return (this->mxx)*(rhs.xx());
}

template<>
SPHERAL_HOST_DEVICE
inline
double
GeomTensor<2>::
doubledot(const GeomTensor<2>& rhs) const {
  return ((this->mxx)*(rhs.xx()) + (this->mxy)*(rhs.yx()) + 
          (this->myx)*(rhs.xy()) + (this->myy)*(rhs.yy()));
}

template<>
SPHERAL_HOST_DEVICE
inline
double
GeomTensor<3>::
doubledot(const GeomTensor<3>& rhs) const {
  return ((this->mxx)*(rhs.xx()) + (this->mxy)*(rhs.yx()) + (this->mxz)*(rhs.zx()) +
          (this->myx)*(rhs.xy()) + (this->myy)*(rhs.yy()) + (this->myz)*(rhs.zy()) +
          (this->mzx)*(rhs.xz()) + (this->mzy)*(rhs.yz()) + (this->mzz)*(rhs.zz()));
}

//------------------------------------------------------------------------------
// Return the doubledot product with a symmetric tensor.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
double
GeomTensor<1>::
doubledot(const GeomSymmetricTensor<1>& rhs) const {
  return (this->mxx)*(rhs.xx());
}

template<>
SPHERAL_HOST_DEVICE
inline
double
GeomTensor<2>::
doubledot(const GeomSymmetricTensor<2>& rhs) const {
  return ((this->mxx)*(rhs.xx()) + (this->mxy)*(rhs.yx()) + 
          (this->myx)*(rhs.xy()) + (this->myy)*(rhs.yy()));
}

template<>
SPHERAL_HOST_DEVICE
inline
double
GeomTensor<3>::
doubledot(const GeomSymmetricTensor<3>& rhs) const {
  return ((this->mxx)*(rhs.xx()) + (this->mxy)*(rhs.yx()) + (this->mxz)*(rhs.zx()) +
          (this->myx)*(rhs.xy()) + (this->myy)*(rhs.yy()) + (this->myz)*(rhs.zy()) +
          (this->mzx)*(rhs.xz()) + (this->mzy)*(rhs.yz()) + (this->mzz)*(rhs.zz()));
}

//------------------------------------------------------------------------------
// Return the doubledot product with ourself.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
double
GeomTensor<1>::
selfDoubledot() const {
  return FastMath::square(this->mxx);
}

template<>
SPHERAL_HOST_DEVICE
inline
double
GeomTensor<2>::
selfDoubledot() const {
  return ((this->mxx)*(this->mxx) + (this->mxy)*(this->myx) + 
          (this->myx)*(this->mxy) + (this->myy)*(this->myy));
}

template<>
SPHERAL_HOST_DEVICE
inline
double
GeomTensor<3>::
selfDoubledot() const {
  return ((this->mxx)*(this->mxx) + (this->mxy)*(this->myx) + (this->mxz)*(this->mzx) +
          (this->myx)*(this->mxy) + (this->myy)*(this->myy) + (this->myz)*(this->mzy) +
          (this->mzx)*(this->mxz) + (this->mzy)*(this->myz) + (this->mzz)*(this->mzz));
}

//------------------------------------------------------------------------------
// Return the square of this tensor.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<1>
GeomTensor<1>::
square() const {
  return GeomTensor<1>((this->mxx)*(this->mxx));
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<2>
GeomTensor<2>::
square() const {
  return GeomTensor<2>((this->mxx)*(this->mxx) + (this->mxy)*(this->myx),
                       (this->mxx)*(this->mxy) + (this->mxy)*(this->myy),
                       (this->myx)*(this->mxx) + (this->myy)*(this->myx),
                       (this->myx)*(this->mxy) + (this->myy)*(this->myy));
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<3>
GeomTensor<3>::
square() const {
  return GeomTensor<3>((this->mxx)*(this->mxx) + (this->mxy)*(this->myx) + (this->mxz)*(this->mzx),
                       (this->mxx)*(this->mxy) + (this->mxy)*(this->myy) + (this->mxz)*(this->mzy),
                       (this->mxx)*(this->mxz) + (this->mxy)*(this->myz) + (this->mxz)*(this->mzz),
                       (this->myx)*(this->mxx) + (this->myy)*(this->myx) + (this->myz)*(this->mzx),
                       (this->myx)*(this->mxy) + (this->myy)*(this->myy) + (this->myz)*(this->mzy),
                       (this->myx)*(this->mxz) + (this->myy)*(this->myz) + (this->myz)*(this->mzz),
                       (this->mzx)*(this->mxx) + (this->mzy)*(this->myx) + (this->mzz)*(this->mzx),
                       (this->mzx)*(this->mxy) + (this->mzy)*(this->myy) + (this->mzz)*(this->mzy),
                       (this->mzx)*(this->mxz) + (this->mzy)*(this->myz) + (this->mzz)*(this->mzz));
}

//------------------------------------------------------------------------------
// Return a new tensor with the elements of this tensor squared.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<1>
GeomTensor<1>::
squareElements() const {
  return GeomTensor<1>((this->mxx)*(this->mxx));
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<2>
GeomTensor<2>::
squareElements() const {
  return GeomTensor<2>((this->mxx)*(this->mxx),
                       (this->mxy)*(this->mxy),
                       (this->myx)*(this->myx),
                       (this->myy)*(this->myy));
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<3>
GeomTensor<3>::
squareElements() const {
  return GeomTensor<3>((this->mxx)*(this->mxx),
                       (this->mxy)*(this->mxy),
                       (this->mxz)*(this->mxz),
                       (this->myx)*(this->myx),
                       (this->myy)*(this->myy),
                       (this->myz)*(this->myz),
                       (this->mzx)*(this->mzx),
                       (this->mzy)*(this->mzy),
                       (this->mzz)*(this->mzz));
}

//------------------------------------------------------------------------------
// Apply a rotational transform to this tensor.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
void
GeomTensor<1>::
rotationalTransform(const GeomTensor<1>& R) {
  CONTRACT_VAR(R);
  REQUIRE(fuzzyEqual(std::abs(R.Determinant()), 1.0, 1.0e-8));
}

template<>
SPHERAL_HOST_DEVICE
inline
void
GeomTensor<2>::
rotationalTransform(const GeomTensor<2>& R) {
  REQUIRE(fuzzyEqual(std::abs(R.Determinant()), 1.0, 1.0e-8));

  const double A0 = this->mxx;
  const double A1 = this->mxy;
  const double A2 = this->myx;
  const double A3 = this->myy;

  const double R0 = R.xx();
  const double R1 = R.xy();
  const double R2 = R.yx();
  const double R3 = R.yy();

  const double T0 = A0*R0 + A2*R1;
  const double T1 = A0*R2 + A2*R3;
  const double T2 = A1*R0 + A3*R1;
  const double T3 = A1*R2 + A3*R3;

  this->mxx = R0*T0 + R1*T2;
  this->mxy = R2*T0 + R3*T2;
  this->myx = R0*T1 + R1*T3;
  this->myy = R2*T1 + R3*T3;
}

template<>
SPHERAL_HOST_DEVICE
inline
void
GeomTensor<3>::
rotationalTransform(const GeomTensor<3>& R) {
  REQUIRE(fuzzyEqual(std::abs(R.Determinant()), 1.0, 1.0e-8));

  const double A0 = this->mxx;
  const double A1 = this->mxy;
  const double A2 = this->mxz;
  const double A3 = this->myx;
  const double A4 = this->myy;
  const double A5 = this->myz;
  const double A6 = this->mzx;
  const double A7 = this->mzy;
  const double A8 = this->mzz;

  const double R0 = R.xx();
  const double R1 = R.xy();
  const double R2 = R.xz();
  const double R3 = R.yx();
  const double R4 = R.yy();
  const double R5 = R.yz();
  const double R6 = R.zx();
  const double R7 = R.zy();
  const double R8 = R.zz();

  const double T0 = A0*R0 + A3*R1 + A6*R2;
  const double T1 = A0*R3 + A3*R4 + A6*R5;
  const double T2 = A0*R6 + A3*R7 + A6*R8;
  const double T3 = A1*R0 + A4*R1 + A7*R2;
  const double T4 = A1*R3 + A4*R4 + A7*R5;
  const double T5 = A1*R6 + A4*R7 + A7*R8;
  const double T6 = A2*R0 + A5*R1 + A8*R2;
  const double T7 = A2*R3 + A5*R4 + A8*R5;
  const double T8 = A2*R6 + A5*R7 + A8*R8;

  this->mxx = R0*T0 + R1*T3 + R2*T6;
  this->mxy = R3*T0 + R4*T3 + R5*T6;
  this->mxz = R6*T0 + R7*T3 + R8*T6;
  this->myx = R0*T1 + R1*T4 + R2*T7;
  this->myy = R3*T1 + R4*T4 + R5*T7;
  this->myz = R6*T1 + R7*T4 + R8*T7;
  this->mzx = R0*T2 + R1*T5 + R2*T8;
  this->mzy = R3*T2 + R4*T5 + R5*T8;
  this->mzz = R6*T2 + R7*T5 + R8*T8;
}

//------------------------------------------------------------------------------
// Return the maximum absolute value of the elements.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
double
GeomTensor<1>::
maxAbsElement() const {
  return std::abs(this->mxx);
}

template<>
SPHERAL_HOST_DEVICE
inline
double
GeomTensor<2>::
maxAbsElement() const {
  return std::max(std::abs(this->mxx),
                  std::max(std::abs(this->mxy), 
                           std::max(std::abs(this->myx),
                                    std::abs(this->myy))));
}

template<>
SPHERAL_HOST_DEVICE
inline
double
GeomTensor<3>::
maxAbsElement() const {
  return std::max(std::abs(this->mxx), 
                  std::max(std::abs(this->mxy), 
                           std::max(std::abs(this->mxz), 
                                    std::max(std::abs(this->myx),
                                             std::max(std::abs(this->myy), 
                                                      std::max(std::abs(this->myz), 
                                                               std::max(std::abs(this->mzx),
                                                                        std::max(std::abs(this->mzy),
                                                                                 std::abs(this->mzz)))))))));
}

//------------------------------------------------------------------------------
// Generate an Eigen Tensor.
//------------------------------------------------------------------------------
template<>
inline
GeomTensor<1>::EigenType
GeomTensor<1>::eigen() const {
  return EigenType(this->mxx);
}

template<>
inline
GeomTensor<2>::EigenType
GeomTensor<2>::eigen() const {
  EigenType result;
  result << this->mxx, this->mxy,
            this->myx, this->myy;
  return result;
}

template<>
inline
GeomTensor<3>::EigenType
GeomTensor<3>::eigen() const {
  EigenType result;
  result << this->mxx, this->mxy, this->mxz,
            this->myx, this->myy, this->myz,
            this->mzx, this->mzy, this->mzz;
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
GeomTensor<nDim>
operator*(double lhs, const GeomTensor<nDim>& rhs) {
  return rhs*lhs;
}

//------------------------------------------------------------------------------
// Input (istream) operator.
//------------------------------------------------------------------------------
template<int nDim>
inline
std::istream&
operator>>(std::istream& is, GeomTensor<nDim>& ten) {
  std::string parenthesis;
  is >> parenthesis;
  for (typename GeomTensor<nDim>::iterator elementItr = ten.begin();
       elementItr != ten.end();
       ++elementItr) {
    is >> *elementItr;
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
operator<<(std::ostream& os, const GeomTensor<nDim>& ten) {
  os << "( ";
  for (typename GeomTensor<nDim>::const_iterator itr = ten.begin();
       itr != ten.end(); ++itr) {
    os << *itr << " ";
  }
  os << ")";
  return os;
}

}

