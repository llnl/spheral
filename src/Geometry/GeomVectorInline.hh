#include "GeomTensor.hh"
#include "GeomSymmetricTensor.hh"
#include "Utilities/DBC.hh"

#include <algorithm>
#include <limits.h>
#include <math.h>
#include <cfloat>
#include <string>

namespace Spheral {

//------------------------------------------------------------------------------
// Construct from an Eigen Vector.
//------------------------------------------------------------------------------
template<>
template<typename Derived>
inline
GeomVector<1>::GeomVector(const Eigen::MatrixBase<Derived>& vec):
  GeomVectorBase<1>(vec(0)) {
}

template<>
template<typename Derived>
inline
GeomVector<2>::GeomVector(const Eigen::MatrixBase<Derived>& vec):
  GeomVectorBase<2>(vec(0), vec(1)) {
}

template<>
template<typename Derived>
inline
GeomVector<3>::GeomVector(const Eigen::MatrixBase<Derived>& vec):
  GeomVectorBase<3>(vec(0), vec(1), vec(2)) {
}

//------------------------------------------------------------------------------
// The assignment operator (Eigen Vector).
//------------------------------------------------------------------------------
template<>
template<typename Derived>
inline
GeomVector<1>&
GeomVector<1>::operator=(const Eigen::MatrixBase<Derived>& vec) {
  this->mx = vec(0);
  return *this;
}

template<>
template<typename Derived>
inline
GeomVector<2>&
GeomVector<2>::operator=(const Eigen::MatrixBase<Derived>& vec) {
  this->mx = vec(0);
  this->my = vec(1);
  return *this;
}

template<>
template<typename Derived>
inline
GeomVector<3>&
GeomVector<3>::operator=(const Eigen::MatrixBase<Derived>& vec) {
  this->mx = vec(0);
  this->my = vec(1);
  this->mz = vec(2);
  return *this;
}

//------------------------------------------------------------------------------
// Set the vector elements to a constant scalar value.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<1>&
GeomVector<1>::operator=(const double val) {
  this->mx = val;
  return *this;
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<2>&
GeomVector<2>::operator=(const double val) {
  this->mx = val;
  this->my = val;
  return *this;
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<3>&
GeomVector<3>::operator=(const double val) {
  this->mx = val;
  this->my = val;
  this->mz = val;
  return *this;
}

//------------------------------------------------------------------------------
// Return the (index) element using the parenthesis operator.
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
double
GeomVector<nDim>::operator()(typename GeomVector<nDim>::size_type index) const {
  REQUIRE(index < nDim);
  return *(begin() + index);
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
double&
GeomVector<nDim>::operator()(typename GeomVector<nDim>::size_type index) {
  REQUIRE(index < nDim);
  return *(begin() + index);
}

//------------------------------------------------------------------------------
// Return the (index) element using the bracket operator.
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
double
GeomVector<nDim>::operator[](typename GeomVector<nDim>::size_type index) const {
  REQUIRE(index < nDim);
  return *(begin() + index);
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
double&
GeomVector<nDim>::operator[](typename GeomVector<nDim>::size_type index) {
  REQUIRE(index < nDim);
  return *(begin() + index);
}

//------------------------------------------------------------------------------
// Return the x (first) element.
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
double
GeomVector<nDim>::x() const {
  return this->mx;
}

//------------------------------------------------------------------------------
// Return the y (second) element
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
double
GeomVector<nDim>::y() const {
  REQUIRE(nDim > 1);
  return this->my;
}

template<>
SPHERAL_HOST_DEVICE
inline
double
GeomVector<1>::y() const {
  return 0.0;
}

//------------------------------------------------------------------------------
// Return the z (third) element
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
double
GeomVector<nDim>::z() const {
  REQUIRE(nDim > 2);
  return this->mz;
}

template<>
SPHERAL_HOST_DEVICE
inline
double
GeomVector<1>::z() const {
  return 0.0;
}

template<>
SPHERAL_HOST_DEVICE
inline
double
GeomVector<2>::z() const {
  return 0.0;
}

//------------------------------------------------------------------------------
// Set the x (first) element.
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
void
GeomVector<nDim>::x(const double val) {
  this->mx = val;
}

//------------------------------------------------------------------------------
// Set the y (second) element
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
void
GeomVector<nDim>::y(const double val) {
  REQUIRE(nDim > 1);
  this->my = val;
}

template<>
SPHERAL_HOST_DEVICE
inline
void
GeomVector<1>::y(const double /*val*/) {
}

//------------------------------------------------------------------------------
// Set the z (third) element
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
void
GeomVector<nDim>::z(const double val) {
  REQUIRE(nDim > 2);
  this->mz = val;
}

template<>
SPHERAL_HOST_DEVICE
inline
void
GeomVector<1>::z(const double /*val*/) {
}

template<>
SPHERAL_HOST_DEVICE
inline
void
GeomVector<2>::z(const double /*val*/) {
}

//------------------------------------------------------------------------------
// Provide begin/end iterators over the elements of the Vector.
//------------------------------------------------------------------------------
// Non-const versions.
template<int nDim>
SPHERAL_HOST_DEVICE
inline
typename GeomVector<nDim>::iterator
GeomVector<nDim>::begin() {
  return &(this->mx);
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
typename GeomVector<nDim>::iterator
GeomVector<nDim>::end() {
  return &(this->mx) + nDim;
}

// Const versions.
template<int nDim>
SPHERAL_HOST_DEVICE
inline
typename GeomVector<nDim>::const_iterator
GeomVector<nDim>::begin() const {
  return &(this->mx);
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
typename GeomVector<nDim>::const_iterator
GeomVector<nDim>::end() const {
  return &(this->mx) + nDim;
}

//------------------------------------------------------------------------------
// Zero out the Vector.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
void
GeomVector<1>::Zero() {
  this->mx = 0.0;
}

template<>
SPHERAL_HOST_DEVICE
inline
void
GeomVector<2>::Zero() {
  this->mx = 0.0;
  this->my = 0.0;
}

template<>
SPHERAL_HOST_DEVICE
inline
void
GeomVector<3>::Zero() {
  this->mx = 0.0;
  this->my = 0.0;
  this->mz = 0.0;
}

//------------------------------------------------------------------------------
// Return the negative of a vector.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<1>
GeomVector<1>::operator-() const {
  return GeomVector<1>(-(this->mx));
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<2>
GeomVector<2>::operator-() const {
  return GeomVector<2>(-(this->mx), -(this->my));
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<3>
GeomVector<3>::operator-() const {
  return GeomVector<3>(-(this->mx), -(this->my), -(this->mz));
}

//------------------------------------------------------------------------------
// Add two vectors.
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
GeomVector<nDim>
GeomVector<nDim>::operator+(const GeomVector<nDim>& vec) const {
  GeomVector<nDim> result(*this);
  result += vec;
  return result;
}

//------------------------------------------------------------------------------
// Subtract a vector from another.
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
GeomVector<nDim>
GeomVector<nDim>::operator-(const GeomVector<nDim>& vec) const {
  GeomVector<nDim> result(*this);
  result -= vec;
  return result;
}

//------------------------------------------------------------------------------
// Mutiply two vectors.  For our purposes this returns the dyad.
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
GeomTensor<nDim>
GeomVector<nDim>::operator*(const GeomVector<nDim>& vec) const {
  return this->dyad(vec);
}

//------------------------------------------------------------------------------
// Multiply a vector by a scalar
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
GeomVector<nDim>
GeomVector<nDim>::operator*(const double val) const {
  GeomVector<nDim> result(*this);
  result *= val;
  return result;
}

//------------------------------------------------------------------------------
// Divide a vector by a scalar
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
GeomVector<nDim>
GeomVector<nDim>::operator/(const double val) const {
  CHECK(val != 0.0);
  GeomVector<nDim> result(*this);
  result /= val;
  return result;
}

//------------------------------------------------------------------------------
// Add two vectors in place.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<1>&
GeomVector<1>::operator+=(const GeomVector<1>& vec) {
  this->mx += vec.mx;
  return *this;
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<2>&
GeomVector<2>::operator+=(const GeomVector<2>& vec) {
  this->mx += vec.mx;
  this->my += vec.my;
  return *this;
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<3>&
GeomVector<3>::operator+=(const GeomVector<3>& vec) {
  this->mx += vec.mx;
  this->my += vec.my;
  this->mz += vec.mz;
  return *this;
}

//------------------------------------------------------------------------------
// Subtract a vector in place.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<1>& 
GeomVector<1>::operator-=(const GeomVector<1>& vec) {
  this->mx -= vec.mx;
  return *this;
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<2>&
GeomVector<2>::operator-=(const GeomVector<2>& vec) {
  this->mx -= vec.mx;
  this->my -= vec.my;
  return *this;
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<3>&
GeomVector<3>::operator-=(const GeomVector<3>& vec) {
  this->mx -= vec.mx;
  this->my -= vec.my;
  this->mz -= vec.mz;
  return *this;
}

//------------------------------------------------------------------------------
// += Eigen vector
//------------------------------------------------------------------------------
template<>
template<typename Derived>
inline
GeomVector<1>&
GeomVector<1>::operator+=(const Eigen::MatrixBase<Derived>& vec) {
  this->mx += vec(0);
  return *this;
}

template<>
template<typename Derived>
inline
GeomVector<2>&
GeomVector<2>::operator+=(const Eigen::MatrixBase<Derived>& vec) {
  this->mx += vec(0);
  this->my += vec(1);
  return *this;
}

template<>
template<typename Derived>
inline
GeomVector<3>&
GeomVector<3>::operator+=(const Eigen::MatrixBase<Derived>& vec) {
  this->mx += vec(0);
  this->my += vec(1);
  this->mz += vec(2);
  return *this;
}

//------------------------------------------------------------------------------
// -= Eigen vector
//------------------------------------------------------------------------------
template<>
template<typename Derived>
inline
GeomVector<1>&
GeomVector<1>::operator-=(const Eigen::MatrixBase<Derived>& vec) {
  this->mx -= vec(0);
  return *this;
}

template<>
template<typename Derived>
inline
GeomVector<2>&
GeomVector<2>::operator-=(const Eigen::MatrixBase<Derived>& vec) {
  this->mx -= vec(0);
  this->my -= vec(1);
  return *this;
}

template<>
template<typename Derived>
inline
GeomVector<3>&
GeomVector<3>::operator-=(const Eigen::MatrixBase<Derived>& vec) {
  this->mx -= vec(0);
  this->my -= vec(1);
  this->mz -= vec(2);
  return *this;
}

//------------------------------------------------------------------------------
// Multiply this vector by a scalar in place.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<1>&
GeomVector<1>::operator*=(const double val) {
  this->mx *= val;
  return *this;
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<2>&
GeomVector<2>::operator*=(const double val) {
  this->mx *= val;
  this->my *= val;
  return *this;
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<3>&
GeomVector<3>::operator*=(const double val) {
  this->mx *= val;
  this->my *= val;
  this->mz *= val;
  return *this;
}

//------------------------------------------------------------------------------
// Divide this vector by a scalar in place.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<1>&
GeomVector<1>::operator/=(const double val) {
  REQUIRE(val != 0.0);
  this->mx /= val;
  return *this;
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<2>&
GeomVector<2>::operator/=(const double val) {
  REQUIRE(val != 0.0);
  const double valInv = 1.0/val;
  this->mx *= valInv;
  this->my *= valInv;
  return *this;
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<3>&
GeomVector<3>::operator/=(const double val) {
  REQUIRE(val != 0.0);
  const double valInv = 1.0/val;
  this->mx *= valInv;
  this->my *= valInv;
  this->mz *= valInv;
  return *this;
}

//------------------------------------------------------------------------------
// Return (-1, 0, 1) if this vector is (less than, equal to, greater than)
// the given vector.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
int
GeomVector<1>::compare(const GeomVector<1>& vec) const {
  return (this->mx < vec.mx ? -1 :
          this->mx > vec.mx ?  1 :
          0);
}

template<>
SPHERAL_HOST_DEVICE
inline
int
GeomVector<2>::compare(const GeomVector<2>& vec) const {
  return (this->my < vec.my ? -1 :
          this->my > vec.my ?  1 :
          this->mx < vec.mx ? -1 :
          this->mx > vec.mx ?  1 :
          0);
}

template<>
SPHERAL_HOST_DEVICE
inline
int
GeomVector<3>::compare(const GeomVector<3>& vec) const {
  return (this->mz < vec.mz ? -1 :
          this->mz > vec.mz ?  1 :
          this->my < vec.my ? -1 :
          this->my > vec.my ?  1 :
          this->mx < vec.mx ? -1 :
          this->mx > vec.mx ?  1 :
          0);
}

//------------------------------------------------------------------------------
// Return (-1, 0, 1) if this vector is (less than, equal to, greater than)
// the given double.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
int
GeomVector<1>::compare(const double val) const {
  return (this->mx < val ? -1 :
          this->mx > val ?  1 :
          0);
}

template<>
SPHERAL_HOST_DEVICE
inline
int
GeomVector<2>::compare(const double val) const {
  return (this->my < val ? -1 :
          this->my > val ?  1 :
          this->mx < val ? -1 :
          this->mx > val ?  1 :
          0);
}

template<>
SPHERAL_HOST_DEVICE
inline
int
GeomVector<3>::compare(const double val) const {
  return (this->mz < val ? -1 :
          this->mz > val ?  1 :
          this->my < val ? -1 :
          this->my > val ?  1 :
          this->mx < val ? -1 :
          this->mx > val ?  1 :
          0);
}

//------------------------------------------------------------------------------
// The equivalence comparator.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
bool
GeomVector<1>::operator==(const GeomVector<1>& vec) const {
  return this->mx == vec.mx;
}

template<>
SPHERAL_HOST_DEVICE
inline
bool
GeomVector<2>::operator==(const GeomVector<2>& vec) const {
  return (this->mx == vec.mx) and (this->my == vec.my);
}

template<>
SPHERAL_HOST_DEVICE
inline
bool
GeomVector<3>::operator==(const GeomVector<3>& vec) const {
  return (this->mx == vec.mx) and (this->my == vec.my) and (this->mz == vec.mz);
}

//------------------------------------------------------------------------------
// The non-equivalence comparator.
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
bool
GeomVector<nDim>::operator!=(const GeomVector<nDim>& vec) const {
  return not (*this == vec);
}

//------------------------------------------------------------------------------
// The less than comparator.
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
bool
GeomVector<nDim>::operator<(const GeomVector<nDim>& rhs) const {
  return this->compare(rhs) == -1;
}

//------------------------------------------------------------------------------
// The greater than comparator.
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
bool
GeomVector<nDim>::operator>(const GeomVector<nDim>& rhs) const {
  return this->compare(rhs) == 1;
}

//------------------------------------------------------------------------------
// The less than or equal comparator.
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
bool
GeomVector<nDim>::operator<=(const GeomVector<nDim>& rhs) const {
  return this->compare(rhs) <= 0;
}

//------------------------------------------------------------------------------
// The greater than or equal comparator.
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
bool
GeomVector<nDim>::operator>=(const GeomVector<nDim>& rhs) const {
  return this->compare(rhs) >= 0;
}

//------------------------------------------------------------------------------
// The equivalence comparator (double).
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
bool
GeomVector<1>::operator==(const double val) const {
  return this->mx == val;
}

template<>
SPHERAL_HOST_DEVICE
inline
bool
GeomVector<2>::operator==(const double val) const {
  return (this->mx == val) and (this->my == val);
}

template<>
SPHERAL_HOST_DEVICE
inline
bool
GeomVector<3>::operator==(const double val) const {
  return (this->mx == val) and (this->my == val) and (this->mz == val);
}

//------------------------------------------------------------------------------
// The non-equivalence comparator (double).
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
bool
GeomVector<nDim>::operator!=(const double val) const {
  return not (*this == val);
}

//------------------------------------------------------------------------------
// The less than comparator (double).
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
bool
GeomVector<nDim>::operator<(const double val) const {
  return this->compare(val) == -1;
}

//------------------------------------------------------------------------------
// The greater than comparator (double).
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
bool
GeomVector<nDim>::operator>(const double val) const {
  return this->compare(val) == 1;
}

//------------------------------------------------------------------------------
// The less than or equal comparator (double).
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
bool
GeomVector<nDim>::operator<=(const double val) const {
  return this->compare(val) <= 0;
}

//------------------------------------------------------------------------------
// The greater than or equal comparator (double).
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
bool
GeomVector<nDim>::operator>=(const double val) const {
  return this->compare(val) >= 0;
}

//------------------------------------------------------------------------------
// Atomic operations.
//------------------------------------------------------------------------------
template<>
template<typename Op>
SPHERAL_HOST_DEVICE
inline
void
GeomVector<1>::atomicOp(const GeomVector<1>& vec) {
  Op::apply(&this->mx, vec.mx);
}

template<>
template<typename Op>
SPHERAL_HOST_DEVICE
inline
void
GeomVector<2>::atomicOp(const GeomVector<2>& vec) {
  Op::apply(&this->mx, vec.mx);
  Op::apply(&this->my, vec.my);
}

template<>
template<typename Op>
SPHERAL_HOST_DEVICE
inline
void
GeomVector<3>::atomicOp(const GeomVector<3>& vec) {
  Op::apply(&this->mx, vec.mx);
  Op::apply(&this->my, vec.my);
  Op::apply(&this->mz, vec.mz);
}

//------------------------------------------------------------------------------
// Dot two vectors.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
double
GeomVector<1>::dot(const GeomVector<1>& vec) const {
  return this->mx*vec.mx;
}

template<>
SPHERAL_HOST_DEVICE
inline
double
GeomVector<2>::dot(const GeomVector<2>& vec) const {
  return this->mx*vec.mx + this->my*vec.my;
}

template<>
SPHERAL_HOST_DEVICE
inline
double
GeomVector<3>::dot(const GeomVector<3>& vec) const {
  return this->mx*vec.mx + this->my*vec.my + this->mz*vec.mz;
}

//------------------------------------------------------------------------------
// Cross two vectors.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<3>
GeomVector<3>::cross(const GeomVector<3>& vec) const {
  return GeomVector<3>(my*vec.mz - this->mz*vec.my,
		       this->mz*vec.mx - this->mx*vec.mz,
		       this->mx*vec.my - this->my*vec.mx);
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<3>
GeomVector<2>::cross(const GeomVector<2>& vec) const {
  return GeomVector<3>(0.0, 0.0,
		       this->mx*vec.my - this->my*vec.mx);
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<3>
GeomVector<1>::cross(const GeomVector<1>&) const {
  return GeomVector<3>(0.0, 0.0, 0.0);
}

//------------------------------------------------------------------------------
// Perform the dyad operation with another vector.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<1>
GeomVector<1>::dyad(const GeomVector<1>& rhs) const {
  return GeomTensor<1>(this->mx*rhs(0));
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<2>
GeomVector<2>::dyad(const GeomVector<2>& rhs) const {
  return GeomTensor<2>(this->mx*rhs(0), this->mx*rhs(1),
                       this->my*rhs(0), this->my*rhs(1));
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomTensor<3>
GeomVector<3>::dyad(const GeomVector<3>& rhs) const {
  return GeomTensor<3>(this->mx*rhs(0), this->mx*rhs(1), this->mx*rhs(2),
                       this->my*rhs(0), this->my*rhs(1), this->my*rhs(2),
                       this->mz*rhs(0), this->mz*rhs(1), this->mz*rhs(2));
}

//------------------------------------------------------------------------------
// Perform the dyad operation with ourself, resulting in a symmetric tensor.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<1>
GeomVector<1>::selfdyad() const {
  return GeomSymmetricTensor<1>((this->mx)*(this->mx));
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<2>
GeomVector<2>::selfdyad() const {
  const double a = (this->mx)*(this->mx);
  const double b = (this->mx)*(this->my);
  const double c = (this->my)*(this->my);
  return GeomSymmetricTensor<2>(a, b,
                                b, c);
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomSymmetricTensor<3>
GeomVector<3>::selfdyad() const {
  const double a = (this->mx)*(this->mx);
  const double b = (this->mx)*(this->my);
  const double c = (this->mx)*(this->mz);
  const double d = (this->my)*(this->my);
  const double e = (this->my)*(this->mz);
  const double f = (this->mz)*(this->mz);
  return GeomSymmetricTensor<3>(a, b, c,
                                b, d, e,
                                c, e, f);
}

//------------------------------------------------------------------------------
// Return a unit vector with the direction of this one.
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
GeomVector<nDim>
GeomVector<nDim>::unitVector() const {
  const double mag = this->magnitude();
  return mag > 1.0e-50 ? (*this)/mag : GeomVector<nDim>(1.0, 0.0, 0.0);
}

//------------------------------------------------------------------------------
// Return the magnitude of the Vector.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
double
GeomVector<1>::magnitude() const {
  return std::abs(this->mx);
}

template<>
SPHERAL_HOST_DEVICE
inline
double
GeomVector<2>::magnitude() const {
  return sqrt((this->mx)*(this->mx) + (this->my)*(this->my));
}

template<>
SPHERAL_HOST_DEVICE
inline
double
GeomVector<3>::magnitude() const {
  return sqrt((this->mx)*(this->mx) + (this->my)*(this->my) + (this->mz)*(this->mz));
}

//------------------------------------------------------------------------------
// Return the square of the magnitude.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
double
GeomVector<1>::magnitude2() const {
  return (this->mx)*(this->mx);
}

template<>
SPHERAL_HOST_DEVICE
inline
double
GeomVector<2>::magnitude2() const {
  return (this->mx)*(this->mx) + (this->my)*(this->my);
}

template<>
SPHERAL_HOST_DEVICE
inline
double
GeomVector<3>::magnitude2() const {
  return (this->mx)*(this->mx) + (this->my)*(this->my) + (this->mz)*(this->mz);
}

//------------------------------------------------------------------------------
// Return the minimum element.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
double
GeomVector<1>::minElement() const {
  return this->mx;
}

template<>
SPHERAL_HOST_DEVICE
inline
double
GeomVector<2>::minElement() const {
  return std::min(this->mx, this->my);
}

template<>
SPHERAL_HOST_DEVICE
inline
double
GeomVector<3>::minElement() const {
  return std::min(this->mx, std::min(this->my, this->mz));
}

//------------------------------------------------------------------------------
// Return the maximum element.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
double
GeomVector<1>::maxElement() const {
  return this->mx;
}

template<>
SPHERAL_HOST_DEVICE
inline
double
GeomVector<2>::maxElement() const {
  return std::max(this->mx, this->my);
}

template<>
SPHERAL_HOST_DEVICE
inline
double
GeomVector<3>::maxElement() const {
  return std::max(this->mx, std::max(this->my, this->mz));
}

//------------------------------------------------------------------------------
// Return the maximum element by absolute value.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
double
GeomVector<1>::maxAbsElement() const {
  return std::abs(this->mx);
}

template<>
SPHERAL_HOST_DEVICE
inline
double
GeomVector<2>::maxAbsElement() const {
  return std::max(std::abs(this->mx), 
                  std::abs(this->my));
}

template<>
SPHERAL_HOST_DEVICE
inline
double
GeomVector<3>::maxAbsElement() const {
  return std::max(std::abs(this->mx),
                  std::max(std::abs(this->my),
                           std::abs(this->mz)));
}

//------------------------------------------------------------------------------
// Return the sum of the elements.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
double
GeomVector<1>::sumElements() const {
  return this->mx;
}

template<>
SPHERAL_HOST_DEVICE
inline
double
GeomVector<2>::sumElements() const {
  return this->mx + this->my;
}

template<>
SPHERAL_HOST_DEVICE
inline
double
GeomVector<3>::sumElements() const {
  return this->mx + this->my + this->mz;
}

//------------------------------------------------------------------------------
// Generate an Eigen Vector.
//------------------------------------------------------------------------------
template<>
inline
GeomVector<1>::EigenType
GeomVector<1>::eigen() const {
  return EigenType(this->mx);
}

template<>
inline
GeomVector<2>::EigenType
GeomVector<2>::eigen() const {
  return EigenType(this->mx, this->my);
}

template<>
inline
GeomVector<3>::EigenType
GeomVector<3>::eigen() const {
  return EigenType(this->mx, this->my, this->mz);
}

//******************************************************************************
// Global functions
//******************************************************************************

//------------------------------------------------------------------------------
// Multiply a scalar by a vector.
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
GeomVector<nDim>
operator*(const double val, const GeomVector<nDim>& vec) {
  return vec*val;
}

//------------------------------------------------------------------------------
// Element wise minimum comparison.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<1>
elementWiseMin(const GeomVector<1>& lhs, const GeomVector<1>& rhs) {
  return GeomVector<1>(std::min(lhs.x(), rhs.x()));
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<2>
elementWiseMin(const GeomVector<2>& lhs, const GeomVector<2>& rhs) {
  return GeomVector<2>(std::min(lhs.x(), rhs.x()),
                       std::min(lhs.y(), rhs.y()));
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<3>
elementWiseMin(const GeomVector<3>& lhs, const GeomVector<3>& rhs) {
  return GeomVector<3>(std::min(lhs.x(), rhs.x()),
                       std::min(lhs.y(), rhs.y()),
                       std::min(lhs.z(), rhs.z()));
}

//------------------------------------------------------------------------------
// Element wise maximum comparison.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<1>
elementWiseMax(const GeomVector<1>& lhs, const GeomVector<1>& rhs) {
  return GeomVector<1>(std::max(lhs.x(), rhs.x()));
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<2>
elementWiseMax(const GeomVector<2>& lhs, const GeomVector<2>& rhs) {
  return GeomVector<2>(std::max(lhs.x(), rhs.x()),
                       std::max(lhs.y(), rhs.y()));
}

template<>
SPHERAL_HOST_DEVICE
inline
GeomVector<3>
elementWiseMax(const GeomVector<3>& lhs, const GeomVector<3>& rhs) {
  return GeomVector<3>(std::max(lhs.x(), rhs.x()),
                       std::max(lhs.y(), rhs.y()),
                       std::max(lhs.z(), rhs.z()));
}

//------------------------------------------------------------------------------
// Input (istream) operator.
//------------------------------------------------------------------------------
template<int nDim>
inline
std::istream&
operator>>(std::istream& is, GeomVector<nDim>& vec) {
  std::string parenthesis;
  is >> parenthesis;
  for (typename GeomVector<nDim>::iterator elementItr = vec.begin();
       elementItr < vec.end();
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
operator<<(std::ostream& os, const GeomVector<nDim>& vec) {
  os << "( ";
  for (typename GeomVector<nDim>::const_iterator elementItr = vec.begin();
       elementItr < vec.end();
       ++elementItr) {
    os << *elementItr << " ";
  }
  os << ")";
  return os;
}

} // namespace Spheral
