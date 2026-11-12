#include "Dimension.hh"
#include "Utilities/SpheralFunctions.hh"
#include "Utilities/pointDistances.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Construct with the given point and normal.
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE 
inline
GeomPlane<Dimension>::GeomPlane(const typename Dimension::Vector& point,
				const typename Dimension::Vector& normal):
  mPoint(point),
  mNormal(normal.unitVector()) {
  ENSURE(valid());
}

//------------------------------------------------------------------------------
// Negative operator (reverse the sign of the normal).
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE 
inline
GeomPlane<Dimension>
GeomPlane<Dimension>::operator-() const {
  GeomPlane<Dimension> result(*this);
  result.normal(-result.normal());
  return result;
}

//------------------------------------------------------------------------------
// Calculate the signed distance from a point to the plane.
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
typename Dimension::Scalar
GeomPlane<Dimension>::
signedDistance(const Vector& point) const {
  REQUIRE(valid());
  return mNormal.dot(point - mPoint);
}

//------------------------------------------------------------------------------
// Calculate the minimum distance from a point to the plane.
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE
inline
typename Dimension::Scalar
GeomPlane<Dimension>::
minimumDistance(const Vector& point) const {
  REQUIRE(valid());
  return std::abs(mNormal.dot(point - mPoint));
}

//------------------------------------------------------------------------------
// Find the closest point on the plane to a given point.
//------------------------------------------------------------------------------
template<>
SPHERAL_HOST_DEVICE 
inline
Dim<1>::Vector
GeomPlane<Dim<1>>::
closestPointOnPlane(const Vector&) const {
  return mPoint;
}

template<>
SPHERAL_HOST_DEVICE 
inline
Dim<2>::Vector
GeomPlane<Dim<2>>::
closestPointOnPlane(const Dim<2>::Vector& point) const {
  const Vector a0p = point - mPoint;
  const Vector direction(-mNormal.y(), mNormal.x());
  const double s = a0p.dot(direction);
  return mPoint + s*direction;
}

template<>
SPHERAL_HOST_DEVICE 
inline
Dim<3>::Vector
GeomPlane<Dim<3>>::
closestPointOnPlane(const Vector& point) const {
  return Spheral::closestPointOnPlane(point, mPoint, mNormal);
}

//------------------------------------------------------------------------------
// Test whether the given plane is parallel to this one.
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE 
inline
bool
GeomPlane<Dimension>::parallel(const GeomPlane<Dimension>& rhs) const {
  REQUIRE(valid());
  REQUIRE(rhs.valid());
  return fuzzyEqual(fabs(mNormal.dot(rhs.normal())), 1.0);
}

//------------------------------------------------------------------------------
// operator ==
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE 
inline
bool
GeomPlane<Dimension>::operator==(const GeomPlane<Dimension>& rhs) const {
  return (fuzzyEqual(mNormal.dot(rhs.normal()), 1.0) &&
	  fuzzyEqual(minimumDistance(rhs.point()), 0.0));
}

//------------------------------------------------------------------------------
// operator !=
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE 
inline
bool
GeomPlane<Dimension>::operator!=(const GeomPlane<Dimension>& rhs) const {
  return !(*this == rhs);
}

//------------------------------------------------------------------------------
// operator <
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE 
inline
bool
GeomPlane<Dimension>::operator<(const GeomPlane<Dimension>& rhs) const {
  return (mPoint < rhs.mPoint   ? true :
          mNormal < rhs.mNormal ? true:
                                  false);
}

//------------------------------------------------------------------------------
// operator == (Point)
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE 
bool
GeomPlane<Dimension>::operator==(const typename Dimension::Vector& rhs) const {
  return compare(rhs) == 0;
}

//------------------------------------------------------------------------------
// operator != (Point)
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE 
bool
GeomPlane<Dimension>::operator!=(const typename Dimension::Vector& rhs) const {
  return compare(rhs) != 0;
}

//------------------------------------------------------------------------------
// operator > (Point)
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE 
bool
GeomPlane<Dimension>::operator>(const typename Dimension::Vector& rhs) const {
  return compare(rhs) == 1;
}

//------------------------------------------------------------------------------
// operator < (Point)
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE 
bool
GeomPlane<Dimension>::operator<(const typename Dimension::Vector& rhs) const {
  return compare(rhs) == -1;
}

//------------------------------------------------------------------------------
// operator >= (Point)
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE 
bool
GeomPlane<Dimension>::operator>=(const typename Dimension::Vector& rhs) const {
  return compare(rhs) != -1;
}

//------------------------------------------------------------------------------
// operator <= (Point)
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE 
bool
GeomPlane<Dimension>::operator<=(const typename Dimension::Vector& rhs) const {
  return compare(rhs) != 1;
}

//------------------------------------------------------------------------------
// Test if the plane is below, equal to or above the point (-1, 0, 1).
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE 
inline
int
GeomPlane<Dimension>::
compare(const typename Dimension::Vector& rhs) const {
  const auto thpt = mNormal.dot(rhs - mPoint);
  return (fuzzyEqual(thpt, 0.0) ?  0 :
          thpt > 0.0            ? -1 :
                                   1);
}


//------------------------------------------------------------------------------
// Valid test.
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE 
bool
GeomPlane<Dimension>::valid() const {
  return fuzzyEqual(mNormal.magnitude2(), 1.0);
}

//------------------------------------------------------------------------------
// Global comparator functions.
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE 
inline
bool
operator==(const typename Dimension::Vector& vec, const GeomPlane<Dimension>& plane) {
  return plane == vec;
}
  
template<typename Dimension>
SPHERAL_HOST_DEVICE 
inline
bool
operator!=(const typename Dimension::Vector& vec, const GeomPlane<Dimension>& plane) {
  return plane != vec;
}
  
template<typename Dimension>
SPHERAL_HOST_DEVICE 
inline
bool
operator<(const typename Dimension::Vector& vec, const GeomPlane<Dimension>& plane) {
  return plane > vec;
}
  
template<typename Dimension>
SPHERAL_HOST_DEVICE 
inline
bool
operator>(const typename Dimension::Vector& vec, const GeomPlane<Dimension>& plane) {
  return plane < vec;
}
  
template<typename Dimension>
SPHERAL_HOST_DEVICE 
inline
bool
operator<=(const typename Dimension::Vector& vec, const GeomPlane<Dimension>& plane) {
  return plane >= vec;
}
  
template<typename Dimension>
SPHERAL_HOST_DEVICE 
inline
bool
operator>=(const typename Dimension::Vector& vec, const GeomPlane<Dimension>& plane) {
  return plane <= vec;
}

//------------------------------------------------------------------------------
// Input (istream) operator.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
std::istream&
operator>>(std::istream& is, GeomPlane<Dimension>& plane) {
  std::string parenthesis;
  typename Dimension::Vector point;
  typename Dimension::Vector normal;
  is >> parenthesis >> point >> normal >> parenthesis;
  plane.point(point);
  plane.normal(normal);
  return is;
}

//------------------------------------------------------------------------------
// Output (std::ostream) operator.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
std::ostream&
operator<<(std::ostream& os, const GeomPlane<Dimension>& plane) {
  os << "( " << plane.point() << plane.normal() << ")";
  return os;
}

}
