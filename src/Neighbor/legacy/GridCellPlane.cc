//---------------------------------Spheral++----------------------------------//
// GridCellPlane -- Geometric Plane Class.
//
// Created by JMO, Thu Feb 24 17:26:32 PST 2000
//----------------------------------------------------------------------------//
#include "GridCellPlane.hh"
#include "GridCellIndex.hh"
#include "Utilities/DBC.hh"
#include "Utilities/SpheralFunctions.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Empty constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
GridCellPlane<Dimension>::GridCellPlane():
  mPoint(),
  mNormal() {
}

//------------------------------------------------------------------------------
// Copy constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
GridCellPlane<Dimension>::GridCellPlane(const GridCellPlane<Dimension>& rhs):
  mPoint(rhs.point()),
  mNormal(rhs.normal()) {
}

//------------------------------------------------------------------------------
// Construct with the given point and normal.
//------------------------------------------------------------------------------
template<typename Dimension>
GridCellPlane<Dimension>::GridCellPlane(const GridCellIndex<Dimension>& point,
                                        const GridCellIndex<Dimension>& normal):
  mPoint(point),
  mNormal(normal) {
  CHECK(mNormal.magnitude2() > 0);
  CHECK(valid());
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
GridCellPlane<Dimension>::~GridCellPlane() {
}

//------------------------------------------------------------------------------
// Assignment.
//------------------------------------------------------------------------------
template<typename Dimension>
GridCellPlane<Dimension>&
GridCellPlane<Dimension>::operator=(const GridCellPlane<Dimension>& rhs) {
  if (this != &rhs) {
    mPoint = rhs.point();
    mNormal = rhs.normal();
  }
  return *this;
}

//------------------------------------------------------------------------------
// Access the point.
//------------------------------------------------------------------------------
template<typename Dimension>
const GridCellIndex<Dimension>&
GridCellPlane<Dimension>::point() const {
  return mPoint;
}

template<typename Dimension>
void
GridCellPlane<Dimension>::setPoint(const GridCellIndex<Dimension>& point) {
  mPoint = point;
}

//------------------------------------------------------------------------------
// Access the normal.
//------------------------------------------------------------------------------
template<typename Dimension>
const GridCellIndex<Dimension>&
GridCellPlane<Dimension>::normal() const {
  return mNormal;
}

template<typename Dimension>
void
GridCellPlane<Dimension>::setNormal(const GridCellIndex<Dimension>& normal) {
  mNormal = normal;
  CHECK(mNormal.magnitude2() > 0);
}

//------------------------------------------------------------------------------
// Calculate the minimum distance from a point to the plane.
//------------------------------------------------------------------------------
template<typename Dimension>
double
GridCellPlane<Dimension>::
minimumDistance(const GridCellIndex<Dimension>& point) const {
  CHECK(valid());
  return abs(normal().dot(point - mPoint))/normal().magnitude();
}

//------------------------------------------------------------------------------
// Test whether the given point is in the plane.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
GridCellPlane<Dimension>::coplanar(const GridCellIndex<Dimension>& otherPoint) const {
  CHECK(valid());
  return normal().dot(otherPoint - point()) == 0;
}

//------------------------------------------------------------------------------
// Test whether the given plane is parallel to this one.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
GridCellPlane<Dimension>::parallel(const GridCellPlane<Dimension>& rhs) const {
  CHECK(valid());
  CHECK(rhs.valid());
  return (abs(mNormal.dot(rhs.normal())) ==
          (int) (mNormal.magnitude()*rhs.normal().magnitude() + 0.5));
}

//------------------------------------------------------------------------------
// operator ==
//------------------------------------------------------------------------------
template<typename Dimension>
bool
GridCellPlane<Dimension>::operator==(const GridCellPlane<Dimension>& rhs) const {
  return (normal() == rhs.normal() &&
          coplanar(rhs.point()));
}

//------------------------------------------------------------------------------
// operator !=
//------------------------------------------------------------------------------
template<typename Dimension>
bool
GridCellPlane<Dimension>::operator!=(const GridCellPlane<Dimension>& rhs) const {
  return !(*this == rhs);
}

//------------------------------------------------------------------------------
// operator > (Point)
//------------------------------------------------------------------------------
template<typename Dimension>
bool
GridCellPlane<Dimension>::operator>(const GridCellIndex<Dimension>& rhs) const {
  return mNormal.dot(rhs - mPoint) < 0;
}

//------------------------------------------------------------------------------
// operator < (Point)
//------------------------------------------------------------------------------
template<typename Dimension>
bool
GridCellPlane<Dimension>::operator<(const GridCellIndex<Dimension>& rhs) const {
  return mNormal.dot(rhs - mPoint) > 0;
}

//------------------------------------------------------------------------------
// operator >= (Point)
//------------------------------------------------------------------------------
template<typename Dimension>
bool
GridCellPlane<Dimension>::operator>=(const GridCellIndex<Dimension>& rhs) const {
  return mNormal.dot(rhs - mPoint) <= 0;
}

//------------------------------------------------------------------------------
// operator <= (Point)
//------------------------------------------------------------------------------
template<typename Dimension>
bool
GridCellPlane<Dimension>::operator<=(const GridCellIndex<Dimension>& rhs) const {
  return mNormal.dot(rhs - mPoint) >= 0;
}

//------------------------------------------------------------------------------
// Valid test.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
GridCellPlane<Dimension>::valid() const {
  return mNormal.magnitude2() > 0;
}

}
