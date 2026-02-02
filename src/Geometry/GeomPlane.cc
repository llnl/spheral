//---------------------------------Spheral++----------------------------------//
// GeomPlane -- Geometric Plane Class.
//
// Created by JMO, Thu Feb 24 17:26:32 PST 2000
//----------------------------------------------------------------------------//
#include "GeomPlane.hh"
#include "Utilities/DBC.hh"
#include "Utilities/SpheralFunctions.hh"
#include "Utilities/pointDistances.hh"
#include "Dimension.hh"

#include <algorithm>
using std::vector;

namespace Spheral {

//------------------------------------------------------------------------------
// Empty constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
GeomPlane<Dimension>::GeomPlane():
  mPoint(),
  mNormal() {
}

//------------------------------------------------------------------------------
// Copy constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
GeomPlane<Dimension>::GeomPlane(const GeomPlane<Dimension>& rhs):
  mPoint(rhs.point()),
  mNormal(rhs.normal()) {
}

//------------------------------------------------------------------------------
// Construct with the given point and normal.
//------------------------------------------------------------------------------
template<typename Dimension>
GeomPlane<Dimension>::GeomPlane(const typename Dimension::Vector& point,
				const typename Dimension::Vector& normal):
  mPoint(point),
  mNormal(normal.unitVector()) {
  CHECK(valid());
}

//------------------------------------------------------------------------------
// Construct as the best fit to a cloud of points.
//------------------------------------------------------------------------------
// 1D
template<>
GeomPlane<Dim<1> >::GeomPlane(const vector<Dim<1>::Vector>& points):
  mPoint(std::accumulate(points.begin(), points.end(), Dim<1>::Vector::zero())/std::max(size_t(1), points.size())),
  mNormal(Dim<1>::Vector::one()) {
  REQUIRE2(points.size() > 0, "GeomPlane(points) ERROR: need at least 1 point.");
}

// 2D
template<>
GeomPlane<Dim<2> >::GeomPlane(const vector<Dim<2>::Vector>& points):
  mPoint(std::accumulate(points.begin(), points.end(), Dim<2>::Vector::zero())/std::max(size_t(1), points.size())),
  mNormal() {
  const unsigned n = points.size();
  CONTRACT_VAR(n);
  REQUIRE2(n >= 2, "GeomPlane(points) ERROR: need at least 2 points.");
  double xy = 0.0, xx = 0.0;
  for (const auto& v: points) {
    const Vector d = v - mPoint;
    xy += d.x()*d.y();
    xx += d.x()*d.x();
  }
  const double m = xy*safeInvVar(xx);
  mNormal = Vector(m, -1.0).unitVector();
}

// 3D
template<>
GeomPlane<Dim<3> >::GeomPlane(const vector<Dim<3>::Vector>& points):
  mPoint(std::accumulate(points.begin(), points.end(), Dim<3>::Vector::zero())/std::max(size_t(1), points.size())),
  mNormal() {
  const unsigned n = points.size();
  CONTRACT_VAR(n);
  REQUIRE2(n >= 3, "GeomPlane(points) ERROR: need at least 3 points.");
  // Calc full 3x3 covariance matrix, excluding symmetries:
  SymTensor A;
  for (const auto& p: points) {
    const auto r = p - mPoint;
    A(0,0) += r.x() * r.x();
    A(0,1) += r.x() * r.y();
    A(0,2) += r.x() * r.z();
    A(1,1) += r.y() * r.y();
    A(1,2) += r.y() * r.z();
    A(2,2) += r.z() * r.z();
  }

  const double det_x = A.yy()*A.zz() - A.yz()*A.yz();
  const double det_y = A.xx()*A.zz() - A.xz()*A.xz();
  const double det_z = A.xx()*A.yy() - A.xy()*A.xy();

  const double det_max = std::max(det_x, std::max(det_y, det_z));
  CHECK2(det_max > 0.0, "The points don't span a plane");

  // Pick path with best conditioning:
  if (det_max == det_x) {
    const double a = (A.xz()*A.yz() - A.xy()*A.zz()) / det_x;
    const double b = (A.xy()*A.yz() - A.xz()*A.yy()) / det_x;
    mNormal = Vector(1.0, a, b).unitVector();
  } else if (det_max == det_y) {
    const double a = (A.yz()*A.xz() - A.xy()*A.zz()) / det_y;
    const double b = (A.xy()*A.xz() - A.yz()*A.xx()) / det_y;
    mNormal = Vector(a, 1.0, b).unitVector();
  } else {
    const double a = (A.yz()*A.xy() - A.xz()*A.yy()) / det_z;
    const double b = (A.xz()*A.xy() - A.yz()*A.xx()) / det_z;
    mNormal = Vector(a, b, 1.0);
  }
}


// // 3D
// template<>
// GeomPlane<Dim<3> >::GeomPlane(const vector<Dim<3>::Vector>& points):
//   mPoint(std::accumulate(points.begin(), points.end(), Dim<3>::Vector::zero())/std::max(size_t(1), points.size())),
//   mNormal() {
//   const unsigned n = points.size();
//   REQUIRE2(n >= 3, "GeomPlane(points) ERROR: need at least 3 points.");
//   SymTensor A(0.0, 0.0, 0.0,
//               0.0, 0.0, 0.0,
//               0.0, 0.0, n);
//   Vector b;
//   for (const auto& v: points) {
//     const Vector d = v - mPoint;
//     A(0,0) += d.x()*d.x();
//     A(0,1) += d.x()*d.y();
//     A(0,2) += d.x();
//     A(1,1) += d.y()*d.y();
//     A(1,2) += d.y();
//     b(0) += d.x()*d.z();
//     b(1) += d.y()*d.z();
//     b(2) += d.z();
//   }
//   cerr << "Points: ";
//   std::copy(points.begin(), points.end(), ostream_iterator<Vector>(cerr, " "));
//   cerr << endl
//        << "mPoint: " << mPoint << endl
//        << "A: " << A << " det " << A.Determinant() << endl
//        << "b: " << b << endl;
//   VERIFY2(abs(A.Determinant()) > 1e-20, "GeomPlane(points) ERROR: input points collinear.");
//   mNormal = (A.Inverse() * b).unitVector();     // I know, I should use a linear solve instead, but this is convenient.
// }

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
GeomPlane<Dimension>::~GeomPlane() {
}

//------------------------------------------------------------------------------
// Assignment.
//------------------------------------------------------------------------------
template<typename Dimension>
GeomPlane<Dimension>&
GeomPlane<Dimension>::operator=(const GeomPlane<Dimension>& rhs) {
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
const typename Dimension::Vector&
GeomPlane<Dimension>::point() const {
  return mPoint;
}

template<typename Dimension>
void
GeomPlane<Dimension>::point(const typename Dimension::Vector& point) {
  mPoint = point;
}

//------------------------------------------------------------------------------
// Access the normal.
//------------------------------------------------------------------------------
template<typename Dimension>
const typename Dimension::Vector&
GeomPlane<Dimension>::normal() const {
  return mNormal;
}

template<typename Dimension>
void
GeomPlane<Dimension>::normal(const typename Dimension::Vector& normal) {
  mNormal = normal.unitVector();
}

//------------------------------------------------------------------------------
// Negative operator (reverse the sign of the normal).
//------------------------------------------------------------------------------
template<typename Dimension>
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
double
GeomPlane<Dimension>::
signedDistance(const typename Dimension::Vector& point) const {
  CHECK(valid());
  return mNormal.dot(point - mPoint);
}

//------------------------------------------------------------------------------
// Calculate the minimum distance from a point to the plane.
//------------------------------------------------------------------------------
template<typename Dimension>
double
GeomPlane<Dimension>::
minimumDistance(const typename Dimension::Vector& point) const {
  CHECK(valid());
  return std::abs(mNormal.dot(point - mPoint));
}

//------------------------------------------------------------------------------
// Find the closest point on the plane to a given point.
//------------------------------------------------------------------------------
template<>
Dim<1>::Vector
GeomPlane<Dim<1> >::
closestPointOnPlane(const Dim<1>::Vector&) const {
  return mPoint;
}

template<>
Dim<2>::Vector
GeomPlane<Dim<2> >::
closestPointOnPlane(const Dim<2>::Vector& point) const {
  const Vector a0p = point - mPoint;
  const Vector direction(-mNormal.y(), mNormal.x());
  const double s = a0p.dot(direction);
  return mPoint + s*direction;
}

template<>
Dim<3>::Vector
GeomPlane<Dim<3> >::
closestPointOnPlane(const Dim<3>::Vector& point) const {
  return Spheral::closestPointOnPlane(point, mPoint, mNormal);
}

//------------------------------------------------------------------------------
// Test whether the given plane is parallel to this one.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
GeomPlane<Dimension>::parallel(const GeomPlane<Dimension>& rhs) const {
  CHECK(valid());
  CHECK(rhs.valid());
  return fuzzyEqual(fabs(mNormal.dot(rhs.normal())), 1.0);
}

//------------------------------------------------------------------------------
// operator ==
//------------------------------------------------------------------------------
template<typename Dimension>
bool
GeomPlane<Dimension>::operator==(const GeomPlane<Dimension>& rhs) const {
  return (fuzzyEqual(mNormal.dot(rhs.normal()), 1.0) &&
	  fuzzyEqual(minimumDistance(rhs.point()), 0.0));
}

//------------------------------------------------------------------------------
// operator !=
//------------------------------------------------------------------------------
template<typename Dimension>
bool
GeomPlane<Dimension>::operator!=(const GeomPlane<Dimension>& rhs) const {
  return !(*this == rhs);
}

//------------------------------------------------------------------------------
// operator <
//------------------------------------------------------------------------------
template<typename Dimension>
bool
GeomPlane<Dimension>::operator<(const GeomPlane<Dimension>& rhs) const {
  return (mPoint < rhs.mPoint   ? true :
          mNormal < rhs.mNormal ? true:
                                  false);
}

//------------------------------------------------------------------------------
// Test if the plane is below, equal to or above the point (-1, 0, 1).
//------------------------------------------------------------------------------
template<typename Dimension>
int
GeomPlane<Dimension>::
compare(const typename Dimension::Vector& rhs) const {
  const double thpt = mNormal.dot(rhs - mPoint);
  return (fuzzyEqual(thpt, 0.0) ?  0 :
          thpt > 0.0            ? -1 :
                                   1);
}


//------------------------------------------------------------------------------
// operator == (Point)
//------------------------------------------------------------------------------
template<typename Dimension>
bool
GeomPlane<Dimension>::operator==(const typename Dimension::Vector& rhs) const {
  return compare(rhs) == 0;
}

//------------------------------------------------------------------------------
// operator != (Point)
//------------------------------------------------------------------------------
template<typename Dimension>
bool
GeomPlane<Dimension>::operator!=(const typename Dimension::Vector& rhs) const {
  return compare(rhs) != 0;
}

//------------------------------------------------------------------------------
// operator > (Point)
//------------------------------------------------------------------------------
template<typename Dimension>
bool
GeomPlane<Dimension>::operator>(const typename Dimension::Vector& rhs) const {
  return compare(rhs) == 1;
}

//------------------------------------------------------------------------------
// operator < (Point)
//------------------------------------------------------------------------------
template<typename Dimension>
bool
GeomPlane<Dimension>::operator<(const typename Dimension::Vector& rhs) const {
  return compare(rhs) == -1;
}

//------------------------------------------------------------------------------
// operator >= (Point)
//------------------------------------------------------------------------------
template<typename Dimension>
bool
GeomPlane<Dimension>::operator>=(const typename Dimension::Vector& rhs) const {
  return compare(rhs) != -1;
}

//------------------------------------------------------------------------------
// operator <= (Point)
//------------------------------------------------------------------------------
template<typename Dimension>
bool
GeomPlane<Dimension>::operator<=(const typename Dimension::Vector& rhs) const {
  return compare(rhs) != 1;
}

//------------------------------------------------------------------------------
// Valid test.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
GeomPlane<Dimension>::valid() const {
  return fuzzyEqual(mNormal.magnitude2(), 1.0);
}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
template class GeomPlane< Dim<1> >;
template class GeomPlane< Dim<2> >;
template class GeomPlane< Dim<3> >;

}
