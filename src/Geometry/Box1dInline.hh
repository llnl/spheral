#include <numeric>
#include "Utilities/DBC.hh"
#include "GeomVector.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructors, destructors.
//------------------------------------------------------------------------------
inline
Box1d::
Box1d():
  mCenter(),
  mExtent(0.0),
  mVertices() {
  mVertices.push_back(Vector());
  mVertices.push_back(Vector());
  REQUIRE(mExtent >= 0.0);
  REQUIRE(mVertices.size() == 2);
}

inline
Box1d::
Box1d(const std::vector<Box1d::Vector>& points):
  mCenter(),
  mExtent(0.0) {
  if (points.size() > 0) {
    double xmin = DBL_MAX;
    double xmax = -DBL_MAX;
    for (std::vector<Vector>::const_iterator itr = points.begin();
         itr != points.end();
         ++itr) {
      xmin = std::min(xmin, itr->x());
      xmax = std::max(xmax, itr->x());
    }
    mCenter = 0.5*(xmin + xmax);
    mExtent = 0.5*(xmax - xmin);
  }
  mVertices.push_back(mCenter - Vector(mExtent));
  mVertices.push_back(mCenter + Vector(mExtent));
}

inline
Box1d::
Box1d(const std::vector<Box1d::Vector>& points,
      const std::vector<std::vector<unsigned> >& facetIndices):
  mCenter(),
  mExtent(0.0) {
  REQUIRE(points.size() >= 2);
  REQUIRE(facetIndices.size() == 2);
  REQUIRE(facetIndices[0].size() == 1);
  REQUIRE(facetIndices[1].size() == 1);
  REQUIRE(facetIndices[0][0] < points.size());
  REQUIRE(facetIndices[1][0] < points.size());

  const unsigned i = facetIndices[0][0],
                 j = facetIndices[1][0];
  const double xmin = std::min(points[i].x(), points[j].x()),
               xmax = std::max(points[i].x(), points[j].x());
  mCenter = 0.5*(xmin + xmax);
  mExtent = 0.5*(xmax - xmin);
  mVertices.push_back(mCenter - Vector(mExtent));
  mVertices.push_back(mCenter + Vector(mExtent));
}

inline
Box1d::
Box1d(const GeomVector<1>& center,
      const double extent):
  mCenter(center),
  mExtent(extent) {
  REQUIRE(mExtent >= 0.0);
  mVertices.push_back(mCenter - Vector(mExtent));
  mVertices.push_back(mCenter + Vector(mExtent));
}

inline
Box1d::
Box1d(const Box1d& rhs):
  mCenter(rhs.mCenter),
  mExtent(rhs.mExtent),
  mVertices(rhs.mVertices) {
  REQUIRE(mExtent >= 0.0);
  REQUIRE(mVertices.size() == 2);
}

inline
Box1d&
Box1d::
operator=(const Box1d& rhs) {
  if (this != &rhs) {
    mCenter = rhs.mCenter;
    mExtent = rhs.mExtent;
    mVertices = rhs.mVertices;
  }
  REQUIRE(mExtent >= 0.0);
  return *this;
}

inline
Box1d::
~Box1d() {
}

//------------------------------------------------------------------------------
// Test if the given point is in the box.
//------------------------------------------------------------------------------
inline
bool
Box1d::
contains(const Box1d::Vector& point,
         const bool countBoundary,
         const double /*tol*/) const {
  if (countBoundary) {
    return std::abs(point.x() - mCenter.x()) <= mExtent;
  } else {
    return std::abs(point.x() - mCenter.x()) <  mExtent;
  }
}

//------------------------------------------------------------------------------
// Test if the given point is in the box.
//------------------------------------------------------------------------------
inline
bool
Box1d::
convexContains(const Box1d::Vector& point,
               const bool countBoundary,
               const double tol) const {
  return this->contains(point, countBoundary, tol);
}

//------------------------------------------------------------------------------
// Test if the boxes intersect.
//------------------------------------------------------------------------------
inline
bool
Box1d::
intersect(const Box1d& rhs) const {
  if (mCenter.x() + mExtent < rhs.mCenter.x() - rhs.mExtent) return false;
  if (mCenter.x() - mExtent > rhs.mCenter.x() + rhs.mExtent) return false;
  return true;
}

//------------------------------------------------------------------------------
// Test if the boxes intersect (convex).
//------------------------------------------------------------------------------
inline
bool
Box1d::
convexIntersect(const Box1d& rhs) const {
  return this->intersect(rhs);
}

//------------------------------------------------------------------------------
// Test if the boxes intersect.
// rhs box reprsented as (xmin, xmax) bounding coordinates.
//------------------------------------------------------------------------------
inline
bool
Box1d::
intersect(const std::pair<Vector, Vector>& rhs) const {
  if (mCenter.x() + mExtent < rhs.first.x()) return false;
  if (mCenter.x() - mExtent > rhs.second.x()) return false;
  return true;
}

//------------------------------------------------------------------------------
// Test if we intersect a line segment (interior counts as intersection).
//------------------------------------------------------------------------------
inline
bool
Box1d::
intersect(const Vector& s0, const Vector& s1) const {
  const auto smin = std::min(s0.x(), s1.x());
  const auto smax = std::max(s0.x(), s1.x());
  if ((smin > mCenter.x() + mExtent) or
      (smax < mCenter.x() - mExtent)) return false;
  return true;
}

//------------------------------------------------------------------------------
// Center attribute.
//------------------------------------------------------------------------------
inline
GeomVector<1>&
Box1d::
center() {
  return mCenter;
}

inline
const GeomVector<1>&
Box1d::
center() const {
  return mCenter;
}

inline
void
Box1d::
center(const GeomVector<1>& val) {
  mCenter = val;
  mVertices.clear();
  mVertices.push_back(mCenter - Vector(mExtent));
  mVertices.push_back(mCenter + Vector(mExtent));
}

//------------------------------------------------------------------------------
// Extent attribute.
//------------------------------------------------------------------------------
inline
double&
Box1d::
extent() {
  return mExtent;
}

inline
double
Box1d::
extent() const {
  return mExtent;
}

inline
void
Box1d::
extent(double val) {
  mExtent = val;
  mVertices.clear();
  mVertices.push_back(mCenter - Vector(mExtent));
  mVertices.push_back(mCenter + Vector(mExtent));
}

//------------------------------------------------------------------------------
// xmin
//------------------------------------------------------------------------------
inline
const GeomVector<1>&
Box1d::
xmin() const {
  return mVertices[0];
}

//------------------------------------------------------------------------------
// xmax
//------------------------------------------------------------------------------
inline
const GeomVector<1>&
Box1d::
xmax() const {
  return mVertices[1];
}

//------------------------------------------------------------------------------
// The vertices.
//------------------------------------------------------------------------------
inline
const std::vector<GeomVector<1> >&
Box1d::
vertices() const {
  return mVertices;
}

//------------------------------------------------------------------------------
// The facets.
//------------------------------------------------------------------------------
inline
const std::vector<GeomFacet1d>&
Box1d::
facets() const {
  mFacets = {Facet(mVertices[0], Vector(-1.0)),
             Facet(mVertices[1], Vector(1.0))};
  return mFacets;
}

//------------------------------------------------------------------------------
// The vertex indices for each "facet".
//------------------------------------------------------------------------------
inline
std::vector<std::vector<unsigned> >
Box1d::
facetVertices() const {
  std::vector<std::vector<unsigned> > result(2);
  result[0] = std::vector<unsigned>(1, 0);
  result[1] = std::vector<unsigned>(1, 1);
  return result;
}

//------------------------------------------------------------------------------
// Useful facet properties.
//------------------------------------------------------------------------------
inline
double
Box1d::
facetArea(const unsigned facetID) const {
  CONTRACT_VAR(facetID);
  REQUIRE(facetID < 2);
  return 1.0;
}

inline
Box1d::Vector
Box1d::
facetAreaNormal(const unsigned facetID) const {
  REQUIRE(facetID < 2);
  if (facetID == 0) {
    return Vector(-1.0);
  } else {
    return Vector( 1.0);
  }
}

//------------------------------------------------------------------------------
// Compute the volume.
//------------------------------------------------------------------------------
inline
double
Box1d::
volume() const {
  return 2.0*mExtent;
}

//------------------------------------------------------------------------------
// Distance to a point.
//------------------------------------------------------------------------------
inline
double
Box1d::
distance(const GeomVector<1>& p) const {
  const double r1 = std::abs(p.x() - mCenter.x());
  return std::abs(r1 - mExtent);
}

//------------------------------------------------------------------------------
// Find the point in the box closest to the given point.
//------------------------------------------------------------------------------
inline
GeomVector<1>
Box1d::
closestPoint(const GeomVector<1>& p) const {
  if (p.x() > mCenter.x()) {
    return GeomVector<1>(p.x() + mExtent);
  } else {
    return GeomVector<1>(p.x() - mExtent);
  }
}

//------------------------------------------------------------------------------
// In 1d, there is no decomposition to be done
//------------------------------------------------------------------------------
inline
void
Box1d::
decompose(std::vector<Box1d>& subcells) const {
  subcells = {*this};
}

//------------------------------------------------------------------------------
// += Vector, shift box in space
//------------------------------------------------------------------------------
inline
Box1d&
Box1d::
operator+=(const Box1d::Vector& rhs) {
  mCenter += rhs;
  mVertices[0] += rhs;
  mVertices[1] += rhs;
  return *this;
}

//------------------------------------------------------------------------------
// -= Vector, shift box in space
//------------------------------------------------------------------------------
inline
Box1d&
Box1d::
operator-=(const Box1d::Vector& rhs) {
  (*this) += -rhs;
  return *this;
}

//------------------------------------------------------------------------------
// + Vector, return shifted box in space
//------------------------------------------------------------------------------
inline
Box1d
Box1d::
operator+(const Box1d::Vector& rhs) const {
  Box1d result(*this);
  result += rhs;
  return result;
}

//------------------------------------------------------------------------------
// - Vector, return shifted box in space
//------------------------------------------------------------------------------
inline
Box1d
Box1d::
operator-(const Box1d::Vector& rhs) const {
  return (*this) + (-rhs);
}

//------------------------------------------------------------------------------
// *= Scalar, scale box
//------------------------------------------------------------------------------
inline
Box1d&
Box1d::
operator*=(const double rhs) {
  mCenter *= rhs;
  mExtent *= rhs;
  mVertices[0] *= rhs;
  mVertices[1] *= rhs;
  return *this;
}

//------------------------------------------------------------------------------
// /= Scalar, scale box
//------------------------------------------------------------------------------
inline
Box1d&
Box1d::
operator/=(const double rhs) {
  (*this) *= 1.0/rhs;
  return *this;
}

//------------------------------------------------------------------------------
// * Scalar, scale box
//------------------------------------------------------------------------------
inline
Box1d
Box1d::
operator*(const double rhs) const {
  Box1d result(*this);
  result *= rhs;
  return result;
}

//------------------------------------------------------------------------------
// / Scalar, scale box
//------------------------------------------------------------------------------
inline
Box1d
Box1d::
operator/(const double rhs) const {
  Box1d result(*this);
  result /= rhs;
  return result;
}

//------------------------------------------------------------------------------
// ==
//------------------------------------------------------------------------------
inline
bool
Box1d::
operator==(const Box1d& rhs) const {
  return mCenter == rhs.mCenter and mExtent == rhs.mExtent;
}

//------------------------------------------------------------------------------
// !=
//------------------------------------------------------------------------------
inline
bool
Box1d::
operator!=(const Box1d& rhs) const {
  return not (*this == rhs);
}

//------------------------------------------------------------------------------
// ostream operator.
//------------------------------------------------------------------------------
inline
std::ostream& operator<<(std::ostream& os, const Box1d& box) {
  os << "Box(" << box.xmin().x() << " " << box.xmax().x() << ")";
  return os;
}

}
