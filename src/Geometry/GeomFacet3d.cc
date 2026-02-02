#include <limits>

#include "GeomVector.hh"
#include "GeomFacet3d.hh"
#include "Utilities/pointDistances.hh"
#include "Utilities/pointInPolygon.hh"
#include "Utilities/lineSegmentIntersections.hh"
#include "Utilities/safeInv.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructors, destructor.
//------------------------------------------------------------------------------
// We really don't want the default constructor, but it's required to have 
// std::vectors of these.
GeomFacet3d::
GeomFacet3d():
  mVerticesPtr(0),
  mPoints(),
  mNormal(1.0, 0.0, 0.0) {
  VERIFY(false);
}

GeomFacet3d::
GeomFacet3d(const std::vector<GeomFacet3d::Vector>& vertices,
            const std::vector<unsigned>& ipoints):
  mVerticesPtr(&vertices),
  mPoints(ipoints),
  mNormal() {
  REQUIRE(mPoints.size() >= 3);
  this->computeNormal();
}

GeomFacet3d::
GeomFacet3d(const GeomFacet3d& rhs):
  mVerticesPtr(rhs.mVerticesPtr),
  mPoints(rhs.mPoints),
  mNormal(rhs.mNormal) {
}

GeomFacet3d&
GeomFacet3d::
operator=(const GeomFacet3d& rhs) {
  if (this != &rhs) {
    mVerticesPtr = rhs.mVerticesPtr;
    mPoints = rhs.mPoints;
    mNormal = rhs.mNormal;
  }
  return *this;
}

GeomFacet3d::
~GeomFacet3d() {
}

//------------------------------------------------------------------------------
// Compare a set of points:
//  1 => all points above.
//  0 => points both above and below (or equal).
// -1 => all points below.
//------------------------------------------------------------------------------
int
GeomFacet3d::
compare(const std::vector<GeomFacet3d::Vector>& points,
        const double tol) const {
  if (points.size() == 0) return 0;
  int result = this->compare(points[0], tol);
  int thpt;
  for (unsigned i = 1; i != points.size(); ++i) {
    thpt = this->compare(points[i], tol);
    if (thpt != result) return 0;
  }
  return result;
}

//------------------------------------------------------------------------------
// Compute the position.
//------------------------------------------------------------------------------
GeomFacet3d::Vector
GeomFacet3d::
position() const {
  const unsigned n = mPoints.size();
  REQUIRE(n >= 3);
  Vector result;
  unsigned i, j;
  double circum = 0.0, dl;
  const auto p0 = point(0);
  for (i = 0; i != n; ++i) {
    j = (i + 1) % n;
    const auto p1 = point(i)-point(0);
    const auto p2 = point(j)-point(0);
    dl = (p1.cross(p2)).magnitude();
    result += (p1+p2) * dl;
    circum += dl;
  }
  result *= safeInvVar(3.0*circum);
  return p0 + result;
}

//------------------------------------------------------------------------------
// Compute the area.
//------------------------------------------------------------------------------
double
GeomFacet3d::
area() const {
  Vector vecsum;
  const Vector cent = this->position();
  unsigned i, j, npts = mPoints.size();
  for (i = 0; i != npts; ++i) {
    j = (i + 1) % npts;
    vecsum += (point(i) - cent).cross(point(j) - cent);
  }
  const double result = 0.5*vecsum.magnitude();
  ENSURE(result >= 0.0);
  return result;
}

//------------------------------------------------------------------------------
// Minimum distance to a point.
//------------------------------------------------------------------------------
double
GeomFacet3d::
distance(const GeomFacet3d::Vector& p) const {
  return (p - this->closestPoint(p)).magnitude();
}

//------------------------------------------------------------------------------
// Find the point on the facet closest to the given point.
//------------------------------------------------------------------------------
GeomFacet3d::Vector
GeomFacet3d::
closestPoint(const GeomFacet3d::Vector& p) const {

  // First check if the closest point on the plane is in the facet.  If so, 
  // that's it!
  const Vector centroid = this->position();
  const Vector nhat = this->normal();
  CHECK(fuzzyEqual(nhat.magnitude2(), 1.0, 1.0e-10));
  Vector result = closestPointOnPlane(p, centroid, nhat);
  if (pointInPolygon(result, *mVerticesPtr, mPoints, nhat)) return result;

  // Nope, so look for the point around the circumference of the facet that is
  // closest.
  unsigned i, j, npts = mPoints.size();
  double r2, minr2 = std::numeric_limits<double>::max();
  Vector potential;
  for (i = 0; i != npts; ++i) {
    j = (i + 1) % npts;
    potential = closestPointOnSegment(p, point(i), point(j));
    r2 = (potential - p).magnitude2();
    if (r2 < minr2) {
      result = potential;
      minr2 = r2;
    }
  }

  // That's it.
  return result;
}

//------------------------------------------------------------------------------
// Return the points for the triangles that decompose this facet.
//------------------------------------------------------------------------------
void
GeomFacet3d::
decompose(std::vector<std::array<Vector, 3>>& subfacets) const {
  const auto numPoints = mPoints.size();
  // switch (numPoints) {
  // We could specialize as below, but that would make the facets inconsistent.
  // case 3:
  //   // Return the input triangle
  //   subfacets = {{point(0), point(1), point(2)}};
  //   break;
  // case 4:
  //   // Split the quadrilateral into two triangles
  //   subfacets = {{point(0), point(1), point(2)},
  //                {point(0), point(2), point(3)}};
  //   break;
  // case 5:
  //   // Split the pentagon into three triangles
  //   subfacets = {{point(0), point(1), point(2)},
  //                {point(0), point(2), point(4)},
  //                {point(2), point(3), point(4)}};
  //   break;
  // case 6:
  //   // Split the hexagon into four triangles
  //   subfacets = {{point(0), point(1), point(2)},
  //                {point(0), point(2), point(3)},
  //                {point(0), point(3), point(5)},
  //                {point(3), point(4), point(5)}};,
  // default:
  const auto centroid = this->position();
  subfacets.resize(numPoints);
  for (auto i = 0u; i < numPoints; ++i) {
    subfacets[i] = {point(i), point((i+1) % numPoints), centroid};
  }
  //   break;
  // }
  
  BEGIN_CONTRACT_SCOPE
  {
    const auto originalArea = this->area();
    auto areasum = 0.;
    CONTRACT_VAR(areasum);
    for (auto& subfacet : subfacets) {
      const auto ab = subfacet[1] - subfacet[0];
      const auto ac = subfacet[2] - subfacet[0];
      const auto subnormal = ab.cross(ac); 
      const auto subarea = 0.5 * subnormal.magnitude();
      CONTRACT_VAR(originalArea);
      CHECK(0 < subarea and subarea < originalArea);
      const auto subnormalUnit = subnormal.unitVector();
      const auto normalUnit = mNormal.unitVector();
      CONTRACT_VAR(subnormalUnit);
      CONTRACT_VAR(normalUnit);
      CHECK2(fuzzyEqual(subnormalUnit(0), normalUnit(0), 1.e-8) &&
             fuzzyEqual(subnormalUnit(1), normalUnit(1), 1.e-8) &&
             fuzzyEqual(subnormalUnit(2), normalUnit(2), 1.e-8), "Normal vector mismatch: " << subnormalUnit << " != " << normalUnit);
      areasum += subarea;
    }
    CHECK(fuzzyEqual(areasum, originalArea));
  }
  END_CONTRACT_SCOPE
}

//------------------------------------------------------------------------------
// Split into triangular sub-facets.
//------------------------------------------------------------------------------
std::vector<GeomFacet3d>
GeomFacet3d::
triangles() const {
  std::vector<GeomFacet3d> result;
  const auto nverts = mPoints.size();
  for (auto k = 1u; k < nverts - 1; ++k) {
    result.emplace_back(*mVerticesPtr, std::vector<unsigned>({mPoints[0], mPoints[k], mPoints[k+1]}));
  }
  return result;
}

//------------------------------------------------------------------------------
// Compute our normal
//------------------------------------------------------------------------------
void
GeomFacet3d::
computeNormal() {
  const auto& vertices = *mVerticesPtr;
  const auto  n = mPoints.size();
  mNormal.Zero();
  for (auto i = 1u; i < n - 1u; ++i) {
    auto j = (i + 1) % n;
    mNormal += (vertices[mPoints[i]] - vertices[mPoints[0]]).cross(vertices[mPoints[j]] - vertices[mPoints[0]]);
  }
  // CHECK2(normal.magnitude2() > 0.0, normal << " " << facetIndices.size());
  mNormal = mNormal.unitVector();
}

//------------------------------------------------------------------------------
// ==
//------------------------------------------------------------------------------
bool
GeomFacet3d::
operator==(const GeomFacet3d& rhs) const {
  return (*mVerticesPtr == *(rhs.mVerticesPtr) and
          mPoints == rhs.mPoints and
          fuzzyEqual(mNormal.dot(rhs.mNormal), 1.0));
}

//------------------------------------------------------------------------------
// !=
//------------------------------------------------------------------------------
bool
GeomFacet3d::
operator!=(const GeomFacet3d& rhs) const {
  return not (*this == rhs);
}

//------------------------------------------------------------------------------
// Output (ostream) operator.
//------------------------------------------------------------------------------
std::ostream&
operator<<(std::ostream& os, const GeomFacet3d& facet) {
  os << "GeomFacet3d( ivertices : ";
  const std::vector<unsigned>& ipoints = facet.ipoints();
  for (unsigned i = 0; i != ipoints.size(); ++i) os << ipoints[i] << " ";
  os << "\n              vertices : ";
  for (unsigned i = 0; i != ipoints.size(); ++i) os << facet.point(i) << " ";
  os << "\n                normal : " << facet.normal() 
     << "\n)";
  return os;
}

}
