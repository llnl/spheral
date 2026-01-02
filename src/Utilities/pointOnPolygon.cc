//------------------------------------------------------------------------------
// pointOnPolygon
//
// Test if a given point is on the boundary of a polygon.
//------------------------------------------------------------------------------
#include "pointOnPolygon.hh"
#include "lineSegmentIntersections.hh"

using std::vector;
using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {

//------------------------------------------------------------------------------
// 2D: Check a polygon
//------------------------------------------------------------------------------
bool pointOnPolygon(const Dim<2>::Vector& p,
                    const Dim<2>::FacetedVolume& polygon,
                    const double tol) {
  const auto& facets = polygon.facets();
  for (const auto& facet: facets) {
    if (between(facet.point1(), facet.point2(), p, tol)) return true;
  }
  return false;
}

//------------------------------------------------------------------------------
// 2D: Check a polygon as a set of vertices
//------------------------------------------------------------------------------
bool pointOnPolygon(const Dim<2>::Vector& p,
                    const vector<Dim<2>::Vector>& vertices,
                    const vector<unsigned>& ipoints,
                    const double tol) {
  const unsigned n = ipoints.size();
  for (unsigned i = 0; i < n; ++i) {
    const unsigned j = (i + 1) % n;
    if (between(vertices[ipoints[i]], vertices[ipoints[j]], p, tol)) return true;
  }
  return false;
}

//------------------------------------------------------------------------------
// 3D: Check a polygon as a set of coplanar vertices
//------------------------------------------------------------------------------
bool pointOnPolygon(const Dim<3>::Vector& p,
                    const vector<Dim<3>::Vector>& vertices,
                    const vector<unsigned>& ipoints,
                    const double tol) {
  const unsigned n = ipoints.size();
  for (unsigned i = 0; i < n; ++i) {
    const unsigned j = (i + 1) % n;
    if (between(vertices[ipoints[i]], vertices[ipoints[j]], p, tol)) return true;
  }
  return false;
}

}
