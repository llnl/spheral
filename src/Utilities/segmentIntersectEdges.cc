//------------------------------------------------------------------------------
// segmentIntersectEdges.
//
// Test if a line segment (a0,a1) intersects any of the edges of a faceted
// volume.
//------------------------------------------------------------------------------
#include "segmentIntersectEdges.hh"
#include "lineSegmentIntersections.hh"

using std::vector;
using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {

//------------------------------------------------------------------------------
// 2-D.
//------------------------------------------------------------------------------
bool segmentIntersectEdges(const Dim<2>::Vector& a0,
                           const Dim<2>::Vector& a1,
                           const Dim<2>::FacetedVolume& poly,
                           const double tol) {
  typedef Dim<2>::Vector Vector;

  // Walk the points.
  const std::vector<Vector>& vertices = poly.vertices();
  unsigned i = 0, j;
  const unsigned n = vertices.size();
  Vector intercept1, intercept2;
  char code;
  bool result = false;
  while (i != n and not result) {
    j = (i + 1) % n;
    code = segmentSegmentIntersection(a0, a1, vertices[i], vertices[j], 
                                      intercept1, intercept2, tol);
    result = (code != '0');
    ++i;
  }
      
  // That's it.
  return result;
}

//------------------------------------------------------------------------------
// 3-D.
//------------------------------------------------------------------------------
bool segmentIntersectEdges(const Dim<3>::Vector& a0,
                           const Dim<3>::Vector& a1,
                           const Dim<3>::FacetedVolume& poly,
                           const double tol) {
  typedef Dim<3>::Vector Vector;
  typedef Dim<3>::FacetedVolume FacetedVolume;
  typedef FacetedVolume::Facet Facet;

  // Determine our normalization.
  const double scale = std::max(tol, (poly.xmax() - poly.xmin()).maxAbsElement());
  const double scaleInv = 1.0/scale;

  // Scale the geometry.
  const Vector aa0 = scaleInv*a0;
  const Vector aa1 = scaleInv*a1;
  vector<Vector> vertices = poly.vertices();
  for (unsigned i = 0; i != vertices.size(); ++i) vertices[i] *= scaleInv;

  // Walk the facets.
  const std::vector<Facet>& facets = poly.facets();
  const unsigned nfacets = facets.size();
  unsigned ifacet = 0, i, j, n;
  bool result = false;
  while (ifacet != nfacets and not result) {
    const std::vector<unsigned>& ipts = facets[ifacet].ipoints();
    n = ipts.size();
    i = 0;
    while (i != n and not result) {
      j = (i + 1) % n;
      result = (segmentSegmentDistance(aa0, aa1, vertices[ipts[i]], vertices[ipts[j]]) < tol);
//       cerr << result << " : " << segmentSegmentDistance(aa0, aa1, vertices[ipts[i]], vertices[ipts[j]]) << " : "
//            << a0 << " " << a1 << " : " << scale*vertices[ipts[i]] << " " << scale*vertices[ipts[j]] << endl;
      ++i;
    }
    ++ifacet;
  }
      
  // That's it.
  return result;
}

}
