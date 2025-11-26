//------------------------------------------------------------------------------
// pointInPolyhedron
//
// Test if a given point is in a polyhedron or not.
//------------------------------------------------------------------------------
#include "pointInPolyhedron.hh"
#include "pointInPolygon.hh"
#include "testBoxIntersection.hh"
#include "pointOnPolyhedron.hh"
#include "segmentIntersectEdges.hh"
#include "lineSegmentIntersections.hh"
#include "rotationMatrix.hh"
#include "Geometry/Dimension.hh"

#include <random>
#include <algorithm>
#include <iostream>
using std::vector;
using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {

//------------------------------------------------------------------------------
// The method in question -- test if a point is on the interior of a polyhedron.
//------------------------------------------------------------------------------
bool pointInPolyhedron(const Dim<3>::Vector& p,
                       const Dim<3>::FacetedVolume& polyhedron,
                       const bool countBoundary,
                       const double tol) {
  typedef Dim<3>::Vector Vector;
  typedef Dim<2>::Vector Vector2d;

  const Vector rayhat(0, 0, 1);    // Project straight up in z.
  const auto pend = p + Vector(0, 0, 1e60);
  const auto& vertices = polyhedron.vertices();
  const auto& facets = polyhedron.facets();
  const auto  nfacets = facets.size();
  const Vector2d pxy(p.x(), p.y()),
                 pzx(p.z(), p.x()),
                 pyz(p.y(), p.z());

  // First reject any points outside the bounding box.
  if (not testPointInBox(p, polyhedron.xmin(), polyhedron.xmax(), tol)) return false;

  // Are we checking the surface?
  if (countBoundary and polyhedron.distance(p) < tol) return true;
  // if (countBoundary and pointOnPolyhedron(p, polyhedron, tol)) return true;

  // Keep count of how many facets we intersect.
  int numxy = 0;

  // Walk the facets.
#pragma omp parallel
  {
    int numxy_thread = 0;
#pragma omp for
    for (auto ifacet = 0u; ifacet < nfacets; ++ifacet) {
      const auto& facet = facets[ifacet];
      const auto& ipts = facet.ipoints();
      const auto  npts = ipts.size();
      const auto& fhat = facet.normal();
      const auto  fdotr = fhat.dot(rayhat);

      // Check if the facet is parallel to ray.
      if (std::abs(fdotr) <= tol) {

        // If they are, we only care if the ray is in the plane of the facet.
        if (std::abs((vertices[ipts[0]] - p).dot(fhat)) <= tol) {

          // The ray and plane are parallel, and the ray is in the plane.
          // Check if the ray crosses any of the facet lines -- if so, counts as an intersection.
          size_t i = 0;
          while (i < npts and
                 segmentSegmentDistance(vertices[ipts[i]], vertices[ipts[(i + 1) % npts]], p, pend) > tol) ++i;
          if (i < npts) ++numxy_thread;
        }

      } else {

        // Find the intersection point of the line corresponding to our ray and the plane of the facet.
        // What we're solving for here is parametric s where the line is represented as
        //   pline = p + s*rayhat
        const auto si = fhat.dot(vertices[ipts[0]] - p)*safeInvVar(fdotr);
        CHECK2(std::abs((p + si*rayhat - vertices[ipts[0]]).dot(fhat)) < tol*std::max(1.0, 1e-3*(p - vertices[ipts[0]]).magnitude()),
               "Point in plane?  " << p << " " << vertices[ipts[0]] << " " << rayhat << " " << si << " : " << std::abs((p + si*rayhat - vertices[ipts[0]]).dot(fhat)) << " !< " << tol);
      
        // We only need to check this facet if si is positive, indicating the the intersection
        // point is in the positive ray direction, i.e., above p in z.
        if (si > 0.0) {

          // Copy the vertices as XY pairs.
          vector<Vector2d> vxy;
          for (size_t i = 0; i < npts; ++i) {
            vxy.push_back(Vector2d(vertices[ipts[i]].x(), vertices[ipts[i]].y()));
          }
          vxy.push_back(vxy[0]);

          // Does a ray straight up the +z direction intersect this projected polygon?
          if (pointInPolygon(pxy, vxy)) {
            // cerr << "Intersect facet: " << si << " " << pxy << endl;
            // for (size_t i = 0; i < npts; ++i) cerr << " " << vxy[i];
            // cerr << endl;
            ++numxy_thread;
          }
        }
      }
    }
#pragma omp critical
    {
      numxy += numxy_thread;
    }
  }

  // If the ray passed through an even number of intersections the point is external to the 
  // polyhedron.
  return ((numxy % 2) == 1);
}

}
