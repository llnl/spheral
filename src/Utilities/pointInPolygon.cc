//------------------------------------------------------------------------------
// pointInPolygon
//
// Test if a given point is in a polygon or not.  This is based on code I found
// at
//
// http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
//
// which I have modified here somewhat to use Spheral structures and such.
//------------------------------------------------------------------------------
#include "pointInPolygon.hh"
#include "testBoxIntersection.hh"
#include "pointOnPolygon.hh"
#include "pointDistances.hh"
#include "lineSegmentIntersections.hh"

#include <limits>
using std::vector;
using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {


//------------------------------------------------------------------------------
// Work with a closed polygon defined by it's vertices.
// We assume the caller has ordered the vertices (clockwise or 
// counter-clockwise), and should have already performed the box exclusion
// test for efficiency.
//------------------------------------------------------------------------------
bool pointInPolygon(const Dim<2>::Vector& p,
                    const std::vector<Dim<2>::Vector>& vertices) {

  // Now we do the test of casting a semi-infinite ray in the x direction from 
  // the point and counting intersections with the polygon.
  // Note this bit of code will reject points lying exactly on the boundary
  // of a polygon.
  unsigned i, j;
  bool result = false;
  const unsigned nvertices = vertices.size();
  const double px = p.x();
  const double py = p.y();
  for (i = 0, j = nvertices - 1; i < nvertices; j = i++) {
    if ( ((vertices[i].y() > py) != (vertices[j].y() > py)) &&
	 (px < (vertices[j].x() - vertices[i].x()) * (py - vertices[i].y()) / (vertices[j].y() - vertices[i].y()) + vertices[i].x()) )
      result = not result;
  }

  // That's it.
  return result;
}

//------------------------------------------------------------------------------
// Test a polygon.
//------------------------------------------------------------------------------
bool pointInPolygon(const Dim<2>::Vector& p,
                    const Dim<2>::FacetedVolume& polygon,
                    const bool countBoundary,
                    const double tol) {

  // Do the quick box rejection test.
  if (not testPointInBox(p, polygon.xmin(), polygon.xmax(), tol)) return false;

  // Check if the point is on the boundary (within tolerance).
  if (pointOnPolygon(p, polygon, tol)) return countBoundary;

  typedef Dim<2>::Vector Vector;
  typedef Dim<2>::FacetedVolume::Facet Facet;
  const vector<Vector>& pverts = polygon.vertices();
  const vector<Facet>& facets = polygon.facets();
  const unsigned nfacets = facets.size();
  CHECK(nfacets >= 3);

  // If there's just one loop it's an easy check.
  if (facets.back().ipoint2() == 0) {
    vector<Vector> verts(pverts.begin(), pverts.end());
    verts.push_back(verts[0]);
    return pointInPolygon(p, verts);
  }

  // If there are multiple loops, we can still do them in one go, but we have to
  // insert (0,0) coordinates between loops.  See the discussion at the above
  // website source for all this for more information.
  vector<Vector> verts(1, Vector::zero());
  verts.reserve(int(1.05*pverts.size()));
  unsigned ifacet = 0;
  while (ifacet < nfacets) {
    const unsigned istart = facets[ifacet].ipoint1();
    verts.push_back(pverts[istart]);
    while (ifacet < nfacets and facets[ifacet].ipoint2() != istart) {
      verts.push_back(facets[ifacet].point2());
      ++ifacet;
    }
    verts.push_back(pverts[istart]);
    verts.push_back(Vector::zero());
    ++ifacet;
  }
  return pointInPolygon(p, verts);
}

//------------------------------------------------------------------------------
// Test a polygon in 3-D.
//------------------------------------------------------------------------------
bool pointInPolygon(const Dim<3>::Vector& p,
                    const vector<Dim<3>::Vector>& vertices,
                    const Dim<3>::Vector& normal) {
  typedef Dim<3>::Vector Vector;

  // Prerequisites.
  const unsigned npts = vertices.size();
  unsigned i, j, k;
  BEGIN_CONTRACT_SCOPE
  {
    REQUIRE(vertices.size() > 2);
    const double normmag = normal.magnitude();
    Vector normi;
    for (i = 1; i != npts; ++i) {
      j = (i + 1) % npts;
      k = (i + 2) % npts;
      normi = (vertices[j] - vertices[i]).cross(vertices[k] - vertices[i]);
      CONTRACT_VAR(normmag);
      REQUIRE(fuzzyEqual(std::abs(normi.dot(normal)), normmag*normi.magnitude(), 1.0e-10));
    }
    REQUIRE(fuzzyEqual(pointPlaneDistance(p, vertices[0], normal.unitVector()), 0.0, 1.0e-10));
  }
  END_CONTRACT_SCOPE

  bool result = false;

  const double px = p.x();
  const double py = p.y();
  const double pz = p.z();

  double fxmin =  1e50, fymin =  1e50, fzmin =  1e50;
  double fxmax = -1e50, fymax = -1e50, fzmax = -1e50;
  for (i = 0; i != npts; ++i) {
    fxmin = std::min(fxmin, vertices[i].x());
    fymin = std::min(fymin, vertices[i].y());
    fzmin = std::min(fzmin, vertices[i].z());
    fxmax = std::max(fxmax, vertices[i].x());
    fymax = std::max(fymax, vertices[i].y());
    fzmax = std::max(fzmax, vertices[i].z());
  }
  const auto fuzz = 1.0e-5*std::max(fxmax - fxmin, std::max(fymax - fymin, fzmax - fzmin));
  fxmin -= fuzz;
  fymin -= fuzz;
  fzmin -= fuzz;
  fxmax += fuzz;
  fymax += fuzz;
  fzmax += fuzz;
  if (px >= fxmin and px <= fxmax and
      py >= fymin and py <= fymax and
      pz >= fzmin and pz <= fzmax) {

    // Check if the point is on the boundary.
    for (i = 0; i != npts; ++i) {
      j = (i + 1) % npts;
      if ((closestPointOnSegment(p, vertices[i], vertices[j]) - p).magnitude2() < 1.0e-10) return true;
    }

    if (std::abs(normal.x()) >= std::abs(normal.y()) and std::abs(normal.x()) >= std::abs(normal.z())) {

      // x plane -- use (y,z) coordinates.
      for (i = 0, j = npts - 1; i < npts; j = i++) {
        if ( ((vertices[i].z() > pz) != (vertices[j].z() > pz)) &&
             (py < (vertices[j].y() - vertices[i].y()) * (pz - vertices[i].z()) / (vertices[j].z() - vertices[i].z()) + vertices[i].y()) )
          result = not result;
      }

    } else if (std::abs(normal.y()) >= std::abs(normal.x()) and std::abs(normal.y()) >= std::abs(normal.z())) {

      // y plane -- use (z,x) coordinates.
      for (i = 0, j = npts - 1; i < npts; j = i++) {
        if ( ((vertices[i].x() > px) != (vertices[j].x() > px)) &&
             (pz < (vertices[j].z() - vertices[i].z()) * (px - vertices[i].x()) / (vertices[j].x() - vertices[i].x()) + vertices[i].z()) )
          result = not result;
      }

    } else {
      CHECK(std::abs(normal.z()) >= std::abs(normal.x()) and std::abs(normal.z()) >= std::abs(normal.y()));

      // z plane -- use (x,y) coordinate.
      for (i = 0, j = npts - 1; i < npts; j = i++) {
        if ( ((vertices[i].y() > py) != (vertices[j].y() > py)) &&
             (px < (vertices[j].x() - vertices[i].x()) * (py - vertices[i].y()) / (vertices[j].y() - vertices[i].y()) + vertices[i].x()) )
          result = not result;
      }

    }
  }

  // That's it.
  return result;
}

//------------------------------------------------------------------------------
// Test a polygon in 3-D.
// This version allows you to pass in an indirect addressing definition of the
// polygon vertices.
// Note -- does not currently work for polygons with more than one loop!
//------------------------------------------------------------------------------
bool pointInPolygon(const Dim<3>::Vector& p,
                    const vector<Dim<3>::Vector>& vertices,
                    const vector<unsigned>& ipoints,
                    const Dim<3>::Vector& normal,
                    const bool countBoundary,
                    const double tol) {

  // Prerequisites.
  const auto npts = ipoints.size();
  unsigned i, j, ik, jk;
  // BEGIN_CONTRACT_SCOPE
  // {
  //   REQUIRE(ipoints.size() > 2);
  //   const double normmag = normal.magnitude();
  //   Vector normi;
  //   for (i = 1; i != npts; ++i) {
  //     j = (i + 1) % npts;
  //     k = (i + 2) % npts;
  //     normi = (vertices[ipoints[j]] - vertices[ipoints[i]]).cross(vertices[ipoints[k]] - vertices[ipoints[i]]);
  //     REQUIRE2(fuzzyEqual(std::abs(normi.dot(normal)), normmag*normi.magnitude(), 1.0e-5), normi << " " << normal << " " << normi.dot(normal) << " " << normmag*normi.magnitude());
  //   }
  //   // REQUIRE2(fuzzyEqual(pointPlaneDistance(p, vertices[ipoints[0]], normal.unitVector()), 0.0, 1.0e-3), pointPlaneDistance(p, vertices[ipoints[0]], normal.unitVector()));
  // }
  // END_CONTRACT_SCOPE

  bool result = false;

  const double px = p.x();
  const double py = p.y();
  const double pz = p.z();

  // Find the bounding box for the polygon.
  double fxmin = std::numeric_limits<double>::max(),    fymin = std::numeric_limits<double>::max(),   fzmin =  std::numeric_limits<double>::max();
  double fxmax = std::numeric_limits<double>::lowest(), fymax = std::numeric_limits<double>::lowest(), fzmax = std::numeric_limits<double>::lowest();
  for (const auto i: ipoints) {
    fxmin = std::min(fxmin, vertices[i].x());
    fymin = std::min(fymin, vertices[i].y());
    fzmin = std::min(fzmin, vertices[i].z());
    fxmax = std::max(fxmax, vertices[i].x());
    fymax = std::max(fymax, vertices[i].y());
    fzmax = std::max(fzmax, vertices[i].z());
  }
  const auto fuzz = std::max(1.0e-5, tol)*std::max(fxmax - fxmin, std::max(fymax - fymin, fzmax - fzmin));
  fxmin -= fuzz;
  fymin -= fuzz;
  fzmin -= fuzz;
  fxmax += fuzz;
  fymax += fuzz;
  fzmax += fuzz;
  if (px >= fxmin and px <= fxmax and
      py >= fymin and py <= fymax and
      pz >= fzmin and pz <= fzmax) {

    // Check if the point is on the boundary (within tolerance).
    if (pointOnPolygon(p, vertices, ipoints, tol)) return countBoundary;

    // Figure out which plane we're going to project to.
    const double nmax = normal.maxAbsElement();

    // x plane -- use (y,z) coordinates.
    if (std::abs(normal.x()) > 0.9*nmax) {
      for (ik = 0, jk = npts - 1; ik < npts; jk = ik++) {
        i = ipoints[ik];
        j = ipoints[jk];
        if ( ((vertices[i].z() > pz) != (vertices[j].z() > pz)) &&
             (py < (vertices[j].y() - vertices[i].y()) * (pz - vertices[i].z()) / (vertices[j].z() - vertices[i].z()) + vertices[i].y()) )
          result = not result;
      }

    // y plane -- use (z,x) coordinates.
    } else if (std::abs(normal.y()) > 0.9*nmax) {
      for (ik = 0, jk = npts - 1; ik < npts; jk = ik++) {
        i = ipoints[ik];
        j = ipoints[jk];
        if ( ((vertices[i].x() > px) != (vertices[j].x() > px)) &&
             (pz < (vertices[j].z() - vertices[i].z()) * (px - vertices[i].x()) / (vertices[j].x() - vertices[i].x()) + vertices[i].z()) )
          result = not result;
      }

    // z plane -- use (x,y) coordinate.
    } else {
      for (ik = 0, jk = npts - 1; ik < npts; jk = ik++) {
        i = ipoints[ik];
        j = ipoints[jk];
        if ( ((vertices[i].y() > py) != (vertices[j].y() > py)) &&
             (px < (vertices[j].x() - vertices[i].x()) * (py - vertices[i].y()) / (vertices[j].y() - vertices[i].y()) + vertices[i].x()) )
          result = not result;
      }
    }

  }

  // That's it.
  return result;
}

}
