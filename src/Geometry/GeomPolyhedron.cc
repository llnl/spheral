//---------------------------------Spheral++----------------------------------//
// GeomPolyhedron -- Geometric polyhedron class.
//
// A 3-D structure representing a polyhedron as a collection of GeomFacets.
//
// Created by JMO, Fri Jan 29 14:36:01 PST 2010
//----------------------------------------------------------------------------//
#include "GeomPolyhedron.hh"
#include "FacetedVolumeUtilities.hh"
#include "Utilities/removeElements.hh"
#include "Utilities/testBoxIntersection.hh"
#include "Utilities/boundingBox.hh"
#include "Utilities/comparisons.hh"
#include "Utilities/pointInPolygon.hh"
#include "Utilities/PairComparisons.hh"
#include "Utilities/CounterClockwiseComparator.hh"
#include "Utilities/pointInPolyhedron.hh"
#include "Utilities/safeInv.hh"

#include "axom/mint.hpp"

#include <algorithm>
#include <iterator>
#include <limits>
#include <map>
#include <vector>

using std::vector;
using std::map;
using std::pair;
using std::make_pair;
using std::min;
using std::max;
using std::cerr;
using std::endl;

extern "C" {
#include "libqhull/qhull_a.h"
}

#include "Utilities/Timer.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Default constructor.
//------------------------------------------------------------------------------
GeomPolyhedron::
GeomPolyhedron():
  mVertices(),
  mFacets(),
  mVertexUnitNorms(),
  mVertexFacetConnectivity(),
  mFacetFacetConnectivity(),
  mXmin(),
  mXmax(),
  mCentroid(),
  mRinterior2(-1.0),
  mConvex(true),
  mSurfaceMeshQueryPtr(nullptr),
  mSignedDistancePtr(nullptr) {
  if (mDevnull == NULL) mDevnull = fopen("/dev/null", "w");
}

//------------------------------------------------------------------------------
// Construct as a convex hull for the given point set.
//------------------------------------------------------------------------------
GeomPolyhedron::
GeomPolyhedron(const vector<GeomPolyhedron::Vector>& points):
  mVertices(),
  mFacets(),
  mVertexUnitNorms(),
  mVertexFacetConnectivity(),
  mFacetFacetConnectivity(),
  mXmin(),
  mXmax(),
  mCentroid(),
  mRinterior2(-1.0),
  mConvex(true),
  mSurfaceMeshQueryPtr(nullptr),
  mSignedDistancePtr(nullptr) {
  TIME_BEGIN("Polyhedron_construct1");
  if (mDevnull == NULL) mDevnull = fopen("/dev/null", "w");

  if (points.size() > 0) {
    REQUIRE(points.size() > 3);

    // Find the appropriate renormalization so that we can give Qhull points
    // in a unit box.  Qhull just seems to work better this way.
    Vector xmin, xmax;
    boundingBox(points, xmin, xmax);
    const double fscale = (xmax - xmin).maxElement();
    CHECK(fscale > 0.0);

    // Copy the point coordinates to a Qhull point array.
    std::vector<coordT> points_qhull;
    points_qhull.reserve(3 * points.size());
    for (const Vector& vec: points) {
      points_qhull.push_back((vec.x() - xmin.x())/fscale);
      points_qhull.push_back((vec.y() - xmin.y())/fscale);
      points_qhull.push_back((vec.z() - xmin.z())/fscale);
    }
    CHECK(points_qhull.size() == 3*points.size());

    // Call Qhull to generate the hull (C interface).
    boolT ismalloc = False;   /* True if qhull should free points in qh_freeqhull() or reallocation */
    char flags[250];          /* option flags for qhull, see qh_opt.htm */
//     FILE *outfile= NULL;      /* output from qh_produce_output()
//                                  use NULL to skip qh_produce_output() */
//     FILE *errfile= stderr;    /* error messages from qhull code */
//     FILE *errfile= fopen("/dev/null", "w");    /* error messages from qhull code */
    //int exitcode;             /* 0 if no error from qhull */
    facetT *facet;            /* set by FORALLfacets */
    int curlong, totlong;     /* memory remaining after qh_memfreeshort */
    sprintf(flags, "qhull s"); // Tcv");
    const int exitcode_qhull = qh_new_qhull(3, points.size(), &points_qhull.front(), ismalloc, flags, mDevnull, mDevnull);

//     if (exitcode_qhull != 0) {
//       // Something didn't work, so spit out the points for diagnostics.
//       vector<pair<double, Vector> > sorted_points;
//       for (vector<Vector>::const_iterator itr = points.begin();
//            itr != points.end();
//            ++itr) sorted_points.push_back(make_pair(itr->magnitude(), *itr));
//       sort(sorted_points.begin(), sorted_points.end(), ComparePairByFirstElement<pair<double, Vector> >());
//       for (size_t k = 0; k != sorted_points.size(); ++k) cerr << "    -----> " << sorted_points[k].first << " " << sorted_points[k].second << endl;

//       // Emit the error message by calling qhull again.
//       FILE *errfile= stderr;    /* error messages from qhull code */
//       const int exitcode_qhull = qh_new_qhull(3, points.size(), &points_qhull.front(), ismalloc, flags, errfile, errfile);
//     }

    // VERIFY2(exitcode_qhull == 0,
    //         "Qhull emitted an error code while generating GeomPolyhedron");

    // If Qhull failed, we fall back on constructing the bounding box.
    if (exitcode_qhull != 0) {
      mVertices.push_back(Vector(xmax.x(), xmin.y(), xmin.z()));
      mVertices.push_back(Vector(xmax.x(), xmax.y(), xmin.z()));
      mVertices.push_back(Vector(xmin.x(), xmax.y(), xmin.z()));
      mVertices.push_back(Vector(xmin.x(), xmin.y(), xmin.z()));
      mVertices.push_back(Vector(xmax.x(), xmin.y(), xmax.z()));
      mVertices.push_back(Vector(xmax.x(), xmax.y(), xmax.z()));
      mVertices.push_back(Vector(xmin.x(), xmax.y(), xmax.z()));
      mVertices.push_back(Vector(xmin.x(), xmin.y(), xmax.z()));
      vector<unsigned> ipoints0, ipoints1, ipoints2, ipoints3, ipoints4, ipoints5;
      ipoints0.push_back(0); ipoints0.push_back(3); ipoints0.push_back(2); ipoints0.push_back(1);
      ipoints1.push_back(0); ipoints1.push_back(1); ipoints1.push_back(5); ipoints1.push_back(4);
      ipoints2.push_back(1); ipoints2.push_back(2); ipoints2.push_back(6); ipoints2.push_back(5);
      ipoints3.push_back(2); ipoints3.push_back(3); ipoints3.push_back(7); ipoints3.push_back(6);
      ipoints4.push_back(0); ipoints4.push_back(4); ipoints4.push_back(7); ipoints4.push_back(3);
      ipoints5.push_back(4); ipoints5.push_back(5); ipoints5.push_back(6); ipoints5.push_back(7);
      mFacets.push_back(Facet(mVertices, ipoints0));
      mFacets.push_back(Facet(mVertices, ipoints1));
      mFacets.push_back(Facet(mVertices, ipoints2));
      mFacets.push_back(Facet(mVertices, ipoints3));
      mFacets.push_back(Facet(mVertices, ipoints4));
      mFacets.push_back(Facet(mVertices, ipoints5));

    } else {

      // Copy Qhull's vertex information.
      vertexT *vertex, **vertexp;
      map<unsigned, unsigned> vertexIDmap;
      {
        unsigned i = 0;
        FORALLvertices {
          mVertices.push_back(fscale*Vector(vertex->point[0], vertex->point[1], vertex->point[2]) + xmin);
          vertexIDmap[vertex->id] = i;
          ++i;
        }
      }
      CHECK(vertexIDmap.size() == mVertices.size());

      // Copy Qhull's facet information to our own format.
      vector<unsigned> vertex_indices;
      FORALLfacets {
        vertex_indices = vector<unsigned>();
        Vector fc;
        FOREACHvertex_(facet->vertices) {
          CHECK(vertexIDmap.find(vertex->id) != vertexIDmap.end());
          const unsigned vertex_i = vertexIDmap[vertex->id];
          CHECK(vertex_i < mVertices.size());
          vertex_indices.push_back(vertex_i);
          fc += mVertices[vertex_indices.back()];
        }
        CHECK(vertex_indices.size() >= 3);
        fc /= vertex_indices.size();
        const Vector normal(facet->normal[0], facet->normal[1], facet->normal[2]);

        // Enforce counter-clockwise ordering (viewed from outside) for the facet vertices.
        CHECK(vertex_indices.size() > 2);
        sort(vertex_indices.begin() + 1, vertex_indices.end(), CounterClockwiseComparator<Vector, vector<Vector> >(mVertices, mVertices[vertex_indices[0]], normal));

//       // If this is a triangular facet, we enforce counter-clockwise ordering (viewed from the outside) for the facet vertices.
//       if (vertex_indices.size() == 3) {
//         if ((mVertices[vertex_indices[1]] - fc).cross(mVertices[vertex_indices[0]] - fc).dot(normal) < 0.0) swap(vertex_indices[0], vertex_indices[1]);
//         if ((mVertices[vertex_indices[2]] - fc).cross(mVertices[vertex_indices[0]] - fc).dot(normal) < 0.0) swap(vertex_indices[0], vertex_indices[2]);
//         if ((mVertices[vertex_indices[2]] - fc).cross(mVertices[vertex_indices[1]] - fc).dot(normal) < 0.0) swap(vertex_indices[1], vertex_indices[2]);
//         CHECK((mVertices[vertex_indices[2]] - mVertices[vertex_indices[0]]).cross(mVertices[vertex_indices[1]] - mVertices[vertex_indices[0]]).dot(normal) > 0.0);
//       }

        mFacets.push_back(Facet(mVertices,
                                vertex_indices));
      }
    }

    // Free Qhull's resources.
    qh_freeqhull(!qh_ALL);
    qh_memfreeshort(&curlong, &totlong);

    // Fill in our bounding box.
    setBoundingBox();

    // Post-conditions.
    BEGIN_CONTRACT_SCOPE
    {
      // All normals should be outward facing.
      for (const auto& facet: mFacets) { 
        CONTRACT_VAR(facet);
        ENSURE2(((facet.position() - mCentroid).dot(facet.normal()) >= 0.0),
               "Inward normal? " << (facet.position() - mCentroid).dot(facet.normal()) << " " << facet.position() << " " << mCentroid << " " << facet.normal());
      }

      // We had better be convex if built from a convex hull.
      ENSURE(convex());

//       // Ensure the seed points are contained.
//       for (vector<Vector>::const_iterator itr = points.begin();
//            itr != points.end();
//            ++itr) ENSURE(convexContains(*itr));
    }
    END_CONTRACT_SCOPE

  }
  TIME_END("Polyhedron_construct1");
}

//------------------------------------------------------------------------------
// Construct given the positions and facet indices.
// We assume here that the nodes for each facet are arranged correctly to
// create outward pointing normals.
//------------------------------------------------------------------------------
GeomPolyhedron::
GeomPolyhedron(const vector<GeomPolyhedron::Vector>& points,
               const vector<vector<unsigned> >& facetIndices):
  mVertices(points),
  mFacets(),
  mVertexUnitNorms(),
  mVertexFacetConnectivity(),
  mFacetFacetConnectivity(),
  mXmin(),
  mXmax(),
  mCentroid(),
  mRinterior2(-1.0),
  mConvex(false),
  mSurfaceMeshQueryPtr(nullptr),
  mSignedDistancePtr(nullptr) {
  TIME_BEGIN("Polyhedron_construct2");

  // Construct the facets.
  mFacets.reserve(facetIndices.size());
  for (const auto& indices: facetIndices) {
    CHECK2(indices.size() > 2, "Need at least two points per facet.");
    CHECK2(*max_element(indices.begin(), indices.end()) < points.size(),
           "Bad vertex index for facet.");
    mFacets.push_back(Facet(mVertices, indices));
  }
  CHECK(mFacets.size() == facetIndices.size());

  // Fill in our bounding box.
  setBoundingBox();
  TIME_END("Polyhedron_construct2");
}

//------------------------------------------------------------------------------
// Copy constructor.
//------------------------------------------------------------------------------
GeomPolyhedron::
GeomPolyhedron(const GeomPolyhedron& rhs):
  mVertices(rhs.mVertices),
  mFacets(rhs.mFacets),
  mVertexUnitNorms(rhs.mVertexUnitNorms),
  mVertexFacetConnectivity(rhs.mVertexFacetConnectivity),
  mFacetFacetConnectivity(rhs.mFacetFacetConnectivity),
  mXmin(rhs.mXmin),
  mXmax(rhs.mXmax),
  mCentroid(rhs.mCentroid),
  mRinterior2(rhs.mRinterior2),
  mConvex(rhs.mConvex),
  mSurfaceMeshQueryPtr(nullptr),
  mSignedDistancePtr(nullptr) {
  for (Facet& facet: mFacets) facet.mVerticesPtr = &mVertices;
}

//------------------------------------------------------------------------------
// Assignment operator.
//------------------------------------------------------------------------------
GeomPolyhedron&
GeomPolyhedron::
operator=(const GeomPolyhedron& rhs) {
  if (this != &rhs) {
    mVertices = rhs.mVertices;
    mFacets = vector<Facet>();
    mFacets.reserve(rhs.mFacets.size());
    for (const Facet& facet: rhs.mFacets) mFacets.push_back(Facet(mVertices,
                                                                  facet.ipoints()));
    mVertexFacetConnectivity = rhs.mVertexFacetConnectivity;
    mFacetFacetConnectivity = rhs.mFacetFacetConnectivity;
    mVertexUnitNorms = rhs.mVertexUnitNorms;
    mXmin = rhs.mXmin;
    mXmax = rhs.mXmax;
    mCentroid = rhs.mCentroid;
    mRinterior2 = rhs.mRinterior2;
    mConvex = rhs.mConvex;
    mSurfaceMeshPtr = nullptr;
    mSurfaceMeshQueryPtr = nullptr;
    mSignedDistancePtr = nullptr;
  }
  ENSURE(mFacets.size() == rhs.mFacets.size());
  return *this;
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
GeomPolyhedron::
~GeomPolyhedron() {
  if (mSurfaceMeshQueryPtr != nullptr) delete mSurfaceMeshQueryPtr;
  if (mSignedDistancePtr != nullptr) delete mSignedDistancePtr;
}

//------------------------------------------------------------------------------
// Test if the given point is on the interior of the polyhedron.
//------------------------------------------------------------------------------
bool
GeomPolyhedron::
contains(const GeomPolyhedron::Vector& point,
         const bool countBoundary,
         const double tol,
         const bool useAxom) const {
  if (useAxom) {

    // Experimental version using Axom
    using AxPoint = axom::quest::InOutOctree<3>::SpacePt;
    if (!mSurfaceMeshPtr) this->buildAxomData();
    const auto inside = mSurfaceMeshQueryPtr->within(AxPoint(&const_cast<Vector&>(point)[0]));
    if (not inside and countBoundary) {
      return this->distance(point) < tol;
    } else {
      return inside;
    }

  } else {

    // Our old internal methods
    if (not testPointInBox(point, mXmin, mXmax, tol)) return false;
    if ((point - mCentroid).magnitude2() < mRinterior2 - tol) return true;
    if (mConvex) {
      return this->convexContains(point, countBoundary, tol);
    } else {
      return pointInPolyhedron(point, *this, countBoundary, tol);
    }

  }
}

//------------------------------------------------------------------------------
// Test if the given point is on the interior of the polyhedron.
// This method only works for convex polyhedra!
//------------------------------------------------------------------------------
bool
GeomPolyhedron::
convexContains(const GeomPolyhedron::Vector& point,
               const bool countBoundary,
               const double tol) const {
  if ((point - mCentroid).magnitude2() < mRinterior2 - tol) return true;
  vector<Facet>::const_iterator facetItr = mFacets.begin();
  bool result = true;
  if (countBoundary) {
    while (facetItr != mFacets.end() and result) {
      result = (facetItr->compare(point, tol) <= 0);
      ++facetItr;
    }
  } else {
    while (facetItr != mFacets.end() and result) {
      result = (facetItr->compare(point, tol) < 0);
      ++facetItr;
    }
  }
  return result;
}

//------------------------------------------------------------------------------
// Test if we intersect the given polyhedron.
//------------------------------------------------------------------------------
bool
GeomPolyhedron::
intersect(const GeomPolyhedron& rhs) const {
  if (not testBoxIntersection(mXmin, mXmax, rhs.mXmin, rhs.mXmax)) return false;
  for (const auto& vec: mVertices) {
    if (rhs.contains(vec)) return true;
  }
  for (const auto& vec: rhs.mVertices) {
    if (this->contains(vec)) return true;
  }
  return false;
}

//------------------------------------------------------------------------------
// Test if we intersect the given polyhedron assuming both are convex.
//------------------------------------------------------------------------------
bool
GeomPolyhedron::
convexIntersect(const GeomPolyhedron& rhs) const {
  REQUIRE(this->convex());
  if (not testBoxIntersection(mXmin, mXmax, rhs.mXmin, rhs.mXmax)) return false;
  
  // Check if we can exclude rhs from us.
  bool outside = true;
  {
    std::vector<Facet>::const_iterator facetItr = mFacets.begin();
    while (outside and facetItr != mFacets.end()) {
      outside = (facetItr->compare(rhs.mVertices) == 1);
      ++facetItr;
    }
    if (outside) return false;
  }

  // Check if we can exclude us from rhs.
  outside = true;
  {
    std::vector<Facet>::const_iterator facetItr = rhs.mFacets.begin();
    while (outside and facetItr != rhs.mFacets.end()) {
      outside = (facetItr->compare(mVertices) == 1);
      ++facetItr;
    }
    if (outside) return false;
  }

  // We can't exclude anybody, so must intersect!
  return true;
}

//------------------------------------------------------------------------------
// Test if we intersect the given box.
//------------------------------------------------------------------------------
bool
GeomPolyhedron::
intersect(const std::pair<Vector, Vector>& rhs) const {
  if (not testBoxIntersection(mXmin, mXmax, rhs.first, rhs.second)) return false;
  
  // Build a GeompPolygon representation of the box and use our generic intersection
  // method.
  // Create the piecewise linear complex representing the box. Note that 
  // the box consists of facets that are defined by their connections to 
  // generating points.
  // Should look like the following:
  //
  //        6--------7            y
  //       /        /|            |
  //      /        / |            |
  //     2--------3  |             ------x
  //     |  .     |  |           /
  //     |  4.....|..5          z
  //     | .      | / 
  //     |.       |/
  //     0--------1             
  //

  // Check if any of our vertices fall in the box.
  const double x1 = rhs.first.x(), y1 = rhs.first.y(), z1 = rhs.first.z(),
               x2 = rhs.second.x(), y2 = rhs.second.y(), z2 = rhs.second.z();
  for (const auto& v: mVertices) {
    if (v.x() >= x1 and v.x() <= x2 and
        v.y() >= y1 and v.y() <= y2 and
        v.z() >= z1 and v.z() <= z2) return true;
  }

  // Check if any of the box vertices are inside this polyhedron
  if (this->contains(Vector(x1, y1, z2))) return true;
  if (this->contains(Vector(x2, y1, z2))) return true;
  if (this->contains(Vector(x1, y2, z2))) return true;
  if (this->contains(Vector(x2, y2, z2))) return true;
  if (this->contains(Vector(x1, y1, z1))) return true;
  if (this->contains(Vector(x2, y1, z1))) return true;
  if (this->contains(Vector(x1, y2, z1))) return true;
  if (this->contains(Vector(x2, y2, z1))) return true;

  // Otherwise nope!
  return false;
}

//------------------------------------------------------------------------------
// Test if we intersect a line segment (interior counts as intersection).
//------------------------------------------------------------------------------
bool
GeomPolyhedron::
intersect(const Vector& s0, const Vector& s1) const {
  if (this->contains(s0) or this->contains(s1)) return true;

  const auto shat = (s1 - s0).unitVector();
  const auto q01 = (s1 - s0).magnitude();

  // Check each facet of the polyhedron.
  Vector inter;
  const auto nf = mFacets.size();
  for (auto k = 0u; k < nf; ++k) {
    const auto& facet = mFacets[k];
    const auto& ipoints = facet.ipoints();
    CHECK(ipoints.size() >= 3);
    const auto& fhat = facet.normal();

    // First does the line segment intersect the facet plane?
    const auto shat_dot_fhat = shat.dot(fhat);
    if (fabs(shat_dot_fhat) > 1.0e-10) {           // Check for segment parallel to plane
      const auto q = (mVertices[ipoints[0]] - s0).dot(fhat)/shat_dot_fhat;
      if (q >= 0.0 and q <= q01) {                // Does the line segment intersect the plane
        return true;
      }
    }
  }

  return false;
}

//------------------------------------------------------------------------------
// Find the Facets and points intersecting a line segment (s0, s1).
//------------------------------------------------------------------------------
void
GeomPolyhedron::
intersections(const Vector& s0, const Vector& s1,
              std::vector<unsigned>& facetIDs,
              std::vector<Vector>& intersections) const {
  facetIDs.clear();
  intersections.clear();

  const auto shat = (s1 - s0).unitVector();
  const auto q01 = (s1 - s0).magnitude();

  // Preserve the distance along the line segment for sorting.
  vector<double> qvals;

  // Check each facet of the polyhedron.
  Vector inter;
  const auto nf = mFacets.size();
  for (auto k = 0u; k < nf; ++k) {
    const auto& facet = mFacets[k];
    const auto& ipoints = facet.ipoints();
    CHECK(ipoints.size() >= 3);
    const auto& fhat = facet.normal();

    // First does the line segment intersect the facet plane?
    const auto shat_dot_fhat = shat.dot(fhat);
    if (fabs(shat_dot_fhat) > 1.0e-10) {           // Check for segment parallel to plane
      const auto q = (mVertices[ipoints[0]] - s0).dot(fhat)/shat_dot_fhat;
      if (q >= 0.0 and q <= q01) {                // Does the line segment intersect the plane
        const auto p = s0 + q*shat;               // point along line in plane
        if (pointInPolygon(p, mVertices, ipoints, fhat, true, 1.0e-10)) {
          auto qitr = qvals.insert(std::upper_bound(qvals.begin(), qvals.end(), q), q);
          auto qpos = std::distance(qvals.begin(), qitr);
          facetIDs.insert(facetIDs.begin() + qpos, k);
          intersections.insert(intersections.begin() + qpos, p);
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
// Compute the centroid.
//------------------------------------------------------------------------------
GeomPolyhedron::Vector
GeomPolyhedron::
centroid() const {
  Vector result;
  double dV, V = 0.0;
  for (const auto& facet: mFacets) {
    const auto& ipoints = facet.ipoints();
    const auto& v0 = mVertices[ipoints[0]];
    const auto  nverts = ipoints.size();
    for (auto i = 1u; i < nverts - 1; ++i) {
      const auto& v1 = mVertices[ipoints[i]];
      const auto& v2 = mVertices[ipoints[i+1]];
      dV = v0.dot(v1.cross(v2));
      V += dV;
      result += dV*(v0 + v1 + v2);
    }
  }
  result *= safeInvVar(4.0*V);
  return result;
}

//------------------------------------------------------------------------------
// Return the edges of the polyhedron as a set of integer pairs for the vertices.
//------------------------------------------------------------------------------
vector<pair<unsigned, unsigned> >
GeomPolyhedron::
edges() const {
  vector<pair<unsigned, unsigned> > result;
  unsigned i, j, k, npoints;
  for (const Facet& facet: mFacets) {
    const vector<unsigned>& ipoints = facet.ipoints();
    npoints = ipoints.size();
    for (k = 0; k != npoints; ++k) {
      i = ipoints[k];
      j = ipoints[(k + 1) % npoints];
      result.push_back(i < j ? make_pair(i, j) : make_pair(j, i));
    }
  }
  sort(result.begin(), result.end(), ComparePairs<pair<unsigned, unsigned> >());
  result.erase(unique(result.begin(), result.end()), result.end());
  return result;
}

//------------------------------------------------------------------------------
// Spit out an encoding of the facets as ordered vertex indices.
//------------------------------------------------------------------------------
vector<vector<unsigned> >
GeomPolyhedron::
facetVertices() const {
  vector<vector<unsigned> > result;
  if (mVertices.size() > 0) {
    for (const Facet& facet: mFacets) {
      vector<unsigned> pts;
      copy(facet.ipoints().begin(), facet.ipoints().end(), back_inserter(pts));
      CHECK(pts.size() == facet.ipoints().size());
      result.push_back(pts);
    }
  }
  return result;
}

//------------------------------------------------------------------------------
// Spit out the facet centroids.
//------------------------------------------------------------------------------
vector<GeomPolyhedron::Vector>
GeomPolyhedron::
facetCentroids() const {
  vector<Vector> result;
  result.reserve(mFacets.size());
  for (const Facet& facet: mFacets) result.push_back(facet.position());
  ENSURE(result.size() == mFacets.size());
  return result;
}

//------------------------------------------------------------------------------
// Spit out the facet normals.
//------------------------------------------------------------------------------
vector<GeomPolyhedron::Vector>
GeomPolyhedron::
facetNormals() const {
  vector<Vector> result;
  result.reserve(mFacets.size());
  for (const Facet& facet: mFacets) result.push_back(facet.normal());
  ENSURE(result.size() == mFacets.size());
  return result;
}

//------------------------------------------------------------------------------
// Spit out the facet normals.
//------------------------------------------------------------------------------
vector<GeomPolyhedron::Vector>
GeomPolyhedron::
facetAreaVectors() const {
  vector<Vector> result;
  result.reserve(mFacets.size());
  for (const Facet& facet: mFacets) result.push_back(facet.normal()*facet.area());
  ENSURE(result.size() == mFacets.size());
  return result;
}

//------------------------------------------------------------------------------
// Reconstruct the internal state given the set of vertices and the enocded 
// facets.
//------------------------------------------------------------------------------
void
GeomPolyhedron::
reconstruct(const vector<GeomPolyhedron::Vector>& vertices,
            const vector<vector<unsigned> >& facetVertices) {
  const unsigned numFacets = facetVertices.size();
  mVertices = vertices;
  mFacets = vector<Facet>();
  mFacets.reserve(numFacets);
  for (unsigned ifacet = 0; ifacet != numFacets; ++ifacet) {
    mFacets.push_back(Facet(mVertices,
                            facetVertices[ifacet]));
  }
  setBoundingBox();
  mConvex = this->convex();
  mVertexFacetConnectivity.clear();
  mFacetFacetConnectivity.clear();
  mVertexUnitNorms.clear();
  // GeometryUtilities::computeAncillaryGeometry(*this, mVertexFacetConnectivity, mFacetFacetConnectivity, mVertexUnitNorms, false);
  ENSURE(mFacets.size() == numFacets);
  ENSURE(mFacetFacetConnectivity.size() == 0); // numFacets);
}

//------------------------------------------------------------------------------
// Compute the volume.
//------------------------------------------------------------------------------
double
GeomPolyhedron::
volume() const {
  double result = 0.0;
  for (const Facet& facet: mFacets) {
    result += facet.area() * facet.normal().dot(facet.point(0));
  }
  ENSURE(result >= 0.0);
  return result/3.0;
}

//------------------------------------------------------------------------------
// Find the facet closest to the given point.
//------------------------------------------------------------------------------
unsigned
GeomPolyhedron::
closestFacet(const GeomPolyhedron::Vector& p) const {
  unsigned result = 0;
  double r2, minr2 = std::numeric_limits<double>::max();
  Vector thpt;
  for (unsigned i = 0; i != mFacets.size(); ++i) {
    thpt = mFacets[i].closestPoint(p);
    r2 = (thpt - p).magnitude2();
    if (r2 < minr2) {
      result = i;
      minr2 = r2;
    }
  }
  ENSURE(result < mFacets.size());
  return result;
}

//------------------------------------------------------------------------------
// Find the minimum distance to a point.
//------------------------------------------------------------------------------
double
GeomPolyhedron::
distance(const GeomPolyhedron::Vector& p,
         const bool useAxom) const {

  if (useAxom) {

    // Experimental version using Axom
    using AxPoint = axom::quest::InOutOctree<3>::SpacePt;
    if (!mSurfaceMeshPtr) this->buildAxomData();
    return std::abs(mSignedDistancePtr->computeDistance(AxPoint(&const_cast<Vector&>(p)[0])));

  } else {
    
    // Not signed currently
    const auto d = p - this->closestPoint(p);
    return d.magnitude();

  }
}

//------------------------------------------------------------------------------
// Find the point in the polyhedron closest to the given point.
//------------------------------------------------------------------------------
GeomPolyhedron::Vector
GeomPolyhedron::
closestPoint(const GeomPolyhedron::Vector& p) const {
  const Facet& f = mFacets[this->closestFacet(p)];
  return f.closestPoint(p);
}

//------------------------------------------------------------------------------
// += Vector, shift polyhedron in space
//------------------------------------------------------------------------------
GeomPolyhedron&
GeomPolyhedron::
operator+=(const GeomPolyhedron::Vector& rhs) {
  for (auto& v: mVertices) v += rhs;
  this->setBoundingBox();
  return *this;
}

//------------------------------------------------------------------------------
// -= Vector, shift polyhedron in space
//------------------------------------------------------------------------------
GeomPolyhedron&
GeomPolyhedron::
operator-=(const GeomPolyhedron::Vector& rhs) {
  (*this) += -rhs;
  return *this;
}

//------------------------------------------------------------------------------
// + Vector, return shifted polyhedron in space
//------------------------------------------------------------------------------
GeomPolyhedron
GeomPolyhedron::
operator+(const GeomPolyhedron::Vector& rhs) const {
  GeomPolyhedron result(*this);
  result += rhs;
  return result;
}

//------------------------------------------------------------------------------
// - Vector, return shifted polyhedron in space
//------------------------------------------------------------------------------
GeomPolyhedron
GeomPolyhedron::
operator-(const GeomPolyhedron::Vector& rhs) const {
  return (*this) + (-rhs);
}

//------------------------------------------------------------------------------
// *= Scalar, scale polyhedron
//------------------------------------------------------------------------------
GeomPolyhedron&
GeomPolyhedron::
operator*=(const double rhs) {
  for (auto& v: mVertices) v *= rhs;
  this->setBoundingBox();
  return *this;
}

//------------------------------------------------------------------------------
// /= Scalar, scale polyhedron
//------------------------------------------------------------------------------
GeomPolyhedron&
GeomPolyhedron::
operator/=(const double rhs) {
  CHECK(rhs != 0.0);
  (*this) *= 1.0/rhs;
  return *this;
}

//------------------------------------------------------------------------------
// * Scalar, scale polyhedron
//------------------------------------------------------------------------------
GeomPolyhedron
GeomPolyhedron::
operator*(const double rhs) const {
  GeomPolyhedron result(*this);
  result *= rhs;
  return result;
}

//------------------------------------------------------------------------------
// / Scalar, scale polyhedron
//------------------------------------------------------------------------------
GeomPolyhedron
GeomPolyhedron::
operator/(const double rhs) const {
  GeomPolyhedron result(*this);
  result /= rhs;
  return result;
}

//------------------------------------------------------------------------------
// transform
//------------------------------------------------------------------------------
GeomPolyhedron&
GeomPolyhedron::
transform(const typename GeomPolyhedron::Tensor& t) {
  for (auto& v: mVertices) v = t*v;
  this->setBoundingBox();
  return *this;
}

//------------------------------------------------------------------------------
// ==
//------------------------------------------------------------------------------
bool
GeomPolyhedron::
operator==(const GeomPolyhedron& rhs) const {
  bool result = (mVertices == rhs.mVertices and
                 mFacets.size() == rhs.mFacets.size());
  size_t i = 0;
  while (result and i != mFacets.size()) {
    result = (mFacets[i] == rhs.mFacets[i]);
    ++i;
  }
  return result;
}

//------------------------------------------------------------------------------
// !=
//------------------------------------------------------------------------------
bool
GeomPolyhedron::
operator!=(const GeomPolyhedron& rhs) const {
  return not (*this == rhs);
}

//------------------------------------------------------------------------------
// Set the minimum and maximum extents of the polyhedron.
// We are also piggybacking updating the Axom data structures here as well.
//------------------------------------------------------------------------------
void
GeomPolyhedron::
setBoundingBox() {
  TIME_BEGIN("Polyhedron_BB");
  boundingBox(mVertices, mXmin, mXmax);

  // Check if we're convex.
  mConvex = this->convex();

  // Compute the ancillary geometry.
  TIME_BEGIN("Polyhedron_BB_ancillary");
  mVertexFacetConnectivity.clear();
  mFacetFacetConnectivity.clear();
  mVertexUnitNorms.clear();
  // GeometryUtilities::computeAncillaryGeometry(*this, mVertexFacetConnectivity, mFacetFacetConnectivity, mVertexUnitNorms, false);
  TIME_END("Polyhedron_BB_ancillary");

  // Have the facets recompute their normals
  for (auto& f: mFacets) f.computeNormal();

  // Stash the centroid and inscribed radius for use in containment.  If the centroid is not contained however,
  // we set this internal radius to zero to disable this accelerated containment checking.
  TIME_BEGIN("Polyhedron_BB_centroid");
  mCentroid = this->centroid();
  TIME_END("Polyhedron_BB_centroid");

  TIME_BEGIN("Polyhedron_BB_R2");
  if (pointInPolyhedron(mCentroid, *this, false, 1.0e-10)) {
    mRinterior2 = std::numeric_limits<double>::max();
    for (const auto& facet: mFacets) mRinterior2 = min(mRinterior2, facet.distance(mCentroid));
    mRinterior2 = FastMath::square(mRinterior2);
  } else {
    mRinterior2 = -1.0;
  }
  mRinterior2 = -1.0; // BLAGO!
  TIME_END("Polyhedron_BB_R2");

  // Clear any existing Axom information, so it's reconstructed if needed
  if (mSurfaceMeshQueryPtr != nullptr) delete mSurfaceMeshQueryPtr;
  if (mSignedDistancePtr != nullptr) delete mSignedDistancePtr;
  mSurfaceMeshPtr.reset();
  mSurfaceMeshQueryPtr = nullptr;
  mSignedDistancePtr = nullptr;
  TIME_END("Polyhedron_BB");
}

//------------------------------------------------------------------------------
// Build the Axom mesh representation
//------------------------------------------------------------------------------
void
GeomPolyhedron::
buildAxomData() const {
  using AxMesh = axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>;
  using AxOctree = axom::quest::InOutOctree<3>;
  using AxBoundingBox = AxOctree::GeometricBoundingBox;
  using AxPoint = AxOctree::SpacePt;
  using AxDistance = axom::quest::SignedDistance<3>;

  // Set the vertex positions
  auto* meshPtr = new AxMesh(3, axom::mint::TRIANGLE);
  for (const auto& v: mVertices) {
    meshPtr->appendNode(v.x(), v.y(), v.z());
  }

  // Set the facets -- note we have to ensure all Axom facets are triangles
  // size_t maxFacet = 0u;
  for (const auto& f: mFacets) {
    const auto& ipts = f.ipoints();
    const auto  n = ipts.size();
    for (auto k = 1u; k < n - 1u; ++k) {  // This only works for convex facets!
      const axom::IndexType tri[] = {axom::IndexType(ipts[0]),
                                     axom::IndexType(ipts[k]),
                                     axom::IndexType(ipts[k + 1])};
      meshPtr->appendCell(tri);
    }
  }

  // Build the query objects
  AxBoundingBox bb;
  Vector xmin = mXmin - 0.01*(mXmax - mXmin);
  Vector xmax = mXmax + 0.01*(mXmax - mXmin);
  bb.addPoint(AxPoint(&xmin[0]));
  bb.addPoint(AxPoint(&xmax[0]));
  axom::mint::write_vtk(meshPtr, "blago.vtk");
  mSurfaceMeshPtr = std::shared_ptr<axom::quest::InOutOctree<3>::SurfaceMesh>(meshPtr);
  mSurfaceMeshQueryPtr = new AxOctree(bb, mSurfaceMeshPtr);
  mSurfaceMeshQueryPtr->generateIndex();
  mSignedDistancePtr = new AxDistance(mSurfaceMeshPtr.get(),
                                      true);               // is_watertight
}

//------------------------------------------------------------------------------
// Test if the polyhedron is convex.
//------------------------------------------------------------------------------
bool
GeomPolyhedron::
convex(const double tol) const {
  TIME_BEGIN("Polyhedron_convex");
  // Do the convex comparison for each vertex.
  bool result = true;
  const double reltol = tol*max(1.0, (mXmax - mXmin).maxAbsElement());
  vector<Vector>::const_iterator vertexItr = mVertices.begin();
  while (vertexItr != mVertices.end() and result) {
    vector<Facet>::const_iterator facetItr = mFacets.begin();
    while (facetItr != mFacets.end() and result) {
      result = (facetItr->compare(*vertexItr, reltol) <= 0);
      ++facetItr;
    }
    ++vertexItr;
  }
  TIME_END("Polyhedron_convex");
  return result;
}

//------------------------------------------------------------------------------
// Decompose the polyhedron into unique pyramids for each facet.
//------------------------------------------------------------------------------
GeomPolyhedron
GeomPolyhedron::
facetSubVolume(const unsigned facetID) const {
  REQUIRE(facetID < mFacets.size());
  const auto& facet = mFacets[facetID];
  std::vector<Vector> points = {this->centroid()};
  const auto& ipoints = facet.ipoints();
  const auto  n = ipoints.size();
  for (auto i = 0u; i < n; ++i) points.push_back(facet.point(i));
  return GeomPolyhedron(points);
}

//------------------------------------------------------------------------------
// Decompose the polyhedron into tetrahedra.
//------------------------------------------------------------------------------
void
GeomPolyhedron::
decompose(std::vector<GeomPolyhedron>& subcells) const {
  const auto originalCentroid = this->centroid();
  const auto numFacets = mFacets.size();
  subcells.clear();
  subcells.reserve(numFacets);
  std::vector<std::array<Vector, 3>> subfacets;
  const std::vector<std::vector<unsigned>> indices = {{0, 1, 2}, {0, 3, 1},
                                                      {1, 3, 2}, {0, 2, 3}};
  for (auto f = 0u; f < numFacets; ++f) {
    const auto& facet = mFacets[f];
    // We don't want to split the facet if we can help it
    if (facet.ipoints().size() == 3) {
      subfacets = {{facet.point(0), facet.point(1), facet.point(2)}};
    }
    else {
      facet.decompose(subfacets);
    }
    
    for (auto& subfacet : subfacets) {
      CHECK(subfacet.size() == 3);
      std::vector<Vector> points = {subfacet[0], subfacet[1],
                                    subfacet[2], originalCentroid};
      subcells.emplace_back(points, indices);
    }
  }

  BEGIN_CONTRACT_SCOPE
  {
    const auto originalVolume = this->volume();
    auto volumesum = 0.;
    CONTRACT_VAR(volumesum);
    for (auto& subcell : subcells) {
      const auto subvolume = subcell.volume();
      CONTRACT_VAR(subvolume);
      CONTRACT_VAR(originalVolume);
      CHECK(0 < subvolume and subvolume < originalVolume);
      volumesum += subcell.volume();
    }
    CHECK(fuzzyEqual(volumesum, originalVolume));
  }
  END_CONTRACT_SCOPE
}

//------------------------------------------------------------------------------
// ostream operator.
//------------------------------------------------------------------------------
std::ostream& operator<<(std::ostream& os, const GeomPolyhedron& polyhedron) {
  typedef GeomPolyhedron::Vector Vector;
  const vector<Vector>& vertices = polyhedron.vertices();
  const vector<vector<unsigned> >& facetIndices = polyhedron.facetVertices();
  os << "Polyhedron( vertices[\n";
  for (size_t i = 0; i != vertices.size(); ++i) os << "                    " << i << " " << vertices[i] << "\n";
  os << "            ]\n           facets[\n";
  for (size_t i = 0; i != facetIndices.size(); ++i) {
    os << "                    " << i << " [";
    std::copy(facetIndices[i].begin(), facetIndices[i].end(), std::ostream_iterator<unsigned>(os, " "));
    os << "]\n";
  }
  os << "])\n";
  return os;
}

//------------------------------------------------------------------------------
// Initialization.
//------------------------------------------------------------------------------
FILE* GeomPolyhedron::mDevnull = NULL;

}


