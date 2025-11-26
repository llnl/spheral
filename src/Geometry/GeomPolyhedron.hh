//---------------------------------Spheral++----------------------------------//
// GeomPolyhedron -- Geometric polyhedron class.
//
// A 3-D structure representing a polyhedron as a collection of GeomFacets.
//
// Created by JMO, Fri Jan 29 14:36:01 PST 2010
//----------------------------------------------------------------------------//
#ifndef __Spheral_GeomPolyhedron__
#define __Spheral_GeomPolyhedron__

#include "GeomVector.hh"
#include "GeomTensor.hh"
#include "GeomFacet3d.hh"

#include "axom/quest.hpp"

#include <vector>
#include <utility>

namespace Spheral {

class GeomPolyhedron {
public:
  //--------------------------- Public Interface ---------------------------//
  using Vector = GeomVector<3>;
  using Tensor = GeomTensor<3>;
  using SymTensor = GeomSymmetricTensor<3>;
  using Facet = GeomFacet3d;

  //----------------------------------------------------------------------------
  // Constructors, assignment, destructor.
  //----------------------------------------------------------------------------
  GeomPolyhedron();

  // Note the following constructor constructs the convex hull of the given points,
  // meaning that the full set of points passed in may not appear in the vertices.
  GeomPolyhedron(const std::vector<Vector>& points);

  // This constructor takes a set of points that define facets of the surface
  // of the polyhedron.
  GeomPolyhedron(const std::vector<Vector>& points,
                 const std::vector<std::vector<unsigned> >& facetIndices);

  GeomPolyhedron(const GeomPolyhedron& rhs);
  GeomPolyhedron& operator=(const GeomPolyhedron& rhs);
  ~GeomPolyhedron() = default;

  // Test if the given point is internal to the polyhedron.
  bool contains(const Vector& point,
                const bool countBoundary = true,
                const double tol = 1.0e-8,
                const bool useAxom = false) const;

  // If we have a convex polyhedron we can do a faster containment test.
  bool convexContains(const Vector& point,
                      const bool countBoundary = true,
                      const double tol = 1.0e-8) const;

  // Test if we intersect another polyhedron.
  bool intersect(const GeomPolyhedron& rhs) const;
  bool convexIntersect(const GeomPolyhedron& rhs) const;
  bool intersect(const std::pair<Vector, Vector>& rhs) const;  // Another way of representing a box.

  // Test if we intersect a line segment (interior counts as intersection).
  bool intersect(const Vector& s0, const Vector& s1) const;

  // Find intersections of this polyhedron with a line segment (s0, s1).  Returns
  // both the intersection points and facet IDs.
  void intersections(const Vector& s0, const Vector& s1,
                     std::vector<unsigned>& facetIDs,
                     std::vector<Vector>& intersections) const;

  // Return the centroid of the vertices.
  Vector centroid() const;

  // Access the internal data.
  const std::vector<Vector>& vertices() const;
  const std::vector<Facet>& facets() const;
  const std::vector<Vector>& vertexUnitNorms() const;
  const std::vector<std::vector<unsigned> >& vertexFacetConnectivity() const;
  const std::vector<std::vector<unsigned> >& facetFacetConnectivity() const;
  const Vector& xmin() const;
  const Vector& xmax() const;
  
  // Get the edges as integer (node) pairs.
  std::vector<std::pair<unsigned, unsigned> > edges() const;

  // Spit out a vector<vector<unsigned> > and vector<Vector> that encode the facets.
  std::vector<std::vector<unsigned> > facetVertices() const;
  std::vector<Vector> facetNormals() const;
  std::vector<Vector> facetAreaVectors() const;
  std::vector<Vector> facetCentroids() const;

  // Useful facet properties.
  double facetArea(const unsigned facetID) const;
  Vector facetAreaNormal(const unsigned facetID) const;
  Vector facetCentroid(const unsigned facetID) const;

  // Reconstruct the internal data given a set of vertices and the vertex
  // indicies that define the facets.
  void reconstruct(const std::vector<Vector>& vertices,
                   const std::vector<std::vector<unsigned> >& facetVertices);

  // Compute the volume.
  double volume() const;

  // Find the facet closest to the given point.
  unsigned closestFacet(const Vector& p) const;

  // Compute the minimum distance to a point.
  double distance(const Vector& p,
                  const bool useAxom = false) const;

  // Find the point in the polyhedron closest to the given point.
  Vector closestPoint(const Vector& p) const;

  // Shift by a Vector delta.
  GeomPolyhedron& operator+=(const Vector& rhs);
  GeomPolyhedron& operator-=(const Vector& rhs);
  GeomPolyhedron operator+(const Vector& rhs) const;
  GeomPolyhedron operator-(const Vector& rhs) const;

  // Scale by a scalar.
  GeomPolyhedron& operator*=(const double rhs);
  GeomPolyhedron& operator/=(const double rhs);
  GeomPolyhedron operator*(const double rhs) const;
  GeomPolyhedron operator/(const double rhs) const;

  // Apply a tensor transformation
  GeomPolyhedron& transform(const Tensor& t);

  // Comparisons.
  bool operator==(const GeomPolyhedron& rhs) const;
  bool operator!=(const GeomPolyhedron& rhs) const;

  // Test if the polyhedron is convex.
  bool convex(const double tol = 1.0e-8) const;

  // Set the bounding box.
  void setBoundingBox();

  // Decompose the polyhedron into pyramids for each facet.
  GeomPolyhedron facetSubVolume(const unsigned facetID) const;

  // Decompose the polyhedron into tetrahedra.
  void decompose(std::vector<GeomPolyhedron>& subcells) const;

private:
  //--------------------------- Private Interface ---------------------------//
  std::vector<Vector> mVertices;
  std::vector<Facet> mFacets;
  std::vector<Vector> mVertexUnitNorms;
  std::vector<std::vector<unsigned> > mVertexFacetConnectivity, mFacetFacetConnectivity;
  Vector mXmin, mXmax, mCentroid;
  double mRinterior2;
  bool mConvex;
  mutable std::shared_ptr<axom::quest::InOutOctree<3>::SurfaceMesh> mSurfaceMeshPtr;
  mutable std::shared_ptr<axom::quest::InOutOctree<3>> mSurfaceMeshQueryPtr;
  mutable std::shared_ptr<axom::quest::SignedDistance<3>> mSignedDistancePtr;

  static FILE* mDevnull;

  // Construct the Axom representation
  void buildAxomData() const;
};

std::ostream& operator<<(std::ostream& os, const GeomPolyhedron& polygon);

}

#include "GeomPolyhedronInline.hh"

#endif

