//---------------------------------Spheral++----------------------------------//
// GeomPolygon -- Geometric polygon class.
//
// A 2-D structure representing a polygon as a collection of GeomFacets.
//
// Created by JMO, Thu Jan 28 11:03:27 PST 2010
//----------------------------------------------------------------------------//
#ifndef __Spheral_GeomPolygon__
#define __Spheral_GeomPolygon__

#include <vector>

#include "GeomVector.hh"
#include "GeomFacet2d.hh"

namespace Spheral {

class GeomPolygon {
public:
  //--------------------------- Public Interface ---------------------------//
  using Vector = GeomVector<2>;
  using Tensor = GeomTensor<2>;
  using SymTensor = GeomSymmetricTensor<2>;
  using Facet = GeomFacet2d;

  //----------------------------------------------------------------------------
  // Constructors, assignment, destructor.
  //----------------------------------------------------------------------------
  GeomPolygon();

  // Note the following constructor constructs the convex hull of the given points,
  // meaning that the full set of points passed in may not appear in the vertices.
  GeomPolygon(const std::vector<Vector>& points);

  // This constructor takes a set of points that define facets of the surface
  // of the polygon.
  GeomPolygon(const std::vector<Vector>& points,
              const std::vector<std::vector<unsigned> >& facetIndices);

  GeomPolygon(const GeomPolygon& rhs);
  GeomPolygon& operator=(const GeomPolygon& rhs);
  ~GeomPolygon();

  // Test if the given point is internal to the polygon.
  // This method works for any polygon.
  bool contains(const Vector& point,
                const bool countBoundary = true,
                const double tol = 1.0e-8) const;

  // If we have a convex polygon we can do a faster containment test.
  bool convexContains(const Vector& point,
                      const bool countBoundary = true,
                      const double tol = 1.0e-8) const;

  // Test if we intersect another polygon.
  bool intersect(const GeomPolygon& rhs) const;
  bool convexIntersect(const GeomPolygon& rhs) const;
  bool intersect(const std::pair<Vector, Vector>& rhs) const;  // Another way of representing a box.

  // Test if we intersect a line segment (interior counts as intersection).
  bool intersect(const Vector& s0, const Vector& s1) const;

  // Find intersections of this polygon with a line segment (s0, s1).  Returns
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

  // Spit out a vector<vector<unsigned> > that encodes the facets.
  std::vector<std::vector<unsigned> > facetVertices() const;

  // Useful facet properties.
  double facetArea(const unsigned facetID) const;
  Vector facetAreaNormal(const unsigned facetID) const;

  // Reconstruct the internal data given a set of vertices and the vertex
  // indicies that define the facets.
  void reconstruct(const std::vector<Vector>& vertices,
                   const std::vector<std::vector<unsigned> >& facetVertices);

  // Compute the volume.
  double volume() const;

  // Find the facet closest to the given point.
  unsigned closestFacet(const Vector& p) const;

  // Compute the minimum distance to a point.
  double distance(const Vector& p) const;

  // Find the point in the polyhedron closest to the given point.
  Vector closestPoint(const Vector& p) const;

  // Shift by a Vector delta.
  GeomPolygon& operator+=(const Vector& rhs);
  GeomPolygon& operator-=(const Vector& rhs);
  GeomPolygon operator+(const Vector& rhs) const;
  GeomPolygon operator-(const Vector& rhs) const;

  // Scale by a scalar.
  GeomPolygon& operator*=(const double rhs);
  GeomPolygon& operator/=(const double rhs);
  GeomPolygon operator*(const double rhs) const;
  GeomPolygon operator/(const double rhs) const;

  // Apply a tensor transformation
  GeomPolygon& transform(const Tensor& t);

  // Comparisons.
  bool operator==(const GeomPolygon& rhs) const;
  bool operator!=(const GeomPolygon& rhs) const;

  // Test if the polygon is convex.
  bool convex(const double tol = 1.0e-8) const;

  // Set the bounding box.
  void setBoundingBox();

  // Decompose the polygon into triangles for each facet.
  GeomPolygon facetSubVolume(const unsigned facetID) const;

  // Decompose the polygon into triangles.
  void decompose(std::vector<GeomPolygon>& subcells) const;
  
private:
  //--------------------------- Private Interface ---------------------------//
  std::vector<Vector> mVertices;
  std::vector<Facet> mFacets;
  std::vector<Vector> mVertexUnitNorms;
  std::vector<std::vector<unsigned> > mVertexFacetConnectivity, mFacetFacetConnectivity;
  Vector mXmin, mXmax;
  bool mConvex;

  static FILE* mDevnull;
};

std::ostream& operator<<(std::ostream& os, const GeomPolygon& polygon);

}

#include "GeomPolygonInline.hh"

#endif

