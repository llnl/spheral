//---------------------------------Spheral++----------------------------------//
// Box1d
//
// The Box is specified by a central position and an extent.  Note that the
// extent here is half the size of the box (i.e., distance from the center to
// the edge of the box).
// This thing is intended to be similar in interface to the Polygon and 
// Polyhedron in 2-D and 3-D.
//
// Created by JMO, Sun Jan 24 21:48:13 PST 2010
//----------------------------------------------------------------------------//
#ifndef __Spheral_Box1d__
#define __Spheral_Box1d__

#include <vector>

#include "GeomFacet1d.hh"

template<int nDim> class GeomVector;

namespace Spheral {

class Box1d {
public:
  //--------------------------- Public Interface ---------------------------//
  using Vector = GeomVector<1>;
  using Tensor = GeomTensor<1>;
  using SymTensor = GeomSymmetricTensor<1>;
  using Facet = GeomFacet1d;

  //----------------------------------------------------------------------------
  // Constructors, assignment, destructor.
  //----------------------------------------------------------------------------
  Box1d();
  Box1d(const std::vector<Vector>& points);
  Box1d(const std::vector<Vector>& points,
        const std::vector<std::vector<unsigned> >& facetIndices);
  Box1d(const Vector& center, const double extent);
  Box1d(const Box1d& rhs);
  Box1d& operator=(const Box1d& rhs);
  ~Box1d();

  // Test if the given point is internal to the box.
  bool contains(const Vector& point,
                const bool countBoundary = true,
                const double tol = 1.0e-8) const;

  bool convexContains(const Vector& point,
                      const bool countBoundary = true,
                      const double tol = 1.0e-8) const;

  // Test if we intersect another box.
  bool intersect(const Box1d& rhs) const;
  bool convexIntersect(const Box1d& rhs) const;
  bool intersect(const std::pair<Vector, Vector>& rhs) const;  // Another way of representing a box.

  // Test if we intersect a line segment (interior counts as intersection).
  bool intersect(const Vector& s0, const Vector& s1) const;

  // Access the attributes.
  Vector& center();
  const Vector& center() const;
  void center(const Vector& val);

  double& extent();
  double extent() const;
  void extent(double val);

  const Vector& xmin() const;
  const Vector& xmax() const;

  Vector centroid() const { return mCenter; }

  const std::vector<Vector>& vertices() const;
  const std::vector<Facet>& facets() const;
  std::vector<std::vector<unsigned> > facetVertices() const;

  // Useful facet properties.
  double facetArea(const unsigned facetID) const;
  Vector facetAreaNormal(const unsigned facetID) const;

  // Compute the minimum distance to a point.
  double distance(const Vector& p) const;

  // Find the point in the box closest to the given point.
  Vector closestPoint(const Vector& p) const;

  // Compute the volume.
  double volume() const;

  // For compatibility with polygon and polyhedron
  void decompose(std::vector<Box1d>& subcells) const;

  // Shift by a Vector delta.
  Box1d& operator+=(const Vector& rhs);
  Box1d& operator-=(const Vector& rhs);
  Box1d operator+(const Vector& rhs) const;
  Box1d operator-(const Vector& rhs) const;

  // Scale by a scalar.
  Box1d& operator*=(const double rhs);
  Box1d& operator/=(const double rhs);
  Box1d operator*(const double rhs) const;
  Box1d operator/(const double rhs) const;

  // Comparisons.
  bool operator==(const Box1d& rhs) const;
  bool operator!=(const Box1d& rhs) const;

private:
  //--------------------------- Private Interface ---------------------------//
  // Our underlying representation.
  Vector mCenter;
  double mExtent;
  std::vector<Vector> mVertices;
  mutable std::vector<Facet> mFacets; // for now, just create this when we need it
};

std::ostream& operator<<(std::ostream& os, const Box1d& box);

}

#include "Box1dInline.hh"

#endif
