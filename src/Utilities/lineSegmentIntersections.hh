//---------------------------------Spheral++----------------------------------//
// lineSegmentIntersections
//
// Methods for computing intersections of line segments with other geometric
// types.
// These are largely adapted from code in "Computational Geometry in C" by
// Joseph O'Rourke.
//
// Created by JMO, Wed Feb  3 10:32:09 PST 2010
//----------------------------------------------------------------------------//
#ifndef __Spheral_lineSegmentIntersections_hh__
#define __Spheral_lineSegmentIntersections_hh__

#include "Geometry/Dimension.hh"
#include "safeInv.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Intersection of two line segments.
// The line segments are characterized by their endpoints: a_seg = (a0, a1)
//                                                         b_seg = (b0, b1)
// Return values are a pair<Vector, char>
//        The Vector is the intersection point (if any)
//        The char is a code characterizing the intersection:
//              "e" -> The segments colinearly overlap (edge)
//              "v" -> The endpoint of one segment lies on the other
//              "1" -> The segments intersect properly
//              "0" -> The segments do not intersect
//------------------------------------------------------------------------------
char
segmentSegmentIntersection(const Dim<1>::Vector& a0,
                           const Dim<1>::Vector& a1,
                           const Dim<1>::Vector& b0,
                           const Dim<1>::Vector& b1,
                           Dim<1>::Vector& result1,
                           Dim<1>::Vector& result2,
                           const double tol = 1.0e-8);

char
segmentSegmentIntersection(const Dim<2>::Vector& a0,
                           const Dim<2>::Vector& a1,
                           const Dim<2>::Vector& b0,
                           const Dim<2>::Vector& b1,
                           Dim<2>::Vector& result1,
                           Dim<2>::Vector& result2,
                           const double tol = 1.0e-8);

char
segmentSegmentIntersection(const Dim<3>::Vector& a0,
                           const Dim<3>::Vector& a1,
                           const Dim<3>::Vector& b0,
                           const Dim<3>::Vector& b1,
                           Dim<3>::Vector& result1,
                           Dim<3>::Vector& result2,
                           const double tol = 1.0e-8);

//------------------------------------------------------------------------------
// Intersection of a line segment with a plane.
// The line segment is characterized by it's endpoints:     seg = (s0, s1)
// The plane is characterized by three points in the plane: plane = (p0, p1, p2)
//
// Return values are a pair<Vector, char>
//        The Vector is the intersection point (if any)
//        The char is a code characterizing the intersection:
//              "p" -> The segment lies in the plane (plane)
//              "d" -> The p points do not define a unique plane (degenerate)
//              "1" -> The segment intersects the plane properly
//              "0" -> The segment does not intersect the plane
//------------------------------------------------------------------------------
// This version uses a (point, normal) plane representation.
char
segmentPlaneIntersection(const Dim<3>::Vector& s0,
                         const Dim<3>::Vector& s1,
                         const Dim<3>::Vector& point,
                         const Dim<3>::Vector& normal,
                         Dim<3>::Vector& result,
                         const double tol = 1.0e-8);

// This is the three plane points form.
char
segmentPlaneIntersection(const Dim<3>::Vector& s0,
                         const Dim<3>::Vector& s1,
                         const Dim<3>::Vector& p0,
                         const Dim<3>::Vector& p1,
                         const Dim<3>::Vector& p2,
                         Dim<3>::Vector& result, 
                         const double tol = 1.0e-8);

//------------------------------------------------------------------------------
// Intersection of a line segment with a polygonal section of a plane.
// The line segment is characterized by it's endpoints:     seg = (s0, s1)
// The polygonal section of the plane is specified by
// a series of points:                                    plane = pverts[0], pverts[1], ...
// Note there must be at least 3 non-collinear points, and they
// must be passed in in order to draw the polygon.
//
// Return values are a pair<Vector, char>
//        The Vector is the intersection point (if any)
//        The char is a code characterizing the intersection:
//              "p" -> The segment lies in the plane (plane)
//              "d" -> The pverts points do not define a unique plane (degenerate)
//              "1" -> The segment intersects the plane properly
//              "0" -> The segment does not intersect the plane
//------------------------------------------------------------------------------
char
segmentPlanarSectionIntersection(const Dim<3>::Vector& s0,
                                 const Dim<3>::Vector& s1,
                                 const std::vector<Dim<3>::Vector>& pverts,
                                 Dim<3>::Vector& result,
                                 const double tol = 1.0e-8);

//------------------------------------------------------------------------------
// Compute the minimum distance between two line segments. (a0,a1) -> (b0,b1)
//------------------------------------------------------------------------------
double
segmentSegmentDistance(const Dim<1>::Vector& a0,
                       const Dim<1>::Vector& a1,
                       const Dim<1>::Vector& b0,
                       const Dim<1>::Vector& b1);

double
segmentSegmentDistance(const Dim<2>::Vector& a0,
                       const Dim<2>::Vector& a1,
                       const Dim<2>::Vector& b0,
                       const Dim<2>::Vector& b1);

double
segmentSegmentDistance(const Dim<3>::Vector& a0,
                       const Dim<3>::Vector& a1,
                       const Dim<3>::Vector& b0,
                       const Dim<3>::Vector& b1);

//------------------------------------------------------------------------------
// Test if points a, b, & c are collinear.
// You should scale the tolerance handed in here!
//------------------------------------------------------------------------------
template<typename Vector>
inline
bool
collinear(const Vector& a, const Vector& b, const Vector& c,
          const double tol = 1.0e-10) {
  const Vector ba = b - a;
  const Vector cb = c - b;
  const double bamag2 = ba.magnitude2();
  const double cbmag2 = cb.magnitude2();
  if (fuzzyEqual(bamag2, 0.0, tol) or fuzzyEqual(cbmag2, 0.0, tol)) return true;
  return fuzzyEqual(std::abs(ba.dot(cb)), sqrt(bamag2*cbmag2), tol);
}

//------------------------------------------------------------------------------
// Test if point c is between a & b.
// You should scale the tolerance handed in here!
//------------------------------------------------------------------------------
template<typename Vector>
inline
bool
between(const Vector& a, const Vector& b, const Vector& c,
        const double tol = 1.0e-10) {

  // If c is equal to either endpoint, we count that as between.
  const Vector ca = c - a;
  const Vector cb = c - b;
  const double ca_mag = ca.magnitude();
  const double cb_mag = cb.magnitude();
  if (fuzzyEqual(ca_mag, 0.0, tol) or fuzzyEqual(cb_mag, 0.0, tol)) return true;

  // If a & b are the same point, but not c it's outside.
  const Vector ba = b - a;
  const double ba_mag = ba.magnitude();
  if (fuzzyEqual(ba_mag, 0.0, tol)) return false;

  // The points are distinct.
  const double thpt = ba.dot(ca)/ba_mag;
  return (fuzzyEqual(thpt, ca_mag, tol) and ca_mag <= ba_mag);
}

template<>
inline
bool
between<Dim<1>::Vector>(const Dim<1>::Vector& a, const Dim<1>::Vector& b, const Dim<1>::Vector& c,
                        const double tol) {

  // If c is equal to either endpoint, we count that as between.
  const double ca = c.x() - a.x();
  const double cb = c.x() - b.x();
  if (fuzzyEqual(std::abs(ca), 0.0, tol) or fuzzyEqual(std::abs(cb), 0.0, tol)) return true;

  // If a & b are the same point, but not c it's outside.
  const double ba = b.x() - a.x();
  if (fuzzyEqual(std::abs(ba), 0.0, tol)) return false;

  // The points are distinct.
  const double thpt = ba*ca/std::abs(ba);
  return (fuzzyEqual(thpt, std::abs(ca), tol) and std::abs(ca) <= std::abs(ba));
}

//------------------------------------------------------------------------------
// Check if two line segments intersect (does not compute intersection point).
// The line segments are characterized by their endpoints: a_seg = (a0, a1)
//                                                         b_seg = (b0, b1)
// The advantage of this method is that it should be faster if you don't need
// the actual intersection.
// Based on stuff from "Computational Geometry in C", Joseph O'Rourke
//------------------------------------------------------------------------------
template<typename Vector>
struct
segmentSegmentIntersectionTester {
  double area2(const Vector& a, const Vector& b, const Vector& c) const                  { return (b - a).cross(c - a).z(); }
  bool left(const Vector& a, const Vector& b, const Vector& c, const double tol) const   { return distinctlyGreaterThan(this->area2(a, b, c), 0.0, tol); }
  bool leftOn(const Vector& a, const Vector& b, const Vector& c, const double /*tol*/) const { return fuzzyGreaterThanOrEqual(this->area2(a, b, c), 0.0); }
  bool Xor(const bool x, const bool y) const                                             { return !x ^ !y; }
  bool intersectProp(const Vector& a, const Vector& b, const Vector& c, const Vector& d, const double tol) const {
    if (collinear(a, b, c, tol) or
        collinear(a, b, d, tol) or
        collinear(c, d, a, tol) or
        collinear(c, d, b, tol)) return false;
    return (this->Xor(this->left(a, b, c, tol), this->left(a, b, d, tol)) and
            this->Xor(this->left(c, d, a, tol), this->left(c, d, b, tol)));
  }
  bool operator()(const Vector& a, const Vector& b, const Vector& c, const Vector& d, const double tol) const {
    if (this->intersectProp(a, b, c, d, tol)) return true;
    if (between(a, b, c, tol) or
        between(a, b, d, tol) or
        between(c, d, a, tol) or
        between(c, d, b, tol)) return true;
    return false;
  }
};

// In 3D we have to pick a positive direction for area2.
// Strictly speaking this no longer returns the area, but all we really care about is the sign.
template<>
inline
double
segmentSegmentIntersectionTester<Dim<3>::Vector>::
area2(const Dim<3>::Vector& a, const Dim<3>::Vector& b, const Dim<3>::Vector& c) const {
  typedef Dim<3>::Vector Vector;
  const Vector ca = c - a;
  const Vector norm = (b - a).cross(c - a);
  const Vector nhat = norm.unitVector();
  if (std::abs(nhat.z()) > 0.5) {
    return ca.dot(Vector(0,0,1));
  } else if (std::abs(nhat.y()) > 0.5) {
    return ca.dot(Vector(0,1,0));
  } else {
    return ca.dot(Vector(1,0,0));
  }
}

// A functional front end for the segment intersection test.
template<typename Vector>
bool 
segmentSegmentIntersection(const Vector& a0,
                           const Vector& a1,
                           const Vector& b0,
                           const Vector& b1,
                           const double tol = 1.0e-10) {
  return segmentSegmentIntersectionTester<Vector>()(a0, a1, b0, b1, tol);
}

//------------------------------------------------------------------------------
// Check if a line segment intersects a polygonal planar section in 3D.
// The line segment is characterized by its endpoints: a_seg = (a0, a1)
// The planar section is characterized the a set of coplanar vertices: vertices.
//
// The advantage of this method is that it should be faster if you don't need
// the actual intersection.
// 
// The second version allows the planar vertices to be represented as a subset
// of the set passed in.
//------------------------------------------------------------------------------
bool
segmentPlaneIntersection(const Dim<3>::Vector& a0,
                         const Dim<3>::Vector& a1,
                         const std::vector<Dim<3>::Vector>& vertices,
                         const Dim<3>::Vector& normal,
                         const double tol = 1.0e-10);

bool
segmentPlaneIntersection(const Dim<3>::Vector& a0,
                         const Dim<3>::Vector& a1,
                         const std::vector<Dim<3>::Vector>& vertices,
                         const std::vector<unsigned>& ipoints,
                         const Dim<3>::Vector& normal,
                         const double tol = 1.0e-10);
}

#endif
