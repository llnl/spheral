//------------------------------------------------------------------------------
// A collection of Spheral wrappers for using R2D/R3D methods.
//------------------------------------------------------------------------------
#include "Geometry/GeomPlane.hh"
#include "Utilities/CounterClockwiseComparator.hh"
#include "Utilities/sort_permutation.hh"

#include "Geometry/PolyClipperUtilities.hh"

#include <algorithm>
#include <set>
#include <iostream>
#include <iterator>
using std::vector;
using std::set;
using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {

//------------------------------------------------------------------------------
// Clip a polygon by a series of planes
//------------------------------------------------------------------------------
Dim<2>::FacetedVolume clipFacetedVolume(const Dim<2>::FacetedVolume& poly,
                                        const std::vector<GeomPlane<Dim<2> > >& planes) {

  typedef Dim<2>::Vector Vector;
  typedef Dim<2>::FacetedVolume FacetedVolume;

  // Is there anything to do?
  const unsigned nplanes = planes.size();
  if (nplanes == 0) return poly;

  // Construct the PolyClipper version of our polygon.
  PolyClipperPolygon poly2d;
  convertToPolyClipper(poly2d, poly);

  // Now the planes.
  vector<PolyClipperPlane2d> planes2d(nplanes);
  for (auto i = 0u; i < nplanes; ++i) {
    const Vector& nhat = planes[i].normal();
    const Vector& p = planes[i].point();
    planes2d[i].normal = nhat;
    planes2d[i].dist = -p.dot(nhat);
  }

  // Sort the planes by distance -- lets us clip more efficiently.
  std::sort(planes2d.begin(), planes2d.end(), [](const PolyClipperPlane2d& lhs, const PolyClipperPlane2d& rhs) { return lhs.dist < rhs.dist; });

  // Do the deed.
  PolyClipper::clipPolygon(poly2d, planes2d);

  // Copy back to a Spheral polygon.
  FacetedVolume result;
  double area;
  Vector cent;
  PolyClipper::moments(area, cent, poly2d);
  const auto tol = 1.0e-10 * area;
  PolyClipper::collapseDegenerates(poly2d, tol);
  convertFromPolyClipper(result, poly2d);
  return result;
}

//------------------------------------------------------------------------------
// Clip a polyhedron by a series of planes
//------------------------------------------------------------------------------
Dim<3>::FacetedVolume clipFacetedVolume(const Dim<3>::FacetedVolume& poly,
                                        const std::vector<GeomPlane<Dim<3> > >& planes) {

  typedef Dim<3>::Vector Vector;
  typedef Dim<3>::FacetedVolume FacetedVolume;

  // Is there anything to do?
  const unsigned nplanes = planes.size();
  if (nplanes == 0) return poly;

  // Construct the PolyClipper version of our polygon.
  PolyClipperPolyhedron poly3d;
  convertToPolyClipper(poly3d, poly);

  // Now the planes.
  vector<PolyClipperPlane3d> planes3d(nplanes);
  for (auto i = 0u; i < nplanes; ++i) {
    const Vector& nhat = planes[i].normal();
    const Vector& p = planes[i].point();
    planes3d[i].normal = nhat;
    planes3d[i].dist = -p.dot(nhat);
  }

  // Sort the planes by distance -- lets us clip more efficiently.
  std::sort(planes3d.begin(), planes3d.end(), [](const PolyClipperPlane3d& lhs, const PolyClipperPlane3d& rhs) { return lhs.dist < rhs.dist; });

  // Do the deed.
  PolyClipper::clipPolyhedron(poly3d, planes3d);

  // Copy back to a Spheral polyhedron.
  FacetedVolume result;
  double area;
  Vector cent;
  PolyClipper::moments(area, cent, poly3d);
  const auto tol = 1.0e-10 * area;
  PolyClipper::collapseDegenerates(poly3d, tol);
  convertFromPolyClipper(result, poly3d);
  return result;
}

//------------------------------------------------------------------------------
// Volume of the cliped polygon
//------------------------------------------------------------------------------
double clippedVolume(const Dim<2>::FacetedVolume& poly,
                     const std::vector<GeomPlane<Dim<2> > >& planes) {

  typedef Dim<2>::Vector Vector;

  // Is there anything to do?
  const unsigned nplanes = planes.size();
  if (nplanes == 0) return poly.volume();

  // Construct the PolyClipper version of our polygon.
  PolyClipperPolygon poly2d;
  convertToPolyClipper(poly2d, poly);

  // Now the planes.
  vector<PolyClipperPlane2d> planes2d(nplanes);
  for (auto i = 0u; i < nplanes; ++i) {
    const Vector& nhat = planes[i].normal();
    const Vector& p = planes[i].point();
    planes2d[i].normal = nhat;
    planes2d[i].dist = -p.dot(nhat);
  }

  // Sort the planes by distance -- lets us clip more efficiently.
  std::sort(planes2d.begin(), planes2d.end(), [](const PolyClipperPlane2d& lhs, const PolyClipperPlane2d& rhs) { return lhs.dist < rhs.dist; });

  // Do the deed.
  PolyClipper::clipPolygon(poly2d, planes2d);

  // Return the volume.
  double area;
  Vector cent;
  PolyClipper::moments(area, cent, poly2d);
  return area;
}

//------------------------------------------------------------------------------
// Volume of the clipped polyhedron
//------------------------------------------------------------------------------
double clippedVolume(const Dim<3>::FacetedVolume& poly,
                     const std::vector<GeomPlane<Dim<3> > >& planes) {

  typedef Dim<3>::Vector Vector;

  // Is there anything to do?
  const unsigned nplanes = planes.size();
  if (nplanes == 0) return poly.volume();

  // Construct the PolyClipper version of our polyhedron.
  PolyClipperPolyhedron poly3d;
  convertToPolyClipper(poly3d, poly);

  // Now the planes.
  vector<PolyClipperPlane3d> planes3d(nplanes);
  for (auto i = 0u; i < nplanes; ++i) {
    const Vector& nhat = planes[i].normal();
    const Vector& p = planes[i].point();
    planes3d[i].normal = nhat.x();
    planes3d[i].dist = -p.dot(nhat);
  }

  // Sort the planes by distance -- lets us clip more efficiently.
  std::sort(planes3d.begin(), planes3d.end(), [](const PolyClipperPlane3d& lhs, const PolyClipperPlane3d& rhs) { return lhs.dist < rhs.dist; });

  // Do the deed.
  PolyClipper::clipPolyhedron(poly3d, planes3d);

  // Return the volume.
  double vol;
  Vector cent;
  PolyClipper::moments(vol, cent, poly3d);
  return vol;
}

}
