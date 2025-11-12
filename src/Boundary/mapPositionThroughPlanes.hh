//---------------------------------Spheral++----------------------------------//
// Function to take a position, and map it through the given enter/exit plane
// pair and return the mapped position.
//
// Created by JMO, Fri Dec  6 13:48:48 PST 2002
//----------------------------------------------------------------------------//
#ifndef __Spheral_mapPositionThroughPlanes__
#define __Spheral_mapPositionThroughPlanes__

#include "Geometry/GeomPlane.hh"
#include "Utilities/SpheralFunctions.hh"
#include "Utilities/DBC.hh"

namespace Spheral {
template<typename Dimension>
typename Dimension::Vector
mapPositionThroughPlanes(const typename Dimension::Vector& position,
                         const GeomPlane<Dimension>& enterPlane,
                         const GeomPlane<Dimension>& exitPlane) {
  REQUIRE(enterPlane.valid() && exitPlane.valid());
  typedef typename Dimension::Vector Vector;
  const Vector p0 = exitPlane.closestPointOnPlane(position);
  const double s = enterPlane.signedDistance(position);
  return p0 - s*exitPlane.normal();
  // const Vector deltaEnter = (position - enterPlane.point()).dot(enterPlane.normal())*enterPlane.normal();
  // const Vector deltaExit = (position - exitPlane.point()).dot(exitPlane.normal())*exitPlane.normal();
  // const double sign = sgn((position - enterPlane.point()).dot(enterPlane.normal()));
  // return position - deltaExit - sign*deltaEnter.magnitude()*exitPlane.normal();
}
}

#endif
