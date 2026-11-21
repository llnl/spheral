//---------------------------------Spheral++----------------------------------//
// pointDistances
//
// Assorted utilities for finding the distances from a point to things.
//
// Created by JMO, Mon Apr 11 21:19:13 PDT 2011
//----------------------------------------------------------------------------//
#ifndef __Spheral_pointDistances_hh__
#define __Spheral_pointDistances_hh__

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "safeInv.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Find the point on a line segment (a0,a1) closest to point (p).
//------------------------------------------------------------------------------
template<typename Vector>
SPHERAL_HOST_DEVICE 
inline
Vector
closestPointOnSegment(const Vector& p,
                      const Vector& a0,
                      const Vector& a1) {
  const auto a01 = a1 - a0;
  const auto a01mag = a01.magnitude();
  const auto ahat = a01*safeInv(a01mag);
  const auto a0p = p - a0;
  const auto ptest = a0p.dot(ahat);
  return a0 + std::max(0.0, std::min(a01mag,ptest))*ahat;
}

// This version returns true if the closest point forms a right angle between
// the line segment and (point, closest) pairs.
template<typename Vector>
SPHERAL_HOST_DEVICE 
inline
bool
closestPointOnSegment(const Vector& p,
                      const Vector& a0,
                      const Vector& a1,
                      Vector& result) {
  const auto a01 = a1 - a0;
  const auto a01mag = a01.magnitude();
  const auto ahat = a01*safeInv(a01mag);
  const auto a0p = p - a0;
  const auto ptest = a0p.dot(ahat);
  result = a0 + std::max(0.0, std::min(a01mag,ptest))*ahat;
  return (ptest >= 0.0 and ptest <= a01mag);
}

//------------------------------------------------------------------------------
// Compute the distance between a point and a plane.
//------------------------------------------------------------------------------
template<typename Vector>
SPHERAL_HOST_DEVICE 
inline
double
pointPlaneDistance(const Vector& point,
                   const Vector& origin,
                   const Vector& unitNormal) {
   REQUIRE(fuzzyEqual(unitNormal.magnitude2(), 1.0, 1.0e-10));
   return std::abs((point - origin).dot(unitNormal));
}

//------------------------------------------------------------------------------
// Find the point on a plane closest to the given point.
//------------------------------------------------------------------------------
SPHERAL_HOST_DEVICE 
inline
Dim<3>::Vector
closestPointOnPlane(const Dim<3>::Vector& p,
                    const Dim<3>::Vector& origin,
                    const Dim<3>::Vector& unitNormal) {
  REQUIRE(fuzzyEqual(unitNormal.magnitude2(), 1.0, 1.0e-10));
  return p + (origin - p).dot(unitNormal)*unitNormal;
}

}

#endif
