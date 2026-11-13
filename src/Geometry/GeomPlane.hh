//---------------------------------Spheral++----------------------------------//
// GeomPlane -- Geometric Plane Class.
//
// Created by JMO, Thu Feb 24 17:26:32 PST 2000
//----------------------------------------------------------------------------//

#ifndef __Spheral_GeomPlane__
#define __Spheral_GeomPlane__

#include "config.hh"

#include <iostream>

namespace Spheral {

template<typename Dimension>
class GeomPlane {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  // Constructors.
  SPHERAL_HOST_DEVICE GeomPlane() = default;
  SPHERAL_HOST_DEVICE GeomPlane(const GeomPlane& rhs) = default;
  SPHERAL_HOST_DEVICE GeomPlane(GeomPlane&& rhs) = default;
  SPHERAL_HOST_DEVICE GeomPlane(const Vector& point, const Vector& normal);
  SPHERAL_HOST_DEVICE GeomPlane& operator=(const GeomPlane& rhs) = default;

  // This special constructor finds a plane fitted to a collection of points.
  template<template<typename ...> typename Container>
  explicit GeomPlane(const Container<Vector>& points);

  // Destructor.
  SPHERAL_HOST_DEVICE ~GeomPlane() = default;

  // Access the point and normal which define the plane.
  SPHERAL_HOST_DEVICE const Vector& point()                const { return mPoint; }
  SPHERAL_HOST_DEVICE void point(const Vector& x)                { mPoint = x; }

  SPHERAL_HOST_DEVICE const Vector& normal()               const { return mNormal; }
  SPHERAL_HOST_DEVICE void normal(const Vector& x)               { mNormal = x.unitVector(); }

  // Negative operator (reverses sign of normal).
  SPHERAL_HOST_DEVICE GeomPlane operator-() const;

  // Calculate the signed distance between a given point and the plane.
  SPHERAL_HOST_DEVICE Scalar signedDistance(const Vector& p) const;

  // Calculate the minimum distance between a given point and the plane.
  SPHERAL_HOST_DEVICE Scalar minimumDistance(const Vector& p) const;

  // Closest point on the plane to a point.
  SPHERAL_HOST_DEVICE Vector closestPointOnPlane(const Vector& p) const;

  // Various tests which can be applied between planes.
  SPHERAL_HOST_DEVICE bool parallel(const GeomPlane& rhs) const;
  SPHERAL_HOST_DEVICE bool operator==(const GeomPlane& rhs) const;
  SPHERAL_HOST_DEVICE bool operator!=(const GeomPlane& rhs) const;
  SPHERAL_HOST_DEVICE bool operator<(const GeomPlane& rhs) const;

  // Some tests which can applied between a plane and a point.
  // These are only meant as tests for whether a point lies "above" or "below" the
  // plane, as represented by the direction of the normal.
  SPHERAL_HOST_DEVICE bool operator==(const Vector& point) const;
  SPHERAL_HOST_DEVICE bool operator!=(const Vector& point) const;
  SPHERAL_HOST_DEVICE bool operator>(const Vector& point) const;
  SPHERAL_HOST_DEVICE bool operator<(const Vector& point) const;
  SPHERAL_HOST_DEVICE bool operator>=(const Vector& point) const;
  SPHERAL_HOST_DEVICE bool operator<=(const Vector& point) const;

  // Test if the plane is below, equal to or above the point (-1, 0, 1).
  SPHERAL_HOST_DEVICE int compare(const Vector& point) const;

  // The usual validity test.
  SPHERAL_HOST_DEVICE bool valid() const;

private:
  //--------------------------- Private Interface ---------------------------//
  Vector mPoint;
  Vector mNormal;
};

// Ostream operators.
template<typename Dimension> std::istream& operator>>(std::istream& is, GeomPlane<Dimension>& vec);
template<typename Dimension> std::ostream& operator<<(std::ostream& os, const GeomPlane<Dimension>& vec);

}

#include "GeomPlaneInline.hh"

// Declare specialized methods.
namespace Spheral {

template<> SPHERAL_HOST_DEVICE Dim<1>::Vector GeomPlane<Dim<1>>::closestPointOnPlane(const Dim<1>::Vector&) const;
template<> SPHERAL_HOST_DEVICE Dim<2>::Vector GeomPlane<Dim<2>>::closestPointOnPlane(const Dim<2>::Vector&) const;
template<> SPHERAL_HOST_DEVICE Dim<3>::Vector GeomPlane<Dim<3>>::closestPointOnPlane(const Dim<3>::Vector&) const;

// template<> GeomPlane<Dim<1> >::GeomPlane(const std::vector<Dim<1>::Vector>& points);
// template<> GeomPlane<Dim<2> >::GeomPlane(const std::vector<Dim<2>::Vector>& points);
// template<> GeomPlane<Dim<3> >::GeomPlane(const std::vector<Dim<3>::Vector>& points);

}

#endif
