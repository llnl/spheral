//---------------------------------Spheral++----------------------------------//
// testBoxIntersection
//
// Some simple methods specialized by dimension to rapidly test if two boxes
// intersect.
//
// Created by JMO, Sun Feb  6 13:44:49 PST 2000
//----------------------------------------------------------------------------//
#ifndef __Spheral_testBoxIntersection__
#define __Spheral_testBoxIntersection__

#include "config.hh"
#include "Geometry/Dimension.hh"

#include <utility>

namespace Spheral {

//------------------------------------------------------------------------------
// These methods specify the boxes by the their min and max coordinates, and
// assume that the boxes are aligned with the coordinate system.
//------------------------------------------------------------------------------
SPHERAL_HOST_DEVICE 
inline
bool
testBoxIntersection(const Dim<1>::Vector& xmin1, const Dim<1>::Vector& xmax1,
                    const Dim<1>::Vector& xmin2, const Dim<1>::Vector& xmax2,
                    const double tol = 1.0e-10) {
  if (xmax1.x() < (xmin2.x() - tol)) return false;
  if (xmax2.x() < (xmin1.x() - tol)) return false;
  return true;
}

SPHERAL_HOST_DEVICE 
inline
bool
testBoxIntersection(const Dim<2>::Vector& xmin1, const Dim<2>::Vector& xmax1,
                    const Dim<2>::Vector& xmin2, const Dim<2>::Vector& xmax2,
                    const double tol = 1.0e-10) {
  if (xmax1.x() < (xmin2.x() - tol)) return false;
  if (xmax2.x() < (xmin1.x() - tol)) return false;

  if (xmax1.y() < (xmin2.y() - tol)) return false;
  if (xmax2.y() < (xmin1.y() - tol)) return false;

  return true;
}

SPHERAL_HOST_DEVICE 
inline
bool
testBoxIntersection(const Dim<3>::Vector& xmin1, const Dim<3>::Vector& xmax1,
                    const Dim<3>::Vector& xmin2, const Dim<3>::Vector& xmax2,
                    const double tol = 1.0e-10) {
  if (xmax1.x() < (xmin2.x() - tol)) return false;
  if (xmax2.x() < (xmin1.x() - tol)) return false;

  if (xmax1.y() < (xmin2.y() - tol)) return false;
  if (xmax2.y() < (xmin1.y() - tol)) return false;

  if (xmax1.z() < (xmin2.z() - tol)) return false;
  if (xmax2.z() < (xmin1.z() - tol)) return false;

  return true;
}

//------------------------------------------------------------------------------
// Same as above using std::pairs to hold the box min/max's.
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE 
inline
bool
testBoxIntersection(const std::pair<typename Dimension::Vector, typename Dimension::Vector>& box1,
                    const std::pair<typename Dimension::Vector, typename Dimension::Vector>& box2,
                    const double tol = 1.0e-10) {
  return testBoxIntersection(box1.first, box1.second, box2.first, box2.second, tol);
}

//------------------------------------------------------------------------------
// Test if a point is contained in a box.
//------------------------------------------------------------------------------
SPHERAL_HOST_DEVICE 
inline
bool
testPointInBox(const Dim<1>::Vector& point,
               const Dim<1>::Vector& xmin, const Dim<1>::Vector& xmax,
               const double tol = 1.0e-10) {
  if (point.x() < (xmin.x() - tol)) return false;
  if (point.x() > (xmax.x() + tol)) return false;
  return true;
}

SPHERAL_HOST_DEVICE 
inline
bool
testPointInBox(const Dim<2>::Vector& point,
               const Dim<2>::Vector& xmin, const Dim<2>::Vector& xmax,
               const double tol = 1.0e-10) {
  if (point.x() < (xmin.x() - tol)) return false;
  if (point.x() > (xmax.x() + tol)) return false;

  if (point.y() < (xmin.y() - tol)) return false;
  if (point.y() > (xmax.y() + tol)) return false;

  return true;
}

SPHERAL_HOST_DEVICE 
inline
bool
testPointInBox(const Dim<3>::Vector& point,
               const Dim<3>::Vector& xmin, const Dim<3>::Vector& xmax,
               const double tol = 1.0e-10) {
  if (point.x() < (xmin.x() - tol)) return false;
  if (point.x() > (xmax.x() + tol)) return false;

  if (point.y() < (xmin.y() - tol)) return false;
  if (point.y() > (xmax.y() + tol)) return false;

  if (point.z() < (xmin.z() - tol)) return false;
  if (point.z() > (xmax.z() + tol)) return false;

  return true;
}

//------------------------------------------------------------------------------
// Same as above but specifying the box min/max using a std::pair.
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE 
inline
bool
testPointInBox(const typename Dimension::Vector& point,
               const std::pair<typename Dimension::Vector, typename Dimension::Vector>& box,
               const double tol = 1.0e-10) {
  return testPointInBox(point, box.first, box.second, tol);
}

}

#endif

