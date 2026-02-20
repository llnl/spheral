//------------------------------------------------------------------------------
// Compute toroidal volume from rotating 2D areal shapes around an axis
//------------------------------------------------------------------------------
#ifndef _Spheral_toroidalVolume_
#define _Spheral_toroidalVolume_

#include "Utilities/DBC.hh"
#include "Utilities/FastMath.hh"
#include <cmath>

namespace Spheral {

//------------------------------------------------------------------------------
// Spin a square around the axis, taking into account of whether the bottom of
// the square touches (or goes below) the axis.
//------------------------------------------------------------------------------
inline
double
cylindricalToroidalVolume(const double d,   // length of one side of the square area we're rotating
                          const double r) { // radius of circle centroid from the axis of rotation
  REQUIRE(d >= 0.0);
  REQUIRE(r >= 0.0);
  return M_PI*(FastMath::square(r + 0.5*d) - FastMath::square(std::max(0.0, r - 0.5*d)))*d;
}

//------------------------------------------------------------------------------
// Compute the volume of a toroid, taking into a account when the circular
// cross-section of the torus contains and is clipped by the axis of rotation
// (i.e., a spindle toroid).
//------------------------------------------------------------------------------
inline
double
circularToroidalVolume(const double R,   // radius of circle we are spinning around the axis
                       const double r) { // radius of circle centroid from the axis of rotation
  REQUIRE(R >= 0.0);
  REQUIRE(r >= 0.0);

  if (r >= R) {
    return 2.0*M_PI*M_PI*R*R*r;
  } else {
    const auto x = std::sqrt(R*R - r*r);
    return 2.0/3.0*M_PI*(2.0*R*R + r*r)*x + M_PI*R*R*r*(M_PI + 2.0*std::atan2(r, x));
  }
}

}

#endif
