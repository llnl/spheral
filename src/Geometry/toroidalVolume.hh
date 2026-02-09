//------------------------------------------------------------------------------
// Compute the volume of a toroid, taking into a account when the circular
// cross-section of the torus contains and is clipped by the axis of rotation
// (i.e., a spindle toroid).
//------------------------------------------------------------------------------
#ifndef _Spheral_toroidalVolume_
#define _Spheral_toroidalVolume_

#include "Utilities/DBC.hh"
#include <cmath>

namespace Spheral {

inline
double
toroidalVolume(const double R,   // radius of circle we are spinning around the axis
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
