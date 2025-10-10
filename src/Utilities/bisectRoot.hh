//---------------------------------Spheral++----------------------------------//
// bisectRoot -- implements the straight-up bisection root-finding search.
//           This is a stripped down version of our newtonRaphson routine.
//
// Created by JMO, Thu Oct 27 10:20:22 PDT 2016
//----------------------------------------------------------------------------//
#ifndef __Spheral_bisectRoot_hh__
#define __Spheral_bisectRoot_hh__

#include "Utilities/DBC.hh"
#include "Utilities/SpheralFunctions.hh"

#include <iostream>

namespace Spheral {

template<typename Function>
inline
double
bisectRoot(const Function& functor,
           double xmin,
           double xmax,
           const double xaccuracy = 1.0e-15,
           const double yaccuracy = 1.0e-10,
           const unsigned maxIterations = 100u,
           const bool verbose = false) {

  // Initialize values for the function and its derivative.
  double fxmin = functor(xmin);
  double fxmax = functor(xmax);

  // Is the root already at the min or max range?
  if (fuzzyEqual(fxmin, 0.0, yaccuracy)) return xmin;
  if (fuzzyEqual(fxmax, 0.0, yaccuracy)) return xmax;

  // Make sure the root is bracketed by the input range.
  VERIFY2(fxmin*fxmax < 0.0,
          "bisectRoot: root must be bracketed by input range:  f(" << xmin << ")=" << fxmin << " ; f(" << xmax << ")=" << fxmax);

  // Initialize the searching parameters.
  double x0, x1;
  if (fxmin < 0.0) {
    CHECK(fxmax > 0.0);
    x0 = xmin;
    x1 = xmax;
  } else {
    CHECK(fxmin > 0.0 and fxmax < 0.0);
    x0 = xmax;
    x1 = xmin;
  }

  // Iterate until we either converge or achieve the desired accuracy.
  double rootSafe, dx, f;
  unsigned iter = 0u;
  while (iter++ < maxIterations) {
    if (verbose) std::cout << "bisectRoot " << iter << ": x in [" << x0 << " " << x1 << "] : F(x) in [" << functor(x0) << " " << functor(x1) << "]" << std::endl;
    dx = 0.5*(x1 - x0);
    rootSafe = x0 + dx;
    if (std::abs(dx) <= xaccuracy) return rootSafe;

    f = functor(rootSafe);
    if (fuzzyEqual(f, 0.0, yaccuracy)) return rootSafe;
    if (f < 0.0) {
      x0 = rootSafe;
    } else {
      x1 = rootSafe;
    }

  }

  VERIFY2(false, "bisectRoot: failed to converge!");

  return 0;
}

}

#endif
