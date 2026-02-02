//---------------------------------Spheral++----------------------------------//
// newtonRaphson -- implements a simple version of the Newton-Raphson root
// finding algorithm, combined with a bisection back-up.
//
// Based on rtsafe from Numerical recipes.
//
// Created by JMO, Thu Sep 23 22:57:33 PDT 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral_newtonRaphson_hh__
#define __Spheral_newtonRaphson_hh__

#include <map>
#include "Utilities/DBC.hh"
#include "Utilities/SpheralFunctions.hh"

namespace Spheral {

template<typename Function>
inline
double
newtonRaphson(const Function& functor,
              double x1,
              double x2,
              const double xaccuracy = 1.0e-15,
              const double yaccuracy = 1.0e-15,
              const unsigned maxIterations = 100) {

  typedef std::pair<double, double> EvalType;

  // Initialize values for the function and its derivative.
  EvalType xminValues = functor(x1);
  EvalType xmaxValues = functor(x2);

  // Is the root already at the min or max range?
  if (fuzzyEqual(xminValues.first, 0.0, yaccuracy)) return x1;
  if (fuzzyEqual(xmaxValues.first, 0.0, yaccuracy)) return x2;

  // Make sure the root is bracketed by the input range.
  VERIFY2(distinctlyLessThan(xminValues.first * xmaxValues.first, 0.0, yaccuracy),
          "newtonRaphson: root must be bracketed by input range:  " << xminValues.first << " " << xmaxValues.first);

  // Initialize the searching parameters.
  double xl, xh;
  if (xminValues.first < 0.0) {
    xl = x1;
    xh = x2;
  } else {
    CHECK(xminValues.first > 0.0 && xmaxValues.first < 0.0);
    xl = x2;
    xh = x1;
  }
  double rootSafe = 0.5*(x1 + x2);
  double dxold = std::abs(x2 - x1);
  double dx = dxold;
  EvalType fdf = functor(rootSafe);
  double f = fdf.first;
  double df = fdf.second;

  // Iterate until we either converge or achieve the desired accuracy.
  unsigned iter = 0;
  while (iter < maxIterations) {
    ++iter;

    // Bisect if Newton out of range or not decreasing fast enough.
    if (((rootSafe - xh)*df - f)*((rootSafe - xl)*df - f) > 0.0 ||
        std::abs(2.0*f) > std::abs(dxold*df)) {
      dxold = dx;
      dx = 0.5*(xh - xl);
      rootSafe = xl + dx;
      if (fuzzyEqual(xl, rootSafe, xaccuracy)) return rootSafe;

    } else {
      // Take a Newton-Raphson step.
      CHECK(!fuzzyEqual(df, 0.0));
      dxold = dx;
      dx = f/df;
      double tmp = rootSafe;
      rootSafe -= dx;
      if (fuzzyEqual(tmp, rootSafe, xaccuracy)) return rootSafe;
    }
    if (std::abs(dx) <= xaccuracy) return rootSafe;

    fdf = functor(rootSafe);
    f = fdf.first;
    df = fdf.second;
    if (f < 0.0) {
      xl = rootSafe;
    } else {
      xh = rootSafe;
    }

  }

  VERIFY2(false, "newtonRaphson: failed to converge!");

  return 0;
}

}

#endif
