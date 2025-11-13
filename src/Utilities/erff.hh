//------------------------------------------------------------------------------
// Compute the error function (erff) and complementary error function (erffc)
// for a real argument.
// 
// Based on the routines in Numerical Recipes in C, pp 220-221.
//------------------------------------------------------------------------------
#ifndef __Spheral_erff__
#define __Spheral_erff__

#include <cmath>

namespace Spheral {

//------------------------------------------------------------------------------
// Complementary error function.
// Note this approximation is only good to a fractional error of 1.2e-7.
//------------------------------------------------------------------------------
inline
double
erffc(double x) {
  const double z = std::abs(x);
  const double t = 1.0/(1.0 + 0.5*z);
  const double ans = t*exp(-z*z - 1.26551223 + 
                           t*(1.00002368 + 
                              t*(0.37409196 + 
                                 t*(0.09678418 + 
                                    t*(-0.18628806 + 
                                       t*(0.27886807 +
                                          t*(-1.13520398 + 
                                             t*(1.48851587 +
                                                t*(-0.82215223 +
                                                   t*(0.17087277))))))))));
  return x >= 0.0 ? ans: 2.0 - ans;
}

//------------------------------------------------------------------------------
// Error function.
//------------------------------------------------------------------------------
inline
double
erff(double x) {
  return 1.0 - erffc(x);
}

}

#endif
