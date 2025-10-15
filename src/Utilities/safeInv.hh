//------------------------------------------------------------------------------
// safeInv
//
// Take the inverse of a number, turning it to zero if the argument is 0.0.
//------------------------------------------------------------------------------
#ifndef __Spheral_safevInv__
#define __Spheral_safevInv__

#include "SpheralFunctions.hh"
#include "config.hh"

namespace Spheral {

template<typename Value>
inline
Value
safeInv(const Value& x,
	const double fuzz = 1.0e-30) {
  return x/(x*x + fuzz);
}

template<typename Value>
inline
Value
safeInvVar(const Value& x,
           const double fuzz = 1.0e-30) {
  return sgn(x)/std::max(fuzz, std::abs(x));
}

}

#endif
