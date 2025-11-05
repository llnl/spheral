//---------------------------------Spheral++----------------------------------//
// SpheralFunctions -- a collection of routines that don't fall in any specific
//   category.
//
// Created by JMO, Sat Jan 22 14:22:33 PST 2000
//----------------------------------------------------------------------------//

#ifndef SpheralFunctions_HH
#define SpheralFunctions_HH

#include "config.hh"

#include <cmath>
#include <algorithm>
#include <limits>

namespace Spheral {

//------------------------------------------------------------------------------
// Fuzzy comparisons.
//------------------------------------------------------------------------------
template<typename DataType>
SPHERAL_HOST_DEVICE
bool fuzzyEqual(const DataType& lhs, const DataType& rhs,
                const double fuzz = 1.0e-15) {
  // return std::abs(lhs - rhs)/std::max(std::numeric_limits<DataType>::epsilon(), std::max(fuzz, std::max(std::abs(lhs), std::abs(rhs)))) < fuzz;
  // const auto absA = std::abs(lhs);
  // const auto absB = std::abs(rhs);
  // const auto diff = std::abs(lhs - rhs);
  // if (lhs == rhs) { // shortcut, handles infinities
  //   return true;
  // } else if (lhs == 0 or rhs == 0 or diff < std::numeric_limits<DataType>::denorm_min()) {
  //   // a or b is zero or both are extremely close to it
  //   // relative error is less meaningful here
  //   return diff < (fuzz * std::numeric_limits<DataType>::denorm_min());
  // } else { // use relative error
  //   return diff / std::min((absA + absB), std::numeric_limits<DataType>::max()) < fuzz;
  // }
  return std::abs(lhs - rhs) <= fuzz * std::max(1.0, std::abs(lhs) + std::abs(rhs));
}

template<typename DataType>
SPHERAL_HOST_DEVICE
bool fuzzyLessThanOrEqual(const DataType& lhs, const DataType& rhs,
                          const double fuzz = 1.0e-15) {
  return lhs < rhs || fuzzyEqual(lhs, rhs, fuzz);
}

template<typename DataType>
SPHERAL_HOST_DEVICE
bool fuzzyGreaterThanOrEqual(const DataType& lhs, const DataType& rhs,
                             const double fuzz = 1.0e-15) {
  return lhs > rhs || fuzzyEqual(lhs, rhs, fuzz);
}

template<typename DataType>
bool distinctlyLessThan(const DataType& lhs, const DataType& rhs,
                        const double fuzz = 1.0e-15) {
  return lhs < rhs && !fuzzyEqual(lhs, rhs, fuzz);
}

template<typename DataType>
bool distinctlyGreaterThan(const DataType& lhs, const DataType& rhs,
                           const double fuzz = 1.0e-15) {
  return lhs > rhs && !fuzzyEqual(lhs, rhs, fuzz);
}

//------------------------------------------------------------------------------
// Return the sign of the argument determined as follows:
//   
//    x >= 0 -> sgn(x) = 1
//    x <  0 -> sgn(x) = -1
//------------------------------------------------------------------------------
SPHERAL_HOST_DEVICE
inline
double
sgn(const double x) {
  return x < 0.0 ? -1.0 : 1.0;
}

inline
int
isgn(const double x) {
  return x < 0.0 ? -1 : 1;
}

inline
int
sgn(const int x) {
  return x < 0 ? -1 : 1;
}

//------------------------------------------------------------------------------
// A varient that allows 0 as an answer:
//   
//    x > 0 -> sgn0(x) =  1
//    x = 0 -> sgn0(x) =  0
//    x < 0 -> sgn0(x) = -1
//------------------------------------------------------------------------------
inline
double
sgn0(const double x) {
  return (x > 0.0 ?  1.0 :
          x < 0.0 ? -1.0 :
          0.0);
}

inline
int
isgn0(const double x) {
  return (x > 0.0 ?  1 :
          x < 0.0 ? -1 :
          0);
}

inline
int
sgn0(const int x) {
  return (x > 0 ?  1 :
          x < 0 ? -1 :
          0);
}

//------------------------------------------------------------------------------
// Generalized versions of min & max that work with Spheral data types.
//------------------------------------------------------------------------------
// template<typename T1>
// T1
// min(const T1& lhs, const T1& rhs) {
//   return std::min(lhs, rhs);
// }

// template<typename T1>
// T1
// max(const T1& lhs, const T1& rhs) {
//   return std::max(lhs, rhs);
// }

}

#endif
