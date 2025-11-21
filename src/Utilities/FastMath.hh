//---------------------------------Spheral++----------------------------------//
// A collection of small numerical routines to speed up common tasks (like cube
// roots or squares).
//
// Some of these methods are open-source things I found online, and include the
// original authors comments and such.
// Note I have modified the CubeRoot method from it's original in that I made
// it sign safe, and handled incorrect answers for zero.
//------------------------------------------------------------------------------
#ifndef __Spheral_FastMath__
#define __Spheral_FastMath__

#include "Utilities/SpheralFunctions.hh"
#include "Utilities/DBC.hh"

namespace Spheral {
namespace FastMath {

//------------------------------------------------------------------------------
// A helpful local class for returning a small but non-zero number for each
// data type we support here (currently just floats and doubles).
//------------------------------------------------------------------------------
template<typename D> D tinyValue();
template<> inline float tinyValue<float>() { return 1.0e-35f; }
template<> inline double tinyValue<double>() { return 1.0e-100; }

//------------------------------------------------------------------------------
// Taking the nth root of a number.
// Note that these methods rely on IEEE 754 floats and doubles, in which case
// the author claims roughly 5-6 digit accuracies for these methods.
// Based on http://metamerist.com/cbrt/cbrt.htm
//------------------------------------------------------------------------------
template<int n>
SPHERAL_HOST_DEVICE
inline 
float 
nth_root(float x) {
   const int ebits = 8;
   const int fbits = 23;

   int& i = (int&) x;
   const int bias = (1 << (ebits-1))-1;
   i = (i - (bias << fbits)) / n + (bias << fbits);

   return x;
}

template<int n>
SPHERAL_HOST_DEVICE
inline 
double
nth_root(double x) {
   const int ebits = 11;
   const int fbits = 52;

//    _int64& i = (_int64&) x;
//    const _int64 bias = (1 << (ebits-1))-1;
   long long& i = (long long&) x;
   const long long bias = (1 << (ebits-1))-1;
   i = (i - (bias << fbits)) / n + (bias << fbits);

   return x;
}

//------------------------------------------------------------------------------
// Quake 3 method for computing a fast 1/sqrt.
// Assumes IEEE floats!
//------------------------------------------------------------------------------
//inline
//float Quake3InvSqrtf(float x) {
//  const float xhalf = 0.5f * x;
//  int i = *(int*)(&x);
//  i = 0x5f3759df - (i >> 1);
//  x = *(float*)(&i);
//  x = x*(1.5f - xhalf*x*x);
//  x = x*(1.5f - xhalf*x*x);
//  return x;
//}
//
//inline
//float Quake3Sqrtf(float x) {
//  return x*Quake3InvSqrtf(x);
//}

//------------------------------------------------------------------------------
// A iterative sqrt estimate.  This appears to start with a very similar trick
// to the Quake 3 idea.
//------------------------------------------------------------------------------
inline
SPHERAL_HOST_DEVICE
double Sqrt(const double y) {
  double x, z, tempf;
  tempf = y;
  unsigned long *tfptr = ((unsigned long *)&tempf) + 1;
  *tfptr = (0xbfcdd90a - *tfptr)>>1; // estimate of 1/sqrt(y)
  x =  tempf;
  z =  y*0.5;                        // hoist out the ~/2~
  x = (1.5*x) - (x*x)*(x*z);         // iteration formula
  x = (1.5*x) - (x*x)*(x*z);
  x = (1.5*x) - (x*x)*(x*z);
  x = (1.5*x) - (x*x)*(x*z);
  x = (1.5*x) - (x*x)*(x*z);
  return x*y;
}

//------------------------------------------------------------------------------
// Cube root approximation using 2 iterations of Halley's method (provided for
// floats and doubles).
// Based on William Kahan's cube root approximation.
// Taken from http://metamerist.com/cbrt/cbrt.htm
//------------------------------------------------------------------------------
// // cube root approximation using bit hack for 32-bit float
// inline
// float
// cbrt_5(float f) {
//   unsigned int* p = (unsigned int *) &f;
//   *p = *p/3 + 709921077;
//   return f;
// }

// // cube root approximation using bit hack for 64-bit float
// // adapted from Kahan's cbrt
// inline
// double
// cbrt_5(double d) {
//   const unsigned int B1 = 715094163;
//   double t = 0.0;
//   unsigned int* pt = (unsigned int*) &t;
//   unsigned int* px = (unsigned int*) &d;
//   pt[1]=px[1]/3+B1;
//   return t;
// }

// iterative cube root approximation using Halley's method
template<typename D>
SPHERAL_HOST_DEVICE
inline
D
cbrta_halley(const D a, const D R) {
  CONTRACT_VAR(R);
#ifdef __INTEL_COMPILER
  const D a3 = a*a*a;
  const D b = a * (a3 + R + R) / (a3 + a3 + R + tinyValue<D>());
  return b;
#else
  return Spheral::sgn(a)*std::pow(std::abs(a), 0.3333333333333333);
#endif
}

// Cube root approximation using 2 iterations of Halley's method
template<typename D>
SPHERAL_HOST_DEVICE
inline
D
CubeRootHalley2(D d) {
#ifdef __INTEL_COMPILER
  const D sgnd = sgn(d);
  d = std::abs(d);
  return sgnd*cbrta_halley(cbrta_halley(nth_root<3>(d), d), d);
#else
  return Spheral::sgn(d)*std::pow(std::abs(d), 0.3333333333333333);
#endif
}

//------------------------------------------------------------------------------
// Square root approximation using 2 iterations of Halley's method (float)
//------------------------------------------------------------------------------
// Single pass Halley iteration for sqrt.
template<typename D>
SPHERAL_HOST_DEVICE
inline
D
sqrta_halley(const D a, const D R) {
  CONTRACT_VAR(R);
#ifdef __INTEL_COMPILER
  const D a2 = a*a;
  return a*(a2 + R + R + R)/(a2 + a2 + a2 + R + tinyValue<D>());
#else
  return sqrt(a);
#endif
}

// Sqrt approximation using 2 iterations of Halley's method.
template<typename D>
SPHERAL_HOST_DEVICE
inline
D
SqrtHalley2(const D x) {
#ifdef __INTEL_COMPILER
  return sqrta_halley(sqrta_halley(nth_root<2>(x), x), x);
#else
  return sqrt(x);
#endif
}

//------------------------------------------------------------------------------
// Square
//------------------------------------------------------------------------------
template<typename T>
SPHERAL_HOST_DEVICE
inline
T
square(const T& x) {
  return x*x;
}

//------------------------------------------------------------------------------
// Cube
//------------------------------------------------------------------------------
template<typename T>
SPHERAL_HOST_DEVICE
inline
T
cube(const T& x) {
  return x*x*x;
}

//------------------------------------------------------------------------------
// powN versions of above.
//------------------------------------------------------------------------------
template<typename T> SPHERAL_HOST_DEVICE inline T pow2(const T& x) { return x*x; }
template<typename T> SPHERAL_HOST_DEVICE inline T pow3(const T& x) { return x*x*x; }
template<typename T> SPHERAL_HOST_DEVICE inline T pow4(const T& x) { return x*x*x*x; }
template<typename T> SPHERAL_HOST_DEVICE inline T pow5(const T& x) { return x*x*x*x*x; }
template<typename T> SPHERAL_HOST_DEVICE inline T pow6(const T& x) { return x*x*x*x*x*x; }
template<typename T> SPHERAL_HOST_DEVICE inline T pow7(const T& x) { return x*x*x*x*x*x*x; }
template<typename T> SPHERAL_HOST_DEVICE inline T pow8(const T& x) { return x*x*x*x*x*x*x*x; }
template<typename T> SPHERAL_HOST_DEVICE inline T pow9(const T& x) { return x*x*x*x*x*x*x*x*x; }

//-----------------------------------------------------------------------------
// Compile time power function
//-----------------------------------------------------------------------------
template<typename T>
SPHERAL_HOST_DEVICE 
constexpr T calcPower(T value, unsigned power) {
  return power == 0 ? 1 : value * calcPower(value, power - 1);
}

// template <int Base, int Exp>
// struct calcPower {
//     static const int value = Base * calcPower<Base, Exp - 1>::value;
// };

// template <int Base>
// struct calcPower<Base, 0u> {
//     static const int value = 1;
// };

//------------------------------------------------------------------------------
// The below CubeRoot method is not recommended for now -- just keeping here in
// case we want to try it again sometime.
//
// CubeRoot
/* Copyright (C) 1997-1999 Ken Turkowski. <turk@computer.org>
 *
 * All rights reserved.
 *
 * Warranty Information
 *  Even though I have reviewed this software, I make no warranty
 *  or representation, either express or implied, with respect to this
 *  software, its quality, accuracy, merchantability, or fitness for a
 *  particular purpose.  As a result, this software is provided "as is,"
 *  and you, its user, are assuming the entire risk as to its quality
 *  and accuracy.
 *
 * This code may be used and freely distributed as long as it includes
 * this copyright notice and the above warranty information.
 */

/*
 * Written by Ken Turkowski.
 */
//------------------------------------------------------------------------------

// inline
// float
// CubeRoot(float x)
// {
//   if (std::abs(x) < 1.0e-60) return 0.0;

//   float fr, r;
//   int ex, shx;

//   const float sgnx = sgn(x);
//   x = std::abs(x);

//   /* Argument reduction */
//   fr = frexp(x, &ex);		/* separate into mantissa and exponent */
//   shx = ex % 3;
//   if (shx > 0)
//     shx -= 3;			/* compute shx such that (ex - shx) is divisible by 3 */
//   ex = (ex - shx) / 3;	/* exponent of cube root */
//   fr = ldexp(fr, shx);										/* 0.125 <= fr < 1.0 */

// #ifdef ITERATE
//   /* Compute seed with a quadratic qpproximation */
//   fr = (-0.46946116F * fr + 1.072302F) * fr + 0.3812513F;		/* 0.5   <= fr < 1.0 */
//   r = ldexp(fr, ex);											/* 6 bits of precision */

//   /* Newton-Raphson iterations */
//   r = (float)(2.0/3.0) * r + (float)(1.0/3.0) * x / (r * r);	/* 12 bits of precision */
//   r = (float)(2.0/3.0) * r + (float)(1.0/3.0) * x / (r * r);	/* 24 bits of precision */
// #else
//   /* Use quadric rational polynomial with error < 2^(-24) */
//   fr = ((((45.2548339756803022511987494 * fr + 192.2798368355061050458134625) * fr + 119.1654824285581628956914143) * fr + 13.43250139086239872172837314) * fr + 0.1636161226585754240958355063) /
//     ((((14.80884093219134573786480845 * fr + 151.9714051044435648658557668) * fr + 168.5254414101568283957668343) * fr + 33.9905941350215598754191872) * fr + 1.0);
//   r = ldexp(fr, ex);											/* 24 bits of precision */
// #endif

//   return sgnx*r;
// }

}
}

#endif
