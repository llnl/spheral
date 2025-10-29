#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Interpolate for the given x value.
//------------------------------------------------------------------------------
SPHERAL_HOST_DEVICE inline
double
CHIView::operator()(const double x) const {
  if (x < mXmin) {
    return mVals[0] + mVals[mN]*(x - mXmin);
  } else if (x > mXmax) {
    return mVals[mN-1u] + mVals[2u*mN-1u]*(x - mXmin);
  } else {
    const auto i0 = lowerBound(x);
    return this->operator()(x, i0);
  }
}

SPHERAL_HOST_DEVICE inline
double
CHIView::operator()(const double x,
                    const size_t i0) const {
  REQUIRE(i0 <= mN - 2u);
  const auto t = std::max(0.0, std::min(1.0, (x - mXmin - i0*mXstep)/mXstep));
  const auto t2 = t*t;
  const auto t3 = t*t2;
  return ((2.0*t3 - 3.0*t2 + 1.0)*mVals[i0] +          // h00
          (-2.0*t3 + 3.0*t2)*mVals[i0 + 1u] +          // h01
          mXstep*((t3 - 2.0*t2 + t)*mVals[mN + i0] +   // h10
                  (t3 - t2)*mVals[mN + i0 + 1u]));     // h11
}

SPHERAL_HOST_DEVICE inline
double
CHIView::operator[](const size_t i) const {
  REQUIRE(size() > 0);
  REQUIRE(i < size());
  return mVals[i];
}

//------------------------------------------------------------------------------
// Interpolate for dy/dx
//------------------------------------------------------------------------------
SPHERAL_HOST_DEVICE inline
double
CHIView::prime(const double x) const {
  if (x < mXmin) {
    return mVals[mN];
  } else if (x > mXmax) {
    return mVals[2u*mN-1u];
  } else {
    const auto i0 = lowerBound(x);
    return this->prime(x, i0);
  }
}

SPHERAL_HOST_DEVICE inline
double
CHIView::prime(const double x,
               const size_t i0) const {
  REQUIRE(i0 <= mN - 2u);
  const auto t = std::max(0.0, std::min(1.0, (x - mXmin - i0*mXstep)/mXstep));
  const auto t2 = t*t;
  return (6.0*(t2 - t)*(mVals[i0] - mVals[i0 + 1u])/mXstep +
          (3.0*t2 - 4.0*t + 1.0)*mVals[mN + i0] +
          (3.0*t2 - 2.0*t)*mVals[mN + i0 + 1u]);
}

//------------------------------------------------------------------------------
// Interpolate for d^2y/dx^2
//------------------------------------------------------------------------------
SPHERAL_HOST_DEVICE inline
double
CHIView::prime2(const double x) const {
  if (x < mXmin or x > mXmax) {
    return 0.0;
  } else {
    const auto i0 = lowerBound(x);
    return prime2(x, i0);
  }
}

SPHERAL_HOST_DEVICE inline
double
CHIView::prime2(const double x,
                const size_t i0) const {
  REQUIRE(i0 <= mN - 2u);
  const auto t = std::max(0.0, std::min(1.0, (x - mXmin - i0*mXstep)/mXstep));
  return 2.0*(3.0*(2.0*t - 1.0)*(mVals[i0] - mVals[i0 + 1u])/mXstep +
              (3.0*t - 2.0)*mVals[mN + i0] +
              (3.0*t - 1.0)*mVals[mN + i0 + 1u])/mXstep;
}

//------------------------------------------------------------------------------
// Return the lower bound entry in the table for the given x coordinate
//------------------------------------------------------------------------------
SPHERAL_HOST_DEVICE inline
size_t
CHIView::lowerBound(const double x) const {
  const auto result = std::min(mN - 2u, size_t(std::max(0.0, x - mXmin)/mXstep));
  ENSURE(result <= mN - 2u);
  return result;
}

}
