#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Interpolate for the given x value.
//------------------------------------------------------------------------------
SPHERAL_HOST_DEVICE inline
double
QuadraticInterpolatorView::operator()(const double x) const {
  const auto i0 = lowerBound(x);
  return mcoeffs[i0] + (mcoeffs[i0 + 1] + mcoeffs[i0 + 2]*x)*x;
}

SPHERAL_HOST_DEVICE inline
double
QuadraticInterpolatorView::operator()(const double x,
                                      const size_t i0) const {
  REQUIRE(i0 <= 3u*mN1);
  return mcoeffs[i0] + (mcoeffs[i0 + 1] + mcoeffs[i0 + 2]*x)*x;
}

SPHERAL_HOST_DEVICE inline
double
QuadraticInterpolatorView::operator[](const size_t i) const {
  REQUIRE(size() > 0);
  REQUIRE(i < size());
  return mcoeffs[i];
}

//------------------------------------------------------------------------------
// Interpolate the first derivative the given x value.
//------------------------------------------------------------------------------
SPHERAL_HOST_DEVICE inline
double
QuadraticInterpolatorView::prime(const double x) const {
  const auto i0 = lowerBound(x);
  return mcoeffs[i0 + 1] + 2.0*mcoeffs[i0 + 2]*x;
}

SPHERAL_HOST_DEVICE inline
double
QuadraticInterpolatorView::prime(const double x,
                                 const size_t i0) const {
  REQUIRE(i0 <= 3u*mN1);
  return mcoeffs[i0 + 1] + 2.0*mcoeffs[i0 + 2]*x;
}

//------------------------------------------------------------------------------
// Interpolate the second derivative for the given x value.
// Just a constant value, so not a great fit.
//------------------------------------------------------------------------------
SPHERAL_HOST_DEVICE inline
double
QuadraticInterpolatorView::prime2(const double x) const {
  const auto i0 = lowerBound(x);
  return 2.0*mcoeffs[i0 + 2];
}

SPHERAL_HOST_DEVICE inline
double
QuadraticInterpolatorView::prime2(const double /*x*/,
                                  const size_t i0) const {
  REQUIRE(i0 <= 3u*mN1);
  return 2.0*mcoeffs[i0 + 2];
}

//------------------------------------------------------------------------------
// Return the lower bound entry in the table for the given x coordinate
//------------------------------------------------------------------------------
SPHERAL_HOST_DEVICE inline
size_t
QuadraticInterpolatorView::lowerBound(const double x) const {
  const auto result = 3u*std::min(mN1, size_t(std::max(0.0, x - mXmin)/mXstep));
  ENSURE(result <= 3u*mN1);
  return result;
}

//------------------------------------------------------------------------------
// Equivalence
//------------------------------------------------------------------------------
SPHERAL_HOST_DEVICE inline
bool
QuadraticInterpolatorView::operator==(const QuadraticInterpolatorView& rhs) const {
  return ((mN1 == rhs.mN1) and
          (mXmin == rhs.mXmin) and
          (mXmax == rhs.mXmax) and
          (mcoeffs == rhs.mcoeffs));
}

}
