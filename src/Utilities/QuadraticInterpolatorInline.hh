#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Construct to fit the given function
//------------------------------------------------------------------------------
template<typename Func>
QuadraticInterpolator::QuadraticInterpolator(double xmin,
                                             double xmax,
                                             size_t n,
                                             const Func& F) {
  initialize(xmin, xmax, n, F);
}

//------------------------------------------------------------------------------
// Initialize to fit the given function
//------------------------------------------------------------------------------
template<typename Func>
void
QuadraticInterpolator::initialize(double xmin,
                                  double xmax,
                                  size_t n,
                                  const Func& F) {
  // Preconditions
  VERIFY2(n > 1, "QuadraticInterpolator requires n > 1 : n=" << n);
  VERIFY2(xmax > xmin, "QuadraticInterpolator requires a positive domain: [" << xmin << " " << xmax << "]");

  // Build up an array of the function values and use the array based initialization.
  if (n % 2 == 0) ++n;  // Need odd number of samples to hit both endpoints of the range
  double xstep = (xmax - xmin)/(n - 1u);
  std::vector<double> yvals(n);
  for (auto i = 0u; i < n; ++i) yvals[i] = F(xmin + i*xstep);
  initialize(xmin, xmax, yvals);
}

//------------------------------------------------------------------------------
// Interpolate for the given x value.
//------------------------------------------------------------------------------
SPHERAL_HOST_DEVICE inline
double
QIView::operator()(const double x) const {
  const auto i0 = lowerBound(x);
  return mcoeffs[i0] + (mcoeffs[i0 + 1] + mcoeffs[i0 + 2]*x)*x;
}

SPHERAL_HOST_DEVICE inline
double
QIView::operator()(const double x,
                   const size_t i0) const {
  REQUIRE(i0 <= 3u*mN1);
  return mcoeffs[i0] + (mcoeffs[i0 + 1] + mcoeffs[i0 + 2]*x)*x;
}

SPHERAL_HOST_DEVICE inline
double
QIView::operator[](const size_t i) const {
  REQUIRE(size() > 0);
  REQUIRE(i < size());
  return mcoeffs[i];
}

//------------------------------------------------------------------------------
// Interpolate the first derivative the given x value.
//------------------------------------------------------------------------------
SPHERAL_HOST_DEVICE inline
double
QIView::prime(const double x) const {
  const auto i0 = lowerBound(x);
  return mcoeffs[i0 + 1] + 2.0*mcoeffs[i0 + 2]*x;
}

SPHERAL_HOST_DEVICE inline
double
QIView::prime(const double x,
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
QIView::prime2(const double x) const {
  const auto i0 = lowerBound(x);
  return 2.0*mcoeffs[i0 + 2];
}

SPHERAL_HOST_DEVICE inline
double
QIView::prime2(const double /*x*/,
               const size_t i0) const {
  REQUIRE(i0 <= 3u*mN1);
  return 2.0*mcoeffs[i0 + 2];
}

//------------------------------------------------------------------------------
// Return the lower bound entry in the table for the given x coordinate
//------------------------------------------------------------------------------
SPHERAL_HOST_DEVICE inline
size_t
QIView::lowerBound(const double x) const {
  const auto result = 3u*std::min(mN1, size_t(std::max(0.0, x - mXmin)/mXstep));
  ENSURE(result <= 3u*mN1);
  return result;
}

//------------------------------------------------------------------------------
// Data accessors
//------------------------------------------------------------------------------
SPHERAL_HOST_DEVICE inline
size_t
QIView::size() const {
  return 3*(mN1 + 1u);
}

SPHERAL_HOST_DEVICE inline
double
QIView::xmin() const {
  return mXmin;
}

SPHERAL_HOST_DEVICE inline
double
QIView::xmax() const {
  return mXmax;
}

SPHERAL_HOST_DEVICE inline
double
QIView::xstep() const {
  return mXstep;
}

//------------------------------------------------------------------------------
// Equivalence
//------------------------------------------------------------------------------
SPHERAL_HOST_DEVICE
inline
bool
QIView::operator==(const QIView& rhs) const {
  return ((mN1 == rhs.mN1) and
          (mXmin == rhs.mXmin) and
          (mXmax == rhs.mXmax) and
          (mcoeffs == rhs.mcoeffs));
}

}
