#include "Utilities/DBC.hh"
#include <algorithm>

#include <Eigen/Dense>

namespace Spheral {

//------------------------------------------------------------------------------
// Construct to fit the given function
//------------------------------------------------------------------------------
template<typename Func>
inline
LinearInterpolator::LinearInterpolator(double xmin,
                                       double xmax,
                                       size_t n,
                                       const Func& F):
  mA(0.0),
  mB(0.0) {
  this->initialize(xmin, xmax, n, F);
}

//------------------------------------------------------------------------------
// Initialize to fit the given function
//------------------------------------------------------------------------------
template<typename Func>
inline
void
LinearInterpolator::initialize(double xmin,
                               double xmax,
                               size_t n,
                               const Func& F) {
  // Preconditions
  VERIFY2(n > 1, "LinearInterpolator requires n > 1 : n=" << n);
  VERIFY2(xmax > xmin, "LinearInterpolator requires a positive domain: [" << xmin << " " << xmax << "]");

  // Build up an array of the function values and use the array based initialization.
  if (n % 2 == 0) ++n;  // Need odd number of samples to hit both endpoints of the range
  const auto dx = (xmax - xmin)/(n - 1u);
  std::vector<double> xvals(n), yvals(n);
  for (auto i = 0u; i < n; ++i) {
    xvals[i] = xmin + i*dx;
    yvals[i] = F(xvals[i]);
  }
  this->initialize(xvals, yvals);
}

//------------------------------------------------------------------------------
// Interpolate for the given x value.
//------------------------------------------------------------------------------
inline
double
LinearInterpolator::operator()(const double x) const {
  return mA*x + mB;
}

//------------------------------------------------------------------------------
// Data accessors
//------------------------------------------------------------------------------
inline
double
LinearInterpolator::slope() const {
  return mA;
}

inline
double
LinearInterpolator::yintercept() const {
  return mB;
}

}
