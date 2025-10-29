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

}
