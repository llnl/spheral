#include "Utilities/DBC.hh"

#include <Eigen/Dense>

namespace Spheral {

//------------------------------------------------------------------------------
// Construct to fit the given function -- we have to estimate the gradient at
// each point.
//------------------------------------------------------------------------------
template<typename Func>
inline
CubicHermiteInterpolator::CubicHermiteInterpolator(const double xmin,
                                                   const double xmax,
                                                   const size_t n,
                                                   const Func& F) {
  initialize(xmin, xmax, n, F);
}

//------------------------------------------------------------------------------
// Construct to fit the given function with it's gradient
//------------------------------------------------------------------------------
template<typename Func, typename GradFunc>
inline
CubicHermiteInterpolator::CubicHermiteInterpolator(const double xmin,
                                                   const double xmax,
                                                   const size_t n,
                                                   const Func& F,
                                                   const GradFunc& Fgrad) {
  initialize(xmin, xmax, n, F, Fgrad);
}

//------------------------------------------------------------------------------
// (Re)initialize from a function
//------------------------------------------------------------------------------
template<typename Func>
inline
void
CubicHermiteInterpolator::initialize(const double xmin,
                                     const double xmax,
                                     const size_t n,
                                     const Func& F) {
  double xstep = (xmax - xmin)/(n - 1u);
  std::vector<double> yvals(n);
  for (auto i = 0u; i < n; ++i) yvals[i] = F(xmin + i*xstep);
  initialize(xmin, xmax, yvals);
}

//------------------------------------------------------------------------------
// (Re)initialize from a function and its gradient
//------------------------------------------------------------------------------
template<typename Func, typename GradFunc>
inline
void
CubicHermiteInterpolator::initialize(const double xmin,
                                     const double xmax,
                                     const size_t n,
                                     const Func& F,
                                     const GradFunc& Fgrad) {
  // Preconditions
  VERIFY2(n > 1u, "CubicHermiteInterpolator requires n >= 2 : n=" << n);
  VERIFY2(xmax > xmin, "CubicHermiteInterpolator requires a positive domain: [" << xmin << " " << xmax << "]");

  mN = n;
  mXmin = xmin;
  mXmax = xmax;
  mXstep = (xmax - xmin)/(n - 1u);
  mVec.resize(2u*n);

  // Compute the function and gradient values
  for (auto i = 0u; i < mN; ++i) {
    const auto xi = xmin + i*mXstep;
    mVec[i] = F(xi);
    mVec[mN + i] = Fgrad(xi);
  }
  initView();
}

}
