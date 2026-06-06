//---------------------------------Spheral++----------------------------------//
// LinearInterpolator
//
// Uses linear regression to compute the best fit line (minimizing vertical
// discrepancy in y values).
//
// Computes linear fit in the form y = A*x + B
//
// Created by JMO, Thu Mar  7 11:28:19 PST 2024
//----------------------------------------------------------------------------//
#include "LinearInterpolator.hh"
#include "Utilities/safeInv.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Default constructor
//------------------------------------------------------------------------------
LinearInterpolator::LinearInterpolator():
  mA(0.0),
  mB(0.0) {
}

//------------------------------------------------------------------------------
// Constructor with sampled values
//------------------------------------------------------------------------------
LinearInterpolator::LinearInterpolator(const std::vector<double>& xvals,
                                       const std::vector<double>& yvals):
  mA(0.0),
  mB(0.0) {
  this->initialize(xvals, yvals);
}

//------------------------------------------------------------------------------
// Initialize the interpolation to fit the given data
//------------------------------------------------------------------------------
void
LinearInterpolator::initialize(const std::vector<double>& xvals,
                               const std::vector<double>& yvals) {

  const auto n = xvals.size();
  VERIFY2(n > 1, "LinearInterpolator::initialize requires at least 3 unique values to fit");
  VERIFY2(yvals.size() == xvals.size(), "LinearInterpolator::initialize requires number of x and y values to be the same");

  // Solve the fit
  double xsum = 0.0;
  double x2sum = 0.0;
  double ysum = 0.0;
  double xysum = 0.0;
  for (auto i = 0u; i < n; ++i) {
    xsum += xvals[i];
    x2sum += xvals[i]*xvals[i];
    ysum += yvals[i];
    xysum += xvals[i]*yvals[i];
  }

  mA = (n*xysum - xsum*ysum)*safeInv(n*x2sum - xsum*xsum);
  mB = (ysum*x2sum - xsum*xysum)*safeInv(n*x2sum - xsum*xsum);
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
LinearInterpolator::~LinearInterpolator() {
}

//------------------------------------------------------------------------------
// Equivalence
//------------------------------------------------------------------------------
bool
LinearInterpolator::
operator==(const LinearInterpolator& rhs) const {
  return ((mA == rhs.mA) and
          (mB == rhs.mB));
}

}
