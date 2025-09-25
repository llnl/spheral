//---------------------------------Spheral++----------------------------------//
// QuadraticInterpolator
//
// Encapsulates the algorithm and data for parabolic interpolation in 1D
// Assumes the results is interpolated as y_interp = a + b*x + c*x^2
//
// Created by JMO, Fri Dec  4 14:28:08 PST 2020
//----------------------------------------------------------------------------//
#include "QuadraticInterpolator.hh"
#include <algorithm>

#include "CHAI_MA_wrapper.hh"
#include <Eigen/Dense>

namespace Spheral {

//------------------------------------------------------------------------------
// Copy constructor
//------------------------------------------------------------------------------
QuadraticInterpolator::QuadraticInterpolator(const QuadraticInterpolator& rhs)
  :
  QIView() {
  mN1 = rhs.mN1;
  mXmin = rhs.mXmin;
  mXmax = rhs.mXmax;
  mXstep = rhs.mXstep;
  mVec = rhs.mVec;
  initMA();
}

//------------------------------------------------------------------------------
// Assignment constructor
//------------------------------------------------------------------------------
QuadraticInterpolator&
QuadraticInterpolator::operator=(const QuadraticInterpolator& rhs) {
  if (this != &rhs) {
    mN1 = rhs.mN1;
    mXmin = rhs.mXmin;
    mXmax = rhs.mXmax;
    mXstep = rhs.mXstep;
    mVec = rhs.mVec;
    initMA();
  }
  return *this;
}

//------------------------------------------------------------------------------
// Constructor with sampled values
//------------------------------------------------------------------------------
QuadraticInterpolator::QuadraticInterpolator(double xmin,
                                             double xmax,
                                             const std::vector<double>& yvals) {
  initialize(xmin, xmax, yvals);
}

//------------------------------------------------------------------------------
// Initialize the interpolation to fit the given data
//------------------------------------------------------------------------------
void
QuadraticInterpolator::initialize(double xmin,
                                  double xmax,
                                  const std::vector<double>& yvals) {
  const auto n = yvals.size();
  VERIFY2(n > 2, "QuadraticInterpolator::initialize requires at least 3 unique values to fit");
  VERIFY2(n % 2 == 1, "QuadraticInterpolator::initialize requires an odd number of tabulated values");
  VERIFY2(xmax > xmin, "QuadraticInterpolator::initialize requires a positive domain: [" << xmin << " " << xmax << "]");

  size_t N1 = (n - 1u)/2u - 1u;  // Maximum index into arrays
  double xstep = (xmax - xmin)/(N1 + 1u);
  size_t N = 3*(N1 + 1u);
  mVec.resize(N);

  typedef Eigen::Matrix<double, 3, 3, Eigen::RowMajor> EMatrix;
  typedef Eigen::Matrix<double, 3, 1> EVector;

  // We use simple least squares fitting for each 3-point interval, giving us an
  // exact parabolic fit for those three points.
  // Find the coefficient fits.
  double x0, x1, x2;
  EMatrix A;
  EVector X, B;
  for (auto i0 = 0u; i0 <= N1; ++i0) {
    x0 = xmin + i0*xstep;
    x1 = x0 + 0.5*xstep;
    x2 = x0 + xstep;
    A << 1.0, x0, x0*x0,
         1.0, x1, x1*x1,
         1.0, x2, x2*x2;
    B << yvals[2u*i0], yvals[2u*i0 + 1u], yvals[2u*i0 + 2u];
    X = A.inverse()*B;
    mVec[3*i0     ] = X(0);
    mVec[3*i0 + 1u] = X(1);
    mVec[3*i0 + 2u] = X(2);
  }
  mN1 = N1;
  mXmin = xmin;
  mXmax = xmax;
  mXstep = xstep;
  initMA();
}

void
QuadraticInterpolator::initMA() {
  initializeManagedArray(mcoeffs, mVec);
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
QuadraticInterpolator::~QuadraticInterpolator() {
  mcoeffs.free();
}

}
