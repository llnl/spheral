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
#ifndef __Spheral_LinearInterpolator__
#define __Spheral_LinearInterpolator__

#include <cstddef>
#include <vector>

namespace Spheral {

class LinearInterpolator {
public:
  //--------------------------- Public Interface ---------------------------//
  // Constructors, destructors
  template<typename Func>
  LinearInterpolator(double xmin, double xmax, size_t n, const Func& F);
  LinearInterpolator(const std::vector<double>& xvals, const std::vector<double>& yvals);
  LinearInterpolator();
  ~LinearInterpolator();

  // Initialize after construction, either with a function or tabulated values
  template<typename Func>
  void initialize(double xmin, double xmax, size_t n, const Func& f);
  void initialize(const std::vector<double>& xvals, const std::vector<double>& yvals);

  // Comparisons
  bool operator==(const LinearInterpolator& rhs) const;

  // Interpolate for the y value
  double operator()(const double x) const;

  // Allow read access the internal data representation
  double slope() const;                  // fitted slope
  double yintercept() const;             // fitted y-intercept
  
private:
  //--------------------------- Private Interface --------------------------//
  // Member data
  double mA, mB;
};

}

#include "LinearInterpolatorInline.hh"

#endif
