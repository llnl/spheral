#include "Geometry/Dimension.hh"
#include "VolumeIntegrationFunctions.hh"
#include "Utilities/DBC.hh"

#include <math.h>
#include <vector>

namespace Spheral {

//------------------------------------------------------------------------------
// Construct with the given order.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
NBSplineKernel<Dimension>::NBSplineKernel(const int order):
  Kernel<Dimension, NBSplineKernel<Dimension> >(),
  mOrder(order) {
  REQUIRE(order > 0);

  // Set the kernels extent and normalization appropriately for the specified
  // order.
  initializeKernel();
  ENSURE(this->kernelExtent() > 0.0);
  ENSURE(this->volumeNormalization() > 0.0);
}

//------------------------------------------------------------------------------
// Return the order of this spline.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
NBSplineKernel<Dimension>::order() const {
  return mOrder;
}

//------------------------------------------------------------------------------
// Set the order.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
NBSplineKernel<Dimension>::setOrder(const int order) {
  REQUIRE(order > 0);
  mOrder = order;

  // Set the kernels extent and normalization appropriately for the specified
  // order.
  initializeKernel();
  ENSURE(this->kernelExtent() > 0.0);
  ENSURE(this->volumeNormalization() > 0.0);
}

//------------------------------------------------------------------------------
// Calculate the factorial of the given integer.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
NBSplineKernel<Dimension>::factorial(const int n) const {
  if (n >= 0) {
    int result = 1;
    for (int i = 1; i < n + 1; ++i) result *= i;
    return result;
  } else {
    return INT_MAX;
  }
}

//------------------------------------------------------------------------------
// Calculate the binomial coefficient of the given integers.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
NBSplineKernel<Dimension>::binomialCoefficient(const int n,
                                               const int m) const {
  REQUIRE(n >= 0);
  REQUIRE(m >= 0);
  REQUIRE(n >= m);
  return factorial(n)/(factorial(n - m)*factorial(m));
}

//------------------------------------------------------------------------------
// Calculate the one sided power function.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
NBSplineKernel<Dimension>::oneSidedPowerFunction(const double x,
                                                 const int exponent) const {
  if (x >= 0.0) {
    return pow(x, exponent);
  } else {
    return 0.0;
  }
}

}
