//---------------------------------Spheral++----------------------------------//
// NBSplineKernel -- The nth order B Spline kernel.
//
// Schoenberg 1969, Journal of Approximation Theory, 2, 167-206.
//
// Created by JMO, Mon Jan 13 22:24:23 PST 2003
//----------------------------------------------------------------------------//

#ifndef __Spheral_NBSplineKernel_hh__
#define __Spheral_NBSplineKernel_hh__

#include "Kernel.hh"

namespace Spheral {

template<typename Dimension>
class NBSplineKernel: 
    public Kernel<Dimension, NBSplineKernel<Dimension> > {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  // Constructor.
  NBSplineKernel(const int order);

  // Return the kernel weight for a given normalized distance or position.
  double kernelValue(double etaij, const double Hdet) const;

  // Return the gradient value for a given normalized distance or position.
  double gradValue(double etaij, const double Hdet) const;

  // Return the second derivative for a given normalized distance or position.
  double grad2Value(double etaij, const double Hdet) const;

  // Access the order of this spline.
  int order() const;
  void setOrder(const int order);

  // Calculate the factorial of the given integer.
  int factorial(const int n) const;

  // Calculate the binomial coefficient for the given integers: (n)
  //                                                            (m)
  int binomialCoefficient(const int n, const int m) const;

  // The one sided power function:  opf(x) = x if x >= 0, 0 otherwise.
  double oneSidedPowerFunction(const double x,
                               const int exponent) const;

private:
  //--------------------------- Private Interface ---------------------------//
  // Order of the kernel.
  int mOrder;

  // Method to handle setting the kernel extent and volume normalization 
  // appropriately for the order of the kernel.
  void initializeKernel();
};

}

#include "NBSplineKernelInline.hh"

#endif
