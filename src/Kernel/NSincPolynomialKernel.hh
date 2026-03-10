//---------------------------------Spheral++----------------------------------//
// NSincPolynomialKernel -- The sinc interpolation kernel: W = sin(pi*eta)/(pi*eta).
//
// Created by JMO, Tue Jan  7 15:01:13 PST 2003
//----------------------------------------------------------------------------//
#ifndef __Spheral_NSincPolynomialKernel_hh__
#define __Spheral_NSincPolynomialKernel_hh__

#include "Kernel.hh"

#include <vector>

namespace Spheral {

template<typename Dimension>
class NSincPolynomialKernel: 
    public Kernel<Dimension, NSincPolynomialKernel<Dimension> > {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  // Constructor.
  NSincPolynomialKernel(const int order);

  // Return the kernel weight for a given normalized distance or position.
  double kernelValue(double etaij, const double Hdet) const;

  // Return the gradient value for a given normalized distance or position.
  double gradValue(double etaij, const double Hdet) const;

  // Return the second derivative for a given normalized distance or position.
  double grad2Value(double etaij, const double Hdet) const;

private:
  //--------------------------- Private Interface ---------------------------//
  // Order of the polynomials.
  int mOrder;

  // The coefficients for each piecewise section.
  std::vector< std::vector<double> > mAij;

  // Private method to fill in the polynomial coefficients.
  void setPolynomialCoefficients(const int order, 
                                 std::vector< std::vector<double> >& Aij) const;

};

}

#include "NSincPolynomialKernelInline.hh"

#endif
