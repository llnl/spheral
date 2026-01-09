//---------------------------------Spheral++----------------------------------//
// PiGaussianKernel -- My own modified, steeper version of the gaussian
// interpolation kernel, intended for use with the artificial viscosity.
//
// The kernel is of the form W(h) = A*exp(-K*h^4)
//
// Owen, Villumsen, Shapiro, & Martel 1998, ApJS, 116, 155
//
// Created by JMO, Thu Dec  2 16:09:27 PST 1999
//----------------------------------------------------------------------------//

#ifndef __Spheral_PiGaussianKernel_hh__
#define __Spheral_PiGaussianKernel_hh__

#include "Kernel.hh"

namespace Spheral {

template<typename Dimension>
class PiGaussianKernel: public Kernel<Dimension, PiGaussianKernel<Dimension> > {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  // Constructor.
  PiGaussianKernel();
  PiGaussianKernel(double K);

  // Return the kernel weight for a given normalized distance or position.
  double kernelValue(double etaij, const double Hdet) const;

  // Return the gradient value for a given normalized distance or position.
  double gradValue(double etaij, const double Hdet) const;

  // Return the second derivative for a given normalized distance or position.
  double grad2Value(double etaij, const double Hdet) const;

  // Get and set the exponential constant K.
  double getK() const;
  void setK(double K);

private:
  //--------------------------- Private Interface ---------------------------//
  double mK;
  double mKV;

};

}

#include "PiGaussianKernelInline.hh"

#endif
