//---------------------------------Spheral++----------------------------------//
// SuperGaussianKernel -- The Super gaussian interpolation kernel, ala Monaghan.
//
// Volume normalization:
// 1-D:  A = 1/sqrt(Pi)
// 2-D:  A = 1/Pi
// 3-D:  A = 1/Pi^(3/2)
//
// Kernel extent:  3.0
//
// Created by JMO, Wed Dec  1 22:38:21 PST 1999
//----------------------------------------------------------------------------//

#ifndef __Spheral_SuperGaussianKernel_hh__
#define __Spheral_SuperGaussianKernel_hh__

#include "Kernel.hh"
#include "Geometry/Dimension.hh"

namespace Spheral {

template<typename Dimension>
class SuperGaussianKernel: public Kernel<Dimension, SuperGaussianKernel<Dimension> > {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  // Constructor.
  SuperGaussianKernel();

  // Return the kernel weight for a given normalized distance or position.
  double kernelValue(double etaij, const double Hdet) const;

  // Return the gradient value for a given normalized distance or position.
  double gradValue(double etaij, const double Hdet) const;

  // Return the gradient value for a given normalized distance or position.
  double grad2Value(double etaij, const double Hdet) const;

private:
  //--------------------------- Private Interface ---------------------------//
  static constexpr double mKW = 0.5*(Dimension::nDim + 2.0);
  static constexpr double mKGW = 0.5*Dimension::nDim;

};

}

#include "SuperGaussianKernelInline.hh"

#endif
