//---------------------------------Spheral++----------------------------------//
// W4 Spline Kernel -- A extrapolation of the Spline kernel as described in
// Monaghan 1983, J. Comp. Phys, 60, 253
// Equation 25
//
// Volume Normalization:
// 1-D:  A = 1
// 2-D:  A = 30/(7*Pi)
// 3-D:  A = 5/(6*Pi)
//
// Kernel extent: 2.0
//
// Created by JMO, Tue Feb  8 10:59:44 PST 2000
//----------------------------------------------------------------------------//
#ifndef __Spheral_W4SplineKernel_hh__
#define __Spheral_W4SplineKernel_hh__

#include "Kernel.hh"

namespace Spheral {

template<typename Dimension>
class W4SplineKernel: public Kernel<Dimension, W4SplineKernel<Dimension> > {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  // Constructors, destructors.
  W4SplineKernel();
  ~W4SplineKernel();

  // Return the kernel weight for a given normalized distance or position.
  double kernelValue(double etaij, const double Hdet) const;

  // Return the gradient value for a given normalized distance or position.
  double gradValue(double etaij, const double Hdet) const;

  // Return the second derivative for a given normalized distance or position.
  double grad2Value(double etaij, const double Hdet) const;

};

}

#include "W4SplineKernelInline.hh"

#endif
