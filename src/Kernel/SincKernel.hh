//---------------------------------Spheral++----------------------------------//
// SincKernel -- The sinc interpolation kernel: W = sin(pi*eta)/(pi*eta).
//
// Created by JMO, Mon Jan  6 22:42:01 PST 2003
//----------------------------------------------------------------------------//
#ifndef __Spheral_SincKernel_hh__
#define __Spheral_SincKernel_hh__

#include "Kernel.hh"

namespace Spheral {

template<typename Dimension>
class SincKernel: public Kernel<Dimension, SincKernel<Dimension> > {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  // Constructors, destructors.
  SincKernel(const double extent);
  ~SincKernel();

  // Return the kernel weight for a given normalized distance or position.
  double kernelValue(double etaij, const double Hdet) const;

  // Return the gradient value for a given normalized distance or position.
  double gradValue(double etaij, const double Hdet) const;

  // Return the second derivative for a given normalized distance or position.
  double grad2Value(double etaij, const double Hdet) const;

};

}

#include "SincKernelInline.hh"

#endif
