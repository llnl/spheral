//---------------------------------Spheral++----------------------------------//
// HatKernel -- The simple linear hat kernel:  W = 2.0 - x.
//
// Volume normalizations:
// 1-D:  A = 1/h
// 2-D:  A = 3/(pi h^2)
// 3-D:  A = 3/(pi h^3)
//
// Kernel extent: 1.0
//
// Created by JMO, Wed Dec 11 17:33:57 PST 2002
//----------------------------------------------------------------------------//
#ifndef __Spheral_HatKernel_hh__
#define __Spheral_HatKernel_hh__

#include "Kernel.hh"

namespace Spheral {

template<typename Dimension>
class HatKernel: public Kernel<Dimension, HatKernel<Dimension> > {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  // Constructors, destructors.
  HatKernel(double eta0, double W0);
  ~HatKernel();

  // Return the kernel weight for a given normalized distance or position.
  double kernelValue(double etaij, const double Hdet) const;

  // Return the gradient value for a given normalized distance or position.
  double gradValue(double etaij, const double Hdet) const;

  // Return the second derivative for a given normalized distance or position.
  double grad2Value(double /*etaij*/, const double /*Hdet*/) const;

  // The x and y kernel intercepts.
  double W0() const;
  double eta0() const;

private:
  //--------------------------- Private Interface ---------------------------//
  // Disallow the default constructor.
  HatKernel();

  // The eta and y intercepts.
  double mEta0;
  double mW0;
  double mSlope;
};

}

#include "HatKernelInline.hh"

#endif
