//---------------------------------Spheral++----------------------------------//
// WendlandC4Kernel -- .
//
//
// Volume normalizations:
// 1-D:  A = 54/37
// 2-D:  A = 9/(5*Pi)
// 3-D:  A = 165/(112*Pi)
//
// Kernel extent: 2.0
//
// Created by CDR, Nov 5 2014
//----------------------------------------------------------------------------//
#ifndef __Spheral_WendlandKernel_hh__
#define __Spheral_WendlandKernel_hh__

#include "Kernel.hh"

namespace Spheral {

template<typename Dimension>
class WendlandC4Kernel: public Kernel<Dimension, WendlandC4Kernel<Dimension> > {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  // Constructors, destructors.
  WendlandC4Kernel();
  ~WendlandC4Kernel();

  // Return the kernel weight for a given normalized distance or position.
  double kernelValue(double etaij, const double Hdet) const;

  // Return the gradient value for a given normalized distance or position.
  double gradValue(double etaij, const double Hdet) const;

  // Return the second derivative value for a given normalized distance or
  // position.
  double grad2Value(double etaij, const double Hdet) const;

};

}

#include "WendlandC4KernelInline.hh"

#endif
