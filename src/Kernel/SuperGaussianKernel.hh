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
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructors, destructors.
  SuperGaussianKernel();
  ~SuperGaussianKernel();

  // Return the kernel weight for a given normalized distance or position.
  double kernelValue(double etaij, const double Hdet) const;

  // Return the gradient value for a given normalized distance or position.
  double gradValue(double etaij, const double Hdet) const;

  // Return the gradient value for a given normalized distance or position.
  double grad2Value(double etaij, const double Hdet) const;

private:
  //--------------------------- Private Interface ---------------------------//
  static const double mKW;
  static const double mKGW;

};

#if defined(SPHERAL1D)
template<> const double SuperGaussianKernel<Dim<1>>::mKW;
template<> const double SuperGaussianKernel<Dim<1>>::mKGW;
#endif

#if defined(SPHERAL2D)
template<> const double SuperGaussianKernel<Dim<2>>::mKW;
template<> const double SuperGaussianKernel<Dim<2>>::mKGW;
#endif

#if defined(SPHERAL3D)
template<> const double SuperGaussianKernel<Dim<3>>::mKW;
template<> const double SuperGaussianKernel<Dim<3>>::mKGW;
#endif

}

#include "SuperGaussianKernelInline.hh"

#endif
