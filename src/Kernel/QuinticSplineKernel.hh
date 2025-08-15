//---------------------------------Spheral++----------------------------------//
// QuinticSplineKernel -- A quintic spline, as described in Dehnen and Aly 
// MNRAS 2012.
//
// Kernel extent: 1.0
//
// Created by JMO, Wed Jul  9 16:24:25 PDT 2008
//----------------------------------------------------------------------------//
#ifndef __Spheral_QuinticSplineKernel_hh__
#define __Spheral_QuinticSplineKernel_hh__

#include "Kernel.hh"
#include "Geometry/Dimension.hh"

namespace Spheral {

template<typename Dimension>
class QuinticSplineKernel: public Kernel<Dimension, QuinticSplineKernel<Dimension> > {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  // Constructors, destructors.
  QuinticSplineKernel();
  ~QuinticSplineKernel();

  // Return the kernel weight for a given normalized distance or position.
  double kernelValue(double etaij, const double Hdet) const;

  // Return the gradient value for a given normalized distance or position.
  double gradValue(double etaij, const double Hdet) const;

  // Return the second derivative value for a given normalized distance or
  // position.
  double grad2Value(double etaij, const double Hdet) const;

};

// Forward declare the specialized constructors.
template<> QuinticSplineKernel<Dim<1> >::QuinticSplineKernel();
template<> QuinticSplineKernel<Dim<2> >::QuinticSplineKernel();
template<> QuinticSplineKernel<Dim<3> >::QuinticSplineKernel();

}

#endif
