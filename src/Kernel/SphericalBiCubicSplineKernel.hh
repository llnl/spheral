//---------------------------------Spheral++----------------------------------//
// SphericalBiCubicSplineKernel
//
// Take a 3D Kernel and build a specialized 1D tabulated version appropriate
// for use with the spherical SPH algorithm described in
// Omang, M., Børve, S., & Trulsen, J. (2006). SPH in spherical and cylindrical coordinates.
// Journal of Computational Physics, 213(1), 391412. https://doi.org/10.1016/j.jcp.2005.08.023
//
// Created by JMO, Wed Dec  2 16:41:20 PST 2020
//----------------------------------------------------------------------------//
#ifndef __Spheral_SphericalBiCubicSplineKernel_hh__
#define __Spheral_SphericalBiCubicSplineKernel_hh__

#include "TableKernel.hh"
#include "Geometry/Dimension.hh"

namespace Spheral {

class SphericalBiCubicSplineKernel {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = Dim<1>::Scalar;
  using Vector = Dim<1>::Vector;
  using Tensor = Dim<1>::Tensor;
  using SymTensor = Dim<1>::SymTensor;

  // Constructor.
  // Takes a normal 3D TableKernel and constructs the integral form appropriate
  // for 1D spherical coordinates.
  SphericalBiCubicSplineKernel(const unsigned numKernel = 200u);
  SphericalBiCubicSplineKernel(const SphericalBiCubicSplineKernel& rhs);

  // Assignment.
  SphericalBiCubicSplineKernel& operator=(const SphericalBiCubicSplineKernel& rhs);

  // Comparisons
  bool operator==(const SphericalBiCubicSplineKernel& rhs) const;

  // These methods taking a Vector eta and Vector position are the special methods
  // allowing this kernel to implement the asymmetric sampling as a function of r.
  // Arguments:
  //  etaj : Vector normalized coordinate: etaj = H*posj
  //  etai : Vector normalized coordinate: etai = H*posi
  //  Hdet  : Determinant of the H tensor used to compute eta
  double operator()(const Vector& etaj, const Vector& etai, const Scalar Hdet) const;
  Vector grad(const Vector& etaj, const Vector& etai, const SymTensor& H) const;
  void kernelAndGrad(const Vector& etaj, const Vector& etai, const SymTensor& H,
                     Scalar& W,
                     Vector& gradW,
                     Scalar& deltaWsum) const;

  // Access our internal data.
  const TableKernel<Dim<3>>& baseKernel3d() const;
  const TableKernel<Dim<1>>& baseKernel1d() const;
  Scalar etamax() const;

private:
  //--------------------------- Private Interface ---------------------------//
  // Data for the kernel tabulation.
  TableKernel<Dim<3>> mBaseKernel3d;
  TableKernel<Dim<1>> mBaseKernel1d;  // Only for use with the IdealH algorithm
  double metamax;
};

}

#include "SphericalBiCubicSplineKernelInline.hh"

#endif

