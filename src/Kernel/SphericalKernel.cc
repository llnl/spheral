//---------------------------------Spheral++----------------------------------//
// SphericalKernel
//
// Take a 3D Kernel and build a specialized 1D tabulated version appropriate
// for use with the spherical SPH algorithm described in
// Omang, M., Børve, S., & Trulsen, J. (2006). SPH in spherical and cylindrical coordinates.
// Journal of Computational Physics, 213(1), 391412. https://doi.org/10.1016/j.jcp.2005.08.023
//
// Created by JMO, Wed Dec  2 16:41:20 PST 2020
//----------------------------------------------------------------------------//

#include "SphericalKernel.hh"
#include "Utilities/simpsonsIntegration.hh"
#include "Utilities/DBC.hh"

#include "Kernel/BSplineKernel.hh"
#include "Kernel/NBSplineKernel.hh"
#include "Kernel/W4SplineKernel.hh"
#include "Kernel/GaussianKernel.hh"
#include "Kernel/SuperGaussianKernel.hh"
#include "Kernel/PiGaussianKernel.hh"
#include "Kernel/HatKernel.hh"
#include "Kernel/SincKernel.hh"
#include "Kernel/NSincPolynomialKernel.hh"
#include "Kernel/QuarticSplineKernel.hh"
#include "Kernel/QuinticSplineKernel.hh"
#include "Kernel/WendlandC2Kernel.hh"
#include "Kernel/WendlandC4Kernel.hh"
#include "Kernel/WendlandC6Kernel.hh"
#include "Kernel/ExpInvKernel.hh"

namespace Spheral {

namespace {  // anonymous

//------------------------------------------------------------------------------
// Functor for doing our volume integration of the kernel  
//------------------------------------------------------------------------------
template<typename KernelType>
struct W3S1Func {
  const KernelType& mW;
  double metaMax;
  size_t mn;
  W3S1Func(const KernelType& W, const size_t n): mW(W), metaMax(mW.kernelExtent()), mn(n) {}

  // Define a local nested functor we'll use to do the volume integral for
  // an (etaj, etai) lookup.
  struct VolFunc {
    const KernelType& mW;
    VolFunc(const KernelType& W): mW(W) {}
    double operator()(const double q) const {
      CHECK(q >= 0.0);
      return q * mW(q, 1.0);
    }
  };

  // Now the lookup based on (etaj, etai) -- the (rprime/h, r/h) from the paper
  double operator()(const double a, const double b) const {
    if (a >= b) return 0.0;
    return simpsonsIntegration<VolFunc, double, double>(VolFunc(mW),
                                                        a,
                                                        b,
                                                        mn);
  }
};

}            // anonymous

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename KernelType>
SphericalKernel::SphericalKernel(const KernelType& kernel,
                                 const unsigned numIntegral,
                                 const unsigned numKernel,
                                 const bool useInterpolation):
  mInterp(0.0, kernel.kernelExtent(),
          0.0, kernel.kernelExtent(),
          numKernel, numKernel, 
          W3S1Func<KernelType>(kernel, numIntegral),
          false, false),
  mBaseKernel3d(kernel, numKernel),
  mBaseKernel1d(kernel, numKernel),
  metamax(kernel.kernelExtent()),
  mNumIntegral(numIntegral),
  mUseInterpolation(useInterpolation) {
}

//------------------------------------------------------------------------------
// Copy constructor
//------------------------------------------------------------------------------
SphericalKernel::SphericalKernel(const SphericalKernel& rhs):
  mInterp(rhs.mInterp),
  mBaseKernel3d(rhs.mBaseKernel3d),
  mBaseKernel1d(rhs.mBaseKernel1d),
  metamax(rhs.metamax),
  mNumIntegral(rhs.mNumIntegral),
  mUseInterpolation(rhs.mUseInterpolation) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
SphericalKernel::~SphericalKernel() {
}

//------------------------------------------------------------------------------
// Assignment
//------------------------------------------------------------------------------
SphericalKernel&
SphericalKernel::operator=(const SphericalKernel& rhs) {
  if (this != &rhs) {
    mInterp = rhs.mInterp;
    mBaseKernel3d = rhs.mBaseKernel3d;
    mBaseKernel1d = rhs.mBaseKernel1d;
    metamax = rhs.metamax;
    mNumIntegral = rhs.mNumIntegral;
    mUseInterpolation = rhs.mUseInterpolation;
  }
  return *this;
}

//------------------------------------------------------------------------------
// Lookup the kernel for (rj/h, ri/h) = (etaj, etai)
//------------------------------------------------------------------------------
double
SphericalKernel::operator()(const Dim<1>::Vector& etaj,
                            const Dim<1>::Vector& etai,
                            const Dim<1>::Scalar  Hdet) const {
  REQUIRE(Hdet >= 0.0);
  const auto ei = std::max(1e-5, etai[0]);
  const auto ej = std::max(1e-5, etaj[0]);
  CHECK(ei > 0.0);
  CHECK(ej > 0.0);
  const auto a = std::abs(ej - ei);            // Lower integration limit
  if (a > metamax) return 0.0;
  const auto b = std::min(metamax, ei + ej);   // Upper integration limit
  return 2.0*M_PI/(ei*ej)*FastMath::cube(Hdet)*integralCorrection(a, b, ei, ej);
}

//------------------------------------------------------------------------------
// Lookup the grad kernel for (rj/h, ri/h) = (etaj, etai)
// Using the Leibniz integral rule to differentiate this integral.
//------------------------------------------------------------------------------
Dim<1>::Vector
SphericalKernel::grad(const Dim<1>::Vector& etaj,
                      const Dim<1>::Vector& etai,
                      const Dim<1>::SymTensor& H) const {
  const auto Hdet = H.Determinant();
  REQUIRE(Hdet >= 0.0);
  const auto ei = std::max(1e-5, etai[0]);
  const auto ej = std::max(1e-5, etaj[0]);
  CHECK(ei > 0.0);
  CHECK(ej > 0.0);
  const auto a = std::abs(ej - ei);            // Lower integration limit
  if (a > metamax) return Vector::zero();
  const auto b = std::min(metamax, ei + ej);   // Upper integration limit
  const auto A = a*mBaseKernel3d.kernelValue(a, 1.0)*sgn0(ei - ej);
  const auto B = (ei + ej >= metamax ?
                  0.0 :
                  b*mBaseKernel3d.kernelValue(b, 1.0));
  return Vector(2.0*M_PI/(ei*ej)*FastMath::pow4(Hdet)*(B - A - integralCorrection(a, b, ei, ej)/ei));
}

//------------------------------------------------------------------------------
// Simultaneously lookup (W,  grad W) for (rj/h, ri/h) = (etaj, etai)
//------------------------------------------------------------------------------
void
SphericalKernel::kernelAndGrad(const Dim<1>::Vector& etaj,
                               const Dim<1>::Vector& etai,
                               const Dim<1>::SymTensor& H,
                               Dim<1>::Scalar& W,
                               Dim<1>::Vector& gradW,
                               Dim<1>::Scalar& deltaWsum) const {
  const auto Hdet = H.Determinant();
  REQUIRE(Hdet >= 0.0);
  const auto ei = std::max(1e-5, etai[0]);
  const auto ej = std::max(1e-5, etaj[0]);
  CHECK(ei > 0.0);
  CHECK(ej > 0.0);
  const auto a = std::abs(ej - ei);            // Lower integration limit
  if (a > metamax) {
    W = 0.0;
    gradW.Zero();
    deltaWsum = 0.0;
  } else {
    const auto b = std::min(metamax, ei + ej); // Upper integration limit
    const auto B = (ei + ej >= metamax ?
                    0.0 :
                    b*mBaseKernel3d.kernelValue(b, 1.0));
    const auto A = a*mBaseKernel3d.kernelValue(a, 1.0)*sgn0(ei - ej);
    const auto pre = 2.0*M_PI/(ei*ej)*FastMath::cube(Hdet);
    const auto interpVal = integralCorrection(a, b, ei, ej);
    W = pre*interpVal;
    gradW = Vector(pre*Hdet*(B - A - interpVal/ei));
    deltaWsum = mBaseKernel1d.gradValue((etaj - etai).magnitude(), Hdet);
  }
}

//------------------------------------------------------------------------------
// Return the integral correction
//------------------------------------------------------------------------------
inline
double
SphericalKernel::integralCorrection(const double a,
                                    const double b,
                                    const double ei,
                                    const double ej) const {
  return (mUseInterpolation and ei > 0.01 and ej > 0.01 ?
          mInterp(a, b) :
          W3S1Func<TableKernel<Dim<3>>>(mBaseKernel3d, mNumIntegral)(a, b));
}

// Constructor instantiations
template SphericalKernel::SphericalKernel(const TableKernel<Dim<3>>&, const unsigned, const unsigned, const bool);
template SphericalKernel::SphericalKernel(const BSplineKernel<Dim<3>>&, const unsigned, const unsigned, const bool);
template SphericalKernel::SphericalKernel(const NBSplineKernel<Dim<3>>&, const unsigned, const unsigned, const bool);
template SphericalKernel::SphericalKernel(const W4SplineKernel<Dim<3>>&, const unsigned, const unsigned, const bool);
template SphericalKernel::SphericalKernel(const GaussianKernel<Dim<3>>&, const unsigned, const unsigned, const bool);
template SphericalKernel::SphericalKernel(const SuperGaussianKernel<Dim<3>>&, const unsigned, const unsigned, const bool);
template SphericalKernel::SphericalKernel(const PiGaussianKernel<Dim<3>>&, const unsigned, const unsigned, const bool);
template SphericalKernel::SphericalKernel(const HatKernel<Dim<3>>&, const unsigned, const unsigned, const bool);
template SphericalKernel::SphericalKernel(const SincKernel<Dim<3>>&, const unsigned, const unsigned, const bool);
template SphericalKernel::SphericalKernel(const NSincPolynomialKernel<Dim<3>>&, const unsigned, const unsigned, const bool);
template SphericalKernel::SphericalKernel(const QuarticSplineKernel<Dim<3>>&, const unsigned, const unsigned, const bool);
template SphericalKernel::SphericalKernel(const QuinticSplineKernel<Dim<3>>&, const unsigned, const unsigned, const bool);
template SphericalKernel::SphericalKernel(const WendlandC2Kernel<Dim<3>>&, const unsigned, const unsigned, const bool);
template SphericalKernel::SphericalKernel(const WendlandC4Kernel<Dim<3>>&, const unsigned, const unsigned, const bool);
template SphericalKernel::SphericalKernel(const WendlandC6Kernel<Dim<3>>&, const unsigned, const unsigned, const bool);

}

// We need to instantiate the special TableKernel constructors we use
#include "TableKernel.cc"
namespace Spheral {
template TableKernel<Dim<1>>::TableKernel(const TableKernel<Dim<3>>&, const unsigned, const double, const double);
template TableKernel<Dim<1>>::TableKernel(const BSplineKernel<Dim<3>>&, const unsigned, const double, const double);
template TableKernel<Dim<1>>::TableKernel(const NBSplineKernel<Dim<3>>&, const unsigned, const double, const double);
template TableKernel<Dim<1>>::TableKernel(const W4SplineKernel<Dim<3>>&, const unsigned, const double, const double);
template TableKernel<Dim<1>>::TableKernel(const GaussianKernel<Dim<3>>&, const unsigned, const double, const double);
template TableKernel<Dim<1>>::TableKernel(const SuperGaussianKernel<Dim<3>>&, const unsigned, const double, const double);
template TableKernel<Dim<1>>::TableKernel(const PiGaussianKernel<Dim<3>>&, const unsigned, const double, const double);
template TableKernel<Dim<1>>::TableKernel(const HatKernel<Dim<3>>&, const unsigned, const double, const double);
template TableKernel<Dim<1>>::TableKernel(const SincKernel<Dim<3>>&, const unsigned, const double, const double);
template TableKernel<Dim<1>>::TableKernel(const NSincPolynomialKernel<Dim<3>>&, const unsigned, const double, const double);
template TableKernel<Dim<1>>::TableKernel(const QuarticSplineKernel<Dim<3>>&, const unsigned, const double, const double);
template TableKernel<Dim<1>>::TableKernel(const QuinticSplineKernel<Dim<3>>&, const unsigned, const double, const double);
template TableKernel<Dim<1>>::TableKernel(const WendlandC2Kernel<Dim<3>>&, const unsigned, const double, const double);
template TableKernel<Dim<1>>::TableKernel(const WendlandC4Kernel<Dim<3>>&, const unsigned, const double, const double);
template TableKernel<Dim<1>>::TableKernel(const WendlandC6Kernel<Dim<3>>&, const unsigned, const double, const double);

template TableKernel<Dim<3>>::TableKernel(const TableKernel<Dim<3>>&, const unsigned, const double, const double);
template TableKernel<Dim<3>>::TableKernel(const BSplineKernel<Dim<3>>&, const unsigned, const double, const double);
template TableKernel<Dim<3>>::TableKernel(const NBSplineKernel<Dim<3>>&, const unsigned, const double, const double);
template TableKernel<Dim<3>>::TableKernel(const W4SplineKernel<Dim<3>>&, const unsigned, const double, const double);
template TableKernel<Dim<3>>::TableKernel(const GaussianKernel<Dim<3>>&, const unsigned, const double, const double);
template TableKernel<Dim<3>>::TableKernel(const SuperGaussianKernel<Dim<3>>&, const unsigned, const double, const double);
template TableKernel<Dim<3>>::TableKernel(const PiGaussianKernel<Dim<3>>&, const unsigned, const double, const double);
template TableKernel<Dim<3>>::TableKernel(const HatKernel<Dim<3>>&, const unsigned, const double, const double);
template TableKernel<Dim<3>>::TableKernel(const SincKernel<Dim<3>>&, const unsigned, const double, const double);
template TableKernel<Dim<3>>::TableKernel(const NSincPolynomialKernel<Dim<3>>&, const unsigned, const double, const double);
template TableKernel<Dim<3>>::TableKernel(const QuarticSplineKernel<Dim<3>>&, const unsigned, const double, const double);
template TableKernel<Dim<3>>::TableKernel(const QuinticSplineKernel<Dim<3>>&, const unsigned, const double, const double);
template TableKernel<Dim<3>>::TableKernel(const WendlandC2Kernel<Dim<3>>&, const unsigned, const double, const double);
template TableKernel<Dim<3>>::TableKernel(const WendlandC4Kernel<Dim<3>>&, const unsigned, const double, const double);
template TableKernel<Dim<3>>::TableKernel(const WendlandC6Kernel<Dim<3>>&, const unsigned, const double, const double);
}

