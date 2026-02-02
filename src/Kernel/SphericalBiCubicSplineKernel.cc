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

#include "SphericalBiCubicSplineKernel.hh"
#include "Kernel/BSplineKernel.hh"
#include "Utilities/FastMath.hh"
#include "Utilities/SpheralFunctions.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

namespace {  // anonymous

//------------------------------------------------------------------------------
// Useful functions from our derivation
//------------------------------------------------------------------------------
inline
double
C(const double q) {
  return q*q + (0.3*q - 0.75)*FastMath::pow4(q);
}

inline
double
D(const double q) {
  return 2.0*(q*q - FastMath::pow3(q)) + (0.75 - 0.1*q)*FastMath::pow4(q);
}

inline
double
gradC(const double q) {
  return 2.0*q + (1.5*q - 3.0)*FastMath::pow3(q);
}

inline
double
gradD(const double q) {
  return 4.0*q - 6.0*q*q + 3.0*FastMath::pow3(q) - 0.5*FastMath::pow4(q);
}

inline
double
W3S1(const double rj,
     const double ri,
     const double h) {
  const auto sigj = rj/h;
  const auto sigi = ri/h;
  const auto sigdiff = std::abs(sigj - sigi);
  const auto sigplus = sigj + sigi;
  auto result = 0.0;
  if (sigplus <= 1.0) {
    result = C(sigplus) - C(sigdiff);
  } else if (sigplus <= 2.0) {
    if (sigdiff < 1.0) {
      result = -0.1 + D(sigplus) - C(sigdiff);
    } else {
      result = D(sigplus) - D(sigdiff);
    }
  } else {
    if (sigdiff < 1.0) {
      result = 0.7 - C(sigdiff);
    } else if (sigdiff < 2.0) {
      result = 0.8 - D(sigdiff);
    }
  }
  return result/(h*rj*ri);
}

inline
double
gradW3S1(const double rj,
         const double ri,
         const double h) {
  const auto sigj = rj/h;
  const auto sigi = ri/h;
  const auto sigdiff = std::abs(sigj - sigi);
  const auto sigplus = sigj + sigi;
  const auto sgndiff = sgn0(sigi - sigj);
  auto result = 0.0;
  if (sigplus <= 1.0) {
    result = -W3S1(rj, ri, h)/ri + (gradC(sigplus) - gradC(sigdiff)*sgndiff)/(h*h*ri*rj);
  } else if (sigplus <= 2.0) {
    if (sigdiff < 1.0) {
      result = -W3S1(rj, ri, h)/ri + (gradD(sigplus) - gradC(sigdiff)*sgndiff)/(h*h*ri*rj);
    } else if (sigdiff < 2.0) {
      result = -W3S1(rj, ri, h)/ri + (gradD(sigplus) - gradD(sigdiff)*sgndiff)/(h*h*ri*rj);
    }
  } else {
    if (sigdiff < 1.0) {
      result = -W3S1(rj, ri, h)/ri - gradC(sigdiff)*sgndiff/(h*h*ri*rj);
    } else if (sigdiff < 2.0) {
      result = -W3S1(rj, ri, h)/ri - gradD(sigdiff)*sgndiff/(h*h*ri*rj);
    }
  }
  return result;
}

}            // anonymous

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
SphericalBiCubicSplineKernel::SphericalBiCubicSplineKernel(const unsigned numKernel):
  mBaseKernel3d(BSplineKernel<Dim<3>>(), numKernel),
  mBaseKernel1d(BSplineKernel<Dim<3>>(), numKernel),
  metamax(2.0) {
}

//------------------------------------------------------------------------------
// Copy constructor
//------------------------------------------------------------------------------
SphericalBiCubicSplineKernel::SphericalBiCubicSplineKernel(const SphericalBiCubicSplineKernel& rhs):
  mBaseKernel3d(rhs.mBaseKernel3d),
  mBaseKernel1d(rhs.mBaseKernel1d),
  metamax(rhs.metamax) {
}

//------------------------------------------------------------------------------
// Assignment
//------------------------------------------------------------------------------
SphericalBiCubicSplineKernel&
SphericalBiCubicSplineKernel::operator=(const SphericalBiCubicSplineKernel& rhs) {
  if (this != &rhs) {
    mBaseKernel3d = rhs.mBaseKernel3d;
    mBaseKernel1d = rhs.mBaseKernel1d;
    metamax = rhs.metamax;
  }
  return *this;
}

//------------------------------------------------------------------------------
// Lookup the kernel for (rj/h, ri/h) = (etaj, etai)
//------------------------------------------------------------------------------
double
SphericalBiCubicSplineKernel::operator()(const Dim<1>::Vector& etaj,
                                         const Dim<1>::Vector& etai,
                                         const Dim<1>::Scalar  Hdet) const {
  REQUIRE(Hdet >= 0.0);
  const auto h = 1.0/Hdet;
  return W3S1(h*etaj.x(), h*etai.x(), h);
}

//------------------------------------------------------------------------------
// Lookup the grad kernel for (rj/h, ri/h) = (etaj, etai)
// Using the Leibniz integral rule to differentiate this integral.
//------------------------------------------------------------------------------
Dim<1>::Vector
SphericalBiCubicSplineKernel::grad(const Dim<1>::Vector& etaj,
                                   const Dim<1>::Vector& etai,
                                   const Dim<1>::SymTensor& H) const {
  const auto Hdet = H.Determinant();
  REQUIRE(Hdet >= 0.0);
  const auto h = 1.0/Hdet;
  return Vector(gradW3S1(h*etaj.x(), h*etai.x(), h));
}

//------------------------------------------------------------------------------
// Simultaneously lookup (W,  grad W) for (rj/h, ri/h) = (etaj, etai)
//------------------------------------------------------------------------------
void
SphericalBiCubicSplineKernel::kernelAndGrad(const Dim<1>::Vector& etaj,
                                            const Dim<1>::Vector& etai,
                                            const Dim<1>::SymTensor& H,
                                            Dim<1>::Scalar& W,
                                            Dim<1>::Vector& gradW,
                                            Dim<1>::Scalar& deltaWsum) const {
  const auto Hdet = H.Determinant();
  REQUIRE(Hdet >= 0.0);
  const auto h = 1.0/Hdet;
  W = W3S1(h*etaj.x(), h*etai.x(), h);
  gradW.x(gradW3S1(h*etaj.x(), h*etai.x(), h));
  deltaWsum = mBaseKernel1d.gradValue((etaj - etai).magnitude(), Hdet);
}

}
