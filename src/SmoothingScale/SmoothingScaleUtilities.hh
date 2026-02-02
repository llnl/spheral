#include "Utilities/GeometricUtilities.hh"
#include "Utilities/rotationMatrix.hh"

namespace Spheral {
namespace SmoothingScaleDetail {

//------------------------------------------------------------------------------
// Convert a given number of neighbors to the equivalent 1D "radius" in nodes.
//------------------------------------------------------------------------------
template<typename Dimension> inline double equivalentRadius(const double n);

// 1D
template<>
inline double
equivalentRadius<Dim<1>>(const double n) {
  return 0.5*n;
}

// 2D
template<>
inline double
equivalentRadius<Dim<2>>(const double n) {
  return std::sqrt(n/M_PI);
}

// 3D
template<>
inline double
equivalentRadius<Dim<3>>(const double n) {
  return Dim<3>::rootnu(3.0*n/(4.0*M_PI));
}

//------------------------------------------------------------------------------
// DH/Dt per dimension
//------------------------------------------------------------------------------
// 1-D case same as SPH.
inline
Dim<1>::SymTensor
smoothingScaleDerivative(const Dim<1>::SymTensor& H,
                         const Dim<1>::Tensor& DvDx) {
  return -H*DvDx.Trace();
}

// 2-D ASPH tensor evolution.
inline
Dim<2>::SymTensor
smoothingScaleDerivative(const Dim<2>::SymTensor& H,
                         const Dim<2>::Tensor& DvDx) {
  REQUIRE(H.Trace() > 0.0);
  const auto thetaDot = (H.xx()*DvDx.xy() - H.yy()*DvDx.yx() - H.yx()*(DvDx.xx() - DvDx.yy()))/H.Trace();
  Dim<2>::SymTensor result;
  result.xx(H.yx()*(thetaDot - DvDx.yx()) - H.xx()*DvDx.xx());
  result.xy(-(H.xx()*thetaDot + H.yx()*DvDx.xx() + H.yy()*DvDx.yx()));
  result.yy(-H.yx()*(thetaDot + DvDx.xy()) - H.yy()*DvDx.yy());
  return result;
}

// 3-D ASPH tensor evolution.
inline
Dim<3>::SymTensor
smoothingScaleDerivative(const Dim<3>::SymTensor& H,
                         const Dim<3>::Tensor& DvDx) {
  REQUIRE(H.Trace() > 0.0);
  const auto AA = H.xx()*DvDx.xy() - H.xy()*(DvDx.xx() - DvDx.yy()) + H.xz()*DvDx.zy() - H.yy()*DvDx.yx() - H.yz()*DvDx.zx();
  const auto BB = H.xx()*DvDx.xz() + H.xy()*DvDx.yz() - H.xz()*(DvDx.xx() - DvDx.zz()) - H.yz()*DvDx.yx() - H.zz()*DvDx.zx();
  const auto CC = H.xy()*DvDx.xz() + H.yy()*DvDx.yz() - H.yz()*(DvDx.yy() - DvDx.zz()) - H.xz()*DvDx.xy() - H.zz()*DvDx.zy();
  const auto thpt = H.yy() + H.zz();
  const auto Ga = (H.xx() + H.yy())*thpt - H.xz()*H.xz();
  const auto Gb = (H.yy() + H.zz())*H.yz() + H.xy()*H.xz();
  const auto Gc = (H.xx() + H.zz())*thpt - H.xy()*H.xy();
  const auto Gd = thpt*AA + H.xz()*CC;
  const auto Ge = thpt*BB - H.xy()*CC;
  const auto ack = 1.0/(Ga*Gc - Gb*Gb);
  const auto Gdot = (Gc*Gd - Gb*Ge)*ack;
  const auto Tdot = (Gb*Gd - Ga*Ge)*ack;
  const auto Phidot = (H.xz()*Gdot + H.xy()*Tdot + CC)/thpt;
  Dim<3>::SymTensor result;
  result.xx(-H.xx()*DvDx.xx() + H.xy()*(Gdot - DvDx.yx()) - H.xz()*(Tdot + DvDx.zx()));
  result.xy(H.yy()*Gdot - H.yz()*Tdot - H.xx()*DvDx.xy() - H.xy()*DvDx.yy() - H.xz()*DvDx.zy());
  result.xz(H.yz()*Gdot - H.zz()*Tdot - H.xx()*DvDx.xz() - H.xy()*DvDx.yz() - H.xz()*DvDx.zz());
  result.yy(H.yz()*(Phidot - DvDx.zy()) - H.xy()*(Gdot + DvDx.xy()) - H.yy()*DvDx.yy());
  result.yz(H.xy()*Tdot - H.yy()*Phidot - H.xz()*DvDx.xy() - H.yz()*DvDx.yy() - H.zz()*DvDx.zy());
  result.zz(H.xz()*(Tdot - DvDx.xz()) - H.yz()*(Phidot + DvDx.yz()) - H.zz()*DvDx.zz());
  return result;
}

//------------------------------------------------------------------------------
// Radial evolution/alignment specialized evolution
//------------------------------------------------------------------------------
// 1-D
inline
Dim<1>::SymTensor
radialEvolution(const Dim<1>::SymTensor& Hi,
                const Dim<1>::Vector& nhat,
                const Dim<1>::Scalar s,
                const Dim<1>::Scalar r0,
                const Dim<1>::Scalar r1) {
  return Hi / s;
}

// 2-D
inline
Dim<2>::SymTensor
radialEvolution(const Dim<2>::SymTensor& Hi,
                const Dim<2>::Vector& nhat,
                const Dim<2>::Scalar s,
                const Dim<2>::Scalar r0,
                const Dim<2>::Scalar r1) {
  const auto T = rotationMatrix(nhat).Transpose();
  const auto hev0 = Hi.eigenVectors();
  Dim<2>::SymTensor result;
  if (abs(hev0.eigenVectors.getColumn(0).dot(nhat)) > abs(hev0.eigenVectors.getColumn(1).dot(nhat))) {
    result(0,0) = hev0.eigenValues(0);
    result(1,1) = hev0.eigenValues(1);
  } else {
    result(0,0) = hev0.eigenValues(1);
    result(1,1) = hev0.eigenValues(0);
  }
  const auto fr = r1*safeInvVar(r0);
  CHECK(fr > 0.0);
  result(0,0) /= s;
  result(1,1) /= fr;
  result.rotationalTransform(T);
  return result;
}

// 3-D
inline
Dim<3>::SymTensor
radialEvolution(const Dim<3>::SymTensor& Hi,
                const Dim<3>::Vector& nhat,
                const Dim<3>::Scalar s,
                const Dim<3>::Scalar r0,
                const Dim<3>::Scalar r1) {
  const auto Tprinciple = rotationMatrix(nhat);
  const auto Tlab = Tprinciple.Transpose();
  auto result = Hi;
  result.rotationalTransform(Tprinciple);
  const auto fr = r1*safeInvVar(r0);
  CHECK(fr > 0.0);
  result(0,0) /= s;
  result(1,1) /= fr;
  result(1,2) /= fr;
  result(2,1) /= fr;
  result(2,2) /= fr;
  result.rotationalTransform(Tlab);
  return result;
}

}
}
