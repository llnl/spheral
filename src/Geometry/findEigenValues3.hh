#ifndef Spheral_findEigenValues3_hh
#define Spheral_findEigenValues3_hh

#include "Utilities/SpheralFunctions.hh"
#include "Utilities/FastMath.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Specialized helper to find the eigen values of a 3x3 tensor.
// This is just solving for the roots of the cubic equation, ala what's in 
// Numerical Recipes.
// Based on the notes by David Eberly at www.geometrictools.com.
//------------------------------------------------------------------------------
template<typename TensorType>
SPHERAL_HOST_DEVICE
inline
GeomVector<3>
findEigenValues3(const TensorType& tin) {

  // It is useful to scale the input so that all elements are in the range [-1,1]
  // to avoid numerical hoo-ha-ness.
  const double fscale = std::max(1.0, tin.maxAbsElement());
  CHECK(fscale >= 1.0);
  const double fscalei = 1.0/fscale;
  const TensorType t = fscalei*tin;

  const double onethird = 1.0/3.0;
  const double onetwentyseventh = 1.0/27.0;
  const double tol = 1.0e-15*std::max(1.0, t.doubledot(t)); // 1.0e-13*max(1.0, std::abs(t.Trace())); // 1.0e-13*sqrt(t.doubledot(t));

  const double c0 = t.Determinant();
  const double c1 = (t.xx()*t.yy() - t.xy()*t.yx() +
                     t.xx()*t.zz() - t.xz()*t.zx() +
                     t.yy()*t.zz() - t.yz()*t.zy());
  const double c2 = t.Trace();

  const double a = c1 - onethird*c2*c2;
  const double b = onetwentyseventh*(-2.0*c2*c2*c2 + 9.0*c1*c2) - c0;
  const double Q = 0.25*b*b + onetwentyseventh*a*a*a;

  const double x = onethird*c2;
  const TensorType diff = t - x*TensorType::one();
  if (fuzzyEqual(sqrt(diff.doubledot(diff)), 0.0, tol)) {

    // Single root case.  The input should just be a scalar multiple of I.
    return fscale*GeomVector<3>(x, x, x);

  } else if (fuzzyGreaterThanOrEqual(Q, 0.0, tol)) {

    // Two distinct roots.
    const double b2 = 0.5*b;
    const double ack = onethird*c2;
    const double thpt = sgn(b2) * FastMath::CubeRootHalley2(std::abs(b2));
    return fscale*GeomVector<3>(ack + thpt, ack + thpt, ack - 2.0*thpt);

  } else {

    // Three distinct roots.
    CHECK(Q < 0.0);
    const double b2 = 0.5*b;
    const double theta = onethird*std::atan2(sqrt(-Q), -b2);
    const double rho = std::sqrt(std::max(0.0, b2*b2 - Q));
    const double rho13 = FastMath::CubeRootHalley2(rho);
    const double ack = onethird*c2;
    const double thpt = rho13*std::cos(theta);
    const double bleah = rho13*std::sqrt(3.0)*std::sin(theta);
    return fscale*GeomVector<3>(ack + 2.0*thpt,
                                ack - thpt - bleah,
                                ack - thpt + bleah);

  }
}

}

#endif
