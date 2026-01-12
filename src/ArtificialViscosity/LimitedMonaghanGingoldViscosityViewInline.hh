//---------------------------------Spheral++----------------------------------//
// Modified form of the standard SPH pair-wise viscosity due to Monaghan &
// Gingold.  This form is modified to use the velocity gradient to limit the
// velocity jump at the mid-point between points.
//
// Created by JMO, Thu Nov 20 14:13:18 PST 2014
//----------------------------------------------------------------------------//
#include "Field/FieldList.hh"

namespace Spheral {

namespace {

//------------------------------------------------------------------------------
// limiter for velocity projection.
//------------------------------------------------------------------------------
//double limiterBJ(const double x) {
//  if (x > 0.0) {
//    return min(1.0, 4.0/(x + 1.0)*min(1.0, x));  // Barth-Jesperson
//  } else {
//    return 0.0;
//  }
//}
//
//double limiterMC(const double x) {
//  if (x > 0.0) {
//    return 2.0/(1.0 + x)*min(2.0*x, min(0.5*(1.0 + x), 2.0));   // monotonized central
//  } else {
//    return 0.0;
//  }
//}

SPHERAL_HOST_DEVICE
double limiterVL(const double x) {
  if (x > 0.0) {
    return 2.0/(1.0 + x)*2.0*x/(1.0 + x);                       // van Leer
  } else {
    return 0.0;
  }
}

//double limiterMM(const double x) {
//  if (x > 0.0) {
//    return 2.0/(1.0 + x)*min(1.0, x);                           // minmod
//  } else {
//    return 0.0;
//  }
//}
//
//double limiterSB(const double x) {
//  if (x > 0.0) {
//    return 2.0/(1.0 + x)*max(min(2.0*x, 1.0), min(x, 2.0));    // superbee
//  } else {
//    return 0.0;
//  }
//}

// template<typename Vector, typename Tensor>
// double limiterConservative(const Vector& vi, const Vector& vj,
//                            const Vector& xi, const Vector& xj,
//                            const Tensor& DvDxi, const Tensor& DvDxj) {
//   const auto xji = xj - xi;
//   const auto vji = vj - vi;
//   const auto di = DvDxi.dot(xji);
//   const auto dj = DvDxj.dot(xji);
//   if (di.dot(dj) <= 0.0 or
//       di.dot(vji) <= 0.0 or
//       dj.dot(vji) <= 0.0) return 0.0;
//   const auto vjimag = vji.magnitude();
//   const auto dimag = di.magnitude();
//   const auto djmag = dj.magnitude();
//   return min(1.0, min(abs(vjimag*safeInv(dimag)), abs(vjimag*safeInv(djmag))));
// }

//double limiterConservative(const double vji, const double deltavi, const double deltavj) {
//  if (deltavi*deltavj <= 0.0 or
//      deltavi*vji <= 0.0 or
//      deltavj*vji <= 0.0) return 0.0;
//  return min(1.0, min(abs(vji*safeInv(deltavi)), abs(vji*safeInv(deltavj))));
//}

}

//------------------------------------------------------------------------------
// Main method -- compute the QPi (P/rho^2) artificial viscosity
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE
void
LimitedMonaghanGingoldViscosityView<Dimension>::
QPiij(Scalar& QPiij, Scalar& QPiji,      // result for QPi (Q/rho^2)
      Scalar& Qij, Scalar& Qji,          // result for viscous pressure
      const size_t nodeListi, const size_t i,
      const size_t nodeListj, const size_t j,
      const Vector& xi,
      const SymTensor& Hi,
      const Vector& etai,
      const Vector& vi,
      const Scalar rhoi,
      const Scalar csi,
      const Vector& xj,
      const SymTensor& Hj,
      const Vector& etaj,
      const Vector& vj,
      const Scalar rhoj,
      const Scalar csj,
      const FieldListView<Dimension, Scalar>& fCl,
      const FieldListView<Dimension, Scalar>& fCq,
      const FieldListView<Dimension, Tensor>& DvDx) const {
  // Preconditions
  REQUIRE(fCl.size() == fCq.size());
  REQUIRE(DvDx.size() > std::max(nodeListi, nodeListj));

  // A few useful constants
  const auto multipliers = fCl.size() > 0u;

  // We need nPerh to figure out our critical folding distance. We assume the first NodeList value for this is
  // correct for all of them...
  const auto nPerhi = DvDx[nodeListi].nodesPerSmoothingScale();
  const auto nPerhj = DvDx[nodeListj].nodesPerSmoothingScale();
  const auto nPerh = 0.5*(nPerhi + nPerhj);
  const auto etaCrit = mEtaCritFrac/nPerh;
  const auto etaFold = mEtaFoldFrac/nPerh;
  CHECK(etaFold > 0.0);

  // Find our linear and quadratic coefficients
  const auto fCli = multipliers ? fCl(nodeListi, i) : 1.0;
  const auto fCqi = multipliers ? fCq(nodeListi, i) : 1.0;
  const auto fClj = multipliers ? fCl(nodeListj, j) : 1.0;
  const auto fCqj = multipliers ? fCq(nodeListj, j) : 1.0;
  const auto& DvDxi = DvDx(nodeListi, i);
  const auto& DvDxj = DvDx(nodeListj, j);
  const auto fshear = (mBalsaraShearCorrection ?
                       0.5*(this->calcBalsaraShearCorrection(DvDxi, Hi, csi) +
                            this->calcBalsaraShearCorrection(DvDxj, Hj, csj)) :
                       1.0);
  const auto Clij = 0.5*(fCli + fClj)*fshear * mClinear;
  const auto Cqij = 0.5*(fCqi + fCqj)*fshear * mCquadratic;

  // Compute our limited velocity gradient along the line connecting these points
  const auto xij = 0.5*(xi - xj);  // midpoint distance
  const auto gradi = (DvDxi.dot(xij)).dot(xij);
  const auto gradj = (DvDxj.dot(xij)).dot(xij);
  const auto ri = gradi/(sgn(gradj)*std::max(1.0e-30, std::abs(gradj)));
  const auto rj = gradj/(sgn(gradi)*std::max(1.0e-30, std::abs(gradi)));
  CHECK(std::min(ri, rj) <= 1.0);
  // const Scalar phi = limiterMM(min(ri, rj));
  auto phi =  limiterVL(std::min(ri, rj));
  
  // const auto xjihat = -xij.unitVector();
  // auto phi = limiterConservative((vj - vi).dot(xjihat), (DvDxi*xjihat).dot(xjihat), (DvDxj*xjihat).dot(xjihat));

  // If the points are getting too close, we let the Q come back full force.
  const auto etaij = std::min(etai.magnitude(), etaj.magnitude());
  // phi *= (etaij2 < etaCrit2 ? 0.0 : 1.0);
  // phi *= min(1.0, etaij2*etaij2/(etaCrit2etaCrit2));
  if (etaij < etaCrit) {
    phi *= std::exp(-FastMath::square((etaij - etaCrit)/etaFold));
  }

  // Compute the corrected velocity difference.
  // "Mike" method.
  const auto vi1 = vi - phi*DvDxi*xij;
  const auto vj1 = vj + phi*DvDxj*xij;

  // const Vector vi1 = vi - phi*DvDxi*xij;
  // const Vector vj1 = vj + phi*DvDxj*xij;
  
  const auto vij = vi1 - vj1;
  
  // Compute mu.
  const auto mui = vij.dot(etai)/(etai.magnitude2() + mEpsilon2);
  const auto muj = vij.dot(etaj)/(etaj.magnitude2() + mEpsilon2);

  // The artificial internal energy.
  const auto ei = -Clij*csi*(mLinearInExpansion    ? mui                : std::min(0.0, mui)) +
    Cqij    *(mQuadraticInExpansion ? -sgn(mui)*mui*mui  : FastMath::square(std::min(0.0, mui)));
  const auto ej = -Clij*csj*(mLinearInExpansion    ? muj                : std::min(0.0, muj)) +
    Cqij    *(mQuadraticInExpansion ? -sgn(muj)*muj*muj  : FastMath::square(std::min(0.0, muj)));
  CHECK2(ei >= 0.0 or (mLinearInExpansion or mQuadraticInExpansion), ei << " " << csi << " " << mui);
  CHECK2(ej >= 0.0 or (mLinearInExpansion or mQuadraticInExpansion), ej << " " << csj << " " << muj);

  // Set the return values
  QPiij = ei/rhoi;
  QPiji = ej/rhoj;
  Qij = rhoi*ei;
  Qji = rhoj*ej;
}

}
