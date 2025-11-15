//---------------------------------Spheral++----------------------------------//
// Modified form of the standard SPH pair-wise viscosity due to Monaghan &
// Gingold.  This form is modified to use the velocity gradient to limit the
// velocity jump at the mid-point between points.
//
// This form specialized for use with the area-weighted RZ formalism.
//
// Created by JMO, Sun May 22 10:45:30 PDT 2016
//----------------------------------------------------------------------------//
#include "LimitedMonaghanGingoldViscosityRZ.hh"
#include "Boundary/Boundary.hh"
#include "DataOutput/Restart.hh"
#include "Field/FieldList.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "NodeList/FluidNodeList.hh"
#include "Neighbor/Neighbor.hh"
#include "Material/EquationOfState.hh"
#include "Boundary/Boundary.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/IncrementState.hh"

namespace Spheral {

using std::min;
using std::max;
using std::abs;
using std::pair;
using std::make_pair;

namespace {

//------------------------------------------------------------------------------
// limiter for velocity projection.
//------------------------------------------------------------------------------
double limiterBJ(const double x) {
  if (x > 0.0) {
    return min(1.0, 4.0/(x + 1.0)*min(1.0, x));  // Barth-Jesperson
  } else {
    return 0.0;
  }
}

double limiterMC(const double x) {
  if (x > 0.0) {
    return 2.0/(1.0 + x)*min(2.0*x, min(0.5*(1.0 + x), 2.0));   // monotonized central
  } else {
    return 0.0;
  }
}

double limiterVL(const double x) {
  if (x > 0.0) {
    return 2.0/(1.0 + x)*2.0*x/(1.0 + x);                       // van Leer
  } else {
    return 0.0;
  }
}

double limiterMM(const double x) {
  if (x > 0.0) {
    return 2.0/(1.0 + x)*min(1.0, x);                           // minmod
  } else {
    return 0.0;
  }
}

double limiterSB(const double x) {
  if (x > 0.0) {
    return 2.0/(1.0 + x)*max(min(2.0*x, 1.0), min(x, 2.0));    // superbee
  } else {
    return 0.0;
  }
}

}

//------------------------------------------------------------------------------
// Construct with the given value for the linear and quadratic coefficients.
//------------------------------------------------------------------------------
LimitedMonaghanGingoldViscosityRZ::
LimitedMonaghanGingoldViscosityRZ(const Scalar Clinear,
                                  const Scalar Cquadratic,
                                  const bool linearInExpansion,
                                  const bool quadraticInExpansion,
                                  const Scalar etaCritFrac,
                                  const Scalar etaFoldFrac):
  LimitedMonaghanGingoldViscosity<Dim<2> >(Clinear,
                                           Cquadratic, 
                                           linearInExpansion,
                                           quadraticInExpansion,
                                           etaCritFrac,
                                           etaFoldFrac) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
LimitedMonaghanGingoldViscosityRZ::
~LimitedMonaghanGingoldViscosityRZ() {
}

//------------------------------------------------------------------------------
// The required method to compute the artificial viscous P/rho^2.
//------------------------------------------------------------------------------
pair<Dim<2>::Tensor,
     Dim<2>::Tensor>
LimitedMonaghanGingoldViscosityRZ::
Piij(const unsigned nodeListi, const unsigned i, 
     const unsigned nodeListj, const unsigned j,
     const Vector& xi,
     const Vector& etai,
     const Vector& vi,
     const Scalar rhoi,
     const Scalar csi,
     const SymTensor& Hi,
     const Vector& xj,
     const Vector& etaj,
     const Vector& vj,
     const Scalar rhoj,
     const Scalar csj,
     const SymTensor& Hj) const {

  double Cl = this->mClinear;
  double Cq = this->mCquadratic;
  const double eps2 = this->mEpsilon2;
  const bool linearInExp = this->linearInExpansion();
  const bool quadInExp = this->quadraticInExpansion();
  const bool balsaraShearCorrection = this->mBalsaraShearCorrection;
  const Tensor& DvDxi = mGradVel(nodeListi, i);
  const Tensor& DvDxj = mGradVel(nodeListj, j);

  // Grab the FieldLists scaling the coefficients.
  // These incorporate things like the Balsara shearing switch or Morris & Monaghan time evolved
  // coefficients.
  const Scalar fCli = this->mClMultiplier(nodeListi, i);
  const Scalar fCqi = this->mCqMultiplier(nodeListi, i);
  const Scalar fClj = this->mClMultiplier(nodeListj, j);
  const Scalar fCqj = this->mCqMultiplier(nodeListj, j);
  Cl *= 0.5*(fCli + fClj);
  Cq *= 0.5*(fCqi + fCqj);

  // Are we applying the shear corrections?
  Scalar fshear = 1.0;
  if (balsaraShearCorrection) {
    const Scalar csneg = this->negligibleSoundSpeed();
    const Scalar hiinv = Hi.Trace()/Dimension::nDim;
    const Scalar hjinv = Hj.Trace()/Dimension::nDim;
    const Scalar ci = max(csneg, csi);
    const Scalar cj = max(csneg, csj);
    const Scalar fi = this->curlVelocityMagnitude(DvDxi)/(this->curlVelocityMagnitude(DvDxi) + abs(DvDxi.Trace()) + eps2*ci*hiinv);
    const Scalar fj = this->curlVelocityMagnitude(DvDxj)/(this->curlVelocityMagnitude(DvDxj) + abs(DvDxj.Trace()) + eps2*cj*hjinv);
    fshear = min(fi, fj);
  }

  // Compute the corrected velocity difference.
  const Vector xij = 0.5*(xi - xj);
  const Scalar gradi = (DvDxi.dot(xij)).dot(xij);
  const Scalar gradj = (DvDxj.dot(xij)).dot(xij);
  const Scalar ri = gradi/(sgn(gradj)*max(1.0e-30, abs(gradj)));
  const Scalar rj = gradj/(sgn(gradi)*max(1.0e-30, abs(gradi)));
  CHECK(min(ri, rj) <= 1.0);
  //const Scalar phi = limiterMM(min(ri, rj));
  Scalar phi = limiterVL(min(ri, rj));

  // If the points are getting too close, we let the Q come back full force.
  const Scalar etaij = min(etai.magnitude(), etaj.magnitude());
  // phi *= (etaij2 < etaCrit2 ? 0.0 : 1.0);
  // phi *= min(1.0, etaij2*etaij2/(etaCrit2etaCrit2));
  if (etaij < mEtaCrit) {
    phi *= exp(-FastMath::square((etaij - mEtaCrit)/mEtaFold));
  }

  // "Mike" method.
  const Vector vi1 = vi - phi*DvDxi*xij;
  const Vector vj1 = vj + phi*DvDxj*xij;
  
  const Vector vij = vi1 - vj1;
  
  // Compute mu.
  const Scalar mui = vij.dot(etai)/(etai.magnitude2() + eps2);
  const Scalar muj = vij.dot(etaj)/(etaj.magnitude2() + eps2);
  const Scalar mui_neg = min(0.0, mui);
  const Scalar muj_neg = min(0.0, muj);

  const Scalar vri = vi1.y();
  const Scalar vrj = vj1.y();
  const Scalar zetai = abs((Hi*xi).y());
  const Scalar zetaj = abs((Hj*xj).y());
  const Scalar muri = vri*safeInvVar(zetai);
  const Scalar murj = vrj*safeInvVar(zetaj);
  Scalar muri_neg = 0.0, murj_neg = 0.0;
  if (vri < 0.0 and vrj < 0.0) {
    const Vector vij0 = vi - vj;
    const Scalar fij = max(0.0, min(1.0, (vij0).dot(vij)*safeInv(vij0.magnitude2(), max(1e-10, eps2*max(csi, csj)))));
    muri_neg = fij*muri;
    murj_neg = fij*murj;
  }

  // The artificial internal energy.
  const Scalar ei = fshear*(-Cl*csi*(mLinearInExpansion    ? mui + muri                                : mui_neg + muri_neg) +
                             Cq    *(mQuadraticInExpansion ? -(sgn(mui)*mui*mui + sgn(muri)*muri*muri) : mui_neg*mui_neg + muri_neg*muri_neg));
  const Scalar ej = fshear*(-Cl*csj*(mLinearInExpansion    ? muj + murj                                : muj_neg + murj_neg) +
                             Cq    *(mQuadraticInExpansion ? -(sgn(muj)*muj*muj + sgn(murj)*murj*murj) : muj_neg*muj_neg + murj_neg*murj_neg));
  CHECK2(ei >= 0.0 or (mLinearInExpansion or mQuadraticInExpansion), ei << " " << csi << " " << mui);
  CHECK2(ej >= 0.0 or (mLinearInExpansion or mQuadraticInExpansion), ej << " " << csj << " " << muj);


  // Now compute the symmetrized artificial viscous pressure.
  return make_pair(ei/rhoi*Tensor::one(),
                   ej/rhoj*Tensor::one());
}

}
