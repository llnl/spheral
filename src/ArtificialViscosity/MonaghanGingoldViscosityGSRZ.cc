//---------------------------------Spheral++----------------------------------//
// A simple form for the artificial viscosity due to Monaghan & Gingold.
//----------------------------------------------------------------------------//
#include "MonaghanGingoldViscosityGSRZ.hh"
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

using std::vector;
using std::pair;
using std::make_pair;

namespace Spheral {

//------------------------------------------------------------------------------
// Construct with the given value for the linear and quadratic coefficients.
//------------------------------------------------------------------------------
MonaghanGingoldViscosityGSRZ::
MonaghanGingoldViscosityGSRZ(const Scalar Clinear,
                             const Scalar Cquadratic,
                             const bool linearInExpansion,
                             const bool quadraticInExpansion):
  MonaghanGingoldViscosity<Dim<2> >(Clinear,
                                    Cquadratic,
                                    linearInExpansion,
                                    quadraticInExpansion) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
MonaghanGingoldViscosityGSRZ::
~MonaghanGingoldViscosityGSRZ() {
}

//------------------------------------------------------------------------------
// The required method to compute the artificial viscous P/rho^2.
//------------------------------------------------------------------------------
pair<Dim<2>::Tensor,
     Dim<2>::Tensor>
MonaghanGingoldViscosityGSRZ::
Piij(const size_t nodeListi, const size_t i, 
     const size_t nodeListj, const size_t j,
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

  // Grab the FieldLists scaling the coefficients.
  // These incorporate things like the Balsara shearing switch or Morris & Monaghan time evolved
  // coefficients.
  const Scalar fCli = this->mClMultiplier(nodeListi, i);
  const Scalar fCqi = this->mCqMultiplier(nodeListi, i);
  const Scalar fClj = this->mClMultiplier(nodeListj, j);
  const Scalar fCqj = this->mCqMultiplier(nodeListj, j);
  const Scalar fshear = std::max(this->mShearCorrection(nodeListi, i), this->mShearCorrection(nodeListj, j));
  Cl *= 0.5*(fCli + fClj)*fshear;
  Cq *= 0.5*(fCqi + fCqj)*fshear;

  // Compute the RZ density and such.
  const Scalar ri = abs(xi.y());
  const Scalar rj = abs(xj.y());
  const Scalar vri = vi.y();
  const Scalar vrj = vj.y();
  const Scalar zetai = abs((Hi*xi).y());
  const Scalar zetaj = abs((Hj*xj).y());
  const Scalar rhoRZi = 2.0*M_PI*ri*rhoi;
  const Scalar rhoRZj = 2.0*M_PI*rj*rhoj;

  // Compute mu.
  const Vector vij = vi - vj;
  const Scalar mui = vij.dot(etai)/(etai.magnitude2() + eps2);
  const Scalar muj = vij.dot(etaj)/(etaj.magnitude2() + eps2);
  const Scalar mui_neg = min(0.0, mui);
  const Scalar muj_neg = min(0.0, muj);

  const Scalar muri = vri*safeInvVar(zetai);
  const Scalar murj = vrj*safeInvVar(zetaj);
  Scalar muri_neg = 0.0, murj_neg = 0.0;
  if (vri < 0.0 and vrj < 0.0) {
    muri_neg = muri;
    murj_neg = murj;
  }

  // The artificial internal energy.
  const Scalar ei = -Cl*csi*(mLinearInExpansion    ? mui + muri                                : mui_neg + muri_neg) +
                     Cq    *(mQuadraticInExpansion ? -(sgn(mui)*mui*mui + sgn(muri)*muri*muri) : mui_neg*mui_neg + muri_neg*muri_neg);
  const Scalar ej = -Cl*csj*(mLinearInExpansion    ? muj + murj                                : muj_neg + murj_neg) +
                     Cq    *(mQuadraticInExpansion ? -(sgn(muj)*muj*muj + sgn(murj)*murj*murj) : muj_neg*muj_neg + murj_neg*murj_neg);
  CHECK2(ei >= 0.0 or (mLinearInExpansion or mQuadraticInExpansion), ei << " " << csi << " " << mui << " " << muri);
  CHECK2(ej >= 0.0 or (mLinearInExpansion or mQuadraticInExpansion), ej << " " << csj << " " << muj << " " << murj);

  // Now compute the symmetrized artificial viscous pressure.
  return make_pair(ei/rhoRZi*Tensor::one(),
                   ej/rhoRZj*Tensor::one());
}

}
