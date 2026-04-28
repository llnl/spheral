//---------------------------------Spheral++----------------------------------//
// A simple form for the artificial viscosity due to Monaghan & Gingold.
//----------------------------------------------------------------------------//
#include "Field/FieldList.hh"
#include "DataBase/DataBase.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Main method -- compute the QPi (P/rho^2) artificial viscosity
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE
void
MonaghanGingoldViscosityView<Dimension>::
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
  REQUIRE((not mBalsaraShearCorrection) or DvDx.size() > std::max(nodeListi, nodeListj));

  // A few useful constants
  const auto multipliers = fCl.size() > 0u;

  // Find our linear and quadratic coefficients
  const auto fCli = multipliers ? fCl(nodeListi, i) : 1.0;
  const auto fCqi = multipliers ? fCq(nodeListi, i) : 1.0;
  const auto fClj = multipliers ? fCl(nodeListj, j) : 1.0;
  const auto fCqj = multipliers ? fCq(nodeListj, j) : 1.0;
  const auto fshear = (mBalsaraShearCorrection ?
                       0.5*(this->calcBalsaraShearCorrection(DvDx(nodeListi, i), Hi, csi) +
                            this->calcBalsaraShearCorrection(DvDx(nodeListj, j), Hj, csj)) :
                       1.0);
  const auto Clij = 0.5*(fCli + fClj)*fshear * mClinear;
  const auto Cqij = 0.5*(fCqi + fCqj)*fshear * mCquadratic;

  // Compute mu.
  const auto vij = vi - vj;
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
