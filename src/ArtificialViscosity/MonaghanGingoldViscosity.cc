//---------------------------------Spheral++----------------------------------//
// A simple form for the artificial viscosity due to Monaghan & Gingold.
//----------------------------------------------------------------------------//
#include "MonaghanGingoldViscosity.hh"
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
#include "Utilities/Timer.hh"

using std::vector;
using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {

//------------------------------------------------------------------------------
// Construct with the given value for the linear and quadratic coefficients.
//------------------------------------------------------------------------------
template<typename Dimension>
MonaghanGingoldViscosity<Dimension>::
MonaghanGingoldViscosity(const Scalar Clinear,
                         const Scalar Cquadratic,
                         const TableKernel<Dimension>& kernel,
                         const bool linearInExpansion,
                         const bool quadraticInExpansion):
  ArtificialViscosity<Dimension>(Clinear, Cquadratic, kernel),
  mLinearInExpansion(linearInExpansion),
  mQuadraticInExpansion(quadraticInExpansion) {
  m_viewPtr = chai::make_managed<MonaghanGingoldViscosityView<Dimension>>(Clinear,
                                                                          Cquadratic,
                                                                          linearInExpansion,
                                                                          quadraticInExpansion);
}

//------------------------------------------------------------------------------
// Main method -- compute the QPi (P/rho^2) artificial viscosity
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE
void
MonaghanGingoldViscosityView<Dimension>::
QPiij(Scalar& QPiij, Scalar& QPiji,      // result for QPi (Q/rho^2)
      Scalar& Qij, Scalar& Qji,          // result for viscous pressure
      const unsigned nodeListi, const unsigned i,
      const unsigned nodeListj, const unsigned j,
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
  const auto ei = -Clij*csi*(mLinearInExpansion    ? mui                : min(0.0, mui)) +
                   Cqij    *(mQuadraticInExpansion ? -sgn(mui)*mui*mui  : FastMath::square(min(0.0, mui)));
  const auto ej = -Clij*csj*(mLinearInExpansion    ? muj                : min(0.0, muj)) +
                   Cqij    *(mQuadraticInExpansion ? -sgn(muj)*muj*muj  : FastMath::square(min(0.0, muj)));
  CHECK2(ei >= 0.0 or (mLinearInExpansion or mQuadraticInExpansion), ei << " " << csi << " " << mui);
  CHECK2(ej >= 0.0 or (mLinearInExpansion or mQuadraticInExpansion), ej << " " << csj << " " << muj);

  // Set the return values
  QPiij = ei/rhoi;
  QPiji = ej/rhoj;
  Qij = rhoi*ei;
  Qji = rhoj*ej;
}

}
