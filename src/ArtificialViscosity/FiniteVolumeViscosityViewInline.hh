//---------------------------------Spheral++----------------------------------//
// A finite-volume based viscosity.  Assumes you have constructred the 
// tessellation in the state.
//
// Created by JMO, Tue Aug 13 09:43:37 PDT 2013
//----------------------------------------------------------------------------//

namespace Spheral {

//------------------------------------------------------------------------------
// Main method -- compute the QPi (P/rho^2) artificial viscosity
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE
void
FiniteVolumeViscosityView<Dimension>::
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

  // Compute the pair QPi
  const auto xji = xj - xi;
  const auto xjihat = xji.unitVector();
  const auto hi = 1.0/(Hi*xjihat).magnitude();
  const auto hj = 1.0/(Hj*xjihat).magnitude();
  const auto DvDxi = std::min(0.0, DvDx(nodeListi, i).Trace());
  const auto DvDxj = std::min(0.0, DvDx(nodeListj, j).Trace());
  QPiij = (-Clij*csi*DvDxi + Cqij*hi*DvDxi*DvDxi)*hi/rhoi;
  QPiji = (-Clij*csj*DvDxj + Cqij*hj*DvDxj*DvDxj)*hj/rhoj;
  Qij = rhoi*rhoi*QPiij;
  Qji = rhoi*rhoi*QPiji;
}

}
