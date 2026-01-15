//---------------------------------Spheral++----------------------------------//
// A modified form of the Monaghan & Gingold viscosity, extended to tensor 
// formalism.
//----------------------------------------------------------------------------//
#include "Utilities/rotationMatrix.hh"
#include "Utilities/GeometricUtilities.hh"

namespace Spheral {

namespace {

//------------------------------------------------------------------------------
// Helper to remove any expansion terms from DvDx
//------------------------------------------------------------------------------
template<typename Tensor>
SPHERAL_HOST_DEVICE
inline
void
removeExpansion(Tensor& DvDx) {
  const auto DvDx_s = DvDx.Symmetric();
  const auto DvDx_a = DvDx.SkewSymmetric();
  const auto eigeni = DvDx_s.eigenVectors();
  DvDx = constructTensorWithMinDiagonal(eigeni.eigenValues, 0.0);
  DvDx.rotationalTransform(eigeni.eigenVectors);
  DvDx += DvDx_a;
}

}

//------------------------------------------------------------------------------
// Main method -- compute the QPi (P/rho^2) artificial viscosity
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE
void
TensorMonaghanGingoldViscosityView<Dimension>::
QPiij(Tensor& QPiij, Tensor& QPiji,      // result for QPi (Q/rho^2)
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
  REQUIRE(DvDx.size() > std::max(nodeListi, nodeListj));

  // A few useful constants
  const auto tiny = 1.0e-20;
  const auto multipliers = fCl.size() > 0u;

  // If the nodes are not closing, then skip the rest and the Q for this pair is zero.
  const auto xij = xi - xj;
  const auto vij = vi - vj;
  if (vij.dot(xij) < 0.0) {

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

    // Some more geometry.
    const auto xij2 = xij.magnitude2();
    const auto xijUnit = xij.unitVector();
    const auto hi2 = xij2/(etai.magnitude2() + tiny);
    const auto hj2 = xij2/(etaj.magnitude2() + tiny);
    const auto hi = sqrt(hi2);
    const auto hj = sqrt(hj2);

    // Build the tensor for grad-v with the pair-wise value sliced in
    Tensor sigmai = DvDxi;
    Tensor sigmaj = DvDxj;
    {
      const auto R = rotationMatrix(xijUnit);
      const auto Rinverse = R.Transpose();
      const auto thpt1 = sqrt(xij2)*(R*vij);
      const auto deltaSigmai = thpt1/(xij2 + mEpsilon2*hi2);
      const auto deltaSigmaj = thpt1/(xij2 + mEpsilon2*hj2);
      sigmai.rotationalTransform(R);
      sigmaj.rotationalTransform(R);
      sigmai.setColumn(0, deltaSigmai);
      sigmaj.setColumn(0, deltaSigmaj);
      sigmai.rotationalTransform(Rinverse);
      sigmaj.rotationalTransform(Rinverse);
    }

    // Remove any expansive components from sigma
    removeExpansion(sigmai);
    removeExpansion(sigmaj);

    // Calculate the tensor viscous internal energy.
    const auto mui = hi*sigmai;
    const auto Qepsi = -Clij*csi*mui.Transpose() + Cqij*mui*mui;

    const auto muj = hj*sigmaj;
    const auto Qepsj = -Clij*csj*muj.Transpose() + Cqij*muj*muj;

    // We now have enough to compute Pi and Q!
    QPiij = Qepsi/rhoi;
    QPiji = Qepsj/rhoj;
    Qij = rhoi*rhoi*(QPiij.diagonalElements().maxAbsElement());
    Qji = rhoj*rhoj*(QPiji.diagonalElements().maxAbsElement());
  } else {
    QPiij = Tensor::zero();
    QPiji = Tensor::zero();
    Qij = 0.0;
    Qji = 0.0;
  }
}

}

