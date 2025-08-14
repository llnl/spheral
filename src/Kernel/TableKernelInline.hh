#include "Geometry/Dimension.hh"
#include "Utilities/DBC.hh"
#include "VolumeIntegrationFunctions.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Return the kernel weight for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE inline
typename Dimension::Scalar
TableKernelView<Dimension>::kernelValue(const Scalar etaij, const Scalar Hdet) const {
  REQUIRE(etaij >= 0.0);
  REQUIRE(Hdet >= 0.0);
  if (etaij < this->mKernelExtent) {
    return Hdet*mInterp(etaij);
  } else {
    return 0.0;
  }
}

//------------------------------------------------------------------------------
// Return the gradient value for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE inline
typename Dimension::Scalar
TableKernelView<Dimension>::gradValue(const Scalar etaij, const Scalar Hdet) const {
  REQUIRE(etaij >= 0.0);
  REQUIRE(Hdet >= 0.0);
  if (etaij < this->mKernelExtent) {
    return Hdet*mGradInterp(etaij);
  } else {
    return 0.0;
  }
}

//------------------------------------------------------------------------------
// Return the second derivative value for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE inline
typename Dimension::Scalar
TableKernelView<Dimension>::grad2Value(const Scalar etaij, const Scalar Hdet) const {
  REQUIRE(etaij >= 0.0);
  REQUIRE(Hdet >= 0.0);
  if (etaij < this->mKernelExtent) {
    return Hdet*mGrad2Interp(etaij);
  } else {
    return 0.0;
  }
}

//------------------------------------------------------------------------------
// Return the kernel and gradient for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE inline
void
TableKernelView<Dimension>::kernelAndGrad(const typename Dimension::Vector& etaj,
                                          const typename Dimension::Vector& etai,
                                          const typename Dimension::SymTensor& H,
                                          typename Dimension::Scalar& W,
                                          typename Dimension::Vector& gradW,
                                          typename Dimension::Scalar& deltaWsum) const {
  const auto etaij = etai - etaj;
  const auto etaijMag = etaij.magnitude();
  const auto Hdet = H.Determinant();
  if (etaijMag < this->mKernelExtent) {
    const auto i0 = mInterp.lowerBound(etaijMag);
    W = Hdet*mInterp(etaijMag, i0);
    deltaWsum = Hdet*mGradInterp(etaijMag, i0);
    gradW = H*etaij.unitVector()*deltaWsum;
  } else {
    W = 0.0;
    deltaWsum = 0.0;
    gradW.Zero();
  }
}

//------------------------------------------------------------------------------
// Return the kernel and gradient value for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE inline
void
TableKernelView<Dimension>::kernelAndGradValue(const Scalar etaij, const Scalar Hdet,
                                               Scalar& Wi, Scalar& gWi) const {
  REQUIRE(etaij >= 0.0);
  REQUIRE(Hdet >= 0.0);
  if (etaij < this->mKernelExtent) {
    const auto i0 = mInterp.lowerBound(etaij);
    Wi = Hdet*mInterp(etaij, i0);
    gWi = Hdet*mGradInterp(etaij, i0);
  } else {
    Wi = 0.0;
    gWi = 0.0;
  }
}

//------------------------------------------------------------------------------
// Return the kernel and gradient values for a set of normalized distances.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
TableKernel<Dimension>::kernelAndGradValues(const std::vector<Scalar>& etaijs,
                                            const std::vector<Scalar>& Hdets,
                                            std::vector<Scalar>& kernelValues,
                                            std::vector<Scalar>& gradValues) const {
  // Preconditions.
  const auto n = etaijs.size();
  BEGIN_CONTRACT_SCOPE
  {
    REQUIRE(Hdets.size() == n);
    for (auto i = 0u; i < n; ++i) {
      REQUIRE(etaijs[i] >= 0.0);
      REQUIRE(Hdets[i] >= 0.0);
    }
  }
  END_CONTRACT_SCOPE

  // Prepare the results.
  kernelValues.resize(n);
  gradValues.resize(n);

  // Fill those suckers in.
  for (auto i = 0u; i < n; ++i) {
    const auto i0 = this->mInterp.lowerBound(etaijs[i]);
    kernelValues[i] = Hdets[i]*this->mInterp(etaijs[i], i0);
    gradValues[i] = Hdets[i]*this->mGradInterp(etaijs[i], i0);
  }
}

}

