#include "Geometry/Dimension.hh"
#include "Utilities/DBC.hh"

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
// Equivalence
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE bool
TableKernelView<Dimension>::
operator==(const TableKernelView<Dimension>& rhs) const {
  return ((mInterp == rhs.mInterp) and
          (mGradInterp == rhs.mGradInterp) and
          (mGrad2Interp == rhs.mGrad2Interp) and
          (mNperhLookup == rhs.mNperhLookup) and
          (mWsumLookup == rhs.mWsumLookup));
}

//------------------------------------------------------------------------------
// Kernel value for SPH smoothing scale nperh lookups
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE typename Dimension::Scalar
TableKernelView<Dimension>::kernelValueSPH(const Scalar etaij) const {
  REQUIRE(etaij >= 0.0);
  if (etaij < this->mKernelExtent) {
    return std::abs(mGradInterp(etaij));
  } else {
    return 0.0;
  }
}

//------------------------------------------------------------------------------
// Kernel value for ASPH smoothing scale nperh lookups
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE typename Dimension::Scalar
TableKernelView<Dimension>::kernelValueASPH(const Scalar etaij, const Scalar nPerh) const {
  REQUIRE(etaij >= 0.0);
  REQUIRE(nPerh > 0.0);
  if (etaij < this->mKernelExtent) {
    const auto deta = 2.0/std::max(2.0, nPerh);
    const auto eta0 = std::max(0.0, 0.5*(this->mKernelExtent - deta));
    const auto eta1 = std::min(this->mKernelExtent, eta0 + deta);
    return (etaij <= eta0 or etaij >= eta1 ?
            0.0 :
            kernelValueSPH((etaij - eta0)/deta));
            // FastMath::square(sin(M_PI*(etaij - eta0)/deta)));
    // return std::abs(mGradInterp(etaij * std::max(1.0, 0.5*nPerh*mKernelExtent))); // * FastMath::square(sin(nPerh*M_PI*etaij));
  } else {
    return 0.0;
  }
}

//------------------------------------------------------------------------------
// Determine the number of nodes per smoothing scale implied by the given
// sum of kernel values (SPH round tensor definition).
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE typename Dimension::Scalar
TableKernelView<Dimension>::
equivalentNodesPerSmoothingScale(const Scalar Wsum) const {
  return std::max(0.0, mNperhLookup(Wsum));
}

//------------------------------------------------------------------------------
// Determine the effective Wsum we would expect for the given n per h.
// (SPH round tensor definition).
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE typename Dimension::Scalar
TableKernelView<Dimension>::
equivalentWsum(const Scalar nPerh) const {
  return std::max(0.0, mWsumLookup(nPerh));
}

}
