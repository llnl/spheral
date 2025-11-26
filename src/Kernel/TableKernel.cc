//---------------------------------Spheral++----------------------------------//
// TableKernel -- Build an interpolation kernel using interpolation between
// tabulated points.
//
// Created by JMO, Mon Jun 19 21:06:28 PDT 2000
//----------------------------------------------------------------------------//
#include "Eigen/Dense"

#include "TableKernel.hh"

#include "Utilities/SpheralFunctions.hh"
#include "Utilities/bisectRoot.hh"
#include "Utilities/simpsonsIntegration.hh"
#include "Utilities/safeInv.hh"

using std::vector;

namespace Spheral {

namespace {  // anonymous

//------------------------------------------------------------------------------
// Sum the Kernel values for the given stepsize (SPH)
//------------------------------------------------------------------------------
inline
double
sumKernelValues(const TableKernel<Dim<1>>& W,
                const double nPerh) {
  REQUIRE(nPerh > 0.0);
  const auto deta = 1.0/nPerh;
  double result = 0.0;
  double etar = deta;
  while (etar < W.kernelExtent()) {
    result += 2.0*W.kernelValueSPH(etar);
    etar += deta;
  }
  return result;
}

inline
double
sumKernelValues(const TableKernel<Dim<2>>& W,
                const double nPerh) {
  REQUIRE(nPerh > 0.0);
  const auto deta = 1.0/nPerh;
  double result = 0.0;
  double etar = deta;
  while (etar < W.kernelExtent()) {
    result += 2.0*M_PI*etar/deta*W.kernelValueSPH(etar);
    etar += deta;
  }
  return sqrt(result);
}

inline
double
sumKernelValues(const TableKernel<Dim<3>>& W,
                const double nPerh) {
  REQUIRE(nPerh > 0.0);
  const auto deta = 1.0/nPerh;
  double result = 0.0;
  double etar = deta;
  while (etar < W.kernelExtent()) {
    result += 4.0*M_PI*FastMath::square(etar/deta)*W.kernelValueSPH(etar);
    etar += deta;
  }
  return pow(result, 1.0/3.0);
}

//------------------------------------------------------------------------------
// Compute the (f1,f2) integrals relation for the given zeta = r/h 
// (RZ corrections).
//------------------------------------------------------------------------------
template<typename KernelType>
class volfunc {
  const KernelType& W;
public:
  volfunc(const KernelType& W): W(W) {}
  double operator()(const double eta) const {
    return W.kernelValue(eta, 1.0);
  }
};

template<typename KernelType>
class f1func {
  const KernelType& W;
  double zeta;
public:
  f1func(const KernelType& W, const double zeta): W(W), zeta(zeta) {}
  double operator()(const double eta) const {
    return abs(safeInvVar(zeta)*eta)*W.kernelValue(abs(zeta - eta), 1.0);
  }
};


template<typename KernelType>
class f2func {
  const KernelType& W;
  double zeta;
public:
  f2func(const KernelType& W, const double zeta): W(W), zeta(zeta) {}
  double operator()(const double eta) const {
    return safeInvVar(zeta*zeta)*eta*abs(eta)*W.kernelValue(abs(zeta - eta), 1.0);
  }
};


template<typename KernelType>
class gradf1func {
  const KernelType& W;
  double zeta;
public:
  gradf1func(const KernelType& W, const double zeta): W(W), zeta(zeta) {}
  double operator()(const double eta) const {
    const double Wu = W.kernelValue(abs(zeta - eta), 1.0);
    const double gWu = W.gradValue(abs(zeta - eta), 1.0);
    const double gf1inv = safeInvVar(zeta)*abs(eta)*gWu - safeInvVar(zeta*zeta)*abs(eta)*Wu;
    if (eta < 0.0) {
      return -gf1inv;
    } else {
      return gf1inv;
    }
  }
};


template<typename KernelType>
double
f1Integral(const KernelType& W,
           const double zeta,
           const unsigned numbins) {
  const double etaMax = W.kernelExtent();
  CHECK(zeta <= etaMax);
  return safeInvVar(simpsonsIntegration<f1func<KernelType>, double, double>(f1func<KernelType>(W, zeta), 
                                                                            zeta - etaMax, 
                                                                            zeta + etaMax,
                                                                            numbins));
}

template<typename KernelType>
double
f2Integral(const KernelType& W,
           const double zeta,
           const unsigned numbins) {
  const double etaMax = W.kernelExtent();
  CHECK(zeta <= etaMax);
  return safeInvVar(simpsonsIntegration<f2func<KernelType>, double, double>(f2func<KernelType>(W, zeta), 
                                                                            zeta - etaMax, 
                                                                            zeta + etaMax,
                                                                            numbins));
}

template<typename KernelType>
double
gradf1Integral(const KernelType& W,
               const double zeta,
               const unsigned numbins) {
  const double etaMax = W.kernelExtent();
  CHECK(zeta <= etaMax);
  return simpsonsIntegration<gradf1func<KernelType>, double, double>(gradf1func<KernelType>(W, zeta), 
                                                                     zeta - etaMax, 
                                                                     zeta + etaMax,
                                                                     numbins);
}

}  // anonymous

//------------------------------------------------------------------------------
// Construct from a kernel.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename KernelType>
TableKernel<Dimension>::TableKernel(const KernelType& kernel,
                                    const unsigned numPoints,
                                    const typename Dimension::Scalar minNperh,
                                    const typename Dimension::Scalar maxNperh):
  Kernel<Dimension, TableKernel<Dimension>>(),
  mNumPoints(numPoints),
  mMinNperh(std::max(minNperh, 1.1/kernel.kernelExtent())),
  mMaxNperh(maxNperh),
  mInterp(0.0, kernel.kernelExtent(), numPoints,      [&](const double x) { return kernel(x, 1.0); }),
  mGradInterp(0.0, kernel.kernelExtent(), numPoints,  [&](const double x) { return kernel.grad(x, 1.0); }),
  mGrad2Interp(0.0, kernel.kernelExtent(), numPoints, [&](const double x) { return kernel.grad2(x, 1.0); }),
  mNperhLookup(),
  mWsumLookup() {

  // Gotta have a minimally reasonable nperh range
  if (mMaxNperh <= mMinNperh) mMaxNperh = 4.0*mMinNperh;

  // Pre-conditions.
  VERIFY2(mNumPoints > 0, "TableKernel ERROR: require numPoints > 0 : " << mNumPoints);
  VERIFY2(mMinNperh > 0.0 and mMaxNperh > mMinNperh, "TableKernel ERROR: Bad (minNperh, maxNperh) range: (" << mMinNperh << ", " << mMaxNperh << ")");

  // Set the volume normalization and kernel extent.
  this->setVolumeNormalization(1.0); // (kernel.volumeNormalization() / Dimension::pownu(hmult));  // We now build this into the tabular kernel values.
  this->setKernelExtent(kernel.kernelExtent());
  this->setInflectionPoint(kernel.inflectionPoint());

  // Set the interpolation methods for looking up nperh (SPH methodology)
  mWsumLookup.initialize(mMinNperh, mMaxNperh, numPoints,
                         [&](const double x) -> double { return sumKernelValues(*this, x); });
  mNperhLookup.initialize(mWsumLookup(mMinNperh), mWsumLookup(mMaxNperh), numPoints,
                          [&](const double Wsum) -> double { return bisectRoot([&](const double nperh) { return mWsumLookup(nperh) - Wsum; }, mMinNperh, mMaxNperh); });

  // Make nperh lookups monotonic
  mWsumLookup.makeMonotonic();
  mNperhLookup.makeMonotonic();
}

//------------------------------------------------------------------------------
// Copy
//------------------------------------------------------------------------------
template<typename Dimension>
TableKernel<Dimension>::
TableKernel(const TableKernel<Dimension>& rhs):
  Kernel<Dimension, TableKernel<Dimension>>(rhs),
  mNumPoints(rhs.mNumPoints),
  mMinNperh(rhs.mMinNperh),
  mMaxNperh(rhs.mMaxNperh),
  mInterp(rhs.mInterp),
  mGradInterp(rhs.mGradInterp),
  mGrad2Interp(rhs.mGrad2Interp),
  mNperhLookup(rhs.mNperhLookup),
  mWsumLookup(rhs.mWsumLookup) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
TableKernel<Dimension>::
~TableKernel() {
}

//------------------------------------------------------------------------------
// Assignment
//------------------------------------------------------------------------------
template<typename Dimension>
TableKernel<Dimension>&
TableKernel<Dimension>::
operator=(const TableKernel<Dimension>& rhs) {
  if (this != &rhs) {
    Kernel<Dimension, TableKernel<Dimension>>::operator=(rhs);
    mNumPoints = rhs.mNumPoints;
    mMinNperh = rhs.mMinNperh;
    mMaxNperh = rhs.mMaxNperh;
    mInterp = rhs.mInterp;
    mGradInterp = rhs.mGradInterp;
    mGrad2Interp = rhs.mGrad2Interp;
    mNperhLookup = rhs.mNperhLookup;
    mWsumLookup = rhs.mWsumLookup;
  }
  return *this;
}

//------------------------------------------------------------------------------
// Equivalence
//------------------------------------------------------------------------------
template<typename Dimension>
bool
TableKernel<Dimension>::
operator==(const TableKernel<Dimension>& rhs) const {
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
typename Dimension::Scalar
TableKernel<Dimension>::kernelValueSPH(const Scalar etaij) const {
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
typename Dimension::Scalar
TableKernel<Dimension>::kernelValueASPH(const Scalar etaij, const Scalar nPerh) const {
  REQUIRE(etaij >= 0.0);
  REQUIRE(nPerh > 0.0);
  if (etaij < mKernelExtent) {
    const auto deta = 2.0/std::max(2.0, nPerh);
    const auto eta0 = std::max(0.0, 0.5*(mKernelExtent - deta));
    const auto eta1 = std::min(mKernelExtent, eta0 + deta);
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
typename Dimension::Scalar
TableKernel<Dimension>::
equivalentNodesPerSmoothingScale(const Scalar Wsum) const {
  return std::max(0.0, mNperhLookup(Wsum));
}

//------------------------------------------------------------------------------
// Determine the effective Wsum we would expect for the given n per h.
// (SPH round tensor definition).
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
TableKernel<Dimension>::
equivalentWsum(const Scalar nPerh) const {
  return std::max(0.0, mWsumLookup(nPerh));
}

}
