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
    return std::abs(safeInvVar(zeta)*eta)*W.kernelValue(std::abs(zeta - eta), 1.0);
  }
};


template<typename KernelType>
class f2func {
  const KernelType& W;
  double zeta;
public:
  f2func(const KernelType& W, const double zeta): W(W), zeta(zeta) {}
  double operator()(const double eta) const {
    return safeInvVar(zeta*zeta)*eta*std::abs(eta)*W.kernelValue(std::abs(zeta - eta), 1.0);
  }
};


template<typename KernelType>
class gradf1func {
  const KernelType& W;
  double zeta;
public:
  gradf1func(const KernelType& W, const double zeta): W(W), zeta(zeta) {}
  double operator()(const double eta) const {
    const double Wu = W.kernelValue(std::abs(zeta - eta), 1.0);
    const double gWu = W.gradValue(std::abs(zeta - eta), 1.0);
    const double gf1inv = safeInvVar(zeta)*std::abs(eta)*gWu - safeInvVar(zeta*zeta)*std::abs(eta)*Wu;
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
  TableKernelView<Dimension>(),
  mInterpVal(0.0, kernel.kernelExtent(), numPoints,      [&](const double x) { return kernel(x, 1.0); }),
  mGradInterpVal(0.0, kernel.kernelExtent(), numPoints,  [&](const double x) { return kernel.grad(x, 1.0); }),
  mGrad2InterpVal(0.0, kernel.kernelExtent(), numPoints, [&](const double x) { return kernel.grad2(x, 1.0); }),
  mNperhLookupVal(),
  mWsumLookupVal() {

  // Gotta have a minimally reasonable nperh range
  this->mMaxNperh = maxNperh;
  this->mMinNperh = std::max(minNperh, 1.1/kernel.kernelExtent());
  this->mNumPoints = numPoints;
  if (this->mMaxNperh <= this->mMinNperh) this->mMaxNperh = 4.0*this->mMinNperh;

  // Pre-conditions.
  VERIFY2(this->mNumPoints > 0, "TableKernel ERROR: require numPoints > 0 : " << this->mNumPoints);
  VERIFY2(this->mMinNperh > 0.0 and this->mMaxNperh > this->mMinNperh, "TableKernel ERROR: Bad (minNperh, maxNperh) range: (" << this->mMinNperh << ", " << this->mMaxNperh << ")");
  this->mInterp = mInterpVal.view();
  this->mGradInterp = mGradInterpVal.view();
  this->mGrad2Interp = mGrad2InterpVal.view();
  // Set the volume normalization and kernel extent.
  this->setVolumeNormalization(1.0); // (kernel.volumeNormalization() / Dimension::pownu(hmult));  // We now build this into the tabular kernel values.
  this->setKernelExtent(kernel.kernelExtent());
  this->setInflectionPoint(kernel.inflectionPoint());

  // Set the interpolation methods for looking up nperh (SPH methodology)
  mWsumLookupVal.initialize(this->mMinNperh, this->mMaxNperh, numPoints,
                            [&](const double x) -> double { return sumKernelValues(*this, x); });
  mNperhLookupVal.initialize(mWsumLookupVal(this->mMinNperh), mWsumLookupVal(this->mMaxNperh), numPoints,
                             [&](const double Wsum) -> double { return bisectRoot([&](const double nperh) { return mWsumLookupVal(nperh) - Wsum; }, this->mMinNperh, this->mMaxNperh); });
  // Make nperh lookups monotonic
  mWsumLookupVal.makeMonotonic();
  mNperhLookupVal.makeMonotonic();
  this->mNperhLookup = mNperhLookupVal.view();
  this->mWsumLookup = mWsumLookupVal.view();
}

//------------------------------------------------------------------------------
// Copy
//------------------------------------------------------------------------------
template<typename Dimension>
TableKernel<Dimension>::
TableKernel(const TableKernel<Dimension>& rhs):
  TableKernelView<Dimension>(rhs),
  mInterpVal(rhs.mInterpVal),
  mGradInterpVal(rhs.mGradInterpVal),
  mGrad2InterpVal(rhs.mGrad2InterpVal),
  mNperhLookupVal(rhs.mNperhLookupVal),
  mWsumLookupVal(rhs.mWsumLookupVal) {
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
    TableKernelView<Dimension>::operator=(rhs);
    mInterpVal = rhs.mInterpVal;
    mGradInterpVal = rhs.mGradInterpVal;
    mGrad2InterpVal = rhs.mGrad2InterpVal;
    mNperhLookupVal = rhs.mNperhLookupVal;
    mWsumLookupVal = rhs.mWsumLookupVal;
  }
  return *this;
}

}
