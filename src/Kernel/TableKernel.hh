//---------------------------------Spheral++----------------------------------//
// TableKernel -- Build an interpolation kernel using interpolation between
// tabulated points.
//
// Created by JMO, Mon Jun 19 21:06:28 PDT 2000
//----------------------------------------------------------------------------//
#ifndef __Spheral_TableKernel_hh__
#define __Spheral_TableKernel_hh__

#include "Kernel.hh"
#include "Utilities/QuadraticInterpolator.hh"
#include "Utilities/CubicHermiteInterpolator.hh"
#include "config.hh"

#include <vector>

namespace Spheral {

template<typename Dimension>
class TableKernelView : public Kernel<Dimension, TableKernelView<Dimension> > {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using InterpolatorType = QuadraticInterpolator;
  using NperhInterpolatorType = CubicHermiteInterpolator;
  using IView = QIView;
  using NperhIView = CHIView;

  SPHERAL_HOST_DEVICE TableKernelView() = default;
  // Equivalence
  SPHERAL_HOST_DEVICE bool operator==(const TableKernelView& rhs) const;

  // Return the kernel weight for a given normalized distance or position.
  SPHERAL_HOST_DEVICE Scalar kernelValue(const Scalar etaij, const Scalar Hdet) const;

  // Return the gradient value for a given normalized distance or position.
  SPHERAL_HOST_DEVICE Scalar gradValue(const Scalar etaij, const Scalar Hdet) const;

  // Return the second derivative value for a given normalized distance or position.
  SPHERAL_HOST_DEVICE Scalar grad2Value(const Scalar etaij, const Scalar Hdet) const;

  // Simultaneously return the kernel value and first derivative.
  SPHERAL_HOST_DEVICE void kernelAndGrad(const Vector& etaj, const Vector& etai, const SymTensor& H,
                                         Scalar& W,
                                         Vector& gradW,
                                         Scalar& deltaWsum) const;
  SPHERAL_HOST_DEVICE void kernelAndGradValue(const Scalar etaij, const Scalar Hdet,
                                              Scalar& W,
                                              Scalar& gW) const;
    // Special kernel values for use in finding smoothing scales (SPH and ASPH versions)
  // ***These are only intended for use adapting smoothing scales***, and are used
  // for the succeeding equivalentNodesPerSmoothingScale lookups!
  SPHERAL_HOST_DEVICE Scalar kernelValueSPH(const Scalar etaij) const;
  SPHERAL_HOST_DEVICE Scalar kernelValueASPH(const Scalar etaij, const Scalar nPerh) const;

  // Return the equivalent number of nodes per smoothing scale implied by the given
  // sum of kernel values, using the zeroth moment SPH algorithm
  SPHERAL_HOST_DEVICE Scalar equivalentNodesPerSmoothingScale(const Scalar Wsum) const;
  SPHERAL_HOST_DEVICE Scalar equivalentWsum(const Scalar nPerh) const;

  // Access the internal data
  SPHERAL_HOST_DEVICE size_t numPoints() const                                    { return mNumPoints; }
  SPHERAL_HOST_DEVICE Scalar minNperhLookup() const                               { return mMinNperh; }
  SPHERAL_HOST_DEVICE Scalar maxNperhLookup() const                               { return mMaxNperh; }
protected:
  //--------------------------- Private Interface ---------------------------//
  // Data for the kernel tabulation.
  size_t mNumPoints = 100u;
  Scalar mMinNperh = 0.25;
  Scalar mMaxNperh = 64.0;
  IView mInterp, mGradInterp, mGrad2Interp; // W, grad W, grad^2 W
  NperhIView mNperhLookup, mWsumLookup;     // SPH nperh lookups

};

template<typename Dimension>
class TableKernel : public TableKernelView<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using InterpolatorType = QuadraticInterpolator;
  using NperhInterpolatorType = CubicHermiteInterpolator;

  // Constructors.
  template<typename KernelType>
  TableKernel(const KernelType& kernel,
              const unsigned numPoints = 100u,
              const Scalar minNperh = 0.25,
              const Scalar maxNperh = 64.0);
  TableKernel(const TableKernel<Dimension>& rhs);

  // Destructor.
  ~TableKernel();

  // Assignment.
  TableKernel& operator=(const TableKernel& rhs);

  // Look up the kernel and first derivative for a set.
  void kernelAndGradValues(const std::vector<Scalar>& etaijs,
                           const std::vector<Scalar>& Hdets,
                           std::vector<Scalar>& kernelValues,
                           std::vector<Scalar>& gradValues) const;

  // Direct access to our interpolators
  const InterpolatorType& Winterpolator() const               { return mInterpVal; }
  const InterpolatorType& gradWinterpolator() const           { return mGradInterpVal; }
  const InterpolatorType& grad2Winterpolator() const          { return mGrad2InterpVal; }
  const NperhInterpolatorType& nPerhInterpolator() const      { return mNperhLookupVal; }
  const NperhInterpolatorType& WsumInterpolator() const       { return mWsumLookupVal; }

  TableKernelView<Dimension> view() {
    return static_cast<TableKernelView<Dimension>>(*this);
  }

private:
  //--------------------------- Private Interface ---------------------------//
  // Data for the kernel tabulation.
  InterpolatorType mInterpVal, mGradInterpVal, mGrad2InterpVal; // W, grad W, grad^2 W
  NperhInterpolatorType mNperhLookupVal, mWsumLookupVal;        // SPH nperh lookups

};

}

#include "TableKernelInline.hh"

#endif

