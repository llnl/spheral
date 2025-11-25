//---------------------------------Spheral++----------------------------------//
// TableKernelView
//
// View for the TableKernel class.
// Created by LDO, Wed Oct 29 16:00:28 PDT 2025
//----------------------------------------------------------------------------//
#ifndef __Spheral_TableKernelView_hh__
#define __Spheral_TableKernelView_hh__

#include "Kernel.hh"
#include "Utilities/QuadraticInterpolatorView.hh"
#include "Utilities/CubicHermiteInterpolatorView.hh"
#include "config.hh"

namespace Spheral {

template<typename Dimension>
class TableKernelView : public Kernel<Dimension, TableKernelView<Dimension> > {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using IView = QuadraticInterpolatorView;
  using NperhIView = CubicHermiteInterpolatorView;

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
  SPHERAL_HOST_DEVICE void kernelAndGrad(const Vector& etaj,
                                         const Vector& etai,
                                         const SymTensor& H,
                                         Scalar& W,
                                         Vector& gradW,
                                         Scalar& deltaWsum) const;

  SPHERAL_HOST_DEVICE void kernelAndGradValue(const Scalar etaij,
                                              const Scalar Hdet,
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
  SPHERAL_HOST_DEVICE size_t numPoints() const      { return mNumPoints; }
  SPHERAL_HOST_DEVICE Scalar minNperhLookup() const { return mMinNperh; }
  SPHERAL_HOST_DEVICE Scalar maxNperhLookup() const { return mMaxNperh; }
protected:
  //--------------------------- Private Interface ---------------------------//
  // Data for the kernel tabulation.
  size_t mNumPoints = 100u;
  Scalar mMinNperh = 0.25;
  Scalar mMaxNperh = 64.0;
  IView mInterp, mGradInterp, mGrad2Interp; // W, grad W, grad^2 W
  NperhIView mNperhLookup, mWsumLookup;     // SPH nperh lookups

};

}

#include "TableKernelViewInline.hh"

#endif

