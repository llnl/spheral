//---------------------------------Spheral++----------------------------------//
// TableKernel -- Build an interpolation kernel using interpolation between
// tabulated points.
//
// Created by JMO, Mon Jun 19 21:06:28 PDT 2000
//----------------------------------------------------------------------------//
#ifndef __Spheral_TableKernel_hh__
#define __Spheral_TableKernel_hh__

#include "Kernel.hh"
#include "TableKernelView.hh"
#include "Utilities/QuadraticInterpolator.hh"
#include "Utilities/CubicHermiteInterpolator.hh"
#include "config.hh"

#include <vector>

namespace Spheral {

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

  // Assignment.
  TableKernel& operator=(const TableKernel& rhs);

  // Look up the kernel and first derivative for a set.
  void kernelAndGradValues(const std::vector<Scalar>& etaijs,
                           const std::vector<Scalar>& Hdets,
                           std::vector<Scalar>& kernelValues,
                           std::vector<Scalar>& gradValues) const;

  // Direct access to our interpolators
  const InterpolatorType& Winterpolator() const          { return mInterpVal; }
  const InterpolatorType& gradWinterpolator() const      { return mGradInterpVal; }
  const InterpolatorType& grad2Winterpolator() const     { return mGrad2InterpVal; }
  const NperhInterpolatorType& nPerhInterpolator() const { return mNperhLookupVal; }
  const NperhInterpolatorType& WsumInterpolator() const  { return mWsumLookupVal; }

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

