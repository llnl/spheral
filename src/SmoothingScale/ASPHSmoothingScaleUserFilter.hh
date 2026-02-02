//---------------------------------Spheral++----------------------------------//
// ASPHSmoothingScaleUserFilter
//
// Provides user-overridable hooks to modify how the ASPH ideal H algorithm
// is applied.
//
// Created by JMO, Mon Sep 23 15:03:26 PDT 2024
//----------------------------------------------------------------------------//
#ifndef __Spheral_ASPHSmooothingScaleUserFilter__
#define __Spheral_ASPHSmooothingScaleUserFilter__

namespace Spheral {

template<typename Dimension>
class ASPHSmoothingScaleUserFilter {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  // Constructors, destructor.
  ASPHSmoothingScaleUserFilter()           {}
  virtual ~ASPHSmoothingScaleUserFilter()  {}

  // Overridable hook called at the start of ASPHSmoothingScale::finalize.
  // Provides the opportunity to prepare for looping over each points new
  // ideal H vote.
  virtual void startFinalize(const Scalar time, 
                             const Scalar dt,
                             DataBase<Dimension>& dataBase, 
                             State<Dimension>& state,
                             StateDerivatives<Dimension>& derivs) {}

  // Overridable hook called for each point with both the old and new
  // H values.  Returns the new H value to use (defaults to the ideal H vote
  // for H1).
  virtual SymTensor __call__(size_t nodeListi,
                             size_t i,
                             const SymTensor& H0,
                             const SymTensor& H1)                 { return H1; }
  virtual SymTensor operator()(size_t nodeListi,
                               size_t i,
                               const SymTensor& H0,
                               const SymTensor& H1)               { return this->__call__(nodeListi, i, H0, H1); }
};

}

#endif
