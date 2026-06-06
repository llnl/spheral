//---------------------------------Spheral++----------------------------------//
// FixedSmoothingScale
//
// Implements the static fixed smoothing scale option.
//
// Created by JMO, Wed Sep 14 13:50:49 PDT 2005
//----------------------------------------------------------------------------//
#ifndef __Spheral_FixedSmooothingScale__
#define __Spheral_FixedSmooothingScale__

#include "SmoothingScale/SmoothingScaleBase.hh"

namespace Spheral {

template<typename Dimension>
class FixedSmoothingScale: public SmoothingScaleBase<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  // Constructors, destructor.
  FixedSmoothingScale(): SmoothingScaleBase<Dimension>(HEvolutionType::FixedH, false, false) {};
  virtual ~FixedSmoothingScale() {};

  // Increment the derivatives.
  virtual 
  void evaluateDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivatives) const override {};

  // It's useful to have labels for Physics packages.  We'll require this to have
  // the same signature as the restart label.
  virtual std::string label() const override { return "FixedSmoothingScale"; }
};

}

#endif
