//---------------------------------Spheral++----------------------------------//
// DeviatoricStressPolicy.
//
// Specialized version of the IncrementFieldList policy, with some criteria for
// zeroing out the deviatoric stress in special cases.
//
// Created by JMO, Mon Feb  6 11:34:57 PST 2006
//----------------------------------------------------------------------------//
#ifndef __Spheral_DeviatoricStress_hh__
#define __Spheral_DeviatoricStress_hh__

#include "DataBase/FieldUpdatePolicy.hh"

#include <string>

namespace Spheral {

template<typename Dimension>
class DeviatoricStressPolicy: public FieldUpdatePolicy<Dimension, typename Dimension::SymTensor> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using KeyType = typename UpdatePolicyBase<Dimension>::KeyType;

  // Constructors, destructor.
  DeviatoricStressPolicy();
  virtual ~DeviatoricStressPolicy() = default;
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt) override;

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const override;

  // Forbidden methods
  DeviatoricStressPolicy(const DeviatoricStressPolicy& rhs) = delete;
  DeviatoricStressPolicy& operator=(const DeviatoricStressPolicy& rhs) = delete;
};

}

#endif
