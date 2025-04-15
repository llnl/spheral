//---------------------------------Spheral++----------------------------------//
// IncrementState -- An implementation of UpdatePolicyBase appropriate for
// when 'ya just want to increment by derivatives:  x1 = x0 + A*dx/dt
//
// Created by JMO, Thu Aug 26 14:28:07 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral_IncrementState_hh__
#define __Spheral_IncrementState_hh__

#include "FieldUpdatePolicy.hh"

namespace Spheral {

template<typename Dimension, typename ValueType>
class IncrementState: public FieldUpdatePolicy<Dimension, ValueType> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using KeyType = typename FieldUpdatePolicy<Dimension, ValueType>::KeyType;

  // Constructors, destructor.
  IncrementState(std::initializer_list<std::string> depends = {},
                 const bool wildCardDerivs = false);
  IncrementState(const bool wildCardDerivs);
  virtual ~IncrementState() = default;
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt) override;

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const override;

  static const std::string prefix() { return "delta "; }

  // Flip whether we try to find multiple registered increment fields.
  bool wildCardDerivs()                                   const { return mWildCardDerivs; }
  void wildCardDerivs(const bool val)                           { mWildCardDerivs = val; }

  // Advance this policy implicitly
  virtual bool independent()                     const override { return true; }
  virtual void serializeDerivatives(std::vector<double>& buf,
                                    const KeyType& key,
                                    const StateDerivatives<Dimension>& derivs) const override;

  // Forbidden methods
  IncrementState(const IncrementState& rhs) = delete;
  IncrementState& operator=(const IncrementState& rhs) = delete;

private:
  //--------------------------- Private Interface ---------------------------//
  // Flag for looking for multiple increment derivatives.
  bool mWildCardDerivs;
};

}

#include "IncrementStateInline.hh"

#endif
