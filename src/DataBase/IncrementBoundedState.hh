//---------------------------------Spheral++----------------------------------//
// IncrementBoundedState -- An implementation of UpdatePolicyBase appropriate for
// when 'ya just want to increment by derivatives:  x1 = x0 + A*dx/dt.
// This version enforces min/max bounds on the result, and therefore is only 
// appropriate for Scalar & SymTensor values.
//
// Created by JMO, Thu Aug 26 14:28:07 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral_IncrementBoundedState_hh__
#define __Spheral_IncrementBoundedState_hh__

#include "IncrementState.hh"

#include <limits>

namespace Spheral {

// Forward declarations.
template<typename Dimension, typename ValueType, typename BoundValueType=ValueType>
class IncrementBoundedState: public IncrementState<Dimension, ValueType> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using KeyType = typename FieldUpdatePolicy<Dimension, ValueType>::KeyType;

  // Constructors, destructor.
  IncrementBoundedState(std::initializer_list<std::string> depends = {},
                        const BoundValueType minValue = BoundValueType(std::numeric_limits<double>::lowest()),
                        const BoundValueType maxValue = BoundValueType(std::numeric_limits<double>::max()),
                        const bool wildCardDerivs = false);
  IncrementBoundedState(const BoundValueType minValue = BoundValueType(std::numeric_limits<double>::lowest()),
                        const BoundValueType maxValue = BoundValueType(std::numeric_limits<double>::max()),
                        const bool wildCardDerivs = false);
  virtual ~IncrementBoundedState() = default;
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt) override;

  // Access the min and max's.
  BoundValueType minValue() const { return mMinValue; }
  BoundValueType maxValue() const { return mMaxValue; }

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const override;

  // Forbidden methods
  IncrementBoundedState(const IncrementBoundedState& rhs) = delete;
  IncrementBoundedState& operator=(const IncrementBoundedState& rhs) = delete;

private:
  //--------------------------- Private Interface ---------------------------//
  BoundValueType mMinValue;
  BoundValueType mMaxValue;
};

}

#include "IncrementBoundedStateInline.hh"

#endif
