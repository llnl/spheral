//---------------------------------Spheral++----------------------------------//
// PureReplaceBoundedState -- An implementation of UpdatePolicyBase appropriate for
// when 'ya just want to replace the state value with the new.
//
// This version enforces min/max bounds on the result, and therefore is only 
// appropriate for Scalar & SymTensor values.
// 
// Created by JMO, Tue Aug 31 14:03:45 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral_PureReplaceBoundedState_hh__
#define __Spheral_PureReplaceBoundedState_hh__

#include "FieldUpdatePolicy.hh"
#include "Utilities/DBC.hh"

#include <limits>

namespace Spheral {

// Forward declarations.
template<typename Dimension> class StateDerivatives;

template<typename Dimension, typename ValueType, typename BoundValueType=ValueType>
class PureReplaceBoundedState: public FieldUpdatePolicy<Dimension, ValueType> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using KeyType = typename FieldUpdatePolicy<Dimension, ValueType>::KeyType;

  // Constructors, destructor.
  PureReplaceBoundedState(const BoundValueType minValue = BoundValueType(std::numeric_limits<double>::lowest()),
                          const BoundValueType maxValue = BoundValueType(std::numeric_limits<double>::max()));
  PureReplaceBoundedState(std::initializer_list<std::string> depends = {},
                          const BoundValueType minValue = BoundValueType(std::numeric_limits<double>::lowest()),
                          const BoundValueType maxValue = BoundValueType(std::numeric_limits<double>::max()));
  virtual ~PureReplaceBoundedState() {}
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt) override;

  // Access the min and max's.
  BoundValueType minValue() const;
  BoundValueType maxValue() const;

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const override;

  static const std::string prefix() { return "new "; }

  // Don't use as evolved state for implicit integration
  virtual bool independent() const override { return false; }

protected:
  //--------------------------- Protected Interface ---------------------------//
  BoundValueType mMinValue;
  BoundValueType mMaxValue;

private:
  //--------------------------- Private Interface ---------------------------//
  PureReplaceBoundedState(const PureReplaceBoundedState& rhs);
  PureReplaceBoundedState& operator=(const PureReplaceBoundedState& rhs);
};

}

#include "PureReplaceBoundedStateInline.hh"

#endif
