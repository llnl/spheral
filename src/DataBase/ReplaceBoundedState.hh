//---------------------------------Spheral++----------------------------------//
// ReplaceBoundedState -- An implementation of UpdatePolicyBase appropriate for
// when 'ya just want to replace the state value with the new.
//
// This version enforces min/max bounds on the result, and therefore is only 
// appropriate for Scalar & SymTensor values.
//
// This version also assumes there is a derivative based update available.
//
// Created by JMO, Tue Aug 31 14:03:45 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral_ReplaceBoundedState_hh__
#define __Spheral_ReplaceBoundedState_hh__

#include "PureReplaceBoundedState.hh"
#include "Utilities/DBC.hh"

#include <limits>

namespace Spheral {

// Forward declarations.
template<typename Dimension> class StateDerivatives;

template<typename Dimension, typename ValueType, typename BoundValueType=ValueType>
class ReplaceBoundedState: public PureReplaceBoundedState<Dimension, ValueType, BoundValueType> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using KeyType = typename FieldUpdatePolicy<Dimension, ValueType>::KeyType;

  // Constructors, destructor.
  ReplaceBoundedState(const BoundValueType minValue = BoundValueType(std::numeric_limits<double>::lowest()),
                      const BoundValueType maxValue = BoundValueType(std::numeric_limits<double>::max()));
  ReplaceBoundedState(std::initializer_list<std::string> depends = {},
                      const BoundValueType minValue = BoundValueType(std::numeric_limits<double>::lowest()),
                      const BoundValueType maxValue = BoundValueType(std::numeric_limits<double>::max()));
  virtual ~ReplaceBoundedState() {}
  
  // An alternate method to be called when you want to specify that the "Replace" information
  // in the derivatives is invalid, and instead the value should be treated as a time advancement
  // algorithm instead.
  virtual void updateAsIncrement(const KeyType& key,
                                 State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivs,
                                 const double multiplier,
                                 const double t,
                                 const double dt) override;

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const override;

  // Don't use as evolved state for implicit integration
  virtual bool independent() const override { return false; }

private:
  //--------------------------- Private Interface ---------------------------//
  using PureReplaceBoundedState<Dimension, ValueType, BoundValueType>::mMinValue;
  using PureReplaceBoundedState<Dimension, ValueType, BoundValueType>::mMaxValue;

  ReplaceBoundedState(const ReplaceBoundedState& rhs);
  ReplaceBoundedState& operator=(const ReplaceBoundedState& rhs);
};

}

#include "ReplaceBoundedStateInline.hh"

#endif
