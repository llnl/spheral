//---------------------------------Spheral++----------------------------------//
// PureReplaceBoundedState -- An implementation of UpdatePolicyBase appropriate for
// when 'ya just want to replace the state value with the new.
//
// This version enforces min/max bounds on the result, and therefore is only 
// appropriate for Scalar & SymTensor values.
// 
// Created by JMO, Tue Aug 31 14:03:45 2004
//----------------------------------------------------------------------------//
#include "State.hh"
#include "StateDerivatives.hh"
#include "Field/Field.hh"
#include "Utilities/DBC.hh"
#include "Utilities/SpheralFunctions.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructors.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType, typename BoundValueType>
inline
PureReplaceBoundedState<Dimension, ValueType, BoundValueType>::
PureReplaceBoundedState(const BoundValueType minValue,
                        const BoundValueType maxValue):
  FieldUpdatePolicy<Dimension, ValueType>({}),
  mMinValue(minValue),
  mMaxValue(maxValue) {
}

template<typename Dimension, typename ValueType, typename BoundValueType>
inline
PureReplaceBoundedState<Dimension, ValueType, BoundValueType>::
PureReplaceBoundedState(std::initializer_list<std::string> depends,
                        const BoundValueType minValue,
                        const BoundValueType maxValue):
  FieldUpdatePolicy<Dimension, ValueType>(depends),
  mMinValue(minValue),
  mMaxValue(maxValue) {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType, typename BoundValueType>
inline
void
PureReplaceBoundedState<Dimension, ValueType, BoundValueType>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double /*multiplier*/,
       const double /*t*/,
       const double /*dt*/) {

  // Find the matching replacement field from the StateDerivatives.
  const auto  replaceKey = prefix() + key;
  auto&       f = state.field(key, ValueType());
  const auto& df = derivs.field(replaceKey, ValueType());

  // Loop over the internal values of the field.
  const auto n = f.nodeList().numInternalNodes();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    f(i) = min(mMaxValue, max(mMinValue, df(i)));
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType, typename BoundValueType>
inline
bool
PureReplaceBoundedState<Dimension, ValueType, BoundValueType>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {
  return dynamic_cast<const PureReplaceBoundedState<Dimension, ValueType, BoundValueType>*>(&rhs) != nullptr;
}

//------------------------------------------------------------------------------
// Min value.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType, typename BoundValueType>
inline
BoundValueType
PureReplaceBoundedState<Dimension, ValueType, BoundValueType>::
minValue() const {
  return mMinValue;
}

//------------------------------------------------------------------------------
// Max value.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType, typename BoundValueType>
inline
BoundValueType
PureReplaceBoundedState<Dimension, ValueType, BoundValueType>::
maxValue() const {
  return mMaxValue;
}

}
