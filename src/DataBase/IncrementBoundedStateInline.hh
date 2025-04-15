//---------------------------------Spheral++----------------------------------//
// IncrementBoundedState -- An implementation of UpdatePolicyBase appropriate for
// when 'ya just want to increment by derivatives:  x1 = x0 + A*dx/dt
//
// Created by JMO, Thu Aug 26 14:28:07 2004
//----------------------------------------------------------------------------//
#include "State.hh"
#include "StateDerivatives.hh"
#include "Field/Field.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructors.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType, typename BoundValueType>
inline
IncrementBoundedState<Dimension, ValueType, BoundValueType>::
IncrementBoundedState(std::initializer_list<std::string> depends,
                      const BoundValueType minValue,
                      const BoundValueType maxValue,
                      const bool wildCardDerivs):
  IncrementState<Dimension, ValueType>(depends, wildCardDerivs),
  mMinValue(minValue),
  mMaxValue(maxValue) {
}

template<typename Dimension, typename ValueType, typename BoundValueType>
inline
IncrementBoundedState<Dimension, ValueType, BoundValueType>::
IncrementBoundedState(const BoundValueType minValue,
                      const BoundValueType maxValue,
                      const bool wildCardDerivs):
  IncrementState<Dimension, ValueType>({}, wildCardDerivs),
  mMinValue(minValue),
  mMaxValue(maxValue) {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType, typename BoundValueType>
inline
void
IncrementBoundedState<Dimension, ValueType, BoundValueType>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double /*t*/,
       const double /*dt*/) {

  // Get the field name portion of the key.
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);

  // Get the state we're updating.
  auto& f = state.field(key, ValueType());

  // Find all the available matching derivative Field keys.
  const auto incrementKey = this->prefix() + fieldKey;
  const auto allkeys = derivs.keys();
  KeyType dfKey, dfNodeListKey;
  auto numDeltaFields = 0u;
  for (const auto& dkey: allkeys) {
    StateBase<Dimension>::splitFieldKey(dkey, dfKey, dfNodeListKey);
    if (dfNodeListKey == nodeListKey and
        dfKey.compare(0, incrementKey.size(), incrementKey) == 0) {
      ++numDeltaFields;

      // This delta field matches the base of increment key, so apply it.
      const auto& df = derivs.field(dkey, ValueType());
      const auto  n = f.numInternalElements();
#pragma omp parallel for
      for (auto i = 0u; i < n; ++i) {
        f(i) = std::min(mMaxValue, std::max(mMinValue, f(i) + multiplier*(df(i))));
      }
    }
  }

  // If we're not allowing wildcard update, there should have only be one match.
  VERIFY2(this->wildCardDerivs() or numDeltaFields == 1,
          "IncrementBoundedState ERROR: unable to find unique match for derivative field key " << incrementKey << " on NodeList " << nodeListKey << " : " << numDeltaFields << " matches");
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType, typename BoundValueType>
inline
bool
IncrementBoundedState<Dimension, ValueType, BoundValueType>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const auto rhsPtr = dynamic_cast<const IncrementBoundedState<Dimension, ValueType, BoundValueType>*>(&rhs);
  if (rhsPtr == nullptr) return false;

  // Ok, now do we agree on min & max?
  return (minValue() == rhsPtr->minValue()) && (maxValue() == rhsPtr->maxValue());
}

}
