//---------------------------------Spheral++----------------------------------//
// ReplaceWithRatioPolicy -- replaces one fieldlist with the ratio of two
// fieldlists from the state.
//
// J.M. Pearl 2022
//----------------------------------------------------------------------------//

#include "DataBase/ReplaceWithRatioPolicy.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/Field.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/DBC.hh"

#include <limits>

namespace Spheral {

//------------------------------------------------------------------------------
// Constructors.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
ReplaceWithRatioPolicy<Dimension, Value>::
ReplaceWithRatioPolicy(const KeyType& numerator,
                       const KeyType& denomenator):
  FieldUpdatePolicy<Dimension, Value>({}),
  mNumerator(numerator),
  mDenomenator(denomenator) {
}

template<typename Dimension, typename Value>
ReplaceWithRatioPolicy<Dimension, Value>::
ReplaceWithRatioPolicy(std::initializer_list<std::string> depends,
                       const KeyType& numerator,
                       const KeyType& denomenator):
  FieldUpdatePolicy<Dimension, Value>(depends),
  mNumerator(numerator),
  mDenomenator(denomenator) {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
void
ReplaceWithRatioPolicy<Dimension, Value>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double /*multiplier*/,
       const double /*t*/,
       const double /*dt*/) {

  const auto tiny = std::numeric_limits<typename Dimension::Scalar>::epsilon();

  // The state we're updating
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  auto& f = state.field(key, Value());

  // Find the matching replacement FieldList from the StateDerivatives.
  const auto  buildKey = [&](const std::string& fkey) { return StateBase<Dimension>::buildFieldKey(fkey, nodeListKey); };
  const auto& numer = state.field(buildKey(mNumerator), Value());
  const auto& denom = state.field(buildKey(mDenomenator), Value());
  CHECK(numer.size() == denom.size());

  // Loop over the internal values of the field.
  const auto n = f.numInternalElements();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    f(i) = numer(i)*safeInvVar(denom(i), tiny);
  }
}


//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
bool
ReplaceWithRatioPolicy<Dimension, Value>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an replace operator.
  const auto* rhsPtr = dynamic_cast<const ReplaceWithRatioPolicy<Dimension, Value>*>(&rhs);
  return rhsPtr != nullptr;
}

}

