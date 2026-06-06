//---------------------------------Spheral++----------------------------------//
// MaxReplaceState -- Replaces the state with the max of the state/deriv value
//
// J.M. Pearl 2024
//----------------------------------------------------------------------------//
#include "State.hh"
#include "StateDerivatives.hh"
#include "Field/Field.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructors.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType>
inline
MaxReplaceState<Dimension, ValueType>::
MaxReplaceState(std::initializer_list<std::string> depends):
  FieldUpdatePolicy<Dimension, ValueType>(depends) {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType>
inline
void
MaxReplaceState<Dimension, ValueType>::
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
    f(i) = std::max(df(i),f(i));
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType>
inline
bool
MaxReplaceState<Dimension, ValueType>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {
  return dynamic_cast<const MaxReplaceState<Dimension, ValueType>*>(&rhs) != nullptr;
}

}

