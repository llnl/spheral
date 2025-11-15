//---------------------------------Spheral++----------------------------------//
// SphericalPositionPolicy
//
// Specializes the IncrementFieldListPolicy for use updating the position in
// spherical coordinates.  This position does not allow points to pass through
// the origin.
//
// Created by JMO, Tue Mar 29 11:12:31 PDT 2022
//----------------------------------------------------------------------------//
#include "SphericalPositionPolicy.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/FieldList.hh"
#include "Utilities/DBC.hh"

#include <regex>
#include <vector>
using std::vector;
using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
SphericalPositionPolicy::
SphericalPositionPolicy():
  UpdatePolicyBase<Dim<1>>() {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
SphericalPositionPolicy::
~SphericalPositionPolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
void
SphericalPositionPolicy::
update(const KeyType& key,
       State<Dim<1>>& state,
       StateDerivatives<Dim<1>>& derivs,
       const double multiplier,
       const double t,
       const double dt) {

  // Get the field name portion of the key.
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);

  // Get the state we're updating.
  auto f = state.field(key, Vector::zero());

  // Find all the available matching derivative Field keys.
  const auto incrementKey = prefix() + fieldKey;
  const auto allkeys = derivs.keys();
  KeyType dfKey, dfNodeListKey;
  for (const auto& key: allkeys) {
    StateBase<Dimension>::splitFieldKey(key, dfKey, dfNodeListKey);
    if (dfNodeListKey == nodeListKey and
        dfKey.compare(0, incrementKey.size(), incrementKey) == 0) {

      // This delta field matches the base of increment key, so apply it.
      const auto& df = derivs.field(key, Vector::zero());
      const auto  n = f.numInternalElements();
#pragma omp parallel for
      for (auto i = 0u; i < n; ++i) {
        // This is where we diverge from the standard IncrementState.  Ensure we cannot cross to
        // negative radius.
        f(i) = std::max(0.5*f(i), f(i) + multiplier*(df(i)));
      }
    }
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
bool
SphericalPositionPolicy::
operator==(const UpdatePolicyBase<Dim<1>>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const auto rhsPtr = dynamic_cast<const SphericalPositionPolicy*>(&rhs);
  return (rhsPtr != nullptr);
}

}

