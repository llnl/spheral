//---------------------------------Spheral++----------------------------------//
// MFVIncrementVelocityPolicy -- specialized policy for hydros that allow for mass
//                      flux between nodes. The momentum time derivative
//                      is used to update the velocity. The "hydro acceleration"
//                      is also added in to be compatible w/ phys packages
//                      that apply a pure acceleration. 
//
// J.M. Pearl 2023
//----------------------------------------------------------------------------//

#include "GSPH/Policies/MFVIncrementVelocityPolicy.hh"
#include "GSPH/GSPHFieldNames.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/Field.hh"
#include "Utilities/DBC.hh"
#include "Hydro/HydroFieldNames.hh"

#include <limits.h>

namespace Spheral {
//------------------------------------------------------------------------------
// Constructors.
//------------------------------------------------------------------------------
template<typename Dimension>
MFVIncrementVelocityPolicy<Dimension>::
MFVIncrementVelocityPolicy(std::initializer_list<std::string> depends):
  FieldUpdatePolicy<Dimension, Vector>(depends){
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MFVIncrementVelocityPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {

  const auto tiny = std::numeric_limits<typename Dimension::Scalar>::epsilon();

  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);

  const auto massKey = StateBase<Dimension>::buildFieldKey(HydroFieldNames::mass, nodeListKey);
  const auto momDerivFieldKey = StateBase<Dimension>::buildFieldKey(prefix() + GSPHFieldNames::momentum, nodeListKey);
  const auto accDerivFieldKey = StateBase<Dimension>::buildFieldKey(HydroFieldNames::hydroAcceleration, nodeListKey);
  
  const auto&  m = state.field(massKey, Scalar());
        auto&  v = state.field(key,     Vector());

  const auto&  DmDt = derivs.field(prefix() + massKey, Scalar());
  const auto&  DpDt = derivs.field(momDerivFieldKey,   Vector());
        auto&  DvDt = derivs.field(accDerivFieldKey,   Vector());

  // set the hydro acceleration
  const auto n = m.numInternalElements();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    const auto m1 = m(i)+DmDt(i)*multiplier;
    const auto DpDti = DpDt(i);
    if (m1 > tiny) DvDt(i) += (DpDti - DmDt(i)*v(i)) * safeInv(m1);
  }

  // now kick to copy/paste of the increment method normally used to update velocity.

  // Find all the available matching derivative Field keys.
  const auto incrementKey = prefix() + fieldKey;
  const auto allkeys = derivs.keys();

  KeyType dfKey, dfNodeListKey;
  auto numDeltaFields = 0u;
  CONTRACT_VAR(numDeltaFields);
  for (const auto& key: allkeys) {
    StateBase<Dimension>::splitFieldKey(key, dfKey, dfNodeListKey);
    if (dfNodeListKey == nodeListKey and
        dfKey.compare(0, incrementKey.size(), incrementKey) == 0) {
      ++numDeltaFields;

      // This delta field matches the base of increment key, so apply it.
      const auto& dv = derivs.field(key, Vector());
      const auto  n = v.numInternalElements();
#pragma omp parallel for
      for (auto i = 0u; i < n; ++i) {
        v(i) += multiplier*(dv(i));
      }
    }
  }

  // If we're not allowing wildcard update, there should have only be one match.
  VERIFY2(numDeltaFields >= 1,
          "IncrementState ERROR: unable to find match for derivative field key " << incrementKey);
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
MFVIncrementVelocityPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an replace operator.
  const auto* rhsPtr = dynamic_cast<const MFVIncrementVelocityPolicy<Dimension>*>(&rhs);
  return rhsPtr != nullptr;
}

}

