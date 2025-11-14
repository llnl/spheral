//---------------------------------Spheral++----------------------------------//
// JohnsonCookDamagePolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the tensor damage variable with the Johnson-Cook model.
//
// Created by JMO, Thu Jul 12 15:33:26 PDT 2018
//----------------------------------------------------------------------------//
#include <algorithm>

#include "JohnsonCookDamagePolicy.hh"
#include "NodeList/SolidNodeList.hh"
#include "Strength/SolidFieldNames.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/UpdatePolicyBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/Field.hh"
#include "Kernel/TableKernel.hh"
#include "Utilities/DBC.hh"

namespace Spheral {


//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
JohnsonCookDamagePolicy<Dimension>::
JohnsonCookDamagePolicy():
  UpdatePolicyBase<Dimension>() {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
JohnsonCookDamagePolicy<Dimension>::
~JohnsonCookDamagePolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
JohnsonCookDamagePolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double /*t*/,
       const double /*dt*/) {
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == SolidFieldNames::tensorDamage);
  auto& stateField = state.field(key, SymTensor::zero);

  // Get the state fields.
  const auto flawsKey = State<Dimension>::buildFieldKey(SolidFieldNames::flaws, nodeListKey);
  const auto psrKey = State<Dimension>::buildFieldKey(SolidFieldNames::plasticStrainRate, nodeListKey);
  CHECK(state.registered(flawsKey));
  CHECK(derivs.registered(psrKey));

  const auto& flaws = state.field(flawsKey, 0.0);
  const auto& plasticStrainRate = derivs.field(psrKey, 0.0);

  // Iterate over the internal nodes.
  const auto ni = stateField.numInternalElements();
#pragma omp parallel for
  for (auto i = 0u; i < ni; ++i) {

    // Update the damage.  We take advantage of the fact this is a scalar update here.
    auto Di = stateField(i).xx();
    Di = std::max(0.0, std::min(1.0, Di + multiplier*plasticStrainRate(i)*safeInvVar(flaws(i))));
    stateField(i) = Di*SymTensor::one;
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
JohnsonCookDamagePolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const JohnsonCookDamagePolicy<Dimension>* rhsPtr = dynamic_cast<const JohnsonCookDamagePolicy<Dimension>*>(&rhs);
  if (rhsPtr == 0) {
    return false;
  } else {
    return true;
  }
}

}

