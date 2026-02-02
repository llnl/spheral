//---------------------------------Spheral++----------------------------------//
// YoungsModulusPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the dependent bulk modulus state.
//
// Created by JMO, Sun Jan 23 14:14:24 PST 2005
//----------------------------------------------------------------------------//
#include "YoungsModulusPolicy.hh"
#include "Strength/SolidFieldNames.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/UpdatePolicyBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/Field.hh"
#include "NodeList/FluidNodeList.hh"
#include "Material/EquationOfState.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

using std::abs;

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
YoungsModulusPolicy<Dimension>::
YoungsModulusPolicy(const SolidNodeList<Dimension>& nodes):
  FieldUpdatePolicy<Dimension, Scalar>({SolidFieldNames::bulkModulus,
                                        SolidFieldNames::shearModulus}),
  mSolidNodeList(nodes) {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
YoungsModulusPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& /*derivs*/,
       const double /*multiplier*/,
       const double /*t*/,
       const double /*dt*/) {
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == SolidFieldNames::YoungsModulus);
  auto& stateField = state.field(key, 0.0);

  // Get the bulk and shear modulus fields from the state.
  const auto  buildKey = [&](const std::string& fkey) { return StateBase<Dimension>::buildFieldKey(fkey, nodeListKey); };
  const auto& K = state.field(buildKey(SolidFieldNames::bulkModulus), 0.0);
  const auto& mu = state.field(buildKey(SolidFieldNames::shearModulus), 0.0);

  // Now set Youngs modulus.
  mSolidNodeList.YoungsModulus(stateField, K, mu);
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
YoungsModulusPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also a Youngs modulus operator.
  return dynamic_cast<const YoungsModulusPolicy<Dimension>*>(&rhs) != nullptr;
}

}

