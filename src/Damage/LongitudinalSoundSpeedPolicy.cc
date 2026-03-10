//---------------------------------Spheral++----------------------------------//
// LongitudinalSoundSpeedPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the dependent longitudinal sound speed.
//
// Created by JMO, Sun Jan 23 14:14:24 PST 2005
//----------------------------------------------------------------------------//
#include "LongitudinalSoundSpeedPolicy.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Strength/SolidFieldNames.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/UpdatePolicyBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/Field.hh"
#include "NodeList/FluidNodeList.hh"
#include "Material/EquationOfState.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/DBC.hh"

using std::vector;
using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
LongitudinalSoundSpeedPolicy<Dimension>::
LongitudinalSoundSpeedPolicy(const SolidNodeList<Dimension>& nodes):
  FieldUpdatePolicy<Dimension, Scalar>({SolidFieldNames::YoungsModulus,
                                        SolidFieldNames::bulkModulus,
                                        SolidFieldNames::shearModulus,
                                        HydroFieldNames::massDensity}),
  mSolidNodeList(nodes) {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
LongitudinalSoundSpeedPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& /*derivs*/,
       const double /*multiplier*/,
       const double /*t*/,
       const double /*dt*/) {
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == SolidFieldNames::longitudinalSoundSpeed);
  auto& stateField = state.field(key, 0.0);

  // Get the bulk, shear, and Youngs modulus, and the mass density fields 
  // from the state.
  const auto  buildKey = [&](const std::string& fkey) { return StateBase<Dimension>::buildFieldKey(fkey, nodeListKey); };
  const auto& rho = state.field(buildKey(HydroFieldNames::massDensity), 0.0);
  const auto& K = state.field(buildKey(SolidFieldNames::bulkModulus), 0.0);
  const auto& mu = state.field(buildKey(SolidFieldNames::shearModulus), 0.0);

  // Now set the longitudinal sound speed.
  mSolidNodeList.longitudinalSoundSpeed(stateField, rho, K, mu);
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
LongitudinalSoundSpeedPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also a Youngs modulus operator.
  return dynamic_cast<const LongitudinalSoundSpeedPolicy<Dimension>*>(&rhs) != nullptr;
}

}

