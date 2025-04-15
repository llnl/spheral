//---------------------------------Spheral++----------------------------------//
// EntropyPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the entropy.
//
// Created by JMO, Sun Mar  6 21:57:47 PST 2016
//----------------------------------------------------------------------------//

#include "EntropyPolicy.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Strength/SolidFieldNames.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/ReplaceState.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/FieldList.hh"
#include "NodeList/FluidNodeList.hh"
#include "Material/EquationOfState.hh"
#include "Utilities/DBC.hh"

namespace Spheral {


//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
EntropyPolicy<Dimension>::
EntropyPolicy():
  FieldUpdatePolicy<Dimension, Scalar>({HydroFieldNames::massDensity,
                                        HydroFieldNames::specificThermalEnergy,
                                        SolidFieldNames::porositySolidDensity}) {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
EntropyPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& /*derivs*/,
       const double multiplier,
       const double t,
       const double dt) {
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == HydroFieldNames::entropy);
  auto& entropy = state.field(key, 0.0);

  // Get the eos.  This cast is ugly, but is a work-around for now.
  const auto* fluidNodeListPtr = dynamic_cast<const FluidNodeList<Dimension>*>(entropy.nodeListPtr());
  VERIFY(fluidNodeListPtr != nullptr);
  const auto& eos = fluidNodeListPtr->equationOfState();

  // Check if we're using porosity
  const auto buildKey = [&](const std::string& fkey) { return StateBase<Dimension>::buildFieldKey(fkey, nodeListKey); };
  const auto usePorosity = state.registered(buildKey(SolidFieldNames::porosityAlpha));

  // Grab the state we need.
  const auto& rhoS = (usePorosity ?
                      state.field(buildKey(SolidFieldNames::porositySolidDensity), 0.0) :
                      state.field(buildKey(HydroFieldNames::massDensity), 0.0));
  const auto& eps = state.field(buildKey(HydroFieldNames::specificThermalEnergy), 0.0);

  // Now set the entropy for this field.
  eos.setEntropy(entropy, rhoS, eps);
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
EntropyPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const auto* rhsPtr = dynamic_cast<const EntropyPolicy<Dimension>*>(&rhs);
  return (rhsPtr != nullptr);
}

}

