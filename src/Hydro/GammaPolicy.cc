//---------------------------------Spheral++----------------------------------//
// GammaPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the dependent gamma state.
//
// Created by JMO, Thu Dec 17 09:36:30 PST 2015
//----------------------------------------------------------------------------//

#include "GammaPolicy.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Strength/SolidFieldNames.hh"
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
GammaPolicy<Dimension>::
GammaPolicy():
  FieldUpdatePolicy<Dimension, Scalar>({HydroFieldNames::massDensity,
                                        HydroFieldNames::specificThermalEnergy,
                                        SolidFieldNames::porositySolidDensity}) {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GammaPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& /*derivs*/,
       const double /*multiplier*/,
       const double /*t*/,
       const double /*dt*/) {
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == HydroFieldNames::gamma);
  auto& gamma = state.field(key, Scalar());

  // Check if this material has porosity and get the state we need
  const auto  buildKey = [&](const std::string& fkey) { return StateBase<Dimension>::buildFieldKey(fkey, nodeListKey); };
  const auto  usePorosity = state.registered(buildKey(SolidFieldNames::porosityAlpha));
  const auto& rhoS = (usePorosity ?
                      state.field(buildKey(SolidFieldNames::porositySolidDensity), 0.0) :
                      state.field(buildKey(HydroFieldNames::massDensity), 0.0));
  const auto& eps = state.field(buildKey(HydroFieldNames::specificThermalEnergy), 0.0);

  // Get the eos.  This cast is ugly, but is a work-around for now.
  const auto* fluidNodeListPtr = dynamic_cast<const FluidNodeList<Dimension>*>(gamma.nodeListPtr());
  CHECK(fluidNodeListPtr != nullptr);
  const auto& eos = fluidNodeListPtr->equationOfState();

  // Now set the gamma for this field.
  eos.setGammaField(gamma, rhoS, eps);
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
GammaPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const auto rhsPtr = dynamic_cast<const GammaPolicy<Dimension>*>(&rhs);
  return (rhsPtr != nullptr);
}

}

