//---------------------------------Spheral++----------------------------------//
// BulkModulusPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the dependent bulk modulus state.
//
// Created by JMO, Tue Oct 5 16:40:54 2004
//----------------------------------------------------------------------------//
#include "BulkModulusPolicy.hh"
#include "SolidFieldNames.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/FieldList.hh"
#include "Material/EquationOfState.hh"
#include "SolidMaterial/StrengthModel.hh"
#include "NodeList/SolidNodeList.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
BulkModulusPolicy<Dimension>::
BulkModulusPolicy():
  FieldUpdatePolicy<Dimension, Scalar>({HydroFieldNames::massDensity,
                                        HydroFieldNames::specificThermalEnergy,
                                        SolidFieldNames::porositySolidDensity}) {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
BulkModulusPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& /*derivs*/,
       const double /*multiplier*/,
       const double /*t*/,
       const double /*dt*/) {
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == SolidFieldNames::bulkModulus);
  auto& K = state.field(key, Scalar());

  // Check if we have a FluidNodeList or SolidNodeList.  Has to be at least fluid!
  const auto* fluidNodeListPtr = dynamic_cast<const FluidNodeList<Dimension>*>(K.nodeListPtr());
  const auto* solidNodeListPtr = dynamic_cast<const SolidNodeList<Dimension>*>(K.nodeListPtr());
  CHECK(fluidNodeListPtr != nullptr);

  // Check if this material has porosity and get the state we need
  const auto  buildKey = [&](const std::string& fkey) { return StateBase<Dimension>::buildFieldKey(fkey, nodeListKey); };
  const auto  usePorosity = state.registered(buildKey(SolidFieldNames::porosityAlpha));
  const auto& rhoS = (usePorosity ?
                      state.field(buildKey(SolidFieldNames::porositySolidDensity), 0.0) :
                      state.field(buildKey(HydroFieldNames::massDensity), 0.0));
  const auto& eps = state.field(buildKey(HydroFieldNames::specificThermalEnergy), 0.0);

  // Is there a strength model that wants to set the bulk modulus?
  if (solidNodeListPtr != nullptr and
      solidNodeListPtr->strengthModel().providesBulkModulus()) {
    solidNodeListPtr->strengthModel().bulkModulus(K, rhoS, eps);
  } else {
    fluidNodeListPtr->equationOfState().setBulkModulus(K, rhoS, eps);
  }

  // If there's porosity, scale K to the bulk value
  if (usePorosity) {
    const auto& alpha = state.field(buildKey(SolidFieldNames::porosityAlpha), 0.0);
    K /= alpha;
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
BulkModulusPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also a bulk modulus operator.
  return dynamic_cast<const BulkModulusPolicy<Dimension>*>(&rhs) != nullptr;
}

}

