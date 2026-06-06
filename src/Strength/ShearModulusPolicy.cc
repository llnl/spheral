//---------------------------------Spheral++----------------------------------//
// ShearModulusPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the dependent shear modulus state.
//
// Created by JMO, Thu Sep 30 23:01:00 PDT 2004
//----------------------------------------------------------------------------//
#include "ShearModulusPolicy.hh"
#include "NodeList/SolidNodeList.hh"
#include "SolidMaterial/StrengthModel.hh"
#include "SolidFieldNames.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Strength/SolidFieldNames.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/FieldList.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
ShearModulusPolicy<Dimension>::
ShearModulusPolicy():
  FieldUpdatePolicy<Dimension, Scalar>({HydroFieldNames::massDensity,
                                        HydroFieldNames::specificThermalEnergy,
                                        HydroFieldNames::pressure,
                                        SolidFieldNames::tensorDamage,
                                        SolidFieldNames::porositySolidDensity,
                                        SolidFieldNames::porosityAlpha}) {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ShearModulusPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& /*derivs*/,
       const double /*multiplier*/,
       const double /*t*/,
       const double /*dt*/) {
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == SolidFieldNames::shearModulus);
  auto& mu = state.field(key, 0.0);

  // Get the strength model.  This cast is ugly, but is a work-around for now.
  const auto* solidNodeListPtr = dynamic_cast<const SolidNodeList<Dimension>*>(mu.nodeListPtr());
  VERIFY(solidNodeListPtr != nullptr);
  const auto& strengthModel = solidNodeListPtr->strengthModel();

  // Check if we're using porosity
  const auto buildKey = [&](const std::string& fkey) { return StateBase<Dimension>::buildFieldKey(fkey, nodeListKey); };
  const auto usePorosity = state.registered(buildKey(SolidFieldNames::porosityAlpha));

  // Get the state fields we need
  const auto& rhoS = (usePorosity ?
                      state.field(buildKey(SolidFieldNames::porositySolidDensity), 0.0) :
                      state.field(buildKey(HydroFieldNames::massDensity), 0.0));
  const auto& eps = state.field(buildKey(HydroFieldNames::specificThermalEnergy), 0.0);
  const auto& P = state.field(buildKey(HydroFieldNames::pressure), 0.0);
  const auto& D = state.field(buildKey(SolidFieldNames::tensorDamage), SymTensor::zero());

  // Setting the shear modulus differs based on the presence of porosity
  if (usePorosity) {

    // We actually need the solid phase pressure
    const auto& alpha = state.field(buildKey(SolidFieldNames::porosityAlpha), 0.0);
    const auto  PS = P*alpha;

    // Set the solid phase shear modulus
    strengthModel.shearModulus(mu, rhoS, eps, PS, D);

    // Scale the result by the distention
    mu /= alpha;

  } else {

    strengthModel.shearModulus(mu, rhoS, eps, P, D);

  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
ShearModulusPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {
  return dynamic_cast<const ShearModulusPolicy<Dimension>*>(&rhs) != nullptr;
}

}

