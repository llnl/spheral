//---------------------------------Spheral++----------------------------------//
// YieldStrengthPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the dependent yield strength state.
//
// Created by JMO, Thu Sep 30 23:01:00 PDT 2004
//----------------------------------------------------------------------------//
#include "YieldStrengthPolicy.hh"
#include "NodeList/SolidNodeList.hh"
#include "SolidMaterial/StrengthModel.hh"
#include "Strength/SolidFieldNames.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/UpdatePolicyBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/IncrementState.hh"
#include "Field/FieldList.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
YieldStrengthPolicy<Dimension>::
YieldStrengthPolicy(const bool scaleWithPorosity):
  FieldUpdatePolicy<Dimension, Scalar>({HydroFieldNames::massDensity,
                                        HydroFieldNames::specificThermalEnergy,
                                        HydroFieldNames::pressure,
                                        SolidFieldNames::plasticStrain,
                                        SolidFieldNames::tensorDamage,
                                        IncrementState<Dimension, Scalar>::prefix() + SolidFieldNames::plasticStrain,
                                        SolidFieldNames::porositySolidDensity,
                                        SolidFieldNames::porosityAlpha}),
  mScaleWithPorosity(scaleWithPorosity) {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
YieldStrengthPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double /*multiplier*/,
       const double /*t*/,
       const double /*dt*/) {

  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == SolidFieldNames::yieldStrength);
  auto& Y = state.field(key, 0.0);

  // Get the strength model.  This cast is ugly, but is a work-around for now.
  const auto* solidNodeListPtr = dynamic_cast<const SolidNodeList<Dimension>*>(Y.nodeListPtr());
  CHECK(solidNodeListPtr != nullptr);
  const auto& strengthModel = solidNodeListPtr->strengthModel();

  // Is there a porosity running around?
  const auto buildKey = [&](const std::string& fkey) { return StateBase<Dimension>::buildFieldKey(fkey, nodeListKey); };
  const auto usePorosity = state.registered(buildKey(SolidFieldNames::porosityAlpha));

  // Get (most of) the state we need
  const auto& eps = state.field(buildKey(HydroFieldNames::specificThermalEnergy), 0.0);
  const auto& P = state.field(buildKey(HydroFieldNames::pressure), 0.0);
  const auto& plasticStrain = state.field(buildKey(SolidFieldNames::plasticStrain), 0.0);
  const auto& plasticStrainRate = derivs.field(buildKey(SolidFieldNames::plasticStrainRate), 0.0);
  const auto& D = state.field(buildKey(SolidFieldNames::tensorDamage), SymTensor::zero);

  // Things change depending on porosity...
  if (usePorosity) {

    // Set the solid phase yield strength
    const auto& rhoS = state.field(buildKey(SolidFieldNames::porositySolidDensity), 0.0);
    const auto& alpha = state.field(buildKey(SolidFieldNames::porosityAlpha), 0.0);
    const auto PS = P*alpha;
    strengthModel.yieldStrength(Y, rhoS, eps, PS, plasticStrain, plasticStrainRate, D);

    // Optionally scale the result by the distention for backwards compatibility
    if (mScaleWithPorosity) Y /= alpha;

  } else {

    // Set the yield strength.
    const auto& rho = state.field(buildKey(HydroFieldNames::massDensity), 0.0);
    const auto& P = state.field(buildKey(HydroFieldNames::pressure), 0.0);
    strengthModel.yieldStrength(Y, rho, eps, P, plasticStrain, plasticStrainRate, D);

  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
YieldStrengthPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const auto* rhsPtr = dynamic_cast<const YieldStrengthPolicy<Dimension>*>(&rhs);
  return (rhsPtr != nullptr);
}

}

