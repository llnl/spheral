//---------------------------------Spheral++----------------------------------//
// SoundSpeedPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the node weight as a dependent quantity.
//
// Created by JMO, Tue Sep 14 22:27:08 2004
//----------------------------------------------------------------------------//

#include "SoundSpeedPolicy.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Strength/SolidFieldNames.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "NodeList/FluidNodeList.hh"
#include "NodeList/SolidNodeList.hh"
#include "Material/EquationOfState.hh"
#include "SolidMaterial/StrengthModel.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
SoundSpeedPolicy<Dimension>::
SoundSpeedPolicy():
  FieldUpdatePolicy<Dimension, Scalar>({HydroFieldNames::massDensity,
                                        HydroFieldNames::specificThermalEnergy,
                                        HydroFieldNames::pressure,
                                        SolidFieldNames::tensorDamage,
                                        SolidFieldNames::porositySolidDensity,
                                        SolidFieldNames::porosityAlpha,
                                        SolidFieldNames::porosityAlpha0,
                                        SolidFieldNames::porosityc0}) {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SoundSpeedPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == HydroFieldNames::soundSpeed);
  auto& cs = state.field(key, 0.0);

  // Get the eos.  This cast is ugly, but is a work-around for now.
  const auto* fluidNodeListPtr = dynamic_cast<const FluidNodeList<Dimension>*>(cs.nodeListPtr());
  CHECK(fluidNodeListPtr != nullptr);
  const auto& eosS = fluidNodeListPtr->equationOfState();

  // Check if there's porosity and get the state we need
  const auto  buildKey = [&](const std::string& fkey) { return StateBase<Dimension>::buildFieldKey(fkey, nodeListKey); };
  const auto  usePorosity = state.registered(buildKey(SolidFieldNames::porosityAlpha));
  const auto& rhoS = (usePorosity ?
                      state.field(buildKey(SolidFieldNames::porositySolidDensity), 0.0) :
                      state.field(buildKey(HydroFieldNames::massDensity), 0.0));
  const auto& eps = state.field(buildKey(HydroFieldNames::specificThermalEnergy), 0.0);

  // Set the starting solid sound speed
  eosS.setSoundSpeed(cs, rhoS, eps);

  // Augment with the strength model if appropriate
  const auto* solidNodeListPtr = dynamic_cast<const SolidNodeList<Dimension>*>(fluidNodeListPtr);
  if (solidNodeListPtr != nullptr) {
    const auto& strengthModelS = solidNodeListPtr->strengthModel();
    if (strengthModelS.providesSoundSpeed()) {
      const auto& P = state.field(buildKey(HydroFieldNames::pressure), 0.0);
      const auto& D = state.field(buildKey(SolidFieldNames::tensorDamage), SymTensor::zero);
      strengthModelS.soundSpeed(cs, rhoS, eps, P, cs, D);
    }
  }

  // Correct for porosity if present
  // Assume the sound speed varies linearly from the initial porous value to
  // the current solid value with distension.
  if (usePorosity) {
    const auto& alpha = state.field(buildKey(SolidFieldNames::porosityAlpha), 0.0);
    const auto& alpha0 = state.field(buildKey(SolidFieldNames::porosityAlpha0), 0.0);
    const auto& cs0 = state.field(buildKey(SolidFieldNames::porosityc0), 0.0);
    const auto n = cs.numInternalElements();
#pragma omp parallel for
    for (auto i = 0u; i < n; ++i) {
      const auto alpha0i = alpha0(i);
      const auto alphai = alpha(i);
      const auto cs0i = cs0(i);
      CHECK(alpha0i >= 1.0 and alphai >= 1.0);
      CHECK(cs0i > 0.0);
      cs(i) += (std::min(alphai, alpha0i) - 1.0)*safeInv(alpha0i - 1.0)*(cs0i - cs(i));
      ENSURE2(cs(i) >= 0.0, "Bad porous sound speed for " << i << " : " << cs(i));
    }
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
SoundSpeedPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {
  return dynamic_cast<const SoundSpeedPolicy<Dimension>*>(&rhs) != nullptr;
}

}

