//---------------------------------Spheral++----------------------------------//
// SpecificThermalEnergyVolumePolicy
//
// Created by JMO, Sun Aug 18 19:21:27 PDT 2013
//----------------------------------------------------------------------------//

#include "SpecificThermalEnergyVolumePolicy.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/Field.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructors.
//------------------------------------------------------------------------------
template<typename Dimension>
SpecificThermalEnergyVolumePolicy<Dimension>::
SpecificThermalEnergyVolumePolicy():
  FieldUpdatePolicy<Dimension, Scalar>({HydroFieldNames::volume}) {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SpecificThermalEnergyVolumePolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double /*t*/,
       const double /*dt*/) {

  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == HydroFieldNames::specificThermalEnergy);

  // Get the state fields.
  const auto  buildKey = [&](const std::string& fkey) { return StateBase<Dimension>::buildFieldKey(fkey, nodeListKey); };
  auto&       eps = state.field(key, 0.0);
  const auto& vol = state.field(buildKey(HydroFieldNames::volume), 0.0);
  const auto& vol0 = state.field(buildKey(HydroFieldNames::volume + "0"), 0.0);
  const auto& P = state.field(buildKey(HydroFieldNames::pressure), 0.0);
  const auto& Q = derivs.field(buildKey(HydroFieldNames::maxViscousPressure), 0.0);

  // Do the deed.
  const auto n = eps.numInternalElements();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    CHECK(vol0(i) > 0.0);
    eps(i) += multiplier*((P(i) + Q(i))*(vol0(i) - vol(i)));
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
SpecificThermalEnergyVolumePolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const auto* rhsPtr = dynamic_cast<const SpecificThermalEnergyVolumePolicy<Dimension>*>(&rhs);
  return rhsPtr != nullptr;
}

}

