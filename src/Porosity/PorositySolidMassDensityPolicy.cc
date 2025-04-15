//---------------------------------Spheral++----------------------------------//
// PorositySolidMassDensityPolicy -- An implementation of UpdatePolicyBase
// specialized for the updating the dependent solid mass density in the presence
// of porosity.
//
// Created by JMO, Wed Oct  4 13:51:43 PDT 2023
//----------------------------------------------------------------------------//

#include "PorositySolidMassDensityPolicy.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Strength/SolidFieldNames.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
PorositySolidMassDensityPolicy<Dimension>::
PorositySolidMassDensityPolicy():
  FieldUpdatePolicy<Dimension, Scalar>({HydroFieldNames::massDensity,
                                        SolidFieldNames::porosityAlpha}) {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PorositySolidMassDensityPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == SolidFieldNames::porositySolidDensity);
  auto& rhoS = state.field(key, Scalar());

  // Get the mass density and alpha fields from the state.
  const auto& rho = state.field(State<Dimension>::buildFieldKey(HydroFieldNames::massDensity, nodeListKey), Scalar());
  const auto& alpha = state.field(State<Dimension>::buildFieldKey(SolidFieldNames::porosityAlpha, nodeListKey), Scalar());

  // Now set the solid density
  const auto n = rhoS.numInternalElements();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    rhoS(i) = alpha(i)*rho(i);
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
PorositySolidMassDensityPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {
  return dynamic_cast<const PorositySolidMassDensityPolicy<Dimension>*>(&rhs) != nullptr;
}

}

