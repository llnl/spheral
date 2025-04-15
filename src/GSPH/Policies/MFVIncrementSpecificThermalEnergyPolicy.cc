//---------------------------------Spheral++----------------------------------//
// MFVIncrementSpecificThermalEnergyPolicy -- This is a specialized increment
//            policy for the specific thermal energy for schemes that allow
//            for flux between nodes. The specific thermal energy is updated
//            based on the time derivative of thermal energy. The mass and 
//            time derivative are needed to got from thermal to specific 
//            thermal. 
//
// J.M. Pearl 2022
//----------------------------------------------------------------------------//
// TODO: the edge case handing for m->0 needs to be improved to robustly
//       handle void when full Eulerian.
//----------------------------------------------------------------------------//

#include "GSPH/Policies/MFVIncrementSpecificThermalEnergyPolicy.hh"
#include "GSPH/GSPHFieldNames.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/Field.hh"
#include "Utilities/DBC.hh"
#include "Hydro/HydroFieldNames.hh"

#include <limits.h>

namespace Spheral {

//------------------------------------------------------------------------------
// Constructors.
//------------------------------------------------------------------------------
template<typename Dimension>
MFVIncrementSpecificThermalEnergyPolicy<Dimension>::
MFVIncrementSpecificThermalEnergyPolicy(std::initializer_list<std::string> depends):
  FieldUpdatePolicy<Dimension, Scalar>(depends){
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MFVIncrementSpecificThermalEnergyPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {

  const auto tiny = std::numeric_limits<typename Dimension::Scalar>::epsilon();

  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);

  const auto massKey = StateBase<Dimension>::buildFieldKey(HydroFieldNames::mass, nodeListKey);
  const auto derivFieldKey = StateBase<Dimension>::buildFieldKey(prefix() + GSPHFieldNames::thermalEnergy, nodeListKey);

  const auto&  m   = state.field(massKey, Scalar());
        auto&  eps = state.field(key,     Scalar());

  const auto&  DmDt = derivs.field(prefix() + massKey, Scalar());
  const auto&  DmepsDt = derivs.field(derivFieldKey,   Scalar());

  const auto n = m.numInternalElements();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    const auto m1 = m(i)+DmDt(i)*multiplier;
    if (m1 > tiny) eps(i) += (DmepsDt(i) - DmDt(i)*eps(i)) * multiplier * safeInv(m1);
  }

}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
MFVIncrementSpecificThermalEnergyPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an replace operator.
  const auto* rhsPtr = dynamic_cast<const MFVIncrementSpecificThermalEnergyPolicy<Dimension>*>(&rhs);
  return rhsPtr != nullptr;
}

}

