//---------------------------------Spheral++----------------------------------//
// ArtificialConductionPolicy -- Override the default energy policy in the
// presence of artificial conduction.
//----------------------------------------------------------------------------//
#include "ArtificialConductionPolicy.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/FieldList.hh"
#include "NodeList/FluidNodeList.hh"
#include "Utilities/DBC.hh"

namespace Spheral {
    
//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
ArtificialConductionPolicy<Dimension>::
ArtificialConductionPolicy(typename State<Dimension>::PolicyPointer& energyPolicy):
  FieldUpdatePolicy<Dimension, Scalar>({}),
  mEnergyPolicy(energyPolicy){
  const auto& dependencies = energyPolicy->dependencies();
  for (const auto& depkey: dependencies) this->addDependency(depkey);
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ArtificialConductionPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {
  mEnergyPolicy->update(key, state, derivs, multiplier, t, dt);
  conduct(key,state,derivs,multiplier,t,dt);
}

template<typename Dimension>
void
ArtificialConductionPolicy<Dimension>::
updateAsIncrement(const KeyType& key,
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs,
                  const double multiplier,
                  const double t,
                  const double dt) {
  mEnergyPolicy->updateAsIncrement(key,state,derivs,multiplier,t,dt);
  conduct(key,state,derivs,multiplier,t,dt);
}

template<typename Dimension>
void
ArtificialConductionPolicy<Dimension>::
conduct(const KeyType& key,
        State<Dimension>& state,
        StateDerivatives<Dimension>& derivs,
        const double multiplier,
        const double /*t*/,
        const double /*dt*/) {

  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == HydroFieldNames::specificThermalEnergy);

  auto&       eps = state.field(key, 0.0);
  const auto& DepsDt = derivs.field(State<Dimension>::buildFieldKey("Artificial Cond DepsDt", nodeListKey), 0.0);

  // Loop over the internal values of the field.
  const auto n = eps.numInternalElements();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    eps(i) += DepsDt(i) * multiplier;
  }
}


//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
ArtificialConductionPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const auto* rhsPtr = dynamic_cast<const ArtificialConductionPolicy<Dimension>*>(&rhs);
  return (rhsPtr != nullptr);
}

}

