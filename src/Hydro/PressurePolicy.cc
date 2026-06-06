//---------------------------------Spheral++----------------------------------//
// PressurePolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the node weight as a dependent quantity.
//
// Created by JMO, Tue Sep 14 22:27:08 2004
//----------------------------------------------------------------------------//

#include "PressurePolicy.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Strength/SolidFieldNames.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/Field.hh"
#include "NodeList/FluidNodeList.hh"
#include "Material/EquationOfState.hh"
#include "SolidMaterial/SolidEquationOfState.hh"
#include "Utilities/DBC.hh"

namespace Spheral {


//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
PressurePolicy<Dimension>::
PressurePolicy():
  FieldUpdatePolicy<Dimension, Scalar>({HydroFieldNames::massDensity,
                                        HydroFieldNames::specificThermalEnergy,
                                        SolidFieldNames::porositySolidDensity,
                                        SolidFieldNames::porosityAlpha,
                                        SolidFieldNames::tensorDamage}) {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PressurePolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE((fieldKey == HydroFieldNames::pressure or 
           fieldKey == SolidFieldNames::damagedPressure));
  auto& P = state.field(key, Scalar());

  // Get the eos.  This cast is ugly, but is a work-around for now.
  const auto* fluidNodeListPtr = dynamic_cast<const FluidNodeList<Dimension>*>(P.nodeListPtr());
  CHECK(fluidNodeListPtr != nullptr);
  const auto& eos = fluidNodeListPtr->equationOfState();

  // Grab some of our required state
  const auto  buildKey = [&](const std::string& fkey) { return StateBase<Dimension>::buildFieldKey(fkey, nodeListKey); };
  const auto& eps = state.field(buildKey(HydroFieldNames::specificThermalEnergy), 0.0);

  // Check if this material has porosity.
  const auto usePorosity = state.registered(buildKey(SolidFieldNames::porosityAlpha));
  if (usePorosity) {

    // Yep, there's porosity
    const auto& rhoS = state.field(buildKey(SolidFieldNames::porositySolidDensity), 0.0);
      
    // Check if we need to update the pressure derivatives by seeing if they're registered for this NodeList
    const auto dPduKey = State<Dimension>::buildFieldKey(HydroFieldNames::partialPpartialEps, nodeListKey);
    const auto dPdrKey = State<Dimension>::buildFieldKey(HydroFieldNames::partialPpartialRho, nodeListKey);
    if (state.registered(dPduKey)) {
      CHECK(state.registered(dPdrKey));
      auto& dPdu = state.field(dPduKey, 0.0);
      auto& dPdr = state.field(dPdrKey, 0.0);
      eos.setPressureAndDerivs(P, dPdu, dPdr, rhoS, eps);
    } else {
      eos.setPressure(P, rhoS, eps);
    }

    // Note we don't scale back to the bulk pressure until we check for any damage effects,
    // which only affect the solid material (not void).
  } else {
      
    // No porosity, so just set the pressure using the straight density
    const auto& rho = state.field(buildKey(HydroFieldNames::massDensity), 0.0);
    eos.setPressure(P, rho, eps);

  }

  // Is someone trying to keep the damaged pressure in an independent Field?
  // (I'm looking at you FSISPH)
  const auto separateDamage = state.registered(buildKey(SolidFieldNames::damagedPressure));
  Field<Dimension, Scalar>* PdPtr = nullptr;
  if (separateDamage) PdPtr = &state.field(buildKey(SolidFieldNames::damagedPressure), 0.0);

  // If there's damage for this material, apply it to the pressure
  // This is complicated by FSISPH, which wants to keep track of the damaged pressure separately,
  // so we check if someone has registered damaged pressure as a field for this NodeList, in
  // which case we apply the damage to a new copy of the pressure.
  if (state.registered(buildKey(SolidFieldNames::tensorDamage))) {
    const auto& D = state.field(buildKey(SolidFieldNames::tensorDamage), SymTensor::zero());

    // Check for the appropriate minimum pressure
    const auto* solidEOSptr = dynamic_cast<const SolidEquationOfState<Dimension>*>(&eos);
    const auto Pmin = eos.minimumPressure();
    const auto PminDamage = (solidEOSptr != nullptr ?
                             solidEOSptr->minimumPressureDamage() :
                             Pmin);
                       
    // Scale by the damage.
    const auto ni = P.numInternalElements();
#pragma omp parallel for
    for (auto i = 0u; i < ni; ++i) {
      const auto Di = std::max(0.0, std::min(1.0, D(i).eigenValues().maxElement()));
      CHECK(Di >= 0.0 and Di <= 1.0);
      const auto Pdi = std::max(Pmin, (1.0 - Di)*P(i)) + std::max(PminDamage, Di*P(i));
      if (separateDamage) {
        (*PdPtr)(i) = Pdi;
      } else {
        P(i) = Pdi;
      }
    }
  }

  // Scale the solid pressure to the bulk value in the case of porosity
  if (usePorosity) {
    const auto& alpha = state.field(buildKey(SolidFieldNames::porosityAlpha), 0.0);
    P /= alpha;
    if (separateDamage) (*PdPtr) /= alpha;
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
PressurePolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const auto rhsPtr = dynamic_cast<const PressurePolicy<Dimension>*>(&rhs);
  return (rhsPtr != nullptr);
}

}

