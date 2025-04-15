//---------------------------------Spheral++----------------------------------//
// CellPressurePolicy
//
// Created by JMO, Sun Aug 25 14:40:42 PDT 2013
//----------------------------------------------------------------------------//
#include "CellPressurePolicy.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/Field.hh"
#include "NodeList/FluidNodeList.hh"
#include "Material/EquationOfState.hh"
#include "Utilities/DBC.hh"

namespace Spheral {


//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
CellPressurePolicy<Dimension>::
CellPressurePolicy():
  FieldUpdatePolicy<Dimension, Scalar>({HydroFieldNames::mass,
                                        HydroFieldNames::volume,
                                        HydroFieldNames::specificThermalEnergy}) {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CellPressurePolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& /*derivs*/,
       const double /*multiplier*/,
       const double /*t*/,
       const double /*dt*/) {

  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == "Cell" + HydroFieldNames::pressure);

  // Get the mass, volume, and specific thermal energy fields from the state.
  const KeyType massKey = StateBase<Dimension>::buildFieldKey(HydroFieldNames::mass, nodeListKey);
  const KeyType volKey = StateBase<Dimension>::buildFieldKey(HydroFieldNames::volume, nodeListKey);
  const KeyType energyKey = StateBase<Dimension>::buildFieldKey(HydroFieldNames::specificThermalEnergy, nodeListKey);
  CHECK(state.registered(massKey));
  CHECK(state.registered(volKey));
  CHECK(state.registered(energyKey));
  Field<Dimension, Scalar>& pressure = state.field(key, Scalar());
  const Field<Dimension, Scalar>& mass = state.field(massKey, Scalar());
  const Field<Dimension, Scalar>& vol = state.field(volKey, Scalar());
  const Field<Dimension, Scalar>& energy = state.field(energyKey, Scalar());

  // Get the eos.  This cast is ugly, but is a work-around for now.
  const FluidNodeList<Dimension>* fluidNodeListPtr = dynamic_cast<const FluidNodeList<Dimension>*>(pressure.nodeListPtr());
  CHECK(fluidNodeListPtr != 0);
  const EquationOfState<Dimension>& eos = fluidNodeListPtr->equationOfState();

  // Compute the cell density.
  const Field<Dimension, Scalar> rho = mass/vol;

  // Now set the pressure.
  eos.setPressure(pressure, rho, energy);
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
CellPressurePolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const auto* rhsPtr = dynamic_cast<const CellPressurePolicy<Dimension>*>(&rhs);
  return (rhsPtr != nullptr);
}

}

