//---------------------------------Spheral++----------------------------------//
// AxisymmetricMassDensityPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the 3D mass density in RZ calculations.
//
// Created by JMO, Thu May  7 15:20:58 PDT 2026
//----------------------------------------------------------------------------//

#include "AxisymmetricMassDensityPolicy.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "SPH/SPHRZUtilities.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
AxisymmetricMassDensityPolicy::
AxisymmetricMassDensityPolicy(const Scalar rhoMin,
                              const Scalar rhoMax):
  FieldUpdatePolicy<Dimension, Scalar>({HydroFieldNames::massDensityRZ,
                                        HydroFieldNames::specificThermalEnergy}),
  mRhoMin(rhoMin),
  mRhoMax(rhoMax) {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
void
AxisymmetricMassDensityPolicy::
update(const KeyType& key,
       State<Dim<2>>& state,
       StateDerivatives<Dim<2>>& derivs,
       const double multiplier,
       const double t,
       const double dt) {
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == HydroFieldNames::massDensity);

  // Grab some of our required state
  const auto   buildKey = [&](const std::string& fkey) { return StateBase<Dimension>::buildFieldKey(fkey, nodeListKey); };
  auto&        massDensity = state.field(buildKey(HydroFieldNames::massDensity), Scalar());
  const auto&  massDensityRZ = state.field(buildKey(HydroFieldNames::massDensityRZ), Scalar());
  const auto&  mass = state.field(buildKey(HydroFieldNames::mass), Scalar());
  const auto&  massRZ = state.field(buildKey(HydroFieldNames::massRZ), Scalar());
  const auto&  position = state.field(buildKey(HydroFieldNames::position), Vector::zero());
  const auto&  H = state.field(buildKey(HydroFieldNames::H), SymTensor::zero());

  // Set the mass density based on the RZ mass density
  SPHRZUtilities::computeMassDensityFromArealMassDensity(massDensity, massDensityRZ, mass, massRZ, position, H, mRhoMin, mRhoMax);
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
bool
AxisymmetricMassDensityPolicy::
operator==(const UpdatePolicyBase<Dim<2>>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const auto rhsPtr = dynamic_cast<const AxisymmetricMassDensityPolicy*>(&rhs);
  return (rhsPtr != nullptr);
}

}

