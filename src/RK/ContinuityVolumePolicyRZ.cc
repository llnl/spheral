//---------------------------------Spheral++----------------------------------//
// ContinuityVolumePolicy -- An implementation of IncrementFieldList
// specialized for time evolving the volume per point using the continuity
// equation.
//
// Specialized for RZ geometry.
//
// Created by JMO, Tue Oct 29 16:14:03 PDT 2019
//----------------------------------------------------------------------------//

#include "ContinuityVolumePolicyRZ.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/UpdatePolicyBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Utilities/DBC.hh"

namespace Spheral {


namespace {
//------------------------------------------------------------------------------
// Compute the dimension dependent volume of the H tensor.
//------------------------------------------------------------------------------
double Hvolume(const Dim<2>::SymTensor& H) {
  return M_PI/H.Determinant();
}
}

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
ContinuityVolumePolicyRZ::
ContinuityVolumePolicyRZ():
  IncrementState<Dim<2>, Dim<2>::Scalar>({HydroFieldNames::mass,
                                          HydroFieldNames::massDensity}) {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
void
ContinuityVolumePolicyRZ::
update(const KeyType& key,
       State<Dim<2>>& state,
       StateDerivatives<Dim<2>>& derivs,
       const double multiplier,
       const double t,
       const double dt) {

  CONTRACT_VAR(t);
  CONTRACT_VAR(dt);

  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == HydroFieldNames::volume);
  auto& volume = state.field(key, 0.0);

  // Get the state we depend on
  const auto  buildKey = [&](const std::string& fkey) { return StateBase<Dimension>::buildFieldKey(fkey, nodeListKey); };
  const auto& pos = state.field(buildKey(HydroFieldNames::position), Vector::zero());
  const auto& mass = state.field(buildKey(HydroFieldNames::mass), 0.0);
  const auto& rho = state.field(buildKey(HydroFieldNames::massDensity), 0.0);
  const auto& H = state.field(buildKey(HydroFieldNames::H), SymTensor::zero());
  const auto& DrhoDt = derivs.field(buildKey(IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity), 0.0);

  // Loop over the internal values of the field.
  const auto n = volume.numInternalElements();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    const auto circi = 2.0*M_PI*abs(pos(i).y());
    CHECK(circi > 0.0);
    const auto volMin = 0.5*mass(i)*safeInvVar(circi*rho(i));
    const auto volMax = Hvolume(H(i));
    const auto dVdt = -mass(i)*safeInvVar(circi*rho(i)*rho(i))*DrhoDt(i);
    volume(i) = std::max(volMin, std::min(volMax, volume(i) + multiplier*dVdt));
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
bool
ContinuityVolumePolicyRZ::
operator==(const UpdatePolicyBase<Dim<2>>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const auto* rhsPtr = dynamic_cast<const ContinuityVolumePolicyRZ*>(&rhs);
  return (rhsPtr != nullptr);
}

}

