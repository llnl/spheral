//---------------------------------Spheral++----------------------------------//
// ContinuityVolumePolicy -- An implementation of IncrementFieldList
// specialized for time evolving the volume per point using the continuity
// equation.
//
// Created by JMO, Tue Sep 20 14:53:32 PDT 2016
//----------------------------------------------------------------------------//

#include "ContinuityVolumePolicy.hh"
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
static inline double Hvolume(const Dim<1>::SymTensor& H) {
  return 2.0/H.xx();
}

static inline double Hvolume(const Dim<2>::SymTensor& H) {
  return M_PI/H.Determinant();
}

static inline double Hvolume(const Dim<3>::SymTensor& H) {
  return 4.0*M_PI/(3.0*H.Determinant());
}

}

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
ContinuityVolumePolicy<Dimension>::
ContinuityVolumePolicy():
  IncrementState<Dimension, typename Dimension::Scalar>({HydroFieldNames::mass,
                                                         HydroFieldNames::massDensity}) {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ContinuityVolumePolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double /*t*/,
       const double /*dt*/) {

  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == HydroFieldNames::volume);
  auto& volume = state.field(key, 0.0);

  // Get the state we depend on
  const auto  buildKey = [&](const std::string& fkey) { return StateBase<Dimension>::buildFieldKey(fkey, nodeListKey); };
  const auto& mass = state.field(buildKey(HydroFieldNames::mass), 0.0);
  const auto& rho = state.field(buildKey(HydroFieldNames::massDensity), 0.0);
  const auto& H = state.field(buildKey(HydroFieldNames::H), SymTensor::zero());
  const auto& DrhoDt = derivs.field(buildKey(IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity), 0.0);

  // Loop over the internal values of the field.
  const auto n = volume.numInternalElements();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    const auto volMin = 0.5*mass(i)*safeInvVar(rho(i));
    const auto volMax = Hvolume(H(i));
    const auto dVdt = -mass(i)*safeInvVar(rho(i)*rho(i))*DrhoDt(i);
    volume(i) = std::max(volMin, std::min(volMax, volume(i) + multiplier*dVdt));
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
ContinuityVolumePolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {
  return dynamic_cast<const ContinuityVolumePolicy<Dimension>*>(&rhs) != nullptr;
}

}

