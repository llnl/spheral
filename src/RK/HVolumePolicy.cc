//---------------------------------Spheral++----------------------------------//
// HVolumePolicy -- An implementation of ReplaceState specialized
// for the updating the volume based on the local hull constructions.
//
// Created by JMO, Wed Aug 13 10:52:16 PDT 2014
//----------------------------------------------------------------------------//

#include "HVolumePolicy.hh"
#include "computeHVolumes.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/FieldList.hh"
#include "Mesh/Mesh.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
HVolumePolicy<Dimension>::
HVolumePolicy(const Scalar kernelExtent):
  UpdatePolicyBase<Dimension>({HydroFieldNames::H}),
  mKernelExtent(kernelExtent) {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
HVolumePolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& /*derivs*/,
       const double /*multiplier*/,
       const double /*t*/,
       const double /*dt*/) {

  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == HydroFieldNames::volume and 
          nodeListKey == UpdatePolicyBase<Dimension>::wildcard());
  FieldList<Dimension, Scalar> volume = state.fields(fieldKey, 0.0);

  // Get the H field from the state, and do the deed.
  const auto H = state.fields(HydroFieldNames::H, SymTensor::zero());
  computeHVolumes(mKernelExtent, H, volume);
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
HVolumePolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {
  return dynamic_cast<const HVolumePolicy<Dimension>*>(&rhs) != nullptr;
}

}

