//---------------------------------Spheral++----------------------------------//
// HullVolumePolicy -- An implementation of ReplaceState specialized
// for the updating the volume based on the local hull constructions.
//
// Created by JMO, Wed Jul  2 22:35:26 PDT 2014
//----------------------------------------------------------------------------//

#include "HullVolumePolicy.hh"
#include "computeHullVolumes.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/UpdatePolicyBase.hh"
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
HullVolumePolicy<Dimension>::
HullVolumePolicy(const ConnectivityMap<Dimension>& connectivityMap):
  mConnectivityMap(connectivityMap),
  ReplaceFieldList<Dimension, typename Dimension::Scalar>(HydroFieldNames::position) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
HullVolumePolicy<Dimension>::
~HullVolumePolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
HullVolumePolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {

  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == HydroFieldNames::volume and 
          nodeListKey == UpdatePolicyBase<Dimension>::wildcard());
  FieldList<Dimension, Scalar> volume = state.fields(fieldKey, Scalar());

  // Get the position from the state, and do the deed.
  const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero());
  computeHullVolumes(mConnectivityMap, position, volume);
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
HullVolumePolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const HullVolumePolicy<Dimension>* rhsPtr = dynamic_cast<const HullVolumePolicy<Dimension>*>(&rhs);
  if (rhsPtr == 0) {
    return false;
  } else {
    return true;
  }
}

}

