//---------------------------------Spheral++----------------------------------//
// MeshIdealHPolicy
//
// Created by JMO, Fri Aug 16 14:47:57 PDT 2013
//----------------------------------------------------------------------------//
#include <vector>

#include "MeshIdealHPolicy.hh"
#include "Hydro/HydroFieldNames.hh"
#include "NodeList/NodeList.hh"
#include "NodeList/FluidNodeList.hh"
#include "Kernel/TableKernel.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Utilities/DBC.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/SpheralFunctions.hh"

namespace Spheral {


//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
MeshIdealHPolicy<Dimension>::
MeshIdealHPolicy(const SmoothingScaleBase<Dimension>& smoothingScaleBase,
                 const Scalar hmin,
                 const Scalar hmax,
                 const Scalar hminratio,
                 const Scalar nPerh):
  UpdatePolicyBase<Dimension>({HydroFieldNames::mesh}),
  mSmoothingScaleBase(smoothingScaleBase),
  mhmin(hmin),
  mhmax(hmax),
  mhminratio(hminratio),
  mnPerh(nPerh) {
  REQUIRE(hmin <= hmax);
  mFired = false;
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
MeshIdealHPolicy<Dimension>::
~MeshIdealHPolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MeshIdealHPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& /*derivs*/,
       const double /*multiplier*/,
       const double /*t*/,
       const double /*dt*/) {
  using Zone = typename Mesh<Dimension>::Zone;

  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == HydroFieldNames::H);

  if (not mFired) {
    mFired = true;

    // Get the state fields.
    FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero());
    const Mesh<Dimension>& mesh = state.mesh();
    const size_t numNodeLists = H.numFields();

    // Walk the NodeLists.
    for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
      const unsigned n = H[nodeListi]->numInternalElements();
      for (unsigned i = 0; i != n; ++i) {
        const Zone& zonei = mesh.zone(nodeListi, i);
        H(nodeListi, i) = mSmoothingScaleBase.idealSmoothingScale(H(nodeListi, i),
                                                                  mesh,
                                                                  zonei,
                                                                  mhmin,
                                                                  mhmax,
                                                                  mhminratio,
                                                                  mnPerh);
      }
    }

    // // Now go through and filter the H tensor locally.
    // FieldList<Dimension, SymTensor> Havg = derivs.fields(ReplaceBoundedState<Dimension, SymTensor, Scalar>::prefix() + HydroFieldNames::H, SymTensor::zero());
    // for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    //   const unsigned n = H[nodeListi]->numInternalElements();
    //   for (unsigned i = 0; i != n; ++i) {
    //     const Zone& zone = mesh.zone(nodeListi, i);
    //     Havg(nodeListi, i).Zero();

    //     // Build the unique set of neighboring zones by nodes.
    //     vector<unsigned> otherZones;
    //     {
    //       const vector<unsigned>& nodeIDs = zone.nodeIDs();
    //       for (vector<unsigned>::const_iterator itr = nodeIDs.begin();
    //            itr != nodeIDs.end();
    //            ++itr) {
    //         const vector<unsigned>& nodeZoneIDs = mesh.node(*itr).zoneIDs();
    //         copy(nodeZoneIDs.begin(), nodeZoneIDs.end(), back_inserter(otherZones));
    //       }
    //       sort(otherZones.begin(), otherZones.end());
    //       otherZones.erase(unique(otherZones.begin(), otherZones.end()), otherZones.end());
    //     }
    //     CHECK(otherZones.size() > 0);

    //     // Build our averaged H tensor.
    //     for (vector<unsigned>::const_iterator itr = otherZones.begin();
    //          itr != otherZones.end();
    //          ++itr) {
    //       unsigned nodeListj = nodeListi, j = i;
    //       if (*itr != Mesh<Dimension>::UNSETID) {
    //         mesh.lookupNodeListID(*itr, nodeListj, j);
    //       }
    //       Havg(nodeListi, i) += H(nodeListj, j).Inverse();
    //     }
    //     Havg(nodeListi, i) = (Havg(nodeListi, i)/otherZones.size()).Inverse();
    //     cerr << " --> " << i << " " << H(nodeListi,i).Inverse() << " " << Havg(nodeListi, i).Inverse() << endl;
    //   }
    // }
    // H.assignFields(Havg);

  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
MeshIdealHPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also a CompatibleSpecificThermalEnergy operator.
  const auto* rhsPtr = dynamic_cast<const MeshIdealHPolicy<Dimension>*>(&rhs);
  return (rhsPtr != nullptr);
}

//------------------------------------------------------------------------------
// Static initializations.
//------------------------------------------------------------------------------
template<typename Dimension> bool MeshIdealHPolicy<Dimension>::mFired = false;

}

