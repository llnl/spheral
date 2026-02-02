//---------------------------------Spheral++----------------------------------//
// MeshPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the Mesh in the state.
//
// Created by JMO, Sat Feb 12 14:37:57 PST 2011
//----------------------------------------------------------------------------//
#include "MeshPolicy.hh"
#include "generateMesh.hh"
#include "Physics/Physics.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Utilities/globalBoundingVolumes.hh"
#include "Utilities/DBC.hh"

#include <vector>
using std::vector;
using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor without specifying bounds.
//------------------------------------------------------------------------------
template<typename Dimension>
MeshPolicy<Dimension>::
MeshPolicy(const Physics<Dimension>& package,
           const double voidThreshold,
           const bool meshGhostNodes,
           const bool generateVoid,
           const bool removeBoundaryZones):
  UpdatePolicyBase<Dimension>({HydroFieldNames::position + 
                               UpdatePolicyBase<Dimension>::wildcard()}),
  mPackage(package),
  mVoidThreshold(voidThreshold),
  mComputeBounds(true),
  mMeshGhostNodes(meshGhostNodes),
  mGenerateVoid(generateVoid),
  mRemoveBoundaryZones(removeBoundaryZones),
  mXmin(),
  mXmax() {
}

//------------------------------------------------------------------------------
// Constructor where we specify bounds.
//------------------------------------------------------------------------------
template<typename Dimension>
MeshPolicy<Dimension>::
MeshPolicy(const Physics<Dimension>& package,
           const Vector& xmin,
           const Vector& xmax,
           const double voidThreshold,
           const bool meshGhostNodes,
           const bool generateVoid,
           const bool removeBoundaryZones):
  UpdatePolicyBase<Dimension>({HydroFieldNames::position + 
                               UpdatePolicyBase<Dimension>::wildcard()}),
  mPackage(package),
  mVoidThreshold(voidThreshold),
  mComputeBounds(false),
  mMeshGhostNodes(meshGhostNodes),
  mGenerateVoid(generateVoid),
  mRemoveBoundaryZones(removeBoundaryZones),
  mXmin(xmin),
  mXmax(xmax) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
MeshPolicy<Dimension>::
~MeshPolicy() {
}

//------------------------------------------------------------------------------
// Update the Mesh.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MeshPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& /*derivs*/,
       const double /*multiplier*/,
       const double /*t*/,
       const double /*dt*/) {
  CONTRACT_VAR(key);
  REQUIRE(key == HydroFieldNames::mesh);

  // Get the state.
  const FieldList<Dimension, Vector> positions = state.fields(HydroFieldNames::position, Vector::zero());
  Mesh<Dimension>& mesh = state.mesh();
  mesh.clear();

  // If required, find the global bounding box.
  if (mComputeBounds) {
    globalBoundingBox<Dimension>(positions, mXmin, mXmax, mMeshGhostNodes);
  }

  // Do the deed.
  NodeList<Dimension> voidNodes("void", 0, 0);
  vector<const NodeList<Dimension>*> nodeLists(positions.nodeListPtrs().begin(),
                                               positions.nodeListPtrs().end());
  nodeLists.push_back(&voidNodes);
  generateMesh<Dimension, 
               typename vector<const NodeList<Dimension>*>::iterator,
               typename Physics<Dimension>::ConstBoundaryIterator>
    (nodeLists.begin(), nodeLists.end(),
     mPackage.boundaryBegin(),
     mPackage.boundaryEnd(),
     mXmin, mXmax,
     mMeshGhostNodes,
     mGenerateVoid,                    // generateVoid
     true,                             // generateParallelConnectivity
     mRemoveBoundaryZones,             // removeBoundaryZones
     2.0,                              // voidThreshold
     mesh,
     voidNodes);
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
MeshPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const MeshPolicy<Dimension>* rhsPtr = dynamic_cast<const MeshPolicy<Dimension>*>(&rhs);
  if (rhsPtr == 0) {
    return false;
  } else {
    return true;
  }
}

}

