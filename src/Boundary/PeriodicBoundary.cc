//---------------------------------Spheral++----------------------------------//
// PeriodicBoundary -- Apply a Periodic boundary condition to Spheral++
// Fields.
//
// Created by JMO, Wed Apr 19 14:19:07 PDT 2000
//----------------------------------------------------------------------------//
#include "PeriodicBoundary.hh"
#include "Geometry/GeomPlane.hh"
#include "Field/Field.hh"

#include "Utilities/DBC.hh"

#include <algorithm>
using std::vector;
using std::map;

namespace Spheral {

//------------------------------------------------------------------------------
// Copy all the boundary nodes from one Boundary to another.
//------------------------------------------------------------------------------
template<typename BNMap>
void
copyAllBoundaryNodes(const BNMap& fromAll,
                     BNMap& toAll) {
  typedef typename BNMap::key_type    NodeListPtr;
  typedef typename BNMap::mapped_type BoundaryNodes;
  for (typename BNMap::const_iterator itr = fromAll.begin();
       itr != fromAll.end();
       ++itr) {
    const NodeListPtr nodeListPtr = itr->first;
    const BoundaryNodes& fromBoundNodes = itr->second;
    BoundaryNodes& toBoundNodes = toAll[nodeListPtr];
    copy(fromBoundNodes.controlNodes.begin(),
         fromBoundNodes.controlNodes.end(),
         back_inserter(toBoundNodes.controlNodes));
    copy(fromBoundNodes.ghostNodes.begin(),
         fromBoundNodes.ghostNodes.end(),
         back_inserter(toBoundNodes.ghostNodes));
  }
}

//------------------------------------------------------------------------------
// Empty constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
PeriodicBoundary<Dimension>::PeriodicBoundary():
  PlanarBoundary<Dimension>(),
  mPlane1Boundary(),
  mPlane2Boundary() {
}

//------------------------------------------------------------------------------
// Construct with the given planes.
//------------------------------------------------------------------------------
template<typename Dimension>
PeriodicBoundary<Dimension>::
PeriodicBoundary(const GeomPlane<Dimension>& plane1,
                 const GeomPlane<Dimension>& plane2):
  PlanarBoundary<Dimension>(plane1, plane2),
  mPlane1Boundary(plane1, plane2),
  mPlane2Boundary(plane2, plane1) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
PeriodicBoundary<Dimension>::~PeriodicBoundary() {
}

//------------------------------------------------------------------------------
// Set the ghost nodes.  We actually want to set the ghost nodes for the 
// member classes, mPlane1Boundary & mPlane2Boundary.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PeriodicBoundary<Dimension>::setGhostNodes(NodeList<Dimension>& nodeList) {

  mPlane1Boundary.setGhostNodes(nodeList);
  mPlane2Boundary.setGhostNodes(nodeList);

  BEGIN_CONTRACT_SCOPE
  {
    const unsigned num = nodeList.numNodes();
    CONTRACT_VAR(num);
    for (unsigned i: mPlane1Boundary.controlNodes(nodeList)) { CONTRACT_VAR(i); CHECK(i < num); }
    for (unsigned i: mPlane1Boundary.ghostNodes(nodeList))   { CONTRACT_VAR(i); CHECK(i < num); }
    for (unsigned i: mPlane2Boundary.controlNodes(nodeList)) { CONTRACT_VAR(i); CHECK(i < num); }
    for (unsigned i: mPlane2Boundary.ghostNodes(nodeList))   { CONTRACT_VAR(i); CHECK(i < num); }
  }
  END_CONTRACT_SCOPE

  // Now add this NodeList to this master Boundary condition.
  this->addNodeList(nodeList);
  auto& boundaryNodes = this->accessBoundaryNodes(nodeList);
  auto& controlNodes = boundaryNodes.controlNodes;
  auto& ghostNodes = boundaryNodes.ghostNodes;

  // Copy the control info from the subclasses to this one, so it can be
  // accessed from the outside world.
  controlNodes.clear();
  controlNodes.reserve(mPlane1Boundary.controlNodes(nodeList).size() +
                       mPlane2Boundary.controlNodes(nodeList).size());
  copy(mPlane1Boundary.controlBegin(nodeList),
       mPlane1Boundary.controlEnd(nodeList),
       back_inserter(controlNodes));
  copy(mPlane2Boundary.controlBegin(nodeList),
       mPlane2Boundary.controlEnd(nodeList),
       back_inserter(controlNodes));

  // Ditto for the ghost nodes.
  ghostNodes.clear();
  ghostNodes.reserve(mPlane1Boundary.ghostNodes(nodeList).size() +
                     mPlane2Boundary.ghostNodes(nodeList).size());
  copy(mPlane1Boundary.ghostBegin(nodeList),
       mPlane1Boundary.ghostEnd(nodeList),
       back_inserter(ghostNodes));
  copy(mPlane2Boundary.ghostBegin(nodeList),
       mPlane2Boundary.ghostEnd(nodeList),
       back_inserter(ghostNodes));
}

//------------------------------------------------------------------------------
// Update the ghost nodes.  We actually want to update the ghost nodes for the 
// member classes, mPlane1Boundary & mPlane2Boundary.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PeriodicBoundary<Dimension>::updateGhostNodes(NodeList<Dimension>& nodeList) {
  mPlane1Boundary.updateGhostNodes(nodeList);
  mPlane2Boundary.updateGhostNodes(nodeList);
}

//------------------------------------------------------------------------------
// Identify the nodes in violation of the boundary.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PeriodicBoundary<Dimension>::setViolationNodes(NodeList<Dimension>& nodeList) {
  mPlane1Boundary.setViolationNodes(nodeList);
  mPlane2Boundary.setViolationNodes(nodeList);

  // Copy the violation node info from the subclasses to this one, so it can be
  // accessed from the outside world.
  this->addNodeList(nodeList);
  auto& boundaryNodes = this->accessBoundaryNodes(nodeList);
  auto& violationNodes = boundaryNodes.violationNodes;
  violationNodes.clear();
  violationNodes.reserve(mPlane1Boundary.violationNodes(nodeList).size() +
			 mPlane2Boundary.violationNodes(nodeList).size());
  copy(mPlane1Boundary.violationBegin(nodeList),
       mPlane1Boundary.violationEnd(nodeList),
       back_inserter(violationNodes));
  copy(mPlane2Boundary.violationBegin(nodeList),
       mPlane2Boundary.violationEnd(nodeList),
       back_inserter(violationNodes));

  // Update the positions and H for the nodes in violation.
  updateViolationNodes(nodeList);
}

//------------------------------------------------------------------------------
// Update  the nodes in violation of the boundary, bringing their nodes and
// H tensors into compliance.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PeriodicBoundary<Dimension>::updateViolationNodes(NodeList<Dimension>& nodeList) {
  // The sub-boundaries have already fired their updateViolationNodes!
  // mPlane1Boundary.updateViolationNodes(nodeList);
  // mPlane2Boundary.updateViolationNodes(nodeList);

  // For periodic boundaries we have the peculiar situation that a point may
  // wrap entirely around the entire volume more than once.  The underlying
  // planar boundaries can't take this into account, so we do it here.
  // Another solution is to limit the timestep by the absolute magnitude of 
  // the velocity, but we hate to do that just for a boundary condition.
  typedef GeomPlane<Dimension> Plane;
  const Plane& plane1 = mPlane1Boundary.enterPlane();
  const Plane& plane2 = mPlane1Boundary.exitPlane();
  const Vector& origin = plane1.point();
  const Vector& nhat = plane1.normal();
  const Scalar L = plane1.minimumDistance(plane2.point());
  const auto& vNodes = this->violationNodes(nodeList);
  auto& position = nodeList.positions();
  for (auto i: vNodes) {
    if (position(i) < plane1 or position(i) < plane2) {
      const Vector p0 = position(i) + (origin - position(i)).dot(nhat)*nhat; // closest point on the plane.
      const Scalar f = max(-1.0, min(1.0, fmod(plane1.signedDistance(position(i))/L, 1.0)));
      position(i) = p0 + f*L*nhat;
      CHECK((position(i) >= plane1) and (position(i) >= plane2));
    }
  }
}

//------------------------------------------------------------------------------
// Access the enter plane.
//------------------------------------------------------------------------------
template<typename Dimension>
const GeomPlane<Dimension>&
PeriodicBoundary<Dimension>::enterPlane() const {
  CHECK(mPlane1Boundary.enterPlane() == mPlane2Boundary.exitPlane());
  return mPlane1Boundary.enterPlane();
}

template<typename Dimension>
void
PeriodicBoundary<Dimension>::
setEnterPlane(const GeomPlane<Dimension>& enterPlane) {
  mPlane1Boundary.setEnterPlane(enterPlane);
  mPlane2Boundary.setExitPlane(enterPlane);
}

//------------------------------------------------------------------------------
// Access the exit plane.
//------------------------------------------------------------------------------
template<typename Dimension>
const GeomPlane<Dimension>&
PeriodicBoundary<Dimension>::exitPlane() const {
  CHECK(mPlane1Boundary.exitPlane() == mPlane2Boundary.enterPlane());
  return mPlane1Boundary.exitPlane();
}

template<typename Dimension>
void
PeriodicBoundary<Dimension>::
setExitPlane(const GeomPlane<Dimension>& exitPlane) {
  mPlane1Boundary.setExitPlane(exitPlane);
  mPlane2Boundary.setEnterPlane(exitPlane);
}

//------------------------------------------------------------------------------
// Dispatch culling of ghost nodes to the individual subboundaries.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PeriodicBoundary<Dimension>::
cullGhostNodes(const FieldList<Dimension, size_t>& flagSet,
               FieldList<Dimension, size_t>& old2newIndexMap,
               vector<size_t>& numNodesRemoved) {
  mPlane1Boundary.cullGhostNodes(flagSet, old2newIndexMap, numNodesRemoved);
  mPlane2Boundary.cullGhostNodes(flagSet, old2newIndexMap, numNodesRemoved);

  // Recreate our set of ghost & control.
  typedef typename Boundary<Dimension>::BoundaryNodes BoundaryNodes;
  typedef map<NodeList<Dimension>*, BoundaryNodes> BNMap;
  BNMap& allBoundaryNodes = this->accessBoundaryNodes();
  allBoundaryNodes = BNMap();
  const BNMap& allBoundaryNodes1 = mPlane1Boundary.accessBoundaryNodes();
  const BNMap& allBoundaryNodes2 = mPlane2Boundary.accessBoundaryNodes();
  copyAllBoundaryNodes(allBoundaryNodes1, allBoundaryNodes);
  copyAllBoundaryNodes(allBoundaryNodes2, allBoundaryNodes);
}

//------------------------------------------------------------------------------
// Apply the boundary condition to fields of different DataTypes.
// For the Periodic boundary, the nested sub classes do all the work.
//------------------------------------------------------------------------------
// FieldBase
template<typename Dimension>
void
PeriodicBoundary<Dimension>::
applyGhostBoundary(FieldBase<Dimension>& fieldBase) const {
  mPlane1Boundary.applyGhostBoundary(fieldBase);
  mPlane2Boundary.applyGhostBoundary(fieldBase);
}

// FacetedVolume Fields.
template<typename Dimension>
void
PeriodicBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::FacetedVolume>& field) const {
  mPlane1Boundary.applyGhostBoundary(field);
  mPlane2Boundary.applyGhostBoundary(field);
}

//------------------------------------------------------------------------------
// Enforce the boundary condition to fields of different DataTypes.
// For the Periodic boundary, the nested sub classes do all the work.
//------------------------------------------------------------------------------
// FacetedVolumeFields.
template<typename Dimension>
void
PeriodicBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::FacetedVolume>& field) const {
  mPlane1Boundary.enforceBoundary(field);
  mPlane2Boundary.enforceBoundary(field);
}

//------------------------------------------------------------------------------
// Clear out any preexisting control/ghost node info.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PeriodicBoundary<Dimension>::reset(const DataBase<Dimension>& dataBase) {
  // We have to clear out not just this boundary, but the two member boundary
  // classes as well.
  Boundary<Dimension>::reset(dataBase);
  mPlane1Boundary.reset(dataBase);
  mPlane2Boundary.reset(dataBase);
}

//------------------------------------------------------------------------------
// Clear out any preexisting control/ghost node info.
//------------------------------------------------------------------------------
template<typename Dimension>
int
PeriodicBoundary<Dimension>::numGhostNodes() const {
  return mPlane1Boundary.numGhostNodes() + mPlane2Boundary.numGhostNodes();
}

//------------------------------------------------------------------------------
// Include methods for the nested class type, PeriodicPlanarBoundary.
//------------------------------------------------------------------------------
#include "PeriodicPlanarBoundary.cc"
}
