//---------------------------------Spheral++----------------------------------//
// SphericalOriginBoundary -- A specialized Reflecting boundary condition for
// use at the origin in Spherical coordinate calculations.
// Note: this boundary is automatically constructred by the Spherical hydro objects,
//       so the user should not explicitly add this boundary.
//
// Created by JMO, Fri Feb 25 15:03:39 PST 2022
//----------------------------------------------------------------------------//
#include "SphericalOriginBoundary.hh"
#include "Geometry/GeomPlane.hh"
#include "Field/Field.hh"
#include "Utilities/DBC.hh"

using std::vector;

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
SphericalOriginBoundary::SphericalOriginBoundary():
  ReflectingBoundary<Dimension>(GeomPlane<Dim<1>>(Vector(0.0, 0.0),
                                                  Vector(0.0, 1.0))) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
SphericalOriginBoundary::~SphericalOriginBoundary() {
}

//------------------------------------------------------------------------------
// Determine the set of ghost nodes for the boundary condition.
//------------------------------------------------------------------------------
void
SphericalOriginBoundary::setGhostNodes(NodeList<Dimension>& nodeList) {
  this->addNodeList(nodeList);
}

//------------------------------------------------------------------------------
// Update the ghost nodes
//------------------------------------------------------------------------------
void
SphericalOriginBoundary::updateGhostNodes(NodeList<Dimension>& nodeList) {
}

//------------------------------------------------------------------------------
// Find the set of nodes in the given NodeList that violate this boundary
// condition.  In this case violation is being "behind" the entrance plane,
// where behind is defined in terms of the plane normal.
//------------------------------------------------------------------------------
void
SphericalOriginBoundary::setViolationNodes(NodeList<Dim<1>>& nodeList) {

  // Get the BoundaryNodes.violationNodes for this NodeList.
  this->addNodeList(nodeList);
  auto& boundaryNodes = this->accessBoundaryNodes(nodeList);
  auto& vNodes = boundaryNodes.violationNodes;
  vNodes.resize(0);

  // Loop over all the internal nodes in the NodeList, and put any that are
  // below r=0
  const auto& pos = nodeList.positions();
  const auto  n = nodeList.numInternalNodes();
  for (auto i = 0u; i < n; ++i) {
    if (pos(i).x() < 0.0) vNodes.push_back(i);
  }

  // Update the positions and H for the nodes in violation.
  updateViolationNodes(nodeList);
}

//------------------------------------------------------------------------------
// Update the nodes in violation of the boundary condition, mapping their
// positions and H tensors.
//------------------------------------------------------------------------------
void
SphericalOriginBoundary::updateViolationNodes(NodeList<Dim<1>>& nodeList) {

  // The effective plane we're reflecting from in eta space.
  GeomPlane<Dim<1>> plane(Vector(0.0), Vector(1.0));

  // Get the set of violation nodes for this NodeList.
  const auto& vNodes = this->violationNodes(nodeList);

  // Loop over these nodes, and reset their positions to valid values.
  auto& pos = nodeList.positions();
  for (auto i: vNodes) {
    auto& posi = pos(i);
    posi.x(-(posi.x()));
  }
}    

}
