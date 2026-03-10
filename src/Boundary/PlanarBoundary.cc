//---------------------------------Spheral++----------------------------------//
// PlanarBoundary -- Abstract base class for the boundary conditions defined
// on planes.  Planar boundaries in general are defined in terms of parallel
// pairs of planes, representing entrance and exit conditions.
//
// Created by JMO, Thu Mar  2 21:34:25 PST 2000
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
#include "findNodesTouchingThroughPlanes.hh"
#include "mapPositionThroughPlanes.hh"
#include "Geometry/GeomPlane.hh"
#include "NodeList/FluidNodeList.hh"
#include "Mesh/Mesh.hh"
#include "Utilities/DBC.hh"

#include "PlanarBoundary.hh"

using std::vector;

namespace Spheral {

//------------------------------------------------------------------------------
// Internal worker method to help with clipping a box range.
//------------------------------------------------------------------------------
// 1D
inline
void
clipBoxWithPlane(const GeomPlane<Dim<1> >& plane,
                 Dim<1>::Vector& point) {
  if (plane.compare(point) == 1) point = plane.point();
}

// 2D
// For now we only check cardinally aligned planes.
inline
void
clipBoxWithPlane(const GeomPlane<Dim<2> >& plane,
                 Dim<2>::Vector& point) {
  if (plane.compare(point) == 1) {
    if (fuzzyEqual(std::abs(plane.normal().x()), 1.0)) {
      point.x(plane.point().x());
    } else if (fuzzyEqual(std::abs(plane.normal().y()), 1.0)) {
      point.y(plane.point().y());
    }
  }
}

// 3D
// For now we only check cardinally aligned planes.
inline
void
clipBoxWithPlane(const GeomPlane<Dim<3> >& plane,
                 Dim<3>::Vector& point) {
  if (plane.compare(point) == 1) {
    if (fuzzyEqual(std::abs(plane.normal().x()), 1.0)) {
      point.x(plane.point().x());
    } else if (fuzzyEqual(std::abs(plane.normal().y()), 1.0)) {
      point.y(plane.point().y());
    } else if (fuzzyEqual(std::abs(plane.normal().z()), 1.0)) {
      point.z(plane.point().z());
    }
  }
}

//------------------------------------------------------------------------------
// Empty constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
PlanarBoundary<Dimension>::PlanarBoundary():
  Boundary<Dimension>(),
  mEnterPlane(),
  mExitPlane(),
  mRestart(registerWithRestart(*this)) {
}

//------------------------------------------------------------------------------
// Construct with the given entrance and exit planes.
//------------------------------------------------------------------------------
template<typename Dimension>
PlanarBoundary<Dimension>::
PlanarBoundary(const GeomPlane<Dimension>& enterPlane,
               const GeomPlane<Dimension>& exitPlane):
  Boundary<Dimension>(),
  mEnterPlane(enterPlane),
  mExitPlane(exitPlane),
  mRestart(registerWithRestart(*this)) {
  ENSURE(valid());
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
PlanarBoundary<Dimension>::~PlanarBoundary() {
}

//------------------------------------------------------------------------------
// Determine the set of ghost nodes for the boundary condition.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PlanarBoundary<Dimension>::setGhostNodes(NodeList<Dimension>& nodeList) {

  // Remember which node list we are setting the ghost nodes for.
  this->addNodeList(nodeList);
  auto& boundaryNodes = this->accessBoundaryNodes(nodeList);
  boundaryNodes.controlNodes = findNodesTouchingThroughPlanes(nodeList, mEnterPlane, mExitPlane);

  // std::sort(controlNodes.begin(), controlNodes.end());
  // controlNodes.erase(std::unique(controlNodes.begin(), controlNodes.end()), controlNodes.end());

  // Set the ghost node indices to correspond to these control nodes.
  setGhostNodeIndices(nodeList);

  // Assign the positions and H's to the new ghost nodes.
  updateGhostNodes(nodeList);
}

//------------------------------------------------------------------------------
// Set the ghost nodes for the given NodeList using the given set of control
// nodes.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PlanarBoundary<Dimension>::
setGhostNodes(NodeList<Dimension>& nodeList, 
              const vector<size_t>& presetControlNodes) {

  typedef typename Boundary<Dimension>::BoundaryNodes BoundaryNodes;

  // Add this NodeList, creating space for control & ghost nodes.
  this->addNodeList(nodeList);

  // Set the list of control nodes.
  BoundaryNodes& boundaryNodes = this->accessBoundaryNodes(nodeList);
  auto& controlNodes = boundaryNodes.controlNodes;
  controlNodes = presetControlNodes;

  // Set the ghost node indices to correspond to these control nodes.
  setGhostNodeIndices(nodeList);

  // Assign the masses and positions to the new ghost nodes.
  updateGhostNodes(nodeList);
}

//------------------------------------------------------------------------------
// Find the set of nodes in the given NodeList that violate this boundary
// condition.  In this case violation is being "behind" the entrance plane,
// where behind is defined in terms of the plane normal.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PlanarBoundary<Dimension>::setViolationNodes(NodeList<Dimension>& nodeList) {

  // Get the BoundaryNodes.violationNodes for this NodeList.
  typedef typename Boundary<Dimension>::BoundaryNodes BoundaryNodes;
  this->addNodeList(nodeList);
  BoundaryNodes& boundaryNodes = this->accessBoundaryNodes(nodeList);
  auto& vNodes = boundaryNodes.violationNodes;
  vNodes.resize(0);

  // Loop over all the internal nodes in the NodeList, and put any that are 
  // below the enter plane in the list of nodes in violation.
  const Field<Dimension, Vector>& positions = nodeList.positions();
  for (auto nodeID = 0u; nodeID < nodeList.numInternalNodes(); ++nodeID) {
    if (positions(nodeID) < mEnterPlane) vNodes.push_back(nodeID);
  }

  // Update the positions and H for the nodes in violation.
  updateViolationNodes(nodeList);
}

//------------------------------------------------------------------------------
// Update the nodes in violation of the boundary condition, mapping their
// positions and H tensors.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PlanarBoundary<Dimension>::updateViolationNodes(NodeList<Dimension>& nodeList) {

  // Get the set of violation nodes for this NodeList.
  const auto& vNodes = this->violationNodes(nodeList);

  // Loop over these nodes, and reset their positions to valid values.
  Field<Dimension, Vector>& positions = nodeList.positions();
  for (auto i: vNodes) {
    CHECK(positions(i) <= enterPlane());
    positions(i) = mapPosition(positions(i), mEnterPlane, mExitPlane);
    // CHECK2((positions(*itr) >= enterPlane()) and
    //        (positions(*itr) >= exitPlane()),
    //        "Bad position mapping: " << *itr << " " << nodeList.firstGhostNode() << " " << positions(*itr));
  }

  // Set the Hfield.
  auto& Hfield = nodeList.Hfield();
  this->enforceBoundary(Hfield);

//   // Update the neighbor information.
//   nodeList.neighbor().updateNodes(); // (vNodes);
}    

//------------------------------------------------------------------------------
// Function to map a position through the enter to the exit plane.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Vector
PlanarBoundary<Dimension>::
mapPosition(const Vector& position,
	    const GeomPlane<Dimension>& enterPlane,
	    const GeomPlane<Dimension>& exitPlane) const {
  REQUIRE(enterPlane.valid() and exitPlane.valid());
  return mapPositionThroughPlanes(position, enterPlane, exitPlane);
//   REQUIRE(enterPlane.parallel(exitPlane));
//   const Vector deltaEnter = (position - enterPlane.point()).dot(enterPlane.normal())*enterPlane.normal();
//   const Vector deltaExit = (position - exitPlane.point()).dot(exitPlane.normal())*exitPlane.normal();
//   double sign = (position - enterPlane.point()).dot(enterPlane.normal());
//   sign /= fabs(sign) + 1.0e-20;
//   return position - deltaExit - sign*deltaEnter.magnitude()*exitPlane.normal();
}

//------------------------------------------------------------------------------
// Dump state.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PlanarBoundary<Dimension>::
dumpState(FileIO& file,
          const std::string& pathName) const {
  file.write(enterPlane(), pathName + "/enterPlane");
  file.write(exitPlane(), pathName + "/exitPlane");
}

//------------------------------------------------------------------------------
// Restore state.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PlanarBoundary<Dimension>::
restoreState(const FileIO& file,
             const std::string& pathName) {
  file.read(mEnterPlane, pathName + "/enterPlane");
  file.read(mExitPlane, pathName + "/exitPlane");
}

//------------------------------------------------------------------------------
// Valid test.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
PlanarBoundary<Dimension>::valid() const {
  return (enterPlane().valid() and
          exitPlane().valid() and
          enterPlane().parallel(exitPlane()));
}

//------------------------------------------------------------------------------
// Internal method to set the ghost node indices once the controls are set.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PlanarBoundary<Dimension>::setGhostNodeIndices(NodeList<Dimension>& nodeList) {

  // Get the sets of control and ghost nodes.
  BoundaryNodes& boundaryNodes = this->accessBoundaryNodes(nodeList);
  auto& controlNodes = boundaryNodes.controlNodes;
  auto& ghostNodes = boundaryNodes.ghostNodes;

  // Create the ghosts of the control nodes.
  auto currentNumGhostNodes = nodeList.numGhostNodes();
  auto firstNewGhostNode = nodeList.numNodes();
  nodeList.numGhostNodes(currentNumGhostNodes + controlNodes.size());
  CHECK(nodeList.numNodes() == firstNewGhostNode + controlNodes.size());

  // Fill in the pointers to the ghost nodes.
  ghostNodes.resize(controlNodes.size());
  for (auto i = 0u; i < controlNodes.size(); ++i) {
    CHECK(i < controlNodes.size());
    CHECK(i < ghostNodes.size());
    ghostNodes[i] = firstNewGhostNode + i;
    CHECK(ghostNodes[i] >= nodeList.numInternalNodes() and
          ghostNodes[i] < nodeList.numNodes());
  }
}

//------------------------------------------------------------------------------
// Internal method to set the minimal field values necessary to make the 
// ghost nodes valid.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PlanarBoundary<Dimension>::updateGhostNodes(NodeList<Dimension>& nodeList) {
  REQUIRE(valid());

  // Get the control and ghost node indices.
  BoundaryNodes& boundaryNodes = this->accessBoundaryNodes(nodeList);
  auto& controlNodes = boundaryNodes.controlNodes;
  auto& ghostNodes = boundaryNodes.ghostNodes;
  CHECK(controlNodes.size() == ghostNodes.size());

  // Loop over the control/ghost node pairs, and set the ghost positions.
  Field<Dimension, Vector>& positions = nodeList.positions();
  auto controlItr = controlNodes.begin();
  auto ghostItr = ghostNodes.begin();
  for (;controlItr < controlNodes.end(); ++controlItr, ++ghostItr) {
    CHECK(*controlItr >= 0 and *controlItr < nodeList.numNodes());
    CHECK2(*ghostItr >= nodeList.firstGhostNode() and *ghostItr < nodeList.numNodes(),
           "Ghost node index out of bounds:  " << *ghostItr << " " << nodeList.firstGhostNode() << " " << nodeList.numNodes());

//     if (!(positions(*controlItr) >= mEnterPlane and
//           positions(*controlItr) >= mExitPlane)) {
//       cerr << "Node out of bounds "
//            << (*controlItr) << " " 
//            << positions(*controlItr) << endl;
//     }
    // CHECK(positions(*controlItr) >= mEnterPlane and
    //       positions(*controlItr) >= mExitPlane);

    positions(*ghostItr) = mapPosition(positions(*controlItr),
                                       mExitPlane,
                                       mEnterPlane);
    // CHECK(positions(*ghostItr) <= mEnterPlane);
  }

  // Set the Hfield.
  Field<Dimension, SymTensor>& Hfield = nodeList.Hfield();
  this->applyGhostBoundary(Hfield);

  // // Update the neighbor information.
  // nodeList.neighbor().updateNodes(); // (ghostNodes);
}

//------------------------------------------------------------------------------
// Clip an input box range.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PlanarBoundary<Dimension>::
clip(typename Dimension::Vector& xmin, typename Dimension::Vector& xmax) const {
  clipBoxWithPlane(mEnterPlane, xmin);
  clipBoxWithPlane(mEnterPlane, xmax);
  clipBoxWithPlane(mExitPlane, xmin);
  clipBoxWithPlane(mExitPlane, xmax);
}

//------------------------------------------------------------------------------
// Find the set of tessellation facets on a plane.
//------------------------------------------------------------------------------
template<typename Dimension>
vector<unsigned>
PlanarBoundary<Dimension>::
facesOnPlane(const Mesh<Dimension>& mesh,
             const GeomPlane<Dimension>& plane,
             const Scalar tol) const {

  typedef typename Mesh<Dimension>::Face Face;

  vector<unsigned> result;

  // Flag all the nodes that are on the plane to the given tolerance.
  const unsigned numNodes = mesh.numNodes();
  vector<unsigned> flagNodes(numNodes, 0);
  for (unsigned inode = 0; inode != numNodes; ++inode) {
    flagNodes[inode] = (plane.minimumDistance(mesh.node(inode).position()) <= tol ? 1 : 0);
  }

  // Look for any faces with all its nodes in the plane.
  const unsigned numFaces = mesh.numFaces();
  for (unsigned iface = 0; iface != numFaces; ++iface) {
    const Face& face = mesh.face(iface);
    const vector<unsigned> nodeIDs = face.nodeIDs();
    unsigned i = 0;
    while (i < nodeIDs.size() and flagNodes[nodeIDs[i]] == 1) ++i;
    if (i == nodeIDs.size()) result.push_back(iface);
  }

  return result;
}

}
