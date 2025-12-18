//---------------------------------Spheral++----------------------------------//
// SphericalBoundary -- Create a 3-D boundary around a line of nodes down
// the x-axis appropriate for enforcing a spherically symmetric system.
//
// Created by JMO, Tue Mar 15 21:39:43 PST 2005
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
#include "Geometry/GeomPlane.hh"
#include "Geometry/innerProduct.hh"
#include "NodeList/FluidNodeList.hh"
#include "Utilities/DBC.hh"

#include "SphericalBoundary.hh"

using std::vector;
using std::string;

namespace Spheral {

//------------------------------------------------------------------------------
// Construct against the given DataBase.
//------------------------------------------------------------------------------
SphericalBoundary::
SphericalBoundary(const DataBase<Dim<3> >& dataBase):
  Boundary<Dim<3> >(),
  mGhostPositions(dataBase.newGlobalFieldList(std::vector<Dim<3>::Vector>(),
                                              "Ghost node positions")),
  mRestart(registerWithRestart(*this)) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
SphericalBoundary::
~SphericalBoundary() {
}

//------------------------------------------------------------------------------
// Set the ghost nodes for the given NodeList.
//------------------------------------------------------------------------------
void
SphericalBoundary::
setGhostNodes(NodeList<Dim<3> >& nodeList) {
  REQUIRE(mGhostPositions.fieldForNodeList(nodeList) < mGhostPositions.end());

  // Add this NodeList, creating space for control & ghost nodes.
  addNodeList(nodeList);
  BoundaryNodes& boundaryNodes = accessBoundaryNodes(nodeList);
  auto& controlNodes = boundaryNodes.controlNodes;
  auto& ghostNodes = boundaryNodes.ghostNodes;
  controlNodes.clear();
  ghostNodes.clear();

  // The field of ghost positions for this NodeList.
  Field<Dim<3>, vector<Vector> >& ghostPositions = **(mGhostPositions.fieldForNodeList(nodeList));

  // The ghost node indices.
  const int currentNumGhostNodes = nodeList.numGhostNodes();
  int ghostNodeIndex = nodeList.numNodes();

  // Get the state for this NodeList.
  Field<Dimension, Vector>& position = nodeList.positions();
  Field<Dimension, SymTensor>& H = nodeList.Hfield();

  // Get the number of nodes per h, and the kernel extent.
  const double nPerh = nodeList.nodesPerSmoothingScale();
  const double kernelExtent = nodeList.neighbor().kernelExtent();
  CHECK(distinctlyGreaterThan(nPerh, 0.0));
  CHECK(distinctlyGreaterThan(kernelExtent, 0.0));

  // The critical phi angle, at which the surface of ghost nodes will 
  // extend for 45 degrees around the full sphere.  At this point we
  // decide it is worth simply filling in the entire spherical surface
  // with ghost nodes.
  const double phicrit = 0.25*M_PI/(nPerh*kernelExtent);

  // Iterate over the internal nodes.
  for (auto i = 0u; i != nodeList.numInternalNodes(); ++i) {

    // Have we computed the set of ghost nodes for this node yet?
    if (ghostPositions(i).size() == 0) {

      // This is a new node we have not computed ghost nodes for yet.
      // Go ahead and create the appropriate set of ghost nodes for
      // it.
      const double ri = position(i).x();
      const SymTensor& Hi = H(i);

      // Determine H extent in the yz plane.
      const double hyz = sqrt((Hi*Vector(1.0, 0.0, 0.0)).magnitude()/Hi.Determinant());

      // The delta expected between nodes, and corresponding angle.
      const double delta = hyz/nPerh;
      const double dphi = 2.0*atan2(delta, 2.0*ri);

      // Do we need to create a full sphere of ghost nodes?
      if (dphi < phicrit) {

        // We're outside the critical point for doing the entire sphere, so 
        // just create a section of the surface with enough ghost nodes that 
        // this one shouldn't see beyond them.
        const int nbins = int(kernelExtent*nPerh + 0.5);
        CHECK(nbins >= 1);
        for (int ibin = 0; ibin != nbins; ++ibin) {
          const double phij = (ibin + 1)*dphi;
          CHECK(phij < 0.5*M_PI);
          const double c = ri*sqrt(2.0*(1.0 - cos(phij)));
          const double ryz = c*sin(0.5*(M_PI - phij));
          const double circ = 2.0*M_PI*ryz;
          const int nTheta = max(1, int(circ/delta + 0.5));
          const double dTheta = 2.0*M_PI/nTheta;
          for (int j = 0; j != nTheta; ++j) {
            const double thetaj = (j + 0.5)*dTheta;
            ghostPositions(i).push_back(Vector(ri*cos(phij),
                                               ri*cos(thetaj)*sin(phij),
                                               ri*sin(thetaj)*sin(phij)));
          }
        }

      } else {

        // We're near enough the origin that we need to create the full sphere
        // of ghost points at this radius.
        const int nbins = max(1, int(0.5*M_PI/dphi + 0.5));
        CHECK(nbins >= 1);
        const double dphieff = 0.5*M_PI/nbins;
        for (int ibin = 0; ibin != nbins; ++ibin) {
          const double phij = (ibin + 1)*dphieff;
          CHECK(phij <= 0.5*M_PI);
          const double c = ri*sqrt(2.0*(1.0 - cos(phij)));
          const double ryz = c*sin(0.5*(M_PI - phij));
          const double circ = 2.0*M_PI*ryz;
          const int nTheta = max(1, int(circ/delta + 0.5));
          const double dTheta = 2.0*M_PI/nTheta;
          for (int j = 0; j != nTheta; ++j) {
            const double thetaj = (j + 0.5)*dTheta;
            ghostPositions(i).push_back(Vector(ri*cos(phij),
                                               ri*cos(thetaj)*sin(phij),
                                               ri*sin(thetaj)*sin(phij)));
            ghostPositions(i).push_back(Vector(ri*cos(phij + M_PI),
                                               ri*cos(thetaj)*sin(phij + M_PI),
                                               ri*sin(thetaj)*sin(phij + M_PI)));
          }
        }

        // Create the mirror image of the control node across the origin.
        ghostPositions(i).push_back(Vector(-ri, 0.0, 0.0));

      }
    }

    // Set the control and ghost node indices for this internal node.
    const int nghosti = ghostPositions(i).size();
    CHECK(nghosti > 0);
    controlNodes.reserve(controlNodes.size() + nghosti);
    ghostNodes.reserve(ghostNodes.size() + nghosti);
    for (int j = 0; j != nghosti; ++j, ++ghostNodeIndex) {
      controlNodes.push_back(i);
      ghostNodes.push_back(ghostNodeIndex);
    }

  }

  // Create the set of ghost nodes on the NodeList.
  const int numNewGhostNodes = ghostNodes.size();
  CHECK((int)controlNodes.size() == numNewGhostNodes);
  nodeList.numGhostNodes(currentNumGhostNodes + numNewGhostNodes);

  // Use the updateGhostNodes method to set the node positions and H's.
  this->updateGhostNodes(nodeList);
}

//------------------------------------------------------------------------------
// Update the state of the current ghost nodes (correcting their positions
// and H tensors.)
//------------------------------------------------------------------------------
void
SphericalBoundary::
updateGhostNodes(NodeList<Dim<3> >& nodeList) {
  REQUIRE(mGhostPositions.fieldForNodeList(nodeList) < mGhostPositions.end());

  // Get the control and ghost nodes for this NodeList.
  const BoundaryNodes& boundaryNodes = accessBoundaryNodes(nodeList);
  const auto& controlNodes = boundaryNodes.controlNodes;
  const auto& ghostNodes = boundaryNodes.ghostNodes;

  // The field of ghost positions for this NodeList.
  Field<Dim<3>, vector<Vector> >& ghostPositions = **(mGhostPositions.fieldForNodeList(nodeList));

  // Get the positions for this NodeList.
  Field<Dimension, Vector>& position = nodeList.positions();

  // Update the positions, forcing the ghost nodes to exist at the
  // same radius as their control nodes.
  CHECK(controlNodes.size() == ghostNodes.size());
  auto controlItr = controlNodes.begin();
  auto ghostItr = ghostNodes.begin();
  for (auto i = 0u; i != nodeList.numInternalNodes(); ++i) {
    const double ri = position(i).magnitude();
    CHECK(ghostPositions(i).size() > 0);
    for (vector<Vector>::iterator ghostPositionItr = ghostPositions(i).begin();
         ghostPositionItr != ghostPositions(i).end();
         ++ghostPositionItr, ++controlItr, ++ghostItr) {
      CHECK(controlItr < controlNodes.end());
      CHECK(ghostItr < ghostNodes.end());
      CHECK(*controlItr == i);
      CHECK(*ghostItr >= nodeList.firstGhostNode() && *ghostItr < nodeList.numNodes());
      const double rj = ghostPositionItr->magnitude();
      *ghostPositionItr *= ri*rj/(rj*rj + 1.0e-50);
      position(*ghostItr) = *ghostPositionItr;
    }
  }
  CHECK(controlItr == controlNodes.end());
  CHECK(ghostItr == ghostNodes.end());

  // We can use the normal ghost boundary enforcement to update the
  // H's of the ghost nodes.
  Field<Dimension, SymTensor>& H = nodeList.Hfield();
  applyGhostBoundary(H);

//   // Update the neighbor information.
//   nodeList.neighbor().updateNodes(); // (ghostNodes);
}

//------------------------------------------------------------------------------
// Apply the ghost boundary condition to fields of different DataTypes.
//------------------------------------------------------------------------------
// Specialization for int fields, just perform a copy.
void
SphericalBoundary::
applyGhostBoundary(Field<Dim<3> , int>& field) const {

  // Apply the boundary condition to all the ghost node values.
  const NodeList<Dimension>& nodeList = field.nodeList();
  CHECK(controlNodes(nodeList).size() == ghostNodes(nodeList).size());
  auto controlItr = controlBegin(nodeList);
  auto ghostItr = ghostBegin(nodeList);
  for (; controlItr < controlEnd(nodeList); ++controlItr, ++ghostItr) {
    CHECK(ghostItr < ghostEnd(nodeList));
    CHECK(*controlItr < nodeList.numNodes());
    CHECK(*ghostItr >= nodeList.firstGhostNode() && *ghostItr < nodeList.numNodes());
    field(*ghostItr) = field(*controlItr);
  }
}

// Specialization for scalar fields, just perform a copy.
void
SphericalBoundary::
applyGhostBoundary(Field<Dim<3> , Dim<3>::Scalar>& field) const {

  // Apply the boundary condition to all the ghost node values.
  const NodeList<Dimension>& nodeList = field.nodeList();
  CHECK(controlNodes(nodeList).size() == ghostNodes(nodeList).size());
  auto controlItr = controlBegin(nodeList);
  auto ghostItr = ghostBegin(nodeList);
  for (; controlItr < controlEnd(nodeList); ++controlItr, ++ghostItr) {
    CHECK(ghostItr < ghostEnd(nodeList));
    CHECK(*controlItr < nodeList.numNodes());
    CHECK(*ghostItr >= nodeList.firstGhostNode() && *ghostItr < nodeList.numNodes());
    field(*ghostItr) = field(*controlItr);
  }
}

// Specialization for Vector fields.
void
SphericalBoundary::
applyGhostBoundary(Field<Dim<3>, Dim<3>::Vector>& field) const {

  // Apply the boundary condition to all the ghost node values.
  const NodeList<Dimension>& nodeList = field.nodeList();
  const Field<Dimension, Vector>& position = nodeList.positions();
  CHECK(controlNodes(nodeList).size() == ghostNodes(nodeList).size());
  auto controlItr = controlBegin(nodeList);
  auto ghostItr = ghostBegin(nodeList);
  for (; controlItr < controlEnd(nodeList); ++controlItr, ++ghostItr) {
    CHECK(ghostItr < ghostEnd(nodeList));
    CHECK(*controlItr < nodeList.numNodes());
    CHECK(*ghostItr >= nodeList.firstGhostNode() && *ghostItr < nodeList.numNodes());
    field(*ghostItr) = reflectOperator(position(*controlItr),
                                       position(*ghostItr))*field(*controlItr);
  }
}

// Specialization for Tensor fields.
void
SphericalBoundary::
applyGhostBoundary(Field<Dim<3>, Dim<3>::Tensor>& field) const {

  // Apply the boundary condition to all the ghost node values.
  const NodeList<Dimension>& nodeList = field.nodeList();
  const Field<Dimension, Vector>& position = nodeList.positions();
  CHECK(controlNodes(nodeList).size() == ghostNodes(nodeList).size());
  auto controlItr = controlBegin(nodeList);
  auto ghostItr = ghostBegin(nodeList);
  for (; controlItr < controlEnd(nodeList); ++controlItr, ++ghostItr) {
    CHECK(ghostItr < ghostEnd(nodeList));
    CHECK(*controlItr < nodeList.numNodes());
    CHECK(*ghostItr >= nodeList.firstGhostNode() && *ghostItr < nodeList.numNodes());
    const Tensor T = reflectOperator(position(*controlItr), position(*ghostItr));
    field(*ghostItr) = T*(field(*controlItr)*T);
  }
}

// Specialization for symmetric tensors.
void
SphericalBoundary::
applyGhostBoundary(Field<Dim<3>, Dim<3>::SymTensor>& field) const {

  // Apply the boundary condition to all the ghost node values.
  const NodeList<Dimension>& nodeList = field.nodeList();
  const Field<Dimension, Vector>& position = nodeList.positions();
  CHECK(controlNodes(nodeList).size() == ghostNodes(nodeList).size());
  auto controlItr = controlBegin(nodeList);
  auto ghostItr = ghostBegin(nodeList);
  for (; controlItr < controlEnd(nodeList); ++controlItr, ++ghostItr) {
    CHECK(ghostItr < ghostEnd(nodeList));
    CHECK(*controlItr < nodeList.numNodes());
    CHECK(*ghostItr >= nodeList.firstGhostNode() && *ghostItr < nodeList.numNodes());
    const Tensor T = reflectOperator(position(*controlItr), position(*ghostItr));
    field(*ghostItr) = T*(field(*controlItr)*T).Symmetric();
  }
}

// Specialization for third rank tensors.
void
SphericalBoundary::
applyGhostBoundary(Field<Dim<3>, Dim<3>::ThirdRankTensor>& field) const {

  // Apply the boundary condition to all the ghost node values.
  const NodeList<Dimension>& nodeList = field.nodeList();
  const Field<Dimension, Vector>& position = nodeList.positions();
  CHECK(controlNodes(nodeList).size() == ghostNodes(nodeList).size());
  auto controlItr = controlBegin(nodeList);
  auto ghostItr = ghostBegin(nodeList);
  ThirdRankTensor val;
  for (; controlItr < controlEnd(nodeList); ++controlItr, ++ghostItr) {
    CHECK(ghostItr < ghostEnd(nodeList));
    CHECK(*controlItr < nodeList.numNodes());
    CHECK(*ghostItr >= nodeList.firstGhostNode() && *ghostItr < nodeList.numNodes());
    const Tensor T = reflectOperator(position(*controlItr), position(*ghostItr));
    val = ThirdRankTensor::zero();
    const ThirdRankTensor& fc = field(*controlItr);
    for (unsigned i = 0; i != Dimension::nDim; ++i) {
      for (unsigned j = 0; j != Dimension::nDim; ++j) {
        for (unsigned k = 0; k != Dimension::nDim; ++k) {
          for (unsigned q = 0; q != Dimension::nDim; ++q) {
            for (unsigned r = 0; r != Dimension::nDim; ++r) {
              for (unsigned s = 0; s != Dimension::nDim; ++s) {
                val(i,j,k) += T(i,q)*T(j,r)*T(k,s)*fc(q,r,s);
              }
            }
          }
        }
      }
    }
    field(*ghostItr) = val;
  }
}

// Specialization for fourth rank tensors.
void
SphericalBoundary::
applyGhostBoundary(Field<Dim<3>, Dim<3>::FourthRankTensor>& field) const {

  // Apply the boundary condition to all the ghost node values.
  const NodeList<Dimension>& nodeList = field.nodeList();
  const Field<Dimension, Vector>& position = nodeList.positions();
  CHECK(controlNodes(nodeList).size() == ghostNodes(nodeList).size());
  auto controlItr = controlBegin(nodeList);
  auto ghostItr = ghostBegin(nodeList);
  FourthRankTensor val;
  for (; controlItr < controlEnd(nodeList); ++controlItr, ++ghostItr) {
    CHECK(ghostItr < ghostEnd(nodeList));
    CHECK(*controlItr < nodeList.numNodes());
    CHECK(*ghostItr >= nodeList.firstGhostNode() && *ghostItr < nodeList.numNodes());
    const Tensor T = reflectOperator(position(*controlItr), position(*ghostItr));
    val = FourthRankTensor::zero();
    const FourthRankTensor& fc = field(*controlItr);
    for (unsigned i = 0; i != Dimension::nDim; ++i) {
      for (unsigned j = 0; j != Dimension::nDim; ++j) {
        for (unsigned k = 0; k != Dimension::nDim; ++k) {
          for (unsigned l = 0; l != Dimension::nDim; ++l) {
            for (unsigned q = 0; q != Dimension::nDim; ++q) {
              for (unsigned r = 0; r != Dimension::nDim; ++r) {
                for (unsigned s = 0; s != Dimension::nDim; ++s) {
                  for (unsigned t = 0; t != Dimension::nDim; ++t) {
                    val(i,j,k,l) += T(i,q)*T(j,r)*T(k,s)*T(l,t)*fc(q,r,s,t);
                  }
                }
              }
            }
          }
        }
      }
    }
    field(*ghostItr) = val;
  }
}

// Specialization for fifth rank tensors.
void
SphericalBoundary::
applyGhostBoundary(Field<Dim<3>, Dim<3>::FifthRankTensor>& field) const {

  // Apply the boundary condition to all the ghost node values.
  const NodeList<Dimension>& nodeList = field.nodeList();
  const Field<Dimension, Vector>& position = nodeList.positions();
  CHECK(controlNodes(nodeList).size() == ghostNodes(nodeList).size());
  auto controlItr = controlBegin(nodeList);
  auto ghostItr = ghostBegin(nodeList);
  FifthRankTensor val;
  for (; controlItr < controlEnd(nodeList); ++controlItr, ++ghostItr) {
    CHECK(ghostItr < ghostEnd(nodeList));
    CHECK(*controlItr < nodeList.numNodes());
    CHECK(*ghostItr >= nodeList.firstGhostNode() && *ghostItr < nodeList.numNodes());
    const Tensor T = reflectOperator(position(*controlItr), position(*ghostItr));
    val = FifthRankTensor::zero();
    const FifthRankTensor& fc = field(*controlItr);
    for (unsigned i = 0; i != Dimension::nDim; ++i) {
      for (unsigned j = 0; j != Dimension::nDim; ++j) {
        for (unsigned k = 0; k != Dimension::nDim; ++k) {
          for (unsigned l = 0; l != Dimension::nDim; ++l) {
            for (unsigned m = 0; m != Dimension::nDim; ++m) {
              for (unsigned q = 0; q != Dimension::nDim; ++q) {
                for (unsigned r = 0; r != Dimension::nDim; ++r) {
                  for (unsigned s = 0; s != Dimension::nDim; ++s) {
                    for (unsigned t = 0; t != Dimension::nDim; ++t) {
                      for (unsigned u = 0; u != Dimension::nDim; ++u) {
                        val(i,j,k,l,u) += T(i,q)*T(j,r)*T(k,s)*T(l,t)*T(m,u)*fc(q,r,s,t,u);
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    field(*ghostItr) = val;
  }
}

// Specialization for FacetedVolume fields.
void
SphericalBoundary::
applyGhostBoundary(Field<Dim<3>, Dim<3>::FacetedVolume>& field) const {

  // Apply the boundary condition to all the ghost node values.
  const auto& nodeList = field.nodeList();
  const auto& position = nodeList.positions();
  CHECK(controlNodes(nodeList).size() == ghostNodes(nodeList).size());
  auto controlItr = controlBegin(nodeList);
  auto ghostItr = ghostBegin(nodeList);
  for (; controlItr < controlEnd(nodeList); ++controlItr, ++ghostItr) {
    CHECK(ghostItr < ghostEnd(nodeList));
    CHECK(*controlItr < nodeList.numNodes());
    CHECK(*ghostItr >= nodeList.firstGhostNode() && *ghostItr < nodeList.numNodes());
    const auto T = reflectOperator(position(*controlItr), position(*ghostItr));
    const auto& poly = field(*controlItr);
    vector<Vector> verts(poly.vertices());
    const auto& facets = poly.facetVertices();
    for (auto& v: verts) v = T*v;
    field(*ghostItr) = FacetedVolume(verts, facets);
  }
}

// Specialization for vector<scalar> fields, just perform a copy.
void
SphericalBoundary::
applyGhostBoundary(Field<Dim<3> , std::vector<Dim<3>::Scalar> >& field) const {

  // Apply the boundary condition to all the ghost node values.
  const NodeList<Dimension>& nodeList = field.nodeList();
  CHECK(controlNodes(nodeList).size() == ghostNodes(nodeList).size());
  auto controlItr = controlBegin(nodeList);
  auto ghostItr = ghostBegin(nodeList);
  for (; controlItr < controlEnd(nodeList); ++controlItr, ++ghostItr) {
    CHECK(ghostItr < ghostEnd(nodeList));
    CHECK(*controlItr < nodeList.numNodes());
    CHECK(*ghostItr >= nodeList.firstGhostNode() && *ghostItr < nodeList.numNodes());
    field(*ghostItr) = field(*controlItr);
  }
}

//------------------------------------------------------------------------------
// Find the set of nodes in the given NodeList that violate this boundary
// condition.  In this case we list all internal nodes as being in violation,
// in order to ensure nodes stay rigorously on the x-axis.
//------------------------------------------------------------------------------
void
SphericalBoundary::
setViolationNodes(NodeList<Dim<3> >& nodeList) {

  // Get the BoundaryNodes.violationNodes for this NodeList.
  BoundaryNodes& boundaryNodes = accessBoundaryNodes(nodeList);
  auto& vNodes = boundaryNodes.violationNodes;
  vNodes.clear();
  vNodes.reserve(nodeList.numInternalNodes());
  for (auto nodeID = 0u; nodeID != nodeList.numInternalNodes(); ++nodeID) 
    vNodes.push_back(nodeID);
  CHECK(vNodes.size() == nodeList.numInternalNodes());

  // Update the positions and H for the nodes in violation.
  updateViolationNodes(nodeList);
}

//------------------------------------------------------------------------------
// Update the nodes in violation of the boundary condition, mapping their
// positions and H tensors.
//------------------------------------------------------------------------------
void
SphericalBoundary::
updateViolationNodes(NodeList<Dim<3> >& nodeList) {

  // Get the set of violation nodes for this NodeList.
  const auto& vNodes = violationNodes(nodeList);

  // Loop over these nodes, and reset their positions to valid values.
  Field<Dim<3>, Vector>& positions = nodeList.positions();
  for (auto i: vNodes) {
    const double ri = positions(i).magnitude();
    positions(i) = Vector(ri, 0.0, 0.0);
  }

  // Set the Hfield.
  Field<Dim<3> , SymTensor>& Hfield = nodeList.Hfield();
  this->enforceBoundary(Hfield);

//   // Update the neighbor information.
//   nodeList.neighbor().updateNodes(); //(vNodes);
}    

//------------------------------------------------------------------------------
// Enforce consistency with the boundary condition on fields of different
// DataTypes.
//------------------------------------------------------------------------------
// Specialization for int fields, no-op.
void
SphericalBoundary::
enforceBoundary(Field<Dim<3> , int>&) const {
}

// Specialization for scalar fields, no-op.
void
SphericalBoundary::
enforceBoundary(Field<Dim<3> , Dim<3>::Scalar>&) const {
}

// Specialization for Vector fields, force radial along the x-axis.
void
SphericalBoundary::
enforceBoundary(Field<Dim<3>, Dim<3>::Vector>& field) const {

  const NodeList<Dimension>& nodeList = field.nodeList();
  for (auto itr = violationBegin(nodeList);
       itr != violationEnd(nodeList);
       ++itr) {
    CHECK(*itr < nodeList.numInternalNodes());
    field(*itr).y(0.0);
    field(*itr).z(0.0);
  }
}

// Specialization for Tensor fields, force all off-diagonal components
// to zero.
void
SphericalBoundary::
enforceBoundary(Field<Dim<3>, Dim<3>::Tensor>& field) const {

  const NodeList<Dimension>& nodeList = field.nodeList();
  for (auto itr = violationBegin(nodeList);
       itr != violationEnd(nodeList);
       ++itr) {
    CHECK(*itr < nodeList.numInternalNodes());
    const double xx = field(*itr).xx();
    const double yy = field(*itr).yy();
    const double zz = field(*itr).zz();
    field(*itr) = Tensor(xx, 0.0, 0.0,
                         0.0, yy, 0.0,
                         0.0, 0.0, zz);
  }
}

// Specialization for SymTensor fields, force all off-diagonal components
// to zero.
void
SphericalBoundary::
enforceBoundary(Field<Dim<3>, Dim<3>::SymTensor>& field) const {

  const NodeList<Dimension>& nodeList = field.nodeList();
  for (auto itr = violationBegin(nodeList);
       itr != violationEnd(nodeList);
       ++itr) {
    CHECK(*itr < nodeList.numInternalNodes());
    const double xx = field(*itr).xx();
    const double yy = field(*itr).yy();
    const double zz = field(*itr).zz();
    field(*itr) = SymTensor(xx, 0.0, 0.0,
                            0.0, yy, 0.0,
                            0.0, 0.0, zz);
  }
}

// Specialization for ThirdRankTensor fields.  Duh?
void
SphericalBoundary::
enforceBoundary(Field<Dim<3>, Dim<3>::ThirdRankTensor>&) const {
}

// Specialization for FourthRankTensor fields.  Duh?
void
SphericalBoundary::
enforceBoundary(Field<Dim<3>, Dim<3>::FourthRankTensor>&) const {
}

// Specialization for FifthRankTensor fields.  Duh?
void
SphericalBoundary::
enforceBoundary(Field<Dim<3>, Dim<3>::FifthRankTensor>&) const {
}

// Specialization for FacetedVolume fields.  Duh?
void
SphericalBoundary::
enforceBoundary(Field<Dim<3>, Dim<3>::FacetedVolume>&) const {
}

//------------------------------------------------------------------------------
// Dump the state to the given file.
//------------------------------------------------------------------------------
void
SphericalBoundary::
dumpState(FileIO& file, const string& pathName) const {
  file.write(mGhostPositions, pathName + "/ghostPositions");
}

//------------------------------------------------------------------------------
// Read the state from the given file.
//------------------------------------------------------------------------------
void
SphericalBoundary::
restoreState(const FileIO& file, const string& pathName) {
  file.read(mGhostPositions, pathName + "/ghostPositions");
}

}
