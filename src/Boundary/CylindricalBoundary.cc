//---------------------------------Spheral++----------------------------------//
// CylindricalBoundary -- Create a 3-D boundary around a set of nodes in the
// xy plane to emulate cylindrical (RZ) coordinates.
//
// Created by JMO, Tue Mar 29 13:48:03 PST 2005
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
#include "Geometry/GeomPlane.hh"
#include "Geometry/innerProduct.hh"
#include "NodeList/FluidNodeList.hh"
#include "Utilities/DBC.hh"
#include "RK/ReproducingKernelMethods.hh"

#include "CylindricalBoundary.hh"

using std::vector;
using std::string;

namespace Spheral {

//------------------------------------------------------------------------------
// Construct against the given DataBase.
//------------------------------------------------------------------------------
CylindricalBoundary::
CylindricalBoundary(const DataBase<Dim<3> >& dataBase):
  Boundary<Dim<3> >(),
  mDeltaPhi(dataBase.newGlobalFieldList(0.0, "Delta angle for generating ghosts")),
  mGhostPositions(dataBase.newGlobalFieldList(Dim<3>::Vector(), "Ghost node positions")),
  mRestart(registerWithRestart(*this)) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
CylindricalBoundary::
~CylindricalBoundary() {
}

//------------------------------------------------------------------------------
// Set the ghost nodes for the given NodeList.
//------------------------------------------------------------------------------
void
CylindricalBoundary::
setGhostNodes(NodeList<Dim<3> >& nodeList) {
  REQUIRE(mDeltaPhi.fieldForNodeList(nodeList) < mDeltaPhi.end());

  // Add this NodeList, creating space for control & ghost nodes.
  addNodeList(nodeList);
  BoundaryNodes& boundaryNodes = accessBoundaryNodes(nodeList);
  auto& controlNodes = boundaryNodes.controlNodes;
  auto& ghostNodes = boundaryNodes.ghostNodes;
  controlNodes.clear();
  ghostNodes.clear();

  // The fields of ghost positions and delta angles for this NodeList.
  Field<Dim<3>, Scalar>& deltaPhi = **(mDeltaPhi.fieldForNodeList(nodeList));

  // The ghost node indices.
  const int currentNumGhostNodes = nodeList.numGhostNodes();
  int ghostNodeIndex = nodeList.numNodes();

  // Get the state for this NodeList.
  Field<Dimension, Vector>& positions = nodeList.positions();
  Field<Dimension, SymTensor>& H = nodeList.Hfield();

  // Get the number of nodes per h, and the kernel extent.
  const double nPerh = nodeList.nodesPerSmoothingScale();
  const double kernelExtent = nodeList.neighbor().kernelExtent();
  CHECK(distinctlyGreaterThan(nPerh, 0.0));
  CHECK(distinctlyGreaterThan(kernelExtent, 0.0));

  // The vector of ghost positions we'll build up, to assign all at once at the
  // end.
  vector<Vector> ghostPositions;

  // Iterate over the internal nodes.
  for (auto i = 0u; i != nodeList.numNodes(); ++i) {

    // The state of the node.
    const double ri = positions(i).y();
    const double zi = positions(i).x();
    const SymTensor& Hi = H(i);
    const Scalar hzi = 1.0/(Hi*Vector(0.0, 0.0, 1.0)).magnitude();
    CHECK(ri >= 0.0);

    // Have we computed the set of ghost nodes for this node yet?
    if (deltaPhi(i) == 0.0) {
      deltaPhi(i) = angularSpacing(ri, hzi, nPerh, kernelExtent);
    }
    const double dphi = deltaPhi(i);
    CHECK(dphi > 0.0);

    // Create the ghost nodes positions.
    const int nbins = max(1, int(min(M_PI - dphi, 1.5*kernelExtent*hzi/(ri + 1.0e-50))/dphi + 0.5));
    if (!(nbins >= 1 && nbins*dphi < M_PI)) {
      cerr << nbins << " "
           << dphi <<  " "
           << hzi << " "
           << ri << " "
           << endl;
    }
    CHECK(nbins >= 1 && nbins*dphi < M_PI);
    for (int ibin = 0; ibin != nbins; ++ibin) {
      const double phij = (ibin + 1)*dphi;
      controlNodes.push_back(i);
      ghostNodes.push_back(ghostNodeIndex++);
      ghostPositions.push_back(Vector(zi,
                                      ri*cos(phij),
                                      ri*sin(phij)));
      controlNodes.push_back(i);
      ghostNodes.push_back(ghostNodeIndex++);
      ghostPositions.push_back(Vector(zi,
                                      ri*cos(phij),
                                      -ri*sin(phij)));
    }

  }

  // Create the set of ghost nodes on the NodeList.
  const int numNewGhostNodes = ghostNodes.size();
  CHECK((int)controlNodes.size() == numNewGhostNodes);
  CHECK((int)ghostPositions.size() == numNewGhostNodes);
  nodeList.numGhostNodes(currentNumGhostNodes + numNewGhostNodes);
  CHECK(ghostNodeIndex == (int)nodeList.numNodes());
  Field<Dimension, Vector>& storedGhostPositions = **(mGhostPositions.fieldForNodeList(nodeList));
  for (int k = 0; k != numNewGhostNodes; ++k) {
    const int i = ghostNodes[k];
    CHECK(i >= (int)nodeList.firstGhostNode() && i < (int)nodeList.numNodes());
    positions(i) = ghostPositions[k];
    storedGhostPositions(i) = ghostPositions[k];
  }

  // Contract checking...
  BEGIN_CONTRACT_SCOPE
  {
    auto controlItr = controlNodes.begin();
    auto ghostItr = ghostNodes.begin();
    for (; controlItr != controlNodes.end(); ++controlItr, ++ghostItr) {
      CHECK(ghostItr < ghostNodes.end());
      const auto i = *controlItr;
      const auto j = *ghostItr;
      const auto ri = positions(i).y();
      const auto zi = positions(i).x();
      CONTRACT_VAR(j);
      CONTRACT_VAR(ri);
      CONTRACT_VAR(zi);
      CHECK(fuzzyEqual(positions(j).x(), zi));
      CHECK(fuzzyEqual(sqrt(FastMath::square(positions(j).y()) + FastMath::square(positions(j).z())), ri));
    }
  }
  END_CONTRACT_SCOPE

  // We can use the normal ghost boundary enforcement to update the
  // H's of the ghost nodes.
  applyGhostBoundary(H);

//   // Update the neighbor information.
//   nodeList.neighbor().updateNodes(); // (ghostNodes);
}

//------------------------------------------------------------------------------
// Update the state of the current ghost nodes (correcting their positions
// and H tensors.)
//------------------------------------------------------------------------------
void
CylindricalBoundary::
updateGhostNodes(NodeList<Dim<3> >& nodeList) {
  REQUIRE(mDeltaPhi.fieldForNodeList(nodeList) < mDeltaPhi.end());

  // Get the control and ghost nodes for this NodeList.
  const BoundaryNodes& boundaryNodes = accessBoundaryNodes(nodeList);
  const auto& controlNodes = boundaryNodes.controlNodes;
  const auto& ghostNodes = boundaryNodes.ghostNodes;

  // Get the positions for this NodeList.
  Field<Dimension, Vector>& positions = nodeList.positions();
  Field<Dimension, Vector>& ghostPositions = **(mGhostPositions.fieldForNodeList(nodeList));

  // Update the positions, forcing the ghost nodes to exist in a ring of 
  // appropriate radius based on the control nodes.
  CHECK(controlNodes.size() == ghostNodes.size());
  auto controlItr = controlNodes.begin();
  auto ghostItr = ghostNodes.begin();
  for (; controlItr != controlNodes.end(); ++controlItr, ++ghostItr) {
    CHECK(ghostItr < ghostNodes.end());
    const int i = *controlItr;
    const int j = *ghostItr;
    const double ri = positions(i).y();
    const double zi = positions(i).x();
    const double f = ri/sqrt(ghostPositions(j).y()*ghostPositions(j).y() + 
                             ghostPositions(j).z()*ghostPositions(j).z() + 1.0e-50);
    positions(j) = Vector(zi,
                          f*ghostPositions(j).y(),
                          f*ghostPositions(j).z());
    ghostPositions(j) = positions(j);
  }
  CHECK(controlItr == controlNodes.end());
  CHECK(ghostItr == ghostNodes.end());

  // Contract checking...
  BEGIN_CONTRACT_SCOPE
  {
    auto controlItr = controlNodes.begin();
    auto ghostItr = ghostNodes.begin();
    for (; controlItr != controlNodes.end(); ++controlItr, ++ghostItr) {
      CHECK(ghostItr < ghostNodes.end());
      const int i = *controlItr;
      const int j = *ghostItr;
      const double ri = positions(i).y();
      const double zi = positions(i).x();
      CONTRACT_VAR(j);
      CONTRACT_VAR(ri);
      CONTRACT_VAR(zi);
      CHECK(fuzzyEqual(positions(j).x(), zi));
      CHECK(fuzzyEqual(sqrt(FastMath::square(positions(j).y()) + FastMath::square(positions(j).z())), ri));
    }
  }
  END_CONTRACT_SCOPE

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
// Specialization for Vector fields.
void
CylindricalBoundary::
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
CylindricalBoundary::
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
CylindricalBoundary::
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
CylindricalBoundary::
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
CylindricalBoundary::
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
CylindricalBoundary::
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
CylindricalBoundary::
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

// Specialization for RKCoefficients fields
void
CylindricalBoundary::
applyGhostBoundary(Field<Dim<3>, RKCoefficients<Dim<3>>>& field) const {
  const auto& nodeList = field.nodeList();
  CHECK(controlNodes(nodeList).size() == ghostNodes(nodeList).size());
  if (not ghostNodes(nodeList).empty()) {

    // Figure out the order of the corrections
    ReproducingKernelMethods<Dim<3>> WR(RKOrder::ZerothOrder);
    switch (field[0].size()) {
    case RKUtilities<Dimension, RKOrder::ZerothOrder>::gradPolynomialSize:
      break;
    case RKUtilities<Dimension, RKOrder::LinearOrder>::gradPolynomialSize:
      WR = ReproducingKernelMethods<Dim<3>>(RKOrder::LinearOrder);
      break;
    case RKUtilities<Dimension, RKOrder::QuadraticOrder>::gradPolynomialSize:
      WR = ReproducingKernelMethods<Dim<3>>(RKOrder::QuadraticOrder);
      break;
    case RKUtilities<Dimension, RKOrder::CubicOrder>::gradPolynomialSize:
      WR = ReproducingKernelMethods<Dim<3>>(RKOrder::CubicOrder);
      break;
    case RKUtilities<Dimension, RKOrder::QuarticOrder>::gradPolynomialSize:
      WR = ReproducingKernelMethods<Dim<3>>(RKOrder::QuarticOrder);
      break;
    case RKUtilities<Dimension, RKOrder::QuinticOrder>::gradPolynomialSize:
      WR = ReproducingKernelMethods<Dim<3>>(RKOrder::QuinticOrder);
      break;
    case RKUtilities<Dimension, RKOrder::SexticOrder>::gradPolynomialSize:
      WR = ReproducingKernelMethods<Dim<3>>(RKOrder::SexticOrder);
      break;
    case RKUtilities<Dimension, RKOrder::SepticOrder>::gradPolynomialSize:
      WR = ReproducingKernelMethods<Dim<3>>(RKOrder::SepticOrder);
      break;
    default:
      VERIFY2(false, "Cylindrical boundary ERROR: unknown order for RKCoefficients");
    };

    const auto& position = nodeList.positions();
    auto controlItr = controlBegin(nodeList);
    auto ghostItr = ghostBegin(nodeList);
    for (; controlItr < controlEnd(nodeList); ++controlItr, ++ghostItr) {
      CHECK(ghostItr < ghostEnd(nodeList));
      CHECK(*controlItr < nodeList.numNodes());
      CHECK(*ghostItr >= nodeList.firstGhostNode() && *ghostItr < nodeList.numNodes());
      field(*ghostItr) = field(*controlItr);
      const auto R = reflectOperator(position(*controlItr), position(*ghostItr));
      const auto T = WR.transformationMatrix(R, false);
      WR.applyTransformation(T, field(*ghostItr));
    }
  }
}

//------------------------------------------------------------------------------
// Find the set of nodes in the given NodeList that violate this boundary
// condition.  In this case we list all internal nodes as being in violation,
// in order to ensure nodes stay rigorously in the xy plane.
//------------------------------------------------------------------------------
void
CylindricalBoundary::
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
CylindricalBoundary::
updateViolationNodes(NodeList<Dim<3> >& nodeList) {

  // Get the set of violation nodes.
  const auto& vNodes = violationNodes(nodeList);

  // State we care about.
  Field<Dim<3>, Vector>& positions = nodeList.positions();
  Field<Dim<3>, Vector>& velocity = nodeList.velocity();
  Field<Dim<3> , SymTensor>& Hfield = nodeList.Hfield();

  // Make sure the positions are sensible.
  for (auto i: vNodes) {

    // Has the node crossed to negative r (y)?
    const double ri = positions(i).y();
    if (ri < 0.0) {
      positions(i).y(-ri);
      velocity(i).y(abs(velocity(i).y()));
    }
    CHECK(positions(i).y() >= 0.0);

    // Make sure there is no component out of the xy plane.
    positions(i).z(0.0);

  }
    
  // The normal Field enforcement policies work H.
  this->enforceBoundary(Hfield);

//   // Update the neighbor information.
//   nodeList.neighbor().updateNodes(); // (vNodes);
}    

//------------------------------------------------------------------------------
// Enforce consistency with the boundary condition on fields of different
// DataTypes.
//------------------------------------------------------------------------------
// Specialization for Vector fields, forcing Vector's to lie in the xy plane.
void
CylindricalBoundary::
enforceBoundary(Field<Dim<3>, Dim<3>::Vector>& field) const {

  const NodeList<Dimension>& nodeList = field.nodeList();
  for (auto itr = violationBegin(nodeList);
       itr != violationEnd(nodeList);
       ++itr) {
    CHECK(*itr < nodeList.numInternalNodes());
    field(*itr).z(0.0);
  }
}

// Specialization for Tensor fields, force all off-diagonal components
// to zero.
void
CylindricalBoundary::
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
CylindricalBoundary::
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


//------------------------------------------------------------------------------
// Compute the target azimuthal angular spacing.
//------------------------------------------------------------------------------
double
CylindricalBoundary::
angularSpacing(const double ri,
               const double hzi,
               const double nodesPerh,
               const double kernelExtent) {
  REQUIRE(ri >= 0.0);
  REQUIRE(hzi >= 0.0);
  REQUIRE(nodesPerh > 0.0);
  REQUIRE(kernelExtent > 0.0);

  // Choose a minimum effective radius such based on whether we see more than
  // 1/4*Pi around the circumference of the ring.
  const double rmin = kernelExtent*hzi/(0.125*M_PI);
  CONTRACT_VAR(rmin);
  CHECK(rmin > 0.0);
  const double reff = ri; // max(rmin, ri);

  // Compute the nominal spacing.
  const double delta = hzi/nodesPerh;
  double dphi = 2.0*atan2(delta, 2.0*reff);
  CHECK(dphi > 0.0);

  // Now find the closest value of dphi that will give us an integer + 1/2
  // number of steps in Pi.
  int twoi = max(2, int(2.0*M_PI/dphi));
  if (twoi % 2 != 0) ++twoi;
  CHECK(twoi % 2 == 0);
  const double result = 2.0*M_PI/(twoi + 1);
  ENSURE(result > 0.0);
  return result;
}

//------------------------------------------------------------------------------
// Dump the state to the given file.
//------------------------------------------------------------------------------
void
CylindricalBoundary::
dumpState(FileIO& file, const string& pathName) const {
  file.write(mDeltaPhi, pathName + "/deltaPhi");
  file.write(mGhostPositions, pathName + "/ghostPositions");
}

//------------------------------------------------------------------------------
// Read the state from the given file.
//------------------------------------------------------------------------------
void
CylindricalBoundary::
restoreState(const FileIO& file, const string& pathName) {
  file.read(mDeltaPhi, pathName + "/deltaPhi");
  file.read(mGhostPositions, pathName + "/ghostPositions");
}

}
