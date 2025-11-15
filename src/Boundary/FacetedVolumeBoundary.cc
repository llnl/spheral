//---------------------------------Spheral++----------------------------------//
// FacetedVolumeBoundary -- Apply a FacetedVolume boundary condition to Spheral++
// Fields.
//----------------------------------------------------------------------------//
#include "findNodesTouchingThroughPlanes.hh"
#include "mapPositionThroughPlanes.hh"
#include "FileIO/FileIO.hh"
#include "Geometry/GeomPlane.hh"
#include "Geometry/innerProduct.hh"
#include "Field/Field.hh"
#include "Utilities/DBC.hh"
#include "Utilities/planarReflectingOperator.hh"
#include "Mesh/Mesh.hh"

#include "FacetedVolumeBoundary.hh"

using std::vector;
using std::pair;
using std::make_pair;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {

namespace {

//------------------------------------------------------------------------------
// Construct the plane associated with a Facet
//------------------------------------------------------------------------------
// 2D
GeomPlane<Dim<2>>
facetPlane(const GeomFacet2d& facet, const bool interiorBoundary) {
  const auto& p1 = facet.point1();
  const auto  nhat = (interiorBoundary ?
                       facet.normal() :
                      -facet.normal());
  return GeomPlane<Dim<2>>(p1, nhat);
}

//..............................................................................
// 3D
GeomPlane<Dim<3>>
facetPlane(const GeomFacet3d& facet, const bool interiorBoundary) {
  const auto& p1 = facet.point(0);
  const auto  nhat = (interiorBoundary ?
                       facet.normal() :
                      -facet.normal());
  return GeomPlane<Dim<3>>(p1, nhat);
}

//------------------------------------------------------------------------------
// Construct the planes bounding the unique facet volume for ghosts
//------------------------------------------------------------------------------
// 2D
std::vector<GeomPlane<Dim<2>>>
facetGhostPlanes(const GeomPolygon& poly,
                 const GeomFacet2d& facet,
                 const bool interiorBoundary) {
  typedef Dim<2>::Vector Vector;
  typedef GeomPlane<Dim<2>> Plane;
  const auto& centroid = poly.centroid();
  const auto& p1 = facet.point1();
  const auto& p2 = facet.point2();
  const auto& vhat1 = (p1 - centroid).unitVector();
  const auto& vhat2 = (p2 - centroid).unitVector();
  return std::vector<Plane>({Plane(p1, Vector(-vhat1.y(),  vhat1.x())),
                             Plane(p2, Vector( vhat2.y(), -vhat2.x())),
                             facetPlane(facet, not interiorBoundary)});
}

// 3D
std::vector<GeomPlane<Dim<3>>>
facetGhostPlanes(const GeomPolyhedron& poly,
                 const GeomFacet3d& facet,
                 const bool interiorBoundary) {
  typedef GeomPlane<Dim<3>> Plane;
  std::vector<Plane> result;
  const auto& centroid = poly.centroid();
  const auto& coords = poly.vertices();
  const auto& iverts = facet.ipoints();
  const auto  nverts = iverts.size();
  for (auto k = 0u; k < nverts; ++k) {
    const auto& p1 = coords[iverts[k]];
    const auto& p2 = coords[iverts[(k + 1) % nverts]];
    result.push_back(Plane(centroid, ((p1 - centroid).cross(p2 - centroid)).unitVector()));
  }
  result.push_back(facetPlane(facet, not interiorBoundary));
  return result;
}

//------------------------------------------------------------------------------
// Check if a point is inside the volume defined by a set of planes
//------------------------------------------------------------------------------
template<typename Vector, typename Plane>
bool
contained(const Vector& p, const std::vector<Plane>& planes) {
  for (const auto& plane: planes) {
    if (plane > p) return false;
  }
  return true;
}

//------------------------------------------------------------------------------
// Find the node touching a given facet
//------------------------------------------------------------------------------
// Polygon
std::vector<int>
nodesTouchingFacet(const NodeList<Dim<2>>& nodes,
                   const GeomFacet2d& facet,
                   const bool interiorBoundary) {
  const auto& p1 = facet.point1();
  const auto& p2 = facet.point2();
  const auto  phat = (p2 - p1).unitVector();
  const auto  smax = (p2 - p1).magnitude();
  const auto  plane = facetPlane(facet, interiorBoundary);
  const auto potentials = findNodesTouchingThroughPlanes(nodes, plane, plane);
  std::vector<int> result;
  const auto& pos = nodes.positions();
  for (const auto i: potentials) {
    const auto s = (pos(i) - p1).dot(phat);
    if (s >= 0.0 and s <= smax) result.push_back(i);
  }
  return result;
}

//..............................................................................
// Polyhedron
std::vector<int>
nodesTouchingFacet(const NodeList<Dim<3>>&,
                   const GeomFacet3d&,
                   const bool /*interiorBoundary*/) {
  std::vector<int> result;
  return result;
}

//------------------------------------------------------------------------------
// Copy control->ghost values
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
void
copyControlValues(Field<Dimension, Value>& field,
                  const std::vector<size_t>& controls,
                  const std::vector<size_t>& ghosts) {
  REQUIRE(controls.size() == ghosts.size());
  auto controlItr = controls.begin();
  auto ghostItr = ghosts.begin();
  for (; controlItr < controls.end(); ++controlItr, ++ghostItr) {
    field(*ghostItr) = field(*controlItr);
  }
}

//------------------------------------------------------------------------------
// The per value type method of mapping values
//------------------------------------------------------------------------------
// Vector
template<typename Dimension>
void
mapValue(typename Dimension::Vector& ghost,
         const typename Dimension::Vector& control,
         const typename Dimension::Tensor& R) {
  ghost = R*control;
}

// Tensor
template<typename Dimension>
void
mapValue(typename Dimension::Tensor& ghost,
         const typename Dimension::Tensor& control,
         const typename Dimension::Tensor& R) {
  ghost = R*control*R;
}

// SymTensor
template<typename Dimension>
void
mapValue(typename Dimension::SymTensor& ghost,
         const typename Dimension::SymTensor& control,
         const typename Dimension::Tensor& R) {
  ghost = (R*control*R).Symmetric();
}

// ThirdRankTensor
template<typename Dimension>
void
mapValue(typename Dimension::ThirdRankTensor& ghost,
         const typename Dimension::ThirdRankTensor& control,
         const typename Dimension::Tensor& R) {
  ghost.Zero();
  for (unsigned i = 0; i != Dimension::nDim; ++i) {
    for (unsigned j = 0; j != Dimension::nDim; ++j) {
      for (unsigned k = 0; k != Dimension::nDim; ++k) {
        for (unsigned q = 0; q != Dimension::nDim; ++q) {
          for (unsigned r = 0; r != Dimension::nDim; ++r) {
            for (unsigned s = 0; s != Dimension::nDim; ++s) {
              ghost(i,j,k) += R(i,q)*R(j,r)*R(k,s)*control(q,r,s);
            }
          }
        }
      }
    }
  }
}

// FourthRankTensor
template<typename Dimension>
void
mapValue(typename Dimension::FourthRankTensor& ghost,
         const typename Dimension::FourthRankTensor& control,
         const typename Dimension::Tensor& R) {
  ghost.Zero();
  for (unsigned i = 0; i != Dimension::nDim; ++i) {
    for (unsigned j = 0; j != Dimension::nDim; ++j) {
      for (unsigned k = 0; k != Dimension::nDim; ++k) {
        for (unsigned l = 0; l != Dimension::nDim; ++l) {
          for (unsigned q = 0; q != Dimension::nDim; ++q) {
            for (unsigned r = 0; r != Dimension::nDim; ++r) {
              for (unsigned s = 0; s != Dimension::nDim; ++s) {
                for (unsigned t = 0; t != Dimension::nDim; ++t) {
                  ghost(i,j,k,l) += R(i,q)*R(j,r)*R(k,s)*R(l,t)*control(q,r,s,t);
                }
              }
            }
          }
        }
      }
    }
  }
}

// FifthRankTensor
template<typename Dimension>
void
mapValue(typename Dimension::FifthRankTensor& ghost,
         const typename Dimension::FifthRankTensor& control,
         const typename Dimension::Tensor& R) {
  ghost.Zero();
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
                      ghost(i,j,k,l,u) += R(i,q)*R(j,r)*R(k,s)*R(l,t)*R(m,u)*control(q,r,s,t,u);
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
}

//------------------------------------------------------------------------------
// Map control->ghost values, depending on the Value type
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
void
mapControlValues(Field<Dimension, Value>& field,
                 const std::vector<typename Dimension::Tensor>& reflectOperators,
                 const std::map<std::string, std::vector<std::vector<int>>>& facetControlNodes,
                 const std::map<std::string, std::vector<std::pair<int,int>>>& facetGhostNodes) {
  const auto& name = field.nodeList().name();
  const auto& controlNodes = facetControlNodes.find(name)->second;
  const auto& ghostNodeRanges = facetGhostNodes.find(name)->second;
  const auto  nfacets = reflectOperators.size();
  REQUIRE(controlNodes.size() == nfacets);
  REQUIRE(ghostNodeRanges.size() == nfacets);
  for (auto ifacet = 0u; ifacet < nfacets; ++ifacet) {
    const auto& R = reflectOperators[ifacet];
    const auto& controls = controlNodes[ifacet];
    auto        ghostID = ghostNodeRanges[ifacet].first;
    CHECK(ghostNodeRanges[ifacet].second - ghostID == (int)controls.size());
    for (const auto i: controls) {
      mapValue<Dimension>(field(ghostID), field(i), R);
      ++ghostID;
    }
    CHECK(ghostID == ghostNodeRanges[ifacet].second);
  }
}

//------------------------------------------------------------------------------
// Map control->ghost values, depending on the Value type
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
void
enforceViolationValues(Field<Dimension, Value>& field,
                       const std::vector<size_t>& violationNodes,
                       const std::vector<typename Dimension::Tensor>& violationOps) {
  const auto n = violationNodes.size();
  REQUIRE(violationOps.size() == n);
  Value newVal;
  for (auto k = 0u; k < n; ++k) {
    const auto i = violationNodes[k];
    mapValue<Dimension>(newVal, field(i), violationOps[k]);
    field(i) = newVal;
  }
}

}

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
FacetedVolumeBoundary<Dimension>::FacetedVolumeBoundary(const FacetedVolume& poly,
                                                        const bool interiorBoundary,
                                                        const bool useGhosts):
  Boundary<Dimension>(),
  mPoly(poly),
  mInteriorBoundary(interiorBoundary),
  mUseGhosts(useGhosts),
  mReflectOperators(),
  mFacetControlNodes(),
  mFacetGhostNodes(),
  mViolationOperators() {

  // Build the reflection operators for each facet of the poly.
  const auto& facets = poly.facets();
  for (const auto& facet: facets) {
    mReflectOperators.push_back(interiorBoundary ?
                                planarReflectingOperator<Dimension>(facet.normal().unitVector()) :
                                planarReflectingOperator<Dimension>(-facet.normal().unitVector()));
  }
  ENSURE(mReflectOperators.size() == facets.size());

}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
FacetedVolumeBoundary<Dimension>::~FacetedVolumeBoundary() {
}

//------------------------------------------------------------------------------
// Create ghost nodes for the NodeList
//------------------------------------------------------------------------------
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::setGhostNodes(NodeList<Dimension>& nodeList) {
  this->addNodeList(nodeList);
  if (mUseGhosts) {
    const auto& facets = mPoly.facets();
    const auto  nfacets = facets.size();
    const auto  name = nodeList.name();
    auto& boundaryNodes = this->accessBoundaryNodes(nodeList);
    boundaryNodes.controlNodes.clear();
    boundaryNodes.ghostNodes.clear();
    auto&       pos = nodeList.positions();
    auto&       H = nodeList.Hfield();

    // Allocate storage for controls and ghosts
    mFacetControlNodes[name] = vector<vector<int>>(nfacets);
    mFacetGhostNodes[name] = vector<pair<int,int>>(nfacets);
    auto& controls = mFacetControlNodes[name];
    auto& ghostRanges = mFacetGhostNodes[name];

    // Build the control points.
    auto firstGhost = nodeList.numNodes();
    vector<Vector> posGhost;
    vector<SymTensor> Hghost;
    for (auto k = 0u; k < nfacets; ++k) {
      const auto& facet = facets[k];
      const auto  boundPlanes = facetGhostPlanes(mPoly, facet, mInteriorBoundary);
      const auto  plane = facetPlane(facet, mInteriorBoundary);
      const auto& R = mReflectOperators[k];
      ghostRanges[k].first = firstGhost;
      const auto  potentials = nodesTouchingFacet(nodeList, facet, mInteriorBoundary);
      for (const auto i: potentials) {
        const auto posj = mapPositionThroughPlanes(pos(i), plane, plane);
        if (contained(posj, boundPlanes)) {
          controls[k].push_back(i);
          firstGhost += 1;
          posGhost.push_back(posj);
          Hghost.push_back((R*H(i)*R).Symmetric());
        }
      }
      ghostRanges[k].second = firstGhost;
      CHECK(posGhost.size() == Hghost.size());
      boundaryNodes.controlNodes.insert(boundaryNodes.controlNodes.end(),
                                        controls[k].begin(),
                                        controls[k].end());
    }
    CHECK(mFacetControlNodes[name].size() == facets.size());
    CHECK(mFacetGhostNodes[name].size() == facets.size());
    CHECK(posGhost.size() == Hghost.size());
    CHECK(posGhost.size() == boundaryNodes.controlNodes.size());

    // Update the ghost node info.
    const auto numNewGhosts = posGhost.size();
    firstGhost = nodeList.numNodes();
    nodeList.numGhostNodes(nodeList.numGhostNodes() + numNewGhosts);
    std::copy(posGhost.begin(), posGhost.end(), &(pos[firstGhost]));
    std::copy(Hghost.begin(), Hghost.end(), &(H[firstGhost]));
    boundaryNodes.ghostNodes.resize(numNewGhosts);
    for (auto k = 0u; k < numNewGhosts; ++k) boundaryNodes.ghostNodes[k] = k + firstGhost;
    CHECK(boundaryNodes.controlNodes.size() == boundaryNodes.ghostNodes.size());
    // std::cerr << "Introduced ghosts in range: " << firstGhost << " " << nodeList.numNodes() << std::endl;
  }
}

//------------------------------------------------------------------------------
// Set the minimal field values necessary to make the ghost nodes valid.
//------------------------------------------------------------------------------
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::updateGhostNodes(NodeList<Dimension>& nodeList) {
  if (mUseGhosts) {

    const auto  name = nodeList.name();
    const auto& controls = mFacetControlNodes[name];
    const auto& ghostRanges = mFacetGhostNodes[name];
    const auto& facets = mPoly.facets();
    const auto  nfacets = facets.size();
    auto&       pos = nodeList.positions();
    auto&       H = nodeList.Hfield();
    CHECK(controls.size() == nfacets);
    CHECK(ghostRanges.size() == nfacets);

    // Walk each of the facets
    for (auto ifacet = 0u; ifacet < nfacets; ++ifacet) {
      const auto  n = controls[ifacet].size();
      const auto  plane = facetPlane(facets[ifacet], mInteriorBoundary);
      const auto& R = mReflectOperators[ifacet];
      CHECK(ghostRanges[ifacet].second - ghostRanges[ifacet].first == (int)n);
      for (auto k = 0u; k < n; ++k) {
        const auto i = controls[ifacet][k];
        const auto j = ghostRanges[ifacet].first + k;
        pos(j) = mapPositionThroughPlanes(pos(i), plane, plane);
        H(j) = (R*H(i)*R).Symmetric();
      }
    }
  }
}

//------------------------------------------------------------------------------
// Apply the ghost boundary condition to fields of different DataTypes.
//------------------------------------------------------------------------------
// Specialization for int fields, just perform a copy.
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, int>& field) const {
  if (mUseGhosts) {
    const auto& nodeList = field.nodeList();
    copyControlValues(field, this->controlNodes(nodeList), this->ghostNodes(nodeList));
  }
}

// Specialization for scalar fields, just perform a copy.
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::Scalar>& field) const {
  if (mUseGhosts) {
    const auto& nodeList = field.nodeList();
    copyControlValues(field, this->controlNodes(nodeList), this->ghostNodes(nodeList));
  }
}

// Specialization for Vector fields.
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::Vector>& field) const {
  if (mUseGhosts) {
    mapControlValues(field, mReflectOperators, mFacetControlNodes, mFacetGhostNodes);
  }
}

// Specialization for Tensor fields.
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::Tensor>& field) const {
  if (mUseGhosts) {
    mapControlValues(field, mReflectOperators, mFacetControlNodes, mFacetGhostNodes);
  }
}

// Specialization for symmetric tensors.
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::SymTensor>& field) const {
  if (mUseGhosts) {
    mapControlValues(field, mReflectOperators, mFacetControlNodes, mFacetGhostNodes);
  }
}

// Specialization for ThirdRankTensor fields.
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::ThirdRankTensor>& field) const {
  if (mUseGhosts) {
    mapControlValues(field, mReflectOperators, mFacetControlNodes, mFacetGhostNodes);
  }
}

// Specialization for FourthRankTensor fields.
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::FourthRankTensor>& field) const {
  if (mUseGhosts) {
    mapControlValues(field, mReflectOperators, mFacetControlNodes, mFacetGhostNodes);
  }
}

// Specialization for FifthRankTensor fields.
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::FifthRankTensor>& field) const {
  if (mUseGhosts) {
    mapControlValues(field, mReflectOperators, mFacetControlNodes, mFacetGhostNodes);
  }
}

// Specialization for FacetedVolumes
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::FacetedVolume>& field) const {
  if (mUseGhosts) {
    const auto& nodeList = field.nodeList();
    copyControlValues(field, this->controlNodes(nodeList), this->ghostNodes(nodeList));   // Punt for FacetedVolumes and just copy them for now
  }
}

//------------------------------------------------------------------------------
// Find the set of nodes in the given NodeList that violate this boundary
// condition.
//------------------------------------------------------------------------------
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::setViolationNodes(NodeList<Dimension>& nodeList) {

  // Get the BoundaryNodes.violationNodes for this NodeList.
  this->addNodeList(nodeList);
  auto& boundaryNodes = this->accessBoundaryNodes(nodeList);
  auto& vNodes = boundaryNodes.violationNodes;
  vNodes.clear();

  // Loop over all the internal nodes in the NodeList and look for any that have
  // wandered into the excluded region.
  const auto& pos = nodeList.positions();
  const auto  n = nodeList.numInternalNodes();
  for (auto i = 0u; i < n; ++i) {
    if ((not mInteriorBoundary) xor mPoly.contains(pos(i), false)) vNodes.push_back(i);
  }

  // std::cerr << "Violation nodes:";
  // for (auto i: vNodes) std::cerr << " " << i;
  // std::cerr << endl;

  // Update the positions and H for the nodes in violation.
  updateViolationNodes(nodeList);
}

//------------------------------------------------------------------------------
// Update the nodes in violation of the boundary condition, mapping their
// positions and H tensors.
//------------------------------------------------------------------------------
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::updateViolationNodes(NodeList<Dimension>& nodeList) {

  // Get the set of violation nodes for this NodeList.
  const auto& vNodes = this->violationNodes(nodeList);
  const auto  numViolation = vNodes.size();
  auto&       violationOps = mViolationOperators[nodeList.name()];
  violationOps = vector<Tensor>(vNodes.size(), Tensor::one());
  const auto& facets = mPoly.facets();

  // Find the longest scale in the FacetedVolume
  const auto chordLength = (mPoly.xmax() - mPoly.xmin()).magnitude();

  // Loop over these nodes, and reset their positions to valid values.
  auto&       pos = nodeList.positions();
  const auto& vel = nodeList.velocity();
  vector<unsigned> potentialFacets;
  vector<Vector> intersections;
  Vector newPos, newVel;
  bool inViolation;
  int iter = 0;
  int minFacet;
  Scalar minDist;
  const int maxIters = 5;
  for (auto k = 0u; k < numViolation; ++k) {
    auto  i = vNodes[k];
    auto& R = violationOps[k];
    newPos = pos(i);
    newVel = vel(i);
    inViolation = true;
    iter = 0;
    while (inViolation and iter++ < maxIters) {
      // Backtrack to which facet we think the point passed through.
      CHECK((not mInteriorBoundary) xor mPoly.contains(newPos));
      const auto backPos = newPos - chordLength*newVel.unitVector();
      mPoly.intersections(backPos, newPos, potentialFacets, intersections);
      CHECK(potentialFacets.size() > 0);
      CHECK(potentialFacets.size() == intersections.size());
      minFacet = potentialFacets[0];
      minDist = (intersections[0] - newPos).magnitude2();
      for (auto kk = 1u; kk < intersections.size(); ++kk) {
        if ((intersections[kk] - newPos).magnitude2() < minDist) {
          minFacet = potentialFacets[kk];
          minDist = (intersections[kk] - newPos).magnitude2();
        }
      }
      const auto plane = facetPlane(facets[minFacet], mInteriorBoundary);
      R *= mReflectOperators[minFacet];
      newPos = mapPositionThroughPlanes(newPos, plane, plane);
      newVel = mReflectOperators[minFacet]*newVel;
      inViolation = ((not mInteriorBoundary) xor mPoly.contains(newPos, false));
    }
    pos(i) = newPos;
  }

  // Set the Hfield.
  auto& Hfield = nodeList.Hfield();
  this->enforceBoundary(Hfield);
}

//------------------------------------------------------------------------------
// Enforce the boundary condition on the set of nodes in violation of the 
// boundary.
//------------------------------------------------------------------------------
// Specialization for int fields.  A no-op.
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::
enforceBoundary(Field<Dimension, int>&) const {
}

// Specialization for scalar fields.  A no-op.
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::Scalar>&) const {
}

// Specialization for vector fields.  Apply the reflection operator.
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::Vector>& field) const {
  enforceViolationValues(field, this->violationNodes(field.nodeList()), mViolationOperators.find(field.nodeList().name())->second);
}

// Specialization for tensor fields.  Apply the reflection operator.
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::Tensor>& field) const {
  enforceViolationValues(field, this->violationNodes(field.nodeList()), mViolationOperators.find(field.nodeList().name())->second);
}

// Specialization for tensor fields.  Apply the reflection operator.
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::SymTensor>& field) const {
  enforceViolationValues(field, this->violationNodes(field.nodeList()), mViolationOperators.find(field.nodeList().name())->second);
}

// Specialization for third rank tensor fields.  Apply the reflection operator.
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::ThirdRankTensor>& field) const {
  enforceViolationValues(field, this->violationNodes(field.nodeList()), mViolationOperators.find(field.nodeList().name())->second);
}

// Specialization for fourth rank tensor fields.  Apply the reflection operator.
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::FourthRankTensor>& field) const {
  enforceViolationValues(field, this->violationNodes(field.nodeList()), mViolationOperators.find(field.nodeList().name())->second);
}

// Specialization for fifth rank tensor fields.  Apply the reflection operator.
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::FifthRankTensor>& field) const {
  enforceViolationValues(field, this->violationNodes(field.nodeList()), mViolationOperators.find(field.nodeList().name())->second);
}

// Specialization for FacetedVolumes
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::FacetedVolume>&) const {
  // Punt on FacetedVolumes for now.
}

//------------------------------------------------------------------------------
// Clear out any NodeList information that is currently present.
//------------------------------------------------------------------------------
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::reset(const DataBase<Dimension>& dataBase) {
  Boundary<Dimension>::reset(dataBase);
  mFacetControlNodes.clear();
  mFacetGhostNodes.clear();
}

//------------------------------------------------------------------------------
// Cull out inactive ghost nodes based on a FieldList of flags.
//------------------------------------------------------------------------------
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::cullGhostNodes(const FieldList<Dimension, size_t>& flagSet,
                                                 FieldList<Dimension, size_t>& old2newIndexMap,
                                                 vector<size_t>& numNodesRemoved) {
  // Base class does some stuff
  Boundary<Dimension>::cullGhostNodes(flagSet, old2newIndexMap, numNodesRemoved);

  // Patch our own node indexing.
  if (mUseGhosts) {
    const auto numNodeLists = flagSet.numFields();
    const auto numFacets = mPoly.facets().size();
    for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
      const auto& nodeList = flagSet[nodeListi]->nodeList();
      const auto  name = nodeList.name();
      auto&       facetControls = mFacetControlNodes[name];
      auto&       facetGhostRanges = mFacetGhostNodes[name];
      CHECK(facetControls.size() == numFacets);
      CHECK(facetGhostRanges.size() == numFacets);
      for (auto ifacet = 0u; ifacet < numFacets; ++ifacet) {
        CHECK((int)facetControls[ifacet].size() == facetGhostRanges[ifacet].second - facetGhostRanges[ifacet].first);

        // Update facet controls
        vector<int> newControls;
        const auto oldFirstGhost = facetGhostRanges[ifacet].first;
        for (auto k = 0u; k < facetControls[ifacet].size(); ++k) {
          const auto i = facetControls[ifacet][k];
          const auto j = oldFirstGhost + k;
          if (flagSet(nodeListi, i) == 1 and flagSet(nodeListi, j) == 1) newControls.push_back(old2newIndexMap(nodeListi, i));
        }
        facetControls[ifacet] = newControls;

        // Update facet ghosts
        int newBegin=-1, numNewGhosts=0;
        for (auto i = facetGhostRanges[ifacet].first; i < facetGhostRanges[ifacet].second; ++i) {
          if (flagSet(nodeListi, i) == 1) {
            if (newBegin == -1) newBegin = old2newIndexMap(nodeListi, i);
            ++numNewGhosts;
          }
        }
        if (newBegin == -1) {
          facetGhostRanges[ifacet] = std::make_pair(nodeList.numNodes(), nodeList.numNodes());
        } else {
          facetGhostRanges[ifacet] = std::make_pair(newBegin, newBegin + numNewGhosts);
        }

        CHECK2((int)facetControls[ifacet].size() == facetGhostRanges[ifacet].second - facetGhostRanges[ifacet].first,
               facetControls[ifacet].size() << " != " << (facetGhostRanges[ifacet].second - facetGhostRanges[ifacet].first)
               << "   ghostRange=(" << facetGhostRanges[ifacet].first << " " << facetGhostRanges[ifacet].second << ")");
      }
    }
  }
}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#ifdef SPHERAL2D
template class FacetedVolumeBoundary<Dim<2>>;
#endif

#ifdef SPHERAL3D
template class FacetedVolumeBoundary<Dim<3>>;
#endif

}
