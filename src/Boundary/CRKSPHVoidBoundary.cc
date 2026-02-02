//---------------------------------Spheral++----------------------------------//
// CRKSPHVoidBoundary -- Special boundary class to create void points off of
// the surface.
//----------------------------------------------------------------------------//
#include "CRKSPHVoidBoundary.hh"
#include "Field/FieldList.hh"
#include "Field/Field.hh"
#include "Field/FieldBase.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Strength/SolidFieldNames.hh"
#include "Utilities/DBC.hh"

#include <limits>

namespace Spheral {

using std::vector;

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
CRKSPHVoidBoundary<Dimension>::
CRKSPHVoidBoundary(const FieldList<Dimension, int>& surfacePoint,
                   const FieldList<Dimension, vector<Vector>>& etaVoidPoints):
  Boundary<Dimension>(),
  mSurfacePoint(surfacePoint),
  mEtaVoidPoints(etaVoidPoints) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
CRKSPHVoidBoundary<Dimension>::~CRKSPHVoidBoundary() {
}

//------------------------------------------------------------------------------
// setGhostNodes
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHVoidBoundary<Dimension>::
setGhostNodes(NodeList<Dimension>& nodeList) {
  this->addNodeList(nodeList);

  typename Boundary<Dimension>::BoundaryNodes& boundaryNodes = this->accessBoundaryNodes(nodeList);
  auto& cNodes = boundaryNodes.controlNodes;
  auto& gNodes = boundaryNodes.ghostNodes;
  cNodes.clear();
  gNodes.clear();
  const unsigned firstNewGhostNode = nodeList.numNodes();
    
  // Line up the surface points as controls for the void points.
  unsigned numVoidPoints = 0;
  unsigned j = firstNewGhostNode;
  if (mSurfacePoint.haveNodeList(nodeList)) {
    const Field<Dimension, int>& surfacePoint = **mSurfacePoint.fieldForNodeList(nodeList);
    const Field<Dimension, vector<Vector>>& etaVoidPoints = **mEtaVoidPoints.fieldForNodeList(nodeList);
    for (unsigned i = 0; i < firstNewGhostNode; ++i) {
      if ((surfacePoint(i) & 1) == 1) {
        const unsigned nv = etaVoidPoints(i).size();
        CHECK(nv > 0);
        numVoidPoints += nv;
        for (unsigned k = 0; k != nv; ++k) {
          cNodes.push_back(i);
          gNodes.push_back(j++);
        }
      }
    }
  }
  CHECK(cNodes.size() == numVoidPoints);
  CHECK(gNodes.size() == numVoidPoints);

  // Create the correct number of new void ghost points.
  if (numVoidPoints > 0) {
    const unsigned numGhostNodes0 = nodeList.numGhostNodes();
    nodeList.numGhostNodes(numGhostNodes0 + numVoidPoints);
  }

  // Update the void point positions and H's.
  this->updateGhostNodes(nodeList);
}

//------------------------------------------------------------------------------
// updateGhostNodes
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHVoidBoundary<Dimension>::
updateGhostNodes(NodeList<Dimension>& nodeList) {
  const auto& cNodes = this->controlNodes(nodeList);
  const auto& gNodes = this->ghostNodes(nodeList);
  const unsigned nvoid = gNodes.size();
  CONTRACT_VAR(cNodes);
  CHECK(cNodes.size() == nvoid);
  
  if (nvoid > 0) {
    Field<Dimension, Vector>& pos = nodeList.positions();
    Field<Dimension, SymTensor>& H = nodeList.Hfield();
    const Field<Dimension, int>& surfacePoint = **mSurfacePoint.fieldForNodeList(nodeList);
    const Field<Dimension, vector<Vector>>& etaVoidPoints = **mEtaVoidPoints.fieldForNodeList(nodeList);

    const unsigned n = pos.numInternalElements();  // Note we assume this is the first BC, and there are no ghost masters.
    unsigned j = gNodes[0];                        // Assuming our ghost nodes are sequential.

    // Update the void points positions based on projecting from the surface point in the normalized eta (H) frame.
    for (unsigned i = 0; i < n; ++i) {
      if ((surfacePoint(i) & 1) == 1) {
        const Vector& xi = pos(i);
        const SymTensor& Hi = H(i);
        const SymTensor Hinv = Hi.Inverse();

        const unsigned nv = etaVoidPoints(i).size();
        CHECK(nv > 0);
        for (unsigned k = 0; k < nv; ++k) {
          Vector& xj = pos(j);
          SymTensor& Hj = H(j);

          xj = xi + Hinv*etaVoidPoints(i)[k];
          Hj = Hi;
          ++j;
        }
      }
    }
    CHECK(j == gNodes[0] + nvoid);
  }
}

//------------------------------------------------------------------------------
// Apply the boundary condition to fields of different DataTypes.
//------------------------------------------------------------------------------
// Specialization for int fields.
template<typename Dimension>
void
CRKSPHVoidBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, int>& field) const {
  const auto& gNodes = this->ghostNodes(field.nodeList());
  const unsigned nvoid = gNodes.size();
  if (field.name() == HydroFieldNames::voidPoint) {
    // voidPoint: flag only ghost void points
    for (unsigned k = 0; k < nvoid; ++k) field(gNodes[k]) = 1;
  } else {
    for (unsigned k = 0; k < nvoid; ++k) field(gNodes[k]) = 0;
  }
}

// Specialization for scalar fields.
template<typename Dimension>
void
CRKSPHVoidBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::Scalar>& field) const {
  const auto& cNodes = this->controlNodes(field.nodeList());
  const auto& gNodes = this->ghostNodes(field.nodeList());
  const unsigned nsurf = cNodes.size();
  CHECK(gNodes.size() == nsurf);
  if (field.name() == HydroFieldNames::volume) {
    // Volume: copy surface->void
    for (unsigned k = 0; k < nsurf; ++k) field(gNodes[k]) = field(cNodes[k]);
  } else if (field.name() == HydroFieldNames::mass or
             field.name() == HydroFieldNames::massDensity) {
    // mass, mass density: negligible but non-zero
    for (unsigned k = 0; k < nsurf; ++k) field(gNodes[k]) = std::numeric_limits<Scalar>::epsilon();
  } else {
    // Default zero.
    for (unsigned k = 0; k < nsurf; ++k) field(gNodes[k]) = 0.0;
  }
}

// Specialization for Vector fields.
template<typename Dimension>
void
CRKSPHVoidBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::Vector>& field) const {
  const auto& cNodes = this->controlNodes(field.nodeList());
  const auto& gNodes = this->ghostNodes(field.nodeList());
  const unsigned nsurf = cNodes.size();
  CHECK(gNodes.size() == nsurf);
  if (field.name() == HydroFieldNames::velocity) {
    // velocity: copy surface->void
    for (unsigned k = 0; k < nsurf; ++k) field(gNodes[k]) = field(cNodes[k]);
  } else {
    // Default zero.
    for (unsigned k = 0; k < nsurf; ++k) field(gNodes[k]) = Vector::zero();
  }
}

// Specialization for Tensor fields.
template<typename Dimension>
void
CRKSPHVoidBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::Tensor>& field) const {
  const auto& gNodes = this->ghostNodes(field.nodeList());
  const unsigned nvoid = gNodes.size();
  for (unsigned k = 0; k < nvoid; ++k) field(gNodes[k]) = Tensor::zero();
}

// Specialization for symmetric tensors.
template<typename Dimension>
void
CRKSPHVoidBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::SymTensor>& field) const {
  const auto& cNodes = this->controlNodes(field.nodeList());
  const auto& gNodes = this->ghostNodes(field.nodeList());
  const unsigned nsurf = cNodes.size();
  CHECK(gNodes.size() == nsurf);
  if (field.name() == SolidFieldNames::deviatoricStress) {
    // deviatoric stress: copy surface->void
    for (unsigned k = 0; k < nsurf; ++k) field(gNodes[k]) = field(cNodes[k]);
  } else {
    // Default zero.
    for (unsigned k = 0; k < nsurf; ++k) field(gNodes[k]) = SymTensor::zero();
  }
}

// Specialization for third rank tensors.
template<typename Dimension>
void
CRKSPHVoidBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::ThirdRankTensor>& field) const {
  const auto& gNodes = this->ghostNodes(field.nodeList());
  const unsigned nvoid = gNodes.size();
  for (unsigned k = 0; k < nvoid; ++k) field(gNodes[k]) = ThirdRankTensor::zero();
}

// Specialization for fourth rank tensors.
template<typename Dimension>
void
CRKSPHVoidBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::FourthRankTensor>& field) const {
  const auto& gNodes = this->ghostNodes(field.nodeList());
  const unsigned nvoid = gNodes.size();
  for (unsigned k = 0; k < nvoid; ++k) field(gNodes[k]) = FourthRankTensor::zero();
}

// Specialization for fifth rank tensors.
template<typename Dimension>
void
CRKSPHVoidBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::FifthRankTensor>& field) const {
  const auto& gNodes = this->ghostNodes(field.nodeList());
  const unsigned nvoid = gNodes.size();
  for (unsigned k = 0; k < nvoid; ++k) field(gNodes[k]) = FifthRankTensor::zero();
}

// Specialization for FacetedVolume fields.
template<typename Dimension>
void
CRKSPHVoidBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::FacetedVolume>& field) const {
  const auto& gNodes = this->ghostNodes(field.nodeList());
  const unsigned nvoid = gNodes.size();
  for (unsigned k = 0; k < nvoid; ++k) field(gNodes[k]) = FacetedVolume();
}

// Specialization for vector<scalar> fields.
template<typename Dimension>
void
CRKSPHVoidBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, std::vector<typename Dimension::Scalar> >& field) const {
  const auto& gNodes = this->ghostNodes(field.nodeList());
  const unsigned nvoid = gNodes.size();
  for (unsigned k = 0; k < nvoid; ++k) field(gNodes[k]) = std::vector<Scalar>();
}

//------------------------------------------------------------------------------
// Provide the setViolationNodes for a NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHVoidBoundary<Dimension>::
setViolationNodes(NodeList<Dimension>& nodeList) {
  this->addNodeList(nodeList);
}

//------------------------------------------------------------------------------
// Provide the updateViolationNodes for a NodeList, correcting positions and 
// H's.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHVoidBoundary<Dimension>::
updateViolationNodes(NodeList<Dimension>&) {
}

}
