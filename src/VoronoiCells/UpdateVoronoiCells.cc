//---------------------------------Spheral++----------------------------------//
// UpdateVoronoiCells
//
// Specialization of UpdatePolicyBase to advance the Vornoi cell geometry
// during a step without actually recomputing the geometry.  Instead we distort
// the cells by the local velocity gradient.
//
// Created by JMO, Mon May 20 16:04:51 PDT 2024
//----------------------------------------------------------------------------//
#include "VoronoiCells/UpdateVoronoiCells.hh"
#include "VoronoiCells/computeVoronoiVolume.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Strength/SolidFieldNames.hh"
#include "Geometry/CellFaceFlag.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Utilities/DBC.hh"

#include <vector>

using std::vector;

namespace Spheral {

//------------------------------------------------------------------------------
// Constructors.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
UpdateVoronoiCells<Dimension>::
UpdateVoronoiCells(FieldList<Dimension, typename Dimension::Scalar>& volume,
                   FieldList<Dimension, typename Dimension::Scalar>& weight,
                   FieldList<Dimension, typename Dimension::Vector>& deltaCentroid,
                   FieldList<Dimension, vector<typename Dimension::Vector>>& etaVoidPoints,
                   const std::vector<Boundary<Dimension>*>& boundaries,
                   const std::vector<FacetedVolume>& facetedBoundaries,
                   const std::vector<std::vector<FacetedVolume>>& facetedHoles):
  UpdatePolicyBase<Dimension>({HydroFieldNames::position,
                               HydroFieldNames::H,
                               HydroFieldNames::mass,
                               HydroFieldNames::massDensity,
                               SolidFieldNames::tensorDamage}),
  mVolume(volume),
  mWeight(weight),
  mDeltaCentroid(deltaCentroid),
  mEtaVoidPoints(etaVoidPoints),
  mBoundaries(boundaries),
  mFacetedBoundaries(facetedBoundaries),
  mFacetedHoles(facetedHoles) {
}

//------------------------------------------------------------------------------
// Update the Voronoi cells
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
UpdateVoronoiCells<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {

  // Check the key 
  BEGIN_CONTRACT_SCOPE;
  {
    KeyType fieldKey, nodeListKey;
    StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
    REQUIRE(fieldKey == HydroFieldNames::cells);
    REQUIRE(nodeListKey == UpdatePolicyBase<Dimension>::wildcard());
  }
  END_CONTRACT_SCOPE;

  // Get the state we're updating.
  auto cells = state.template fields<FacetedVolume>(HydroFieldNames::cells);
  auto surfacePoint = state.fields(HydroFieldNames::surfacePoint, 0);
  auto cellFaceFlags = state.template fields<vector<CellFaceFlag>>(HydroFieldNames::cellFaceFlags);

  // State we need to compute the Voronoi cells
  const auto& cm = state.connectivityMap();
  const auto  pos = state.fields(HydroFieldNames::position, Vector::zero());
  const auto  H = state.fields(HydroFieldNames::H, SymTensor::zero());
  const auto  mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto  rho = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto  D = state.fields(SolidFieldNames::tensorDamage, SymTensor::zero());

  // Use m/rho to estimate our weighting to roughly match cell volumes
  const auto numNodeLists = cells.numFields();
  for (auto k = 0u; k < numNodeLists; ++k) {
    const auto n = mass[k]->numInternalElements();
#pragma omp parallel for
    for (auto i = 0u; i < n; ++i) {
      CHECK(rho(k,i) > 0.0);
      mVolume(k,i) = mass(k,i)/rho(k,i);
      // mWeight(k,i) = 1.0/Dimension::rootnu(mVolume(k,i));
    }
  }
  
  for (auto* bcPtr: mBoundaries) {
    bcPtr->applyFieldListGhostBoundary(mVolume);
    // bcPtr->applyFieldListGhostBoundary(mWeight);
  }
  for (auto* bcPtr: mBoundaries) bcPtr->finalizeGhostBoundary();

  // Compute the cell data.  Note we are using the fact the state versions of the things
  // we're updating (mSurfacePoint, mCells, etc.) are just pointing at our internal fields.
  computeVoronoiVolume(pos, H, cm, D, mFacetedBoundaries, mFacetedHoles, mBoundaries, mWeight,
                       surfacePoint, mVolume, mDeltaCentroid, mEtaVoidPoints, cells, cellFaceFlags);
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
UpdateVoronoiCells<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {
  return dynamic_cast<const UpdateVoronoiCells<Dimension>*>(&rhs) != nullptr;
}

}

