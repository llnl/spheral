//---------------------------------Spheral++----------------------------------//
// IncrementVoronoiCells
//
// Specialization of UpdatePolicyBase to advance the Vornoi cell geometry
// during a step without actually recomputing the geometry.  Instead we distort
// the cells by the local velocity gradient.
//
// Created by JMO, Mon May 20 16:04:51 PDT 2024
//----------------------------------------------------------------------------//
#include "VoronoiCells/IncrementVoronoiCells.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Utilities/DBC.hh"

#include <vector>

using std::vector;

namespace Spheral {

namespace {   // anonymous

//------------------------------------------------------------------------------
// Function to apply a velocity distortion to a polytope.  Need to specialize
// for 1D...
//------------------------------------------------------------------------------
template<typename Poly>
inline
Poly
distortPolytope(const Poly& poly,
                const typename Poly::Vector& ri,
                const typename Poly::Vector& vi,
                const typename Poly::Tensor& DvDxi,
                const double dt) {
  using Vector = typename Poly::Vector;
  const auto& facetIndices = poly.facetVertices();
  vector<Vector> verts(poly.vertices());
  for (auto& v: verts) {
    const auto dr = v - ri;
    v += (vi + DvDxi.dot(dr))*dt;
  }
  return Poly(verts, facetIndices);
}

//..............................................................................
// 1D
template<>
inline
Dim<1>::FacetedVolume
distortPolytope<Dim<1>::FacetedVolume>(const Dim<1>::FacetedVolume& poly,
                                       const Dim<1>::FacetedVolume::Vector& ri,
                                       const Dim<1>::FacetedVolume::Vector& vi,
                                       const Dim<1>::FacetedVolume::Tensor& DvDxi,
                                       const double dt) {
  using Vector = Dim<1>::Vector;
  auto updateCoords = [&](const Vector& pos) { return pos + (vi + DvDxi.dot(pos - ri))*dt; };
  const auto v1 = updateCoords(poly.xmin());
  const auto v2 = updateCoords(poly.xmax());
  return Dim<1>::FacetedVolume(0.5*(v1 + v2), 0.5*abs(v2.x() - v1.x()));
}

}             // anonymous

//------------------------------------------------------------------------------
// Constructors.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
IncrementVoronoiCells<Dimension>::
IncrementVoronoiCells():
  UpdatePolicyBase<Dimension>({HydroFieldNames::velocity}) {
}

//------------------------------------------------------------------------------
// Update the Voronoi cells
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
IncrementVoronoiCells<Dimension>::
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

  // We depend on velocity information.  If there is no velocity there's nothing to do.
  if (state.fieldNameRegistered(HydroFieldNames::velocity)) {
    const auto pos = state.fields(HydroFieldNames::position, Vector::zero());
    const auto vel = state.fields(HydroFieldNames::velocity, Vector::zero());
    // const auto DvDt = derivs.fields(HydroFieldNames::hydroAcceleration, Vector::zero());
    const auto DvDx = derivs.fields(HydroFieldNames::internalVelocityGradient, Tensor::zero());

    const auto numNodeLists = cells.numFields();
    for (auto k = 0u; k < numNodeLists; ++k) {
      const auto n = cells[k]->numInternalElements();
#pragma omp parallel for
      for (auto i = 0u; i < n; ++i) {
        const auto& ri = pos(k,i);
        const auto& vi = vel(k,i);
        // const auto  vi = vel(k,i) - multiplier*DvDt(k,i);   // velocity was already advanced, so back it up
        const auto& DvDxi = DvDx(k,i);
        auto& celli = cells(k,i);
        celli = distortPolytope(celli, ri, vi, DvDxi, multiplier);
      }
    }
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
IncrementVoronoiCells<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {
  return dynamic_cast<const IncrementVoronoiCells<Dimension>*>(&rhs) != nullptr;
}

}

