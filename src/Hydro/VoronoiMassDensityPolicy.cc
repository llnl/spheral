//---------------------------------Spheral++----------------------------------//
// VoronoiMassDensityPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the mass density according to the specific volume from the
// Voronoi tesselation.
//
// Created by JMO, Tue Sep 14 22:27:08 2004
//----------------------------------------------------------------------------//
#include "VoronoiMassDensityPolicy.hh"
#include "HydroFieldNames.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/FieldList.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/DBC.hh"

using std::vector;
using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
VoronoiMassDensityPolicy<Dimension>::
VoronoiMassDensityPolicy(const double rhoMin, const double rhoMax):
  UpdatePolicyBase<Dimension>({HydroFieldNames::mass,
                               HydroFieldNames::volume}),
  mRhoMin(rhoMin),
  mRhoMax(rhoMax) {
  REQUIRE(rhoMin <= rhoMax);
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
VoronoiMassDensityPolicy<Dimension>::
~VoronoiMassDensityPolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
VoronoiMassDensityPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& /*derivs*/,
       const double /*multiplier*/,
       const double /*t*/,
       const double /*dt*/) {

  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == HydroFieldNames::massDensity and 
          nodeListKey == UpdatePolicyBase<Dimension>::wildcard());
  auto       massDensity = state.fields(fieldKey, Scalar());
  const auto numFields = massDensity.numFields();

  // Get the mass and volume from the state.
  const auto mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto volume = state.fields(HydroFieldNames::volume, 0.0);

  // Set the mass density.
  for (auto i = 0u; i < numFields; ++i) {
    const auto n = massDensity[i]->numInternalElements();
#pragma omp parallel for
    for (auto j = 0u; j < n; ++j) {
      massDensity(i,j) = max(mRhoMin, min(mRhoMax, mass(i,j) * safeInv(volume(i,j))));
    }
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
VoronoiMassDensityPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const auto rhsPtr = dynamic_cast<const VoronoiMassDensityPolicy<Dimension>*>(&rhs);
  return (rhsPtr != nullptr);
}

}

