//---------------------------------Spheral++----------------------------------//
// IncrementCullenMultipliers
//
// Created by JMO, Mon Dec 28 22:18:44 PST 2015
//----------------------------------------------------------------------------//
#include "IncrementCullenMultipliers.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/FieldList.hh"
#include "Utilities/DBC.hh"
#include "Hydro/HydroFieldNames.hh"

using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {

//------------------------------------------------------------------------------
// Constructors.
//------------------------------------------------------------------------------
template<typename Dimension>
IncrementCullenMultipliers<Dimension>::
IncrementCullenMultipliers(const typename Dimension::Scalar minValue,
                           const typename Dimension::Scalar maxValue,
                           const bool hopkinsForm):
  UpdatePolicyBase<Dimension>(),
  mHopkinsForm(hopkinsForm),
  mAlphaMin(minValue),
  mAlphaMax(maxValue) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
IncrementCullenMultipliers<Dimension>::
~IncrementCullenMultipliers() {
}

//------------------------------------------------------------------------------
// Update the Cullen Q multiplier factors.
//------------------------------------------------------------------------------
template<typename Dimension>
void
IncrementCullenMultipliers<Dimension>::
update(const KeyType&,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double /*t*/,
       const double /*dt*/) {

  // Get the state we're advancing and the needed derivatives.
  auto rvQ = state.fields(HydroFieldNames::ArtificialViscousCqMultiplier, 0.0);
  auto rvL = state.fields(HydroFieldNames::ArtificialViscousClMultiplier, 0.0);
  auto alpha0 = state.fields("mCullAlpha", 0.0);
  const auto alpha_local = derivs.fields("Cullen alpha local", 0.0);
  const auto DalphaDt = derivs.fields("Cullen alpha delta", 0.0);
  const auto alpha_tmp = derivs.fields("mCullAlpha2", 0.0);

  const auto numNodeLists = rvQ.size();
  CHECK(rvL.size() == numNodeLists);
  CHECK(alpha0.size() == numNodeLists);
  CHECK(alpha_local.size() == numNodeLists);
  CHECK(DalphaDt.size() == numNodeLists);

  // Set the new values.
  for (auto k = 0u; k < numNodeLists; ++k) {
    const auto n = rvQ[k]->numInternalElements();
#pragma omp parallel for
    for (auto i = 0u; i < n; ++i) {
      if (mHopkinsForm) {
        // Hopkins 2014
        const auto alphai = std::max(mAlphaMin, std::min(mAlphaMax, alpha_local(k, i)));
        rvQ(k, i) = alphai;
        rvL(k, i) = alphai;
        alpha0(k, i) = DalphaDt(k, i);
      } else {
        // Cullen & Dehnen 2010
        const auto alphai = std::max(mAlphaMin, std::min(mAlphaMax, max(alpha_local(k, i), rvQ(k, i)) + multiplier*DalphaDt(k, i)));
        rvQ(k, i) = alphai;
        rvL(k, i) = alphai;
      }
    }
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
IncrementCullenMultipliers<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const auto* rhsPtr = dynamic_cast<const IncrementCullenMultipliers<Dimension>*>(&rhs);
  if (rhsPtr == nullptr) return false;

  // Ok, now do we agree on min & max?
  return (mAlphaMin == rhsPtr->mAlphaMin) && (mAlphaMax == rhsPtr->mAlphaMax);
}

}

