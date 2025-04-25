//---------------------------------Spheral++----------------------------------//
// DeviatoricStressPolity -- An implementation of UpdatePolicyBase appropriate for
// when 'ya just want to increment by derivatives:  x1 = x0 + A*dx/dt
//
// Created by JMO, Thu Aug 26 14:28:07 2004
//----------------------------------------------------------------------------//
#include "DeviatoricStressPolicy.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/IncrementBoundedState.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Strength/SolidFieldNames.hh"
#include "Kernel/TableKernel.hh"
#include "Material/EquationOfState.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Geometry/GeometryRegistrar.hh"
#include "Utilities/DBC.hh"

#include <algorithm>

namespace Spheral {

//------------------------------------------------------------------------------
// Function to enforce zeroing the trace of the deviatoric stress in 3D only
//------------------------------------------------------------------------------
namespace {  // anonymous

template<typename Dimension>
inline
void
zeroTrace(typename Dimension::SymTensor& Si) {
}

template<>
inline
void
zeroTrace<Dim<3>>(Dim<3>::SymTensor& Si) {
  const auto dS = Si.Trace()/3.0;
  Si[0] -= dS;
  Si[3] -= dS;
  Si[5] -= dS;
}

}            // anonymous

//------------------------------------------------------------------------------
// Constructors.
//------------------------------------------------------------------------------
template<typename Dimension>
DeviatoricStressPolicy<Dimension>::
DeviatoricStressPolicy():
  FieldUpdatePolicy<Dimension, SymTensor>({}) {
}

//------------------------------------------------------------------------------
// Update the FieldList.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DeviatoricStressPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {

  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == SolidFieldNames::deviatoricStress);

  // Alias for shorter call building State Field keys
  auto buildKey = [&](const std::string& fkey) -> std::string { return StateBase<Dimension>::buildFieldKey(fkey, nodeListKey); };

  // Get the state we're advancing.
  auto&       S = state.field(key, SymTensor::zero);
  const auto& DSDt = derivs.field(IncrementState<Dimension, SymTensor>::prefix() + key, SymTensor::zero);

  // Check if a porosity model has registered a modifier for the deviatoric stress.
  // They should have added it as a dependency of this policy if so.
  const auto porosityScaling = state.registered(buildKey(SolidFieldNames::fDSjutzi));
  const Field<Dimension, Scalar>* fDSptr = nullptr;
  const Field<Dimension, Scalar>* alphaPtr = nullptr;
  const Field<Dimension, Scalar>* DalphaDtPtr = nullptr;
  if (porosityScaling) {
    fDSptr = &state.field(buildKey(SolidFieldNames::fDSjutzi), 0.0);
    alphaPtr = &state.field(buildKey(SolidFieldNames::porosityAlpha), 0.0);
    DalphaDtPtr = &derivs.field(buildKey(IncrementBoundedState<Dimension, Scalar>::prefix() + SolidFieldNames::porosityAlpha), 0.0);
  }

  // Iterate over the internal nodes.
  const auto n = S.numInternalElements();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {

    // Get the deviatoric time derivative, and if necessary scale by the porosity factor
    auto DSDti = DSDt(i);
    if (porosityScaling) {
      const auto fDSi = (*fDSptr)(i);
      const auto alphai = (*alphaPtr)(i);
      const auto DalphaDti = (*DalphaDtPtr)(i);
      CHECK(alphai >= 1.0);
      DSDti = (fDSi*DSDti - S(i)*DalphaDti/alphai)/alphai;
    }

    // Update S
    // Note -- purely elastic flow.  The plastic yielding is accounted for when we update the plastic strain.
    S(i) += multiplier*DSDti;                                          // Elastic prediction for the new deviatoric stress

    // We only want to enforce zeroing the trace in 3D Cartesian coordinates.  In lower
    // dimensions we assume the missing components on the diagonal sum to -Trace(S).
    zeroTrace<Dimension>(S(i));
  }

//     // Finally apply the pressure limits to the allowed deviatoric stress.
//     S(i) = max(Pmin, min(Pmax, S(i)));
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
DeviatoricStressPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is a DeviatoricStress operator, and has
  // the same internal parameters
  const auto rhsPtr = dynamic_cast<const DeviatoricStressPolicy<Dimension>*>(&rhs);
  return (rhsPtr != nullptr);
}

}

