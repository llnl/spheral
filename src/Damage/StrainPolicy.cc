//---------------------------------Spheral++----------------------------------//
// StrainPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the strain.
//
// Created by JMO, Tue Oct 12 16:43:23 PDT 2004
//----------------------------------------------------------------------------//
#include "StrainPolicy.hh"
#include "NodeList/NodeList.hh"
#include "Strength/SolidFieldNames.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/UpdatePolicyBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/Field.hh"
#include "Utilities/DBC.hh"
#include "Kernel/HatKernel.hh"

#include <vector>
using std::vector;
using std::string;
using std::pair;
using std::make_pair;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
StrainPolicy<Dimension>::
StrainPolicy():
  FieldUpdatePolicy<Dimension, SymTensor>({HydroFieldNames::position,
                                           HydroFieldNames::H,
                                           SolidFieldNames::YoungsModulus,
                                           HydroFieldNames::pressure,
                                           SolidFieldNames::deviatoricStress}) {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
StrainPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& /*derivs*/,
       const double /*multiplier*/,
       const double /*t*/,
       const double /*dt*/) {
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == SolidFieldNames::strain);
  auto& stateField = state.field(key, 0.0);

  const auto tiny = 1.0e-30;

  // Get the state fields.
  auto buildKey = [&](const std::string& fkey) -> std::string { return StateBase<Dimension>::buildFieldKey(fkey, nodeListKey); };
  const auto& E = state.field(buildKey(SolidFieldNames::YoungsModulus), 0.0);
  const auto& P = state.field(buildKey(HydroFieldNames::pressure), 0.0);
  const auto& S = state.field(buildKey(SolidFieldNames::deviatoricStress), SymTensor::zero);

  // Iterate over the internal nodes.
  const auto n = stateField.numInternalElements();
#pragma omp parallel for
  for (auto i = 0u; i < n ; ++i) {

    // Compute the maximum tensile stress.
//     Scalar Pi = P(i);
//     if (Pi < 0.0) Pi *= 1.0 - D(i);
    const auto sigmai = S(i) - P(i)*SymTensor::one;
    const auto sigmati = max(0.0, sigmai.eigenValues().maxElement());
    CHECK(sigmati >= 0.0);

    // Compute the initial strain for this node.
    CHECK(E(i) >= 0.0);
    stateField(i) = sigmati/(E(i) + tiny*max(1.0, sigmati));
    CHECK(stateField(i) >= 0.0);
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
StrainPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {
  return dynamic_cast<const StrainPolicy<Dimension>*>(&rhs) != nullptr;
}

}

