//---------------------------------Spheral++----------------------------------//
// PairMaxDamageNodeCoupling
//
// A functor class encapsulating how we couple solid nodes in the presence of
// multiple materials and damage.
//
// Created by JMO, Fri Jul 31 14:46:25 PDT 2015
//----------------------------------------------------------------------------//
#include "Damage/PairMaxDamageNodeCoupling.hh"
#include "Utilities/DBC.hh"
#include "Field/FieldList.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Strength/SolidFieldNames.hh"

namespace Spheral {

namespace {

//------------------------------------------------------------------------------
// Extract the effective scalar damage on a node.
//------------------------------------------------------------------------------
template<typename FieldListType>
inline
double scalarDamage(const FieldListType& damage,
                    const unsigned nodeListi, const unsigned i) {
  auto sDi = std::max(0.0, std::min(1.0, damage(nodeListi, i).Trace()/FieldListType::FieldDimension::nDim));
  if (sDi > 1.0 - 1.0e-3) sDi = 1.0;
  return sDi;
}

// // Construct a unit vector of the argument, going to zero as the magnitude falls below a given "fuzz".
// Vector unitVectorWithZero(const Vector& x,
//                           const double fuzz = 0.01) const {
//   if (x.magnitude2() < fuzz) {
//     return Vector::zero();
//   } else {
//     return x.unitVector();
//   }
// }

// // Damage direction based on the gradient.
// Vector damageDirection(const unsigned nodeListi, const unsigned i) const {
//   return unitVectorWithZero(mDamageGradient(nodeListi, i)*Dimension::nDim/mH(nodeListi, i).Trace());
// }

}

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename Dimension>
PairMaxDamageNodeCoupling<Dimension>::
PairMaxDamageNodeCoupling(const State<Dimension>& state,
                          NodePairList& pairs):
  NodeCoupling() {
  const auto n = pairs.size();
  const auto pos = state.fields(HydroFieldNames::position, Vector::zero());
  const auto D = state.fields(SolidFieldNames::tensorDamage, SymTensor::zero());
  const auto fragIDs = state.fields(SolidFieldNames::fragmentIDs, int(0));
#pragma omp for
  for (auto k = 0u; k < n; ++k) {
    auto& pair = pairs[k];
    if ((pair.i_list != pair.j_list) or
        (fragIDs(pair.i_list, pair.i_node) != fragIDs(pair.j_list, pair.j_node))) {
      pair.f_couple = 0.0;
    } else {
      const auto xijhat = (pos(pair.i_list, pair.i_node) - pos(pair.j_list, pair.j_node)).unitVector();
      const auto sDi = (D(pair.i_list, pair.i_node)*xijhat).magnitude();
      const auto sDj = (D(pair.j_list, pair.j_node)*xijhat).magnitude();
      const auto sDmin = std::min(sDi, sDj);
      const auto sDmax = std::max(sDi, sDj);
      pair.f_couple = std::max(0.0, (1.0 - sDmax)*(1.0 - sDmin) + sDmin);  // sDmin weighting to return to fully coupled for fully damaged nodes
    }
  }
}

}
