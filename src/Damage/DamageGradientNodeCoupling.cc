//---------------------------------Spheral++----------------------------------//
// DamageGradientNodeCoupling
//
// A functor class encapsulating how we couple solid nodes in the presence of
// multiple materials and damage.
//
// This one attempts to mock up the shielding effect of ThreePointDamagedNodeCoupling
// by using local damage gradient to estimate when nodes are separated by
// regions of greater damage (or fractures).
//
// Created by JMO, Fri Jul 31 14:46:25 PDT 2015
//----------------------------------------------------------------------------//
#include "Damage/DamageGradientNodeCoupling.hh"
#include "Utilities/pointDistances.hh"
#include "Utilities/DBC.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Strength/SolidFieldNames.hh"
#include "Field/FieldList.hh"
#include "FieldOperations/FieldListFunctions.hh"
#include "Boundary/Boundary.hh"
#include "Utilities/Timer.hh"

#include <vector>

using std::vector;

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename Dimension>
DamageGradientNodeCoupling<Dimension>::
DamageGradientNodeCoupling(const State<Dimension>& state,
                           const TableKernel<Dimension>& W,
                           ConstBoundaryIterator boundaryBegin,
                           ConstBoundaryIterator boundaryEnd,
                           NodePairList& pairs):
  NodeCoupling() {

  TIME_BEGIN("DamageGradientCoupling");

  const auto npairs = pairs.size();
  const auto position = state.fields(HydroFieldNames::position, Vector::zero());
  const auto H = state.fields(HydroFieldNames::H, SymTensor::zero());
  const auto D = state.fields(SolidFieldNames::tensorDamage, SymTensor::zero());
  const auto mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto rho = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto numNodeLists = position.numFields();

  // Compute a scalar damage estimate to take the gradient of.
  TIME_BEGIN("DamageGradientCoupling_grad");
  FieldList<Dimension, Scalar> Dtrace(FieldStorageType::CopyFields);
  for (auto k = 0u; k < numNodeLists; ++k) {
    Dtrace.appendNewField("grad D", D[k]->nodeList(), 0.0);
    const auto n = D[k]->numElements();
#pragma omp parallel for
    for (auto i = 0u; i < n; ++i) {
      Dtrace(k,i) = D(k,i).Trace();
    }
  }

  // Compute the damage gradient.
  const auto weight = mass/rho;
  auto gradD = gradient(Dtrace, position, weight, mass, rho, H, W);  // third rank tensor
  for (auto boundaryItr = boundaryBegin; boundaryItr < boundaryEnd; ++boundaryItr) {
    (*boundaryItr)->applyFieldListGhostBoundary(gradD);
    (*boundaryItr)->finalizeGhostBoundary();
  }
  TIME_END("DamageGradientCoupling_grad");

  // For each interacting pair we need to compute the effective damage shielding, expressed
  // as the f_couple parameter in the NodePairIdxType.
  TIME_BEGIN("DamageGradientCoupling_pairs");
  size_t i, il, j, jl;
  Scalar hij;
  Vector xjihat;
#pragma omp parallel for private(i, il, j, jl, hij, xjihat)
  for (auto k = 0u; k < npairs; ++k) {
    i = pairs[k].i_node;
    il = pairs[k].i_list;
    j = pairs[k].j_node;
    jl = pairs[k].j_list;
    xjihat = (position(jl, j) - position(il, i)).unitVector();
    // di = (D(il, i)*xijhat).magnitude();
    // dj = (D(jl, j)*xijhat).magnitude();
    // CHECK(di >= 0.0 and di <= 1.0);
    // CHECK(dj >= 0.0 and dj <= 1.0);
    hij = 2.0/((H(il,i)*xjihat).magnitude() + (H(jl,j)*xjihat).magnitude());
    pairs[k].f_couple = std::max(0.0, 1.0 + std::min(0.0, gradD(il, i).dot(gradD(jl, j))*4.0*hij));
    // if (i == 47 and j == 49) {
    //   std::cerr << " --> " << pairs[k] << " " << D(il,i) << " " << D(jl,j) << " " << gradD(il,i) << " " << gradD(jl,j) << " " << gradD(il,i).dot(gradD(jl,j)) << " " << hij << " : " << pairs[k].f_couple << std::endl;
    // }
    CHECK(pairs[k].f_couple >= 0.0 and pairs[k].f_couple <= 1.0);
  }
  TIME_END("DamageGradientCoupling_pairs");
  TIME_END("DamageGradientCoupling");
}

}
