//---------------------------------Spheral++----------------------------------//
// RZCompatibleDifferenceSpecificThermalEnergyPolicy -- An implementation of 
// UpdatePolicyBase specialized for the updating the specific thermal energy 
// as a dependent quantity.
// 
// This version is specialized for materials with different properties. A 
// compatible energy discretization in which pairwise work allows for opposite
// sign pair-wise work. DepsDti and  DepsDtj are used as weights and the 
// difference between the conservative and consistent formulations is added 
// back in.
//
// This version is for use with RZ axisymmetric symmetry.
//
//----------------------------------------------------------------------------//
#include "Hydro/RZCompatibleDifferenceSpecificThermalEnergyPolicy.hh"
#include "Hydro/HydroFieldNames.hh"
#include "NodeList/NodeList.hh"
#include "NodeList/FluidNodeList.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/IncrementState.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Neighbor/PairwiseField.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Utilities/DBC.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/SpheralFunctions.hh"

#include <vector>
#include <limits>
using std::vector;
using std::numeric_limits;
using std::abs;
using std::min;
using std::max;

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
RZCompatibleDifferenceSpecificThermalEnergyPolicy::
RZCompatibleDifferenceSpecificThermalEnergyPolicy(const DataBase<Dim<2>>& dataBase):
  UpdatePolicyBase<Dim<2>>(),
  mDataBasePtr(&dataBase) {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
void
RZCompatibleDifferenceSpecificThermalEnergyPolicy::
update(const KeyType& key,
       State<Dim<2>>& state,
       StateDerivatives<Dim<2>>& derivs,
       const double multiplier,
       const double /*t*/,
       const double /*dt*/) {

  using PairAccelerationsType = PairwiseField<Dim<2>, Vector>;
  using PairWorkType = PairwiseField<Dim<2>, Scalar, 2u>;

  KeyType fieldKey, nodeListKey;
  StateBase<Dim<2>>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == HydroFieldNames::specificThermalEnergy and 
          nodeListKey == UpdatePolicyBase<Dim<2>>::wildcard());
  auto eps = state.fields(fieldKey, Scalar());
  const auto numFields = eps.numFields();

  // constant we'll need for the weighting scheme
  const auto tiny = numeric_limits<Scalar>::epsilon();

  // Get the state fields.
  const auto  pos = state.fields(HydroFieldNames::position, Vector::zero);
  const auto  mass = state.fields(HydroFieldNames::mass, Scalar());
  const auto  velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const auto  acceleration = derivs.fields(HydroFieldNames::hydroAcceleration, Vector::zero);
  const auto& pairAccelerations = derivs.template get<PairAccelerationsType>(HydroFieldNames::pairAccelerations);
  const auto& pairDepsDt = derivs.template get<PairWorkType>(HydroFieldNames::pairWork);
  const auto& connectivityMap = mDataBasePtr->connectivityMap();
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();

  CHECK(pairAccelerations.size() == npairs);
  CHECK(pairDepsDt.size() == npairs);

  auto  DepsDt = derivs.fields(IncrementState<Dim<2>, Field<Dim<2>, Scalar> >::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
  DepsDt.Zero();

  const auto hdt = 0.5*multiplier;
  
  // Walk all pairs and figure out the discrete work for each point
#pragma omp parallel
  {
    // Thread private variables
    auto DepsDt_thread = DepsDt.threadCopy();

#pragma omp for
    for (auto kk = 0u; kk < npairs; ++kk) {
      const auto i = pairs[kk].i_node;
      const auto j = pairs[kk].j_node;
      const auto nodeListi = pairs[kk].i_list;
      const auto nodeListj = pairs[kk].j_list;

      const auto& paccij = pairAccelerations[kk];
      const auto& DepsDt0i = pairDepsDt[kk][0];
      const auto& DepsDt0j = pairDepsDt[kk][1];

      const auto  ri = abs(pos(nodeListi, i).y());
      const auto  mi = mass(nodeListi, i)/(2.0*M_PI*ri);
      const auto& vi = velocity(nodeListi, i);
      const auto& ai = acceleration(nodeListi, i);

      const auto  rj = abs(pos(nodeListj, j).y());
      const auto  mj = mass(nodeListj, j)/(2.0*M_PI*rj);
      const auto& vj = velocity(nodeListj, j);
      const auto& aj = acceleration(nodeListj, j);

      // half-step velocity
      const auto vi12 = vi + ai*hdt;
      const auto vj12 = vj + aj*hdt;
      const auto vij = vi12 - vj12;

      // weighting scheme
      const auto weighti = abs(DepsDt0i) + tiny;
      const auto weightj = abs(DepsDt0j) + tiny;
      const auto wi = weighti/(weighti+weightj);

      // difference between assessed derivs and conserative ones
      const Scalar delta_duij = vij.dot(paccij)-DepsDt0i-DepsDt0j;

      CHECK(wi >= 0.0 and wi <= 1.0);

      // make conservative
      DepsDt_thread(nodeListi, i) += mj*(wi*delta_duij+DepsDt0i);
      DepsDt_thread(nodeListj, j) += mi*((1.0-wi)*delta_duij+DepsDt0j);

    }

#pragma omp critical
    {
      DepsDt_thread.threadReduce();
    }
  }

  // Now we can update the energy.
  for (auto nodeListi = 0u; nodeListi < numFields; ++nodeListi) {
    const auto n = eps[nodeListi]->numInternalElements();
#pragma omp parallel for
    for (auto i = 0u; i < n; ++i) {
      eps(nodeListi, i) += DepsDt(nodeListi, i)*multiplier;
    }
  }
}

//------------------------------------------------------------------------------
// Update the field using increments
//------------------------------------------------------------------------------
void
RZCompatibleDifferenceSpecificThermalEnergyPolicy::
updateAsIncrement(const KeyType& key,
                  State<Dim<2>>& state,
                  StateDerivatives<Dim<2>>& derivs,
                  const double multiplier,
                  const double t,
                  const double dt) {

  KeyType fieldKey, nodeListKey;
  StateBase<Dim<2>>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == HydroFieldNames::specificThermalEnergy and 
          nodeListKey == UpdatePolicyBase<Dim<2>>::wildcard());
  auto eps = state.fields(fieldKey, Scalar());

  // Build an increment policy to use.
  IncrementState<Dim<2>, Scalar> fpolicy;

  // Do the deed for each of our Fields.
  for (auto fptr: eps) {
    fpolicy.updateAsIncrement(State<Dim<2>>::key(*fptr),
                              state, derivs, multiplier, t, dt);
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
bool
RZCompatibleDifferenceSpecificThermalEnergyPolicy::
operator==(const UpdatePolicyBase<Dim<2>>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const auto* rhsPtr = dynamic_cast<const RZCompatibleDifferenceSpecificThermalEnergyPolicy*>(&rhs);
  return (rhsPtr != nullptr);
}

}

