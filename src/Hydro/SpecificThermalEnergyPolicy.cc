//---------------------------------Spheral++----------------------------------//
// SpecificThermalEnergyPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the specific thermal energy as a dependent quantity.
// 
// This version is specialized for the compatible energy discretization 
// method.
//
// Created by JMO, Tue Sep 14 22:27:08 2004
//----------------------------------------------------------------------------//
#include "SpecificThermalEnergyPolicy.hh"
#include "HydroFieldNames.hh"
#include "entropyWeightingFunction.hh"
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
template<typename Dimension>
SpecificThermalEnergyPolicy<Dimension>::
SpecificThermalEnergyPolicy(const DataBase<Dimension>& dataBase):
  UpdatePolicyBase<Dimension>(),
  mDataBasePtr(&dataBase) {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SpecificThermalEnergyPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {

//   // HACK!
//   std::cerr.setf(std::ios::scientific, std::ios::floatfield);
//   std::cerr.precision(15);

  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == HydroFieldNames::specificThermalEnergy and 
          nodeListKey == UpdatePolicyBase<Dimension>::wildcard());
  auto eps = state.fields(fieldKey, Scalar());
  const auto numFields = eps.numFields();

  // Get the state field lists
  const auto  mass = state.fields(HydroFieldNames::mass, Scalar());
  const auto  velocity = state.fields(HydroFieldNames::velocity, Vector::zero());
  const auto  DvDt = derivs.fields(HydroFieldNames::hydroAcceleration, Vector::zero());
  const auto& pairAccelerations = derivs.template get<PairwiseField<Dimension, Vector>>(HydroFieldNames::pairAccelerations);
  const auto  selfAccelerations = derivs.fields(HydroFieldNames::selfAccelerations, Vector::zero(), true);
  const auto  DepsDt0 = derivs.fields(IncrementState<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
  const auto& connectivityMap = mDataBasePtr->connectivityMap();
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();
  CHECK(pairAccelerations.size() == npairs);
  CHECK(selfAccelerations.numFields() == 0 or selfAccelerations.numFields() == numFields);
  const bool selfInteraction = selfAccelerations.numFields() == numFields;

  const auto hdt = 0.5*multiplier;
  auto DepsDt = mDataBasePtr->newFluidFieldList(0.0, "delta E");

  // Check that the partial accelerations sum to the total hydro acceleration, or this isn't going to conserve
  BEGIN_CONTRACT_SCOPE {
    auto DvDt_check = mDataBasePtr->newFluidFieldList(Vector::zero(), "hydro acceleration check");
    for (auto kk = 0u; kk < npairs; ++kk) {
      const auto i = pairs[kk].i_node;
      const auto j = pairs[kk].j_node;
      const auto nodeListi = pairs[kk].i_list;
      const auto nodeListj = pairs[kk].j_list;
      const auto& paccij = pairAccelerations[kk];
      const auto  mi = mass(nodeListi, i);
      const auto  mj = mass(nodeListj, j);
      DvDt_check(nodeListi, i) += paccij;
      DvDt_check(nodeListj, j) -= paccij * mi/mj;
    }
    const auto numNodeLists = mDataBasePtr->numFluidNodeLists();
    for (auto k = 0u; k < numNodeLists; ++k) {
      const auto n = DvDt_check[k]->numInternalElements();
      for (auto i = 0u; i < n; ++i) {
        CHECK2(fuzzyEqual(DvDt_check(k, i).dot(DvDt(k, i)), DvDt(k, i).magnitude2(), 1.0e-8),
               DvDt_check(k, i) << " != " << DvDt(k, i) << " for (NodeList,i) = " << k << " " << i);
      }
    }
  }
  END_CONTRACT_SCOPE;

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

      // State for node i.
      const auto  mi = mass(nodeListi, i);
      const auto& vi = velocity(nodeListi, i);
      const auto& ai = DvDt(nodeListi, i);
      const auto  vi12 = vi + ai*hdt;
      const auto& paccij = pairAccelerations[kk];

      // State for node j.
      const auto  mj = mass(nodeListj, j);
      const auto& vj = velocity(nodeListj, j);
      const auto& aj = DvDt(nodeListj, j);
      const auto  vj12 = vj + aj*hdt;

      const auto  vji12 = vj12 - vi12;
      const Scalar duij = vji12.dot(paccij);
      const auto weighti = std::max(numeric_limits<Scalar>::epsilon(), DepsDt0(nodeListi, i)*sgn(duij));
      const auto weightj = std::max(numeric_limits<Scalar>::epsilon(), DepsDt0(nodeListj, j)*sgn(duij));
      const Scalar wi = weighti/(weighti + weightj);         // Du/Dt weighting
      // const Scalar wi = entropyWeighting(si, sj, duij);   // entropy weighting
      CHECK(wi >= 0.0 and wi <= 1.0);

      DepsDt_thread(nodeListi, i) += wi*duij;
      DepsDt_thread(nodeListj, j) += (1.0 - wi)*duij*mi/mj;
    }

#pragma omp critical
    {
      DepsDt_thread.threadReduce();
    }
  }

  // Now we can update the energy.
  auto offset = npairs;
  for (auto nodeListi = 0u; nodeListi < numFields; ++nodeListi) {
    const auto n = eps[nodeListi]->numInternalElements();
#pragma omp parallel for
    for (auto i = 0u; i < n; ++i) {

      // Add the self-contribution if any
      if (selfInteraction) {
        const auto& vi = velocity(nodeListi, i);
        const auto& ai = DvDt(nodeListi, i);
        const auto  vi12 = vi + ai*hdt;
        const auto duii = -vi12.dot(selfAccelerations(nodeListi, i));
        DepsDt(nodeListi, i) += duii;
      }

      eps(nodeListi, i) += DepsDt(nodeListi, i)*multiplier;
    }
    offset += n;
  }
}

//------------------------------------------------------------------------------
// Update the field using increments
//------------------------------------------------------------------------------
template<typename Dimension>
void
SpecificThermalEnergyPolicy<Dimension>::
updateAsIncrement(const KeyType& key,
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs,
                  const double multiplier,
                  const double t,
                  const double dt) {

  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == HydroFieldNames::specificThermalEnergy and 
          nodeListKey == UpdatePolicyBase<Dimension>::wildcard());
  auto eps = state.fields(fieldKey, Scalar());

  // Build an increment policy to use.
  IncrementState<Dimension, Scalar> fpolicy;

  // Do the deed for each of our Fields.
  for (auto fptr: eps) {
    fpolicy.updateAsIncrement(State<Dimension>::key(*fptr),
                              state, derivs, multiplier, t, dt);
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
SpecificThermalEnergyPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const auto rhsPtr = dynamic_cast<const SpecificThermalEnergyPolicy<Dimension>*>(&rhs);
  return (rhsPtr != nullptr);
}

}

