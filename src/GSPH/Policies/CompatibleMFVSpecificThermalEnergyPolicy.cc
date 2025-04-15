//---------------------------------Spheral++----------------------------------//
// CompatibleMFVSpecificThermalEnergyPolicy -- This is a generalization of the 
//     Lagrangian compatible energy scheme to ALE-based scheme with mass flux
//     between nodes. 
//
// J.M. Pearl 2023
//----------------------------------------------------------------------------//

#include "GSPH/Policies/CompatibleMFVSpecificThermalEnergyPolicy.hh"
#include "GSPH/GSPHFieldNames.hh"

#include "Hydro/HydroFieldNames.hh"

#include "NodeList/NodeList.hh"
#include "NodeList/FluidNodeList.hh"

#include "DataBase/IncrementState.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"

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
CompatibleMFVSpecificThermalEnergyPolicy<Dimension>::
CompatibleMFVSpecificThermalEnergyPolicy(const DataBase<Dimension>& dataBase):
  UpdatePolicyBase<Dimension>(),
  mDataBasePtr(&dataBase){
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CompatibleMFVSpecificThermalEnergyPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double /*t*/,
       const double /*dt*/) {

  using PairAccelerationsType = PairwiseField<Dimension, Vector>;
  using PairWorkType = PairwiseField<Dimension, Scalar, 2u>;
  using PairMassFluxType = PairwiseField<Dimension, Scalar, 1u>;

  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == HydroFieldNames::specificThermalEnergy and 
          nodeListKey == UpdatePolicyBase<Dimension>::wildcard());
  auto eps = state.fields(fieldKey, Scalar());
  const auto numFields = eps.numFields();

  // constant we'll need for the weighting scheme
  const auto tiny = numeric_limits<Scalar>::epsilon();

  // Get the state fields.
  const auto  mass = state.fields(HydroFieldNames::mass, Scalar());
  const auto  velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const auto  DmassDt = derivs.fields(IncrementState<Dimension, Vector>::prefix() + HydroFieldNames::mass, 0.0);
  const auto  DmomentumDt = derivs.fields(IncrementState<Dimension, Vector>::prefix() + GSPHFieldNames::momentum, Vector::zero);
  const auto& pairAccelerations = derivs.template get<PairAccelerationsType>(HydroFieldNames::pairAccelerations);
  const auto& pairDepsDt = derivs.template get<PairWorkType>(HydroFieldNames::pairWork);
  const auto& pairMassFlux = derivs.template get<PairMassFluxType>(GSPHFieldNames::pairMassFlux);

  const auto& connectivityMap = mDataBasePtr->connectivityMap();
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();

  CHECK(pairAccelerations.size() == npairs);
  CHECK(pairMassFlux.size() == npairs);
  CHECK(pairDepsDt.size() == npairs);

  auto  DepsDt = derivs.fields(IncrementState<Dimension, Scalar >::prefix() + GSPHFieldNames::thermalEnergy, 0.0);
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
      const auto& massFlux = pairMassFlux[kk];

      const auto  mi = mass(nodeListi, i);
      const auto  pi = mi*velocity(nodeListi, i);
      const auto& DPDti = DmomentumDt(nodeListi, i);
      const auto& DmDti = DmassDt(nodeListi,i);

      const auto  mj = mass(nodeListj, j);
      const auto  pj = mj*velocity(nodeListj, j);
      const auto& DPDtj = DmomentumDt(nodeListj, j);
      const auto& DmDtj = DmassDt(nodeListj,j);

      // half-step momenta
      const auto pi12 = pi + DPDti*hdt;
      const auto pj12 = pj + DPDtj*hdt;
      //const auto pij = pi12 - pj12;
      
      // weighting scheme
      const auto weighti = abs(DepsDt0i) + tiny;
      const auto weightj = abs(DepsDt0j) + tiny;
      const auto wi = weighti/(weighti+weightj);

      // safeInv
      const auto invmi0 = safeInv(mi);
      const auto invmj0 = safeInv(mj);
      const auto invmi1 = safeInv(mi+DmDti*multiplier);
      const auto invmj1 = safeInv(mj+DmDtj*multiplier);

      const Scalar delta_duij = (pi12*invmi1 - pj12*invmj1).dot(paccij)
                              + (pj.dot(pj)*invmj0*invmj1 - pi.dot(pi)*invmi0*invmi1) * massFlux*0.5
                              - DepsDt0i-DepsDt0j;

      CHECK(wi >= 0.0 and wi <= 1.0);
      CHECK(invmi0 >= 0.0);
      CHECK(invmj0 >= 1.0);

      const auto depsi =      (wi *delta_duij+DepsDt0i);
      const auto depsj = ((1.0-wi)*delta_duij+DepsDt0j);

      // make conservative
      DepsDt_thread(nodeListi, i) += depsi;
      DepsDt_thread(nodeListj, j) += depsj;

    }

#pragma omp critical
    {
      DepsDt_thread.threadReduce();
    }
  }

//   // Now we can update the energy.
  for (auto nodeListi = 0u; nodeListi < numFields; ++nodeListi) {
    const auto n = eps[nodeListi]->numInternalElements();
#pragma omp parallel for
    for (auto i = 0u; i < n; ++i) {
      const auto m1 = mass(nodeListi,i)+DmassDt(nodeListi,i)*multiplier;
      if (m1 > tiny) eps(nodeListi, i) += (DepsDt(nodeListi, i) - DmassDt(nodeListi, i)*eps(nodeListi, i)) * multiplier * safeInv(m1);
    }
  }

}

//------------------------------------------------------------------------------
// Update the field using increments
//------------------------------------------------------------------------------
template<typename Dimension>
void
CompatibleMFVSpecificThermalEnergyPolicy<Dimension>::
updateAsIncrement(const KeyType& key,
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs,
                  const double multiplier,
                  const double t,
                  const double dt) {

  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE2(fieldKey == HydroFieldNames::specificThermalEnergy and 
           nodeListKey == UpdatePolicyBase<Dimension>::wildcard(),
           "Bad key choice: " << key << " " << fieldKey << " " << nodeListKey);
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
CompatibleMFVSpecificThermalEnergyPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const auto* rhsPtr = dynamic_cast<const CompatibleMFVSpecificThermalEnergyPolicy<Dimension>*>(&rhs);
  return rhsPtr != nullptr;
}

}

