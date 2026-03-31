//---------------------------------Spheral++----------------------------------//
// RZNonSymmetricSpecificThermalEnergyPolicy -- An implementation of 
// UpdatePolicyBase specialized for the updating the specific thermal energy
// as a dependent quantity.
// 
// This version is for use with RZ axisymmetric symmetry.
//
// Created by JMO, Mon Feb  3 21:23:11 PST 2020
//----------------------------------------------------------------------------//
#include "RZNonSymmetricSpecificThermalEnergyPolicy.hh"
#include "HydroFieldNames.hh"
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
#include "Geometry/toroidalVolume.hh"

#include <vector>
#include <limits>
#include <cmath>
using std::vector;
using std::numeric_limits;
using std::cerr;
using std::endl;

namespace Spheral {

// namespace {

// inline
// double
// integrate_vr_over_r(const double vr,
//                     const double r,
//                     const double ar,
//                     const double dt,
//                     const bool barf = false) {
//   const double tiny = 1.0e-10;
//   VERIFY(r > 0.0);
//   VERIFY(r + vr*dt > 0.0);
//   if (fuzzyEqual(vr, 0.0, tiny) and fuzzyEqual(ar, 0.0, tiny)) return 0.0;
//   const auto q = 4.0*ar*r - vr*vr;
//   if (barf) cerr << "q: " << q << endl;
//   if (fuzzyEqual(q, 0.0, 1.0e-10)) {
//     const auto a = vr*safeInv(2.0*ar);
//     if (fuzzyEqual(a, 0.0, tiny)) return 0.0;
//     // auto Xt =  [&](const double t) { return ar*FastMath::square(a + t); };
//     auto F0t = [&](const double t) { return -ar*safeInv(a + t); };
//     auto F1t = [&](const double t) { return ar*(log(std::abs(a + t)) + a*safeInv(a + t)); };
//     auto Ft =  [&](const double t) { return vr*F0t(t) + ar*F1t(t); };
//     return Ft(dt) - Ft(0.0);
//   } else {
//     auto Xt =  [&](const double t) { return r + vr*t + ar*t*t; };
//     auto F0t = [&](const double t) {
//                  const auto thpt = std::sqrt(std::abs(q));
//                  if (q > 0.0) {
//                    return 2.0/thpt*atan2(2.0*ar*t + vr, thpt);
//                  } else {
//                    CHECK(q < 0.0);
//                    const auto ack = 2.0*ar*t + vr;
//                    if (fuzzyEqual(ack, thpt, tiny)) return 0.0;
//                    return 1.0/thpt*log(std::abs((thpt - ack)*safeInvVar(thpt + ack)));
//                  }
//                };
//     auto F1t = [&](const double t) { return (log(Xt(t)) - vr*F0t(t))*safeInv(2.0*ar); };
//     auto Ft =  [&](const double t) { return vr*F0t(t) + ar*F1t(t); };
//     return Ft(dt) - Ft(0.0);
//   }  
// }

// }

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
RZNonSymmetricSpecificThermalEnergyPolicy::
RZNonSymmetricSpecificThermalEnergyPolicy(const DataBase<Dim<2>>& dataBase):
  UpdatePolicyBase<Dimension>(),
  mDataBasePtr(&dataBase) {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
void
RZNonSymmetricSpecificThermalEnergyPolicy::
update(const KeyType& key,
       State<Dim<2>>& state,
       StateDerivatives<Dim<2>>& derivs,
       const double multiplier,
       const double /*t*/,
       const double /*dt*/) {

//   // HACK!
//   std::cerr.setf(std::ios::scientific, std::ios::floatfield);
//   std::cerr.precision(15);

  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == HydroFieldNames::specificThermalEnergy and 
          nodeListKey == UpdatePolicyBase<Dimension>::wildcard());
  auto eps = state.fields(fieldKey, Scalar());
  const auto numFields = eps.numFields();
  const auto tiny = 1.0e-20;
  const auto nodeListPtrs = eps.nodeListPtrs();

  // Get the state fields.
  const auto  mass = state.fields(HydroFieldNames::mass, Scalar());
  const auto  velocity = state.fields(HydroFieldNames::velocity, Vector::zero());
  const auto  acceleration = derivs.fields(HydroFieldNames::hydroAcceleration, Vector::zero());
  const auto& pairAccelerations = derivs.template get<PairwiseField<Dimension, Vector, 2u>>(HydroFieldNames::pairAccelerations);
  const auto& pairWork = derivs.template get<PairwiseField<Dimension, Scalar, 4u>>(HydroFieldNames::pairWork);
  const auto  selfAccelerations = derivs.fields(HydroFieldNames::selfAccelerations, Vector::zero(), true);
  const auto& connectivityMap = mDataBasePtr->connectivityMap();
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();
  CHECK(pairAccelerations.size() == npairs);
  CHECK(pairWork.size() == npairs);
  CHECK(selfAccelerations.numFields() == 0 or selfAccelerations.numFields() == numFields);
  const bool selfInteraction = selfAccelerations.numFields() == numFields;

  const auto hdt = 0.5*multiplier;
  auto DepsDt = mDataBasePtr->newFluidFieldList(0.0, "delta E");

  // Walk all pairs and figure out the discrete work for each point
#pragma omp parallel
  {
    // Thread private variables
    auto DepsDt_thread = DepsDt.threadCopy();
    // auto dEtot_thread = 0.0;

#pragma omp for
    for (auto kk = 0u; kk < npairs; ++kk) {
      const auto i = pairs[kk].i_node;
      const auto j = pairs[kk].j_node;
      const auto nodeListi = pairs[kk].i_list;
      const auto nodeListj = pairs[kk].j_list;

      // State for node i.
      const auto  mi = mass(nodeListi, i);
      const auto& vi = velocity(nodeListi, i);
      const auto& ai = acceleration(nodeListi, i);
      const auto  vi12 = vi + ai*hdt;
      const auto& pacci = pairAccelerations[kk][0];
      const auto  pworkzi = pairWork[kk][0];
      const auto  pworkri = pairWork[kk][1];

      // State for node j.
      const auto  mj = mass(nodeListj, j);
      const auto& vj = velocity(nodeListj, j);
      const auto& aj = acceleration(nodeListj, j);
      const auto  vj12 = vj + aj*hdt;
      const auto& paccj = pairAccelerations[kk][1];
      const auto  pworkzj = pairWork[kk][2];
      const auto  pworkrj = pairWork[kk][3];

      // Correct just the r-component error (multiplicative form)
      const auto dErij = -(mi*vi12[1]*pacci[1] + mj*vj12[1]*paccj[1]);
      const auto chi = dErij*safeInv(mi*pworkri + mj*pworkrj, tiny);
      DepsDt_thread(nodeListi, i) += pworkzi + chi*pworkri;
      DepsDt_thread(nodeListj, j) += pworkzj + chi*pworkrj;

      // // Correct just the r-component error (additive form)
      // const auto durij = (-(mi*vi12[1]*pacci[1] + mj*vj12[1]*paccj[1]) - (mi*pworkri + mj*pworkrj))/(mi + mj);
      // DepsDt_thread(nodeListi, i) += pworkzi + pworkri + durij;
      // DepsDt_thread(nodeListj, j) += pworkzj + pworkrj + durij;
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

      // Add the self-contribution if any (RZ with strength does this for instance).
      if (selfInteraction) {
        const auto& vi = velocity(nodeListi, i);
        const auto& ai = acceleration(nodeListi, i);
        const auto  vi12 = vi + ai*hdt;
        DepsDt(nodeListi, i) -= vi12.dot(selfAccelerations(nodeListi, i));
      }

      eps(nodeListi, i) += DepsDt(nodeListi, i)*multiplier;
    }
  }
}

//------------------------------------------------------------------------------
// Update the field using increments
//------------------------------------------------------------------------------
void
RZNonSymmetricSpecificThermalEnergyPolicy::
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
bool
RZNonSymmetricSpecificThermalEnergyPolicy::
operator==(const UpdatePolicyBase<Dim<2>>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const auto* rhsPtr = dynamic_cast<const RZNonSymmetricSpecificThermalEnergyPolicy*>(&rhs);
  return (rhsPtr != nullptr);
}

}

