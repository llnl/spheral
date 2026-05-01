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

#ifdef SPHERAL_ENABLE_MPI
#include "Distributed/TreeDistributedBoundary.hh"
#endif
#include "Utilities/globalNodeIDs.hh"

#include <vector>
#include <limits>
#include <cmath>
using std::vector;
using std::numeric_limits;
using std::cerr;
using std::endl;

namespace Spheral {

namespace {

template<typename Dimension>
inline
bool
isCommunicatedNode(const NodeList<Dimension>& nodes,
                   const size_t i) {
#ifdef SPHERAL_ENABLE_MPI
  const auto& bound = TreeDistributedBoundary<Dimension>::instance();
  return bound.isGhostNode(nodes, i);
#else
  return false;
#endif
}

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

}

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

  auto globalIDs = globalNodeIDs(*mDataBasePtr);
  const auto& bound = TreeDistributedBoundary<Dimension>::instance();
  bound.applyFieldListGhostBoundary(globalIDs);
  bound.finalizeGhostBoundary();

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
  const auto  massRZ = state.fields(HydroFieldNames::massRZ, Scalar());
  const auto  position = state.fields(HydroFieldNames::position, Vector::zero());
  const auto  H = state.fields(HydroFieldNames::H, SymTensor::zero());
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
  auto DepsDt = mDataBasePtr->newFluidFieldList(0.0, "DepsDt");
  auto deltaE0 = 0.0;
  auto deltaEconserve = 0.0;
  auto wsum = 0.0;
  auto pairCount = 0.0;

  // const auto numInternalNodesPerNodeList = mDataBasePtr->numInternalNodesPerFluidNodeList();
  // CHECK(numInternalNodesPerNodeList.size() == numFields);
  // auto isGhost = [&](const size_t k, const size_t i) { return false; }; //return i >= numInternalNodesPerNodeList[k]; };

  // Walk all pairs to find the total energy changes
#pragma omp parallel
  {
    // Thread private variables
    auto deltaEconserve_thread = 0.0;
    auto deltaE0_thread = 0.0;
    auto wsum_thread = 0.0;
    auto pairCount_thread = 0.0;

#pragma omp for
    for (auto kk = 0u; kk < npairs; ++kk) {
      const auto i = pairs[kk].i_node;
      const auto j = pairs[kk].j_node;
      const auto nodeListi = pairs[kk].i_list;
      const auto nodeListj = pairs[kk].j_list;

      // const auto ifactor = isGhost(nodeListi, i) ? 0.0 : 1.0;
      // const auto jfactor = isGhost(nodeListj, j) ? 0.0 : 1.0;
      const auto ifactor = isCommunicatedNode(*nodeListPtrs[nodeListi], i) ? 0.0 : 1.0;
      const auto jfactor = isCommunicatedNode(*nodeListPtrs[nodeListj], j) ? 0.0 : 1.0;
      const auto pairScale = 0.5*(ifactor + jfactor);
      CHECK(fuzzyEqual(pairScale, 0.5) or fuzzyEqual(pairScale, 1.0));
      pairCount_thread += pairScale;

      // const auto gIDi = globalIDs(nodeListi, i);
      // const auto gIDj = globalIDs(nodeListj, j);
      // // if (ifactor + jfactor < 2.0) std::cerr << "--> (" << globalIDs(nodeListi, i) << " " << globalIDs(nodeListj, j) << ") (" << ifactor << " " << jfactor << ")\n";
      // if ((gIDi == 297 or gIDi == 1000) and (gIDj == 297 or gIDj == 1000)) std::cerr << "--> (" << gIDi << " " << gIDj << ") (" << ifactor << " " << jfactor << ")\n";
      
      // State for node i.
      const auto  mi = mass(nodeListi, i);
      const auto& vi = velocity(nodeListi, i);
      const auto& ai = acceleration(nodeListi, i);
      const auto  vi12 = vi + ai*hdt;
      const auto& pacci = pairAccelerations[kk][0];
      const auto  pworkzi = pairWork[kk][0];
      const auto  pworkri = pairWork[kk][1];
      const auto  epsi = eps(nodeListi, i);

      // State for node j.
      const auto  mj = mass(nodeListj, j);
      const auto& vj = velocity(nodeListj, j);
      const auto& aj = acceleration(nodeListj, j);
      const auto  vj12 = vj + aj*hdt;
      const auto& paccj = pairAccelerations[kk][1];
      const auto  pworkzj = pairWork[kk][2];
      const auto  pworkrj = pairWork[kk][3];
      const auto  epsj = eps(nodeListj, j);

      // Exact conservation energy change
      const auto dEij = -(mi*vi12.dot(pacci) + mj*vj12.dot(paccj));
      const auto dE0ij = (mi*(pworkzi + pworkri) + mj*(pworkzj + pworkrj));
      deltaEconserve_thread += dEij*pairScale;
      deltaE0_thread += dE0ij*pairScale;
      wsum_thread += (std::abs(epsi) + std::abs(epsj) + tiny)*pairScale;
    }

#pragma omp critical
    {
      deltaEconserve += deltaEconserve_thread;
      deltaE0 += deltaE0_thread;
      wsum += wsum_thread;
      pairCount += pairCount_thread;
    }
  }

  // Find the global delta for conservation
  deltaEconserve = allReduce(deltaEconserve, SPHERAL_OP_SUM);
  deltaE0 = allReduce(deltaE0, SPHERAL_OP_SUM);
  wsum = allReduce(wsum, SPHERAL_OP_SUM);
  CHECK(wsum > 0.0);
  const auto dEtot = (deltaEconserve - deltaE0)/wsum;

  pairCount = allReduce(pairCount, SPHERAL_OP_SUM);
  if (Process::getRank() == 0) std::cerr << "PAIR COUNT: " << pairCount << " " << wsum << " " << deltaEconserve << " " << deltaE0 << std::endl;

  // Walk all pairs again and update the energy derivative
  auto dEcheck = 0.0;
  auto wsumCheck = 0.0;
#pragma omp parallel
  {
    // Thread private variables
    auto DepsDt_thread = DepsDt.threadCopy();
    auto dEcheck_thread = 0.0;
    auto wsumCheck_thread = 0.0;

#pragma omp for
    for (auto kk = 0u; kk < npairs; ++kk) {
      const auto i = pairs[kk].i_node;
      const auto j = pairs[kk].j_node;
      const auto nodeListi = pairs[kk].i_list;
      const auto nodeListj = pairs[kk].j_list;

      // const auto ifactor = isGhost(nodeListi, i) ? 0.0 : 1.0;
      // const auto jfactor = isGhost(nodeListj, j) ? 0.0 : 1.0;
      const auto ifactor = isCommunicatedNode(*nodeListPtrs[nodeListi], i) ? 0.0 : 1.0;
      const auto jfactor = isCommunicatedNode(*nodeListPtrs[nodeListj], j) ? 0.0 : 1.0;
      const auto pairScale = 0.5*(ifactor + jfactor);
      CHECK(fuzzyEqual(pairScale, 0.5) or fuzzyEqual(pairScale, 1.0));

      // State for node i.
      const auto  mi = mass(nodeListi, i);
      const auto  pworkzi = pairWork[kk][0];
      const auto  pworkri = pairWork[kk][1];
      const auto  epsi = eps(nodeListi, i);

      // State for node j.
      const auto  mj = mass(nodeListj, j);
      const auto  pworkzj = pairWork[kk][2];
      const auto  pworkrj = pairWork[kk][3];
      const auto  epsj = eps(nodeListj, j);

      // Update the energy derivative
      const auto dE0ij = mi*(pworkzi + pworkri) + mj*(pworkzj + pworkrj);
      const auto deltaij = dEtot*(std::abs(dE0ij) + tiny);
      dEcheck_thread += deltaij*pairScale;
      wsumCheck_thread += (std::abs(epsi) + std::abs(epsj) + tiny)*pairScale;
      // wsumCheck_thread += (std::abs(dE0ij) + tiny)*pairScale;
      auto wi = 0.5; //mi*std::abs(pworkzi + pworkri);
      auto wj = 0.5; //mj*std::abs(pworkzj + pworkrj);
      const auto thpt = safeInv(wi + wj, tiny);
      wi *= thpt;
      wj *= thpt;
      if (not fuzzyEqual(wi + wj, 1.0, 1.0e-10)) {
        wi = 0.5;
        wj = 0.5;
      }
      DepsDt_thread(nodeListi, i) += pworkzi + pworkri + wi*deltaij/mi;
      DepsDt_thread(nodeListj, j) += pworkzj + pworkrj + wj*deltaij/mj;

      // const auto duij = deltaij/(mi + mj);
      // DepsDt_thread(nodeListi, i) += pworkzi + pworkri + duij;
      // DepsDt_thread(nodeListj, j) += pworkzj + pworkrj + duij;
    }

#pragma omp critical
    {
      DepsDt_thread.threadReduce();
      dEcheck += dEcheck_thread;
      wsumCheck += wsumCheck_thread;
    }
  }
  dEcheck = allReduce(dEcheck, SPHERAL_OP_SUM);
  wsumCheck = allReduce(wsumCheck, SPHERAL_OP_SUM);
  // VERIFY2(fuzzyEqual(dEcheck, deltaEconserve - deltaE0, 1.0e-10),
  //         "Bad energy correction: " << dEcheck << " " << (deltaEconserve - deltaE0) << " : " << wsum << " " << wsumCheck);
  // VERIFY2(fuzzyEqual(wsumCheck, wsum, 1.0e-10),
  //         "Bad wsum correction: " << dEcheck << " " << (deltaEconserve - deltaE0) << " : " << wsum << " " << wsumCheck);

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

