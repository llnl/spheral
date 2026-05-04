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
#include "Utilities/range.hh"
#include "Geometry/toroidalVolume.hh"

#ifdef SPHERAL_ENABLE_MPI
#include "Distributed/TreeDistributedBoundary.hh"
#endif

#include <vector>
#include <limits>
#include <cmath>
using std::vector;
using std::numeric_limits;
using std::cerr;
using std::endl;

namespace Spheral {

namespace {

//------------------------------------------------------------------------------
// Functor class to check if a point is a communicated ghost point from another
// domain.
// This implementation assumes the ghost point indices are contiguous, and that
// we're using TreeDistributedBoundary.
//------------------------------------------------------------------------------
template<typename Dimension>
class CheckCommunicatedNodes {
public:
  CheckCommunicatedNodes(const std::vector<NodeList<Dimension>*>& nodeListPtrs):
    mGhostRanges(nodeListPtrs.size(), std::pair<size_t, size_t>(0u, 0u)) {
#ifdef SPHERAL_ENABLE_MPI
    if (Process::getTotalNumberOfProcesses() > 1) {
      const auto& bound = TreeDistributedBoundary<Dimension>::instance();
      for (const auto [k, nptr]: enumerate(nodeListPtrs)) {
        const auto& ghosts = bound.ghostNodes(*nptr);
        const auto [imin, imax] = std::minmax_element(ghosts.begin(), ghosts.end());
        CHECK(imin < ghosts.end());
        CHECK(imax < ghosts.end());
        mGhostRanges[k] = std::make_pair(*imin, *imax);
      }
    }
#endif
  }
  ~CheckCommunicatedNodes() = default;
  bool operator()(const size_t nodeListID,
                  const size_t i) const {
    REQUIRE(nodeListID < mGhostRanges.size());
    const auto [imin, imax] = mGhostRanges[nodeListID];
    return (i >=  imin and i < imax);
  }

private:
  std::vector<std::pair<size_t, size_t>> mGhostRanges;
};

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

  // // HACK!
  // std::cerr.setf(std::ios::scientific, std::ios::floatfield);
  // std::cerr.precision(15);

  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == HydroFieldNames::specificThermalEnergy and 
          nodeListKey == UpdatePolicyBase<Dimension>::wildcard());
  auto eps = state.fields(fieldKey, Scalar());
  const auto numFields = eps.numFields();
  const auto nodeListPtrs = eps.nodeListPtrs();

  // Function to check if the given node is a ghost node
  const auto numInternalNodesPerNodeList = mDataBasePtr->numInternalNodesPerFluidNodeList();
  CHECK(numInternalNodesPerNodeList.size() == numFields);
  auto isGhost = [&](const size_t k, const size_t i) { return i >= numInternalNodesPerNodeList[k]; };

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
  auto weight = mDataBasePtr->newFluidFieldList(0.0, "weight");
  auto deltaEconserve = 0.0;
  auto wsum = 0.0;

  // Walk all pairs to find the total energy changes
#pragma omp parallel
  {
    // Thread private variables
    typename SpheralThreads<Dim<2>>::FieldListStack threadStack;
    auto DepsDt_thread = DepsDt.threadCopy(threadStack);
    auto weight_thread = weight.threadCopy(threadStack);
    auto deltaEconserve_thread = 0.0;
    auto wsum_thread = 0.0;

#pragma omp for
    for (auto kk = 0u; kk < npairs; ++kk) {
      const auto i = pairs[kk].i_node;
      const auto j = pairs[kk].j_node;
      const auto nodeListi = pairs[kk].i_list;
      const auto nodeListj = pairs[kk].j_list;

      const auto ifactor = isGhost(nodeListi, i) ? 0.0 : 1.0;
      const auto jfactor = isGhost(nodeListj, j) ? 0.0 : 1.0;

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

      // Exact conservation energy change
      const auto dEij = -(mi*vi12.dot(pacci) + mj*vj12.dot(paccj));
      const auto dE0ij = (mi*(pworkzi + pworkri) + mj*(pworkzj + pworkrj));
      const auto duij = (dEij - dE0ij)/(mi + mj);
      const auto wi = mi*(std::abs(pworkzi) + std::abs(pworkri))*ifactor;
      const auto wj = mj*(std::abs(pworkzj) + std::abs(pworkrj))*jfactor;
      deltaEconserve_thread += mi*duij*ifactor + mj*duij*jfactor;
      wsum_thread += wi + wj;

      // Add the uncorrected energy change vote
      DepsDt_thread(nodeListi, i) += pworkzi + pworkri;
      DepsDt_thread(nodeListj, j) += pworkzj + pworkrj;
      weight_thread(nodeListi, i) += wi;
      weight_thread(nodeListj, j) += wj;
    }

#pragma omp critical
    {
      deltaEconserve += deltaEconserve_thread;
      wsum += wsum_thread;
    }
    threadReduceFieldLists<Dimension>(threadStack);
  }

  // Find the global delta for conservation
  deltaEconserve = allReduce(deltaEconserve, SPHERAL_OP_SUM);
  wsum = allReduce(wsum, SPHERAL_OP_SUM);
  CHECK(wsum >= 0.0);
  const auto dEtot = wsum > 0.0 ? deltaEconserve/wsum : 0.0;

  // Go through the internal points, add the conservative energy fix and update the final energy
  auto dEcheck = 0.0;
  auto wsumCheck = 0.0;
  for (auto nodeListi = 0u; nodeListi < numFields; ++nodeListi) {
    const auto n = eps[nodeListi]->numInternalElements();
#pragma omp parallel for
    for (auto i = 0u; i < n; ++i) {
      auto& DepsDti = DepsDt(nodeListi, i);

      // Add the conservative fix to this point
      const auto mi = mass(nodeListi, i);
      const auto wi = weight(nodeListi, i);
      DepsDti += wi*dEtot/mi;

      // Add the self-contribution if any (RZ with strength does this for instance).
      if (selfInteraction) {
        const auto& vi = velocity(nodeListi, i);
        const auto& ai = acceleration(nodeListi, i);
        const auto  vi12 = vi + ai*hdt;
        DepsDti -= vi12.dot(selfAccelerations(nodeListi, i));
      }

      // Update the specific thermal energy
      eps(nodeListi, i) += DepsDti*multiplier;

#pragma omp critical
      {
        wsumCheck += wi;
        dEcheck += wi*dEtot;
      }
    }
  }

  dEcheck = allReduce(dEcheck, SPHERAL_OP_SUM);
  wsumCheck = allReduce(wsumCheck, SPHERAL_OP_SUM);
  // if (Process::getRank() == 0) printf("wsum=%f %f      dE=%f %f\n", wsum, wsumCheck, deltaEconserve, dEcheck);
  CHECK2(wsum == 0.0 or fuzzyEqual(dEcheck, deltaEconserve, 1.0e-10),
         "Bad energy correction: " << dEcheck << " " << deltaEconserve << " : " << wsum << " " << wsumCheck);
  CHECK2(fuzzyEqual(wsumCheck, wsum, 1.0e-10),
         "Bad wsum correction: " << dEcheck << " " << deltaEconserve << " : " << wsum << " " << wsumCheck);
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

