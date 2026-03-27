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
  const auto tiny = 1.0e-10;
  const auto nodeListPtrs = eps.nodeListPtrs();

  // Get the state fields.
  const auto  pos = state.fields(HydroFieldNames::position, Vector::zero());
  const auto  mass = state.fields(HydroFieldNames::mass, Scalar());
  const auto  massRZ = state.fields(HydroFieldNames::massRZ, Scalar());
  const auto  velocity = state.fields(HydroFieldNames::velocity, Vector::zero());
  // const auto  rho = state.fields(HydroFieldNames::massDensity, Scalar());
  // const auto  pressure = state.fields(HydroFieldNames::pressure, Scalar());
  const auto  acceleration = derivs.fields(HydroFieldNames::hydroAcceleration, Vector::zero());
  const auto& pairAccelerations = derivs.template get<PairwiseField<Dimension, Vector, 2u>>(HydroFieldNames::pairAccelerations);
  const auto& pairWork = derivs.template get<PairwiseField<Dimension, Scalar, 4u>>(HydroFieldNames::pairWork);
  const auto  selfAccelerations = derivs.fields(HydroFieldNames::selfAccelerations, Vector::zero(), true);
  const auto  DepsDt0 = derivs.fields(IncrementState<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
  const auto& connectivityMap = mDataBasePtr->connectivityMap();
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();
  CHECK(pairAccelerations.size() == npairs);
  CHECK(pairWork.size() == npairs);
  CHECK(selfAccelerations.numFields() == 0 or selfAccelerations.numFields() == numFields);
  const bool selfInteraction = selfAccelerations.numFields() == numFields;

  const auto hdt = 0.5*multiplier;
  auto DepsDt = mDataBasePtr->newFluidFieldList(0.0, "delta E");
  // auto dEtot = 0.0;

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
      const auto  ri = std::abs(pos(nodeListi, i).y());
      const auto  mRZi = massRZ(nodeListi, i);
      const auto  mi = mass(nodeListi, i);
      const auto  mri = 2.0*M_PI*ri*mRZi;
      const auto& vi = velocity(nodeListi, i);
      const auto& ai = acceleration(nodeListi, i);
      const auto  vi12 = vi + ai*hdt;
      const auto& pacci = pairAccelerations[kk][0];
      const auto  pworkzi = pairWork[kk][0];
      const auto  pworkri = pairWork[kk][1];
      const auto  pworki = pworkzi + pworkri;
      // const auto  vri = vi12.y();

      // State for node j.
      const auto  rj = std::abs(pos(nodeListj, j).y());
      const auto  mRZj = massRZ(nodeListj, j);
      const auto  mj = mass(nodeListj, j);
      const auto  mrj = 2.0*M_PI*rj*mRZj;
      // const auto  rhoj = rho(nodeListj, j);
      // const auto  Pj = pressure(nodeListj, j);
      const auto& vj = velocity(nodeListj, j);
      const auto& aj = acceleration(nodeListj, j);
      const auto  vj12 = vj + aj*hdt;
      const auto& paccj = pairAccelerations[kk][1];
      const auto  pworkzj = pairWork[kk][2];
      const auto  pworkrj = pairWork[kk][3];
      const auto  pworkj = pworkzj + pworkrj;
      // const auto  vrj = vj12.y();

      // Conservative energy definition
      const auto dEij = -(mi*vi12.dot(pacci) + mj*vj12.dot(paccj));
      // const auto duij = (dEij - (mi*pworki + mj*pworkj))/(mi + mj);

      // // Additive correction for RZ frame
      // const auto dERZij = -(mRZi*vi12.dot(pacci) + mRZj*vj12.dot(paccj));
      // const auto duij = (dERZij - (mRZi*pworki + mRZj*pworkj))/(mRZi + mRZj);

      // // Hybrid correction
      // const auto duZij = (-(mRZi*vi12[0]*pacci[0] + mRZj*vj12[0]*paccj[0]) - mRZi*pworkzi - mRZj*pworkzj)/(mRZi + mRZj);
      // const auto duRij = (-(mri* vi12[1]*pacci[1] + mrj* vj12[1]*paccj[1]) - mri* pworkri - mrj* pworkrj)/(mri + mrj);
      // // const auto duij = duRij + duZij + (dEij - mi*pworki - mj*pworkj - (mi + mj)*(duRij + duZij))*safeInv(mi + mj);
      // // VERIFY2(fuzzyEqual(mi*pworki + mj*pworkj + (mi + mj)*duij, dEij, tiny), "Energy balance error: (" << i << " " << j << "): " << dEij << " " << (mi*pworki + mj*pworkj + (mi + mj)*duij) << " " << duij << " " << duZij << " " << duRij);
      // const auto chi = (duij0 - duZij)*safeInv(duRij, tiny);
      // const auto duij = duZij + chi*duRij;

      // // Sum the conservative energy change
      // if (i < nodeListPtrs[nodeListi]->firstGhostNode()) dEtot -= mi*vi12.dot(pacci);
      // if (j < nodeListPtrs[nodeListj]->firstGhostNode()) dEtot -= mj*vj12.dot(paccj);

      // Now figure out the remaining total energy error, and apply that assuming
      // it's due to the r evolution term
      // const auto dEij = -(mi*vi12.dot(pacci) + mj*vj12.dot(paccj));
      // const auto dEerr = dEij - mi*(pworki + duRZij) - mj*(pworkj + duRZij);

      // const auto rpworki = -Pi*vri*safeInv(rhoi*ri, tiny);
      // const auto rpworkj = -Pj*vrj*safeInv(rhoj*rj, tiny);
      // const auto chi = std::max(0.0, std::min(1.0, Eerr*safeInv(mi*rpworki + mj*rpworkj, tiny)));
      // const auto duri = chi*rpworki;
      // const auto durj = chi*rpworkj;

      // const auto wi = std::abs(rpworki)*std::max(0.0, sgn(rpworki*Eerr)) + tiny;
      // const auto wj = std::abs(rpworkj)*std::max(0.0, sgn(rpworkj*Eerr)) + tiny;
      // VERIFY(wi > 0.0);
      // VERIFY(wj > 0.0);
      // const auto duri = wi/(wi + wj)*Eerr/mi;
      // const auto durj = wj/(wi + wj)*Eerr/mj;
      // // const auto duri = (Eerr - mi*(pworki + duRZij) - mj*(pworkj + duRZij))*wi*safeInv((wi + wj)*(mi + mj), tiny);
      // // const auto durj = duri*wj*safeInv(wi, tiny);

      // const auto duri = Eerr/(mi + mj);
      // const auto durj = Eerr/(mi + mj);
      // // const auto duij = (dEij - mi*(pworki + duRZij + rpworki) - mj*(pworkj + duRZij + rpworkj))/(mi + mj);

      // // Apply additive corrections
      // DepsDt_thread(nodeListi, i) += pworki + duij;
      // DepsDt_thread(nodeListj, j) += pworkj + duij;

      // Multiplicative correction
      const auto chi = dEij*safeInv(mi*pworki + mj*pworkj);
      DepsDt_thread(nodeListi, i) += chi*pworki;
      DepsDt_thread(nodeListj, j) += chi*pworkj;
    }

#pragma omp critical
    {
      DepsDt_thread.threadReduce();
      // dEtot += dEtot_thread;
    }
  }

  // // Correct for total energy conservation
  // auto dEsum0 = 0.0;
  // for (auto k = 0u; k < numFields; ++k) {
  //   const auto n = eps[k]->numInternalElements();
  //   for (auto i = 0u; i < n; ++i) {
  //     dEsum0 += mass(k,i) * DepsDt(k,i);
  //   }
  // }
  // dEsum0 = allReduce(dEsum0, SPHERAL_OP_SUM);
  // dEtot = allReduce(dEtot, SPHERAL_OP_SUM);
  // DepsDt *= dEtot*safeInv(dEsum0);
  // cerr << "BLAGO: " << dEtot*safeInv(dEsum0) << endl;

  // Now we can update the energy.
  for (auto nodeListi = 0u; nodeListi < numFields; ++nodeListi) {
    const auto n = eps[nodeListi]->numInternalElements();
#pragma omp parallel for
    for (auto i = 0u; i < n; ++i) {
      auto duii = 0.0;

      // // Add self-interaction contributions
      // const auto  mi = mass(nodeListi, i);
      // const auto  Pi = pressure(nodeListi, i);
      // const auto  Qi = effViscousPressure(nodeListi, i);
      // const auto  rhoi = rho(nodeListi, i);
      // const auto  ri = pos(nodeListi, i).y();
      // const auto& vi = velocity(nodeListi, i);
      // const auto& ai = acceleration(nodeListi, i);
      // const auto  vi12 = vi + ai*hdt;
      // const auto  vri12 = vi12.y();
      // const auto  ri12 = ri + vri12*hdt;
      // VERIFY2(ri > 0.0, "BLAGO: " << i << " " << ri << " " << vi);
      // duii = -(Pi + Qi)/rhoi*integrate_vr_over_r(vi.y(), ri, ai.y(), multiplier)*safeInv(multiplier);
      // duii = -(Pi + Qi)/rhoi*vri12*safeInv(ri12);
      // // duii *= std::max(0.0, std::min(1.0, Eerror(nodeListi, i)*safeInv(mi*duii)));

      // Add the self-contribution if any (RZ with strength does this for instance).
      if (selfInteraction) {
        const auto& vi = velocity(nodeListi, i);
        const auto& ai = acceleration(nodeListi, i);
        const auto  vi12 = vi + ai*hdt;
        duii -= vi12.dot(selfAccelerations(nodeListi, i));
        // duii -= 2.0*vi12.dot(selfAccelerations(nodeListi, i));
      }

      DepsDt(nodeListi, i) += duii;
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

