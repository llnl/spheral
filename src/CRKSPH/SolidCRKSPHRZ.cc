//---------------------------------Spheral++----------------------------------//
// SolidCRKSPH -- The CRKSPH/ACRKSPH solid material hydrodynamic
// package for Spheral++.
//
// This is the area-weighted RZ specialization.
//
// Created by JMO, Fri May 13 10:50:36 PDT 2016
//----------------------------------------------------------------------------//
#include "CRKSPH/SolidCRKSPHRZ.hh"
#include "FileIO/FileIO.hh"
#include "RK/ReproducingKernel.hh"
#include "RK/RKFieldNames.hh"
#include "CRKSPH/CRKSPH.hh"
#include "CRKSPH/computeCRKSPHSumMassDensity.hh"
#include "Physics/GenericHydro.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Hydro/RZNonSymmetricSpecificThermalEnergyPolicy.hh"
#include "Hydro/SpecificFromTotalThermalEnergyPolicy.hh"
#include "Strength/SolidFieldNames.hh"
#include "NodeList/SolidNodeList.hh"
#include "Strength/BulkModulusPolicy.hh"
#include "Strength/ShearModulusPolicy.hh"
#include "Strength/YieldStrengthPolicy.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/IncrementBoundedState.hh"
#include "DataBase/ReplaceState.hh"
#include "DataBase/ReplaceBoundedState.hh"
#include "RK/ContinuityVolumePolicyRZ.hh"
#include "ArtificialViscosity/ArtificialViscosity.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "Boundary/Boundary.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Neighbor/PairwiseField.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/range.hh"
#include "Utilities/NodeCoupling.hh"
#include "SolidMaterial/SolidEquationOfState.hh"
#include "Geometry/GeometryRegistrar.hh"

#include <algorithm>
#include <fstream>
#include <vector>
using std::vector;
using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {

//------------------------------------------------------------------------------
// Compute the artificial tensile stress correction tensor for the given 
// stress tensor
//------------------------------------------------------------------------------
inline
Dim<1>::SymTensor
tensileStressCorrection(const Dim<1>::SymTensor& sigma) {
  if (sigma.xx() > 0.0) {
    return -sigma;
  } else {
    return Dim<1>::SymTensor();
  }
}

inline
Dim<2>::SymTensor
tensileStressCorrection(const Dim<2>::SymTensor& sigma) {
  const EigenStruct<2> eigen = sigma.eigenVectors();
  const double lambdax = eigen.eigenValues.x();
  const double lambday = eigen.eigenValues.y();
  Dim<2>::SymTensor result((lambdax > 0.0 ? -lambdax : 0.0), 0.0,
                           0.0, (lambday > 0.0 ? -lambday : 0.0));
  result.rotationalTransform(eigen.eigenVectors);
  return result;
}

inline
Dim<3>::SymTensor
tensileStressCorrection(const Dim<3>::SymTensor& sigma) {
  const EigenStruct<3> eigen = sigma.eigenVectors();
  const double lambdax = eigen.eigenValues.x();
  const double lambday = eigen.eigenValues.y();
  const double lambdaz = eigen.eigenValues.z();
  Dim<3>::SymTensor result((lambdax > 0.0 ? -lambdax : 0.0), 0.0, 0.0,
                           0.0, (lambday > 0.0 ? -lambday : 0.0), 0.0,
                           0.0, 0.0, (lambdaz > 0.0 ? -lambdaz : 0.0));
  result.rotationalTransform(eigen.eigenVectors);
  return result;
}

//------------------------------------------------------------------------------
// Construct with the given artificial viscosity and kernels.
//------------------------------------------------------------------------------
SolidCRKSPHRZ::
SolidCRKSPHRZ(DataBase<Dimension>& dataBase,
              ArtificialViscosityHandle<Dimension>& Q,
              const RKOrder order,
              const double cfl,
              const bool useVelocityMagnitudeForDt,
              const bool compatibleEnergyEvolution,
              const bool evolveTotalEnergy,
              const bool XSPH,
              const MassDensityType densityUpdate,
              const double epsTensile,
              const double nTensile,
              const bool damageRelieveRubble):
  SolidCRKSPH<Dimension>(dataBase,
                         Q,
                         order,
                         cfl,
                         useVelocityMagnitudeForDt,
                         compatibleEnergyEvolution,
                         evolveTotalEnergy,
                         XSPH,
                         densityUpdate,
                         epsTensile,
                         nTensile,
                         damageRelieveRubble),
  mPairAccelerationsPtr(),
  mSelfAccelerations(FieldStorageType::CopyFields) {
}

//------------------------------------------------------------------------------
// On problem start up, we need to initialize our internal data.
//------------------------------------------------------------------------------
void
SolidCRKSPHRZ::
initializeProblemStartupDependencies(DataBase<Dimension>& dataBase,
                                     State<Dimension>& state,
                                     StateDerivatives<Dimension>& derivs) {

  // Correct the mass to mass/r.
  auto mass = dataBase.fluidMass();
  const auto pos = dataBase.fluidPosition();
  const unsigned numNodeLists = mass.numFields();
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned n = mass[nodeListi]->numElements();
    for (unsigned i = 0; i != n; ++i) {
      const Scalar circi = 2.0*M_PI*abs(pos(nodeListi, i).y());
      mass(nodeListi, i) /= circi;
    }
  }

  // Call the ancestor.
  SolidCRKSPH<Dimension>::initializeProblemStartupDependencies(dataBase, state, derivs);

  // Convert back to mass.
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned n = mass[nodeListi]->numElements();
    for (unsigned i = 0; i != n; ++i) {
      const Scalar circi = 2.0*M_PI*abs(pos(nodeListi, i).y());
      mass(nodeListi, i) *= circi;
    }
  }
}

//------------------------------------------------------------------------------
// Register the state we need/are going to evolve.
//------------------------------------------------------------------------------
void
SolidCRKSPHRZ::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {

  // Call the ancestor.
  SolidCRKSPH<Dimension>::registerState(dataBase, state);

  // // Reregister the deviatoric stress and plastic strain policies to the RZ specialized versions
  // // that account for the theta-theta component of the stress.
  // auto ps = state.fields(SolidFieldNames::plasticStrain, 0.0);
  // PolicyPointer plasticStrainPolicy(new RZPlasticStrainPolicy());
  // state.enroll(ps, plasticStrainPolicy);

  // Reregister the volume update
  auto vol = state.fields(HydroFieldNames::volume, 0.0);
  state.enroll(vol, make_policy<ContinuityVolumePolicyRZ>());

  // We have to choose either compatible or total energy evolution.
  const auto compatibleEnergy = this->compatibleEnergyEvolution();
  const auto evolveTotalEnergy = this->evolveTotalEnergy();
  VERIFY2(not (compatibleEnergy and evolveTotalEnergy),
          "SPH error : you cannot simultaneously use both compatibleEnergyEvolution and evolveTotalEnergy");

  // Register the specific thermal energy.
  // Note in RZ we require the specific thermal energy go before the position so we can use the r position
  // during update.  This is why we make position update dependent on the thermal energy in SPHBase.
  auto specificThermalEnergy = dataBase.fluidSpecificThermalEnergy();
  if (compatibleEnergy) {
    state.enroll(specificThermalEnergy, make_policy<RZNonSymmetricSpecificThermalEnergyPolicy>(dataBase));

  } else if (evolveTotalEnergy) {
    // If we're doing total energy, we register the specific energy to advance with the
    // total energy policy.
    state.enroll(specificThermalEnergy, make_policy<SpecificFromTotalThermalEnergyPolicy<Dimension>>());

  } else {
    // Otherwise we're just time-evolving the specific energy.
    state.enroll(specificThermalEnergy, make_policy<IncrementState<Dimension, Scalar>>());
  }
}

//------------------------------------------------------------------------------
// Register the state derivative fields.
//------------------------------------------------------------------------------
void
SolidCRKSPHRZ::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {
  SolidCRKSPH<Dimension>::registerDerivatives(dataBase, derivs);
  const auto compatibleEnergy = this->compatibleEnergyEvolution();
  if (compatibleEnergy) {
    const auto& connectivityMap = dataBase.connectivityMap();
    mPairAccelerationsPtr = std::make_unique<PairAccelerationsType>(connectivityMap);
    dataBase.resizeFluidFieldList(mSelfAccelerations, Vector::zero(), HydroFieldNames::selfAccelerations, false);
    derivs.enroll(HydroFieldNames::pairAccelerations, *mPairAccelerationsPtr);
    derivs.enroll(HydroFieldNames::selfAccelerations, mSelfAccelerations);
  }
}

//------------------------------------------------------------------------------
// This method is called once at the beginning of a timestep, after all state registration.
//------------------------------------------------------------------------------
void
SolidCRKSPHRZ::
preStepInitialize(const DataBase<Dimension>& dataBase, 
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {

  // Convert the mass to mass per unit length first.
  auto mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto pos = state.fields(HydroFieldNames::position, Vector::zero());
  const unsigned numNodeLists = mass.numFields();
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned n = mass[nodeListi]->numElements();
    for (unsigned i = 0; i != n; ++i) {
      const auto circi = 2.0*M_PI*abs(pos(nodeListi, i).y());
#ifdef WIN32
      if (circi > 0.0) mass(nodeListi, i) /= circi;
#else
      mass(nodeListi, i) /= circi;
#endif
    }
  }

  // Base class finalization does most of the work.
  CRKSPH<Dimension>::preStepInitialize(dataBase, state, derivs);

  // Now convert back to true masses.
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned n = mass[nodeListi]->numElements();
    for (unsigned i = 0; i != n; ++i) {
      const auto& xi = pos(nodeListi, i);
      const auto circi = 2.0*M_PI*abs(xi.y());
#ifdef WIN32
      if (circi > 0.0) mass(nodeListi, i) *= circi;
#else
      mass(nodeListi, i) *= circi;
#endif
    }
  }
}

//------------------------------------------------------------------------------
// Determine the principle derivatives.
//------------------------------------------------------------------------------
void
SolidCRKSPHRZ::
evaluateDerivatives(const Dimension::Scalar time,
                    const Dimension::Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivatives) const {

  // Depending on the type of the ArtificialViscosity, dispatch the call to
  // the secondDerivativesLoop
  auto& Qhandle = this->artificialViscosity();
  if (Qhandle.QPiTypeIndex() == std::type_index(typeid(Scalar))) {
      const auto& Q = dynamic_cast<const ArtificialViscosity<Dimension, Scalar>&>(Qhandle);
      this->evaluateDerivativesImpl(time, dt, dataBase, state, derivatives, Q);
  } else {
    CHECK(Qhandle.QPiTypeIndex() == std::type_index(typeid(Tensor)));
    const auto& Q = dynamic_cast<const ArtificialViscosity<Dimension, Tensor>&>(Qhandle);
    this->evaluateDerivativesImpl(time, dt, dataBase, state, derivatives, Q);
  }
}
  
//------------------------------------------------------------------------------
// Determine the principle derivatives.
//------------------------------------------------------------------------------
template<typename QType>
void
SolidCRKSPHRZ::
evaluateDerivativesImpl(const Dimension::Scalar /*time*/,
                        const Dimension::Scalar dt,
                        const DataBase<Dimension>& dataBase,
                        const State<Dimension>& state,
                        StateDerivatives<Dimension>& derivs,
                        const QType& Q) const {

  using QPiType = typename QType::ReturnType;

  // The kernels and such.
  const auto  order = this->correctionOrder();
  const auto& WR = state.template get<ReproducingKernel<Dimension>>(RKFieldNames::reproducingKernel(order));

  // A few useful constants we'll use in the following loop.
  //const double tiny = 1.0e-30;
  const auto  compatibleEnergy = this->compatibleEnergyEvolution();
  const auto  XSPH = this->XSPH();
  const auto  damageRelieveRubble = this->damageRelieveRubble();

  // The connectivity.
  const auto& connectivityMap = dataBase.connectivityMap();
  const auto& nodeLists = connectivityMap.nodeLists();
  const auto  numNodeLists = nodeLists.size();
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();

  // Get the state and derivative FieldLists.
  // State FieldLists.
  const auto mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto volume = state.fields(HydroFieldNames::volume, 0.0);
  const auto position = state.fields(HydroFieldNames::position, Vector::zero());
  const auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero());
  const auto massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  const auto H = state.fields(HydroFieldNames::H, SymTensor::zero());
  const auto pressure = state.fields(HydroFieldNames::pressure, 0.0);
  const auto soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  const auto S = state.fields(SolidFieldNames::deviatoricStress, SymTensor::zero());
  const auto mu = state.fields(SolidFieldNames::shearModulus, 0.0);
  const auto damage = state.fields(SolidFieldNames::tensorDamage, SymTensor::zero());
  const auto fragIDs = state.fields(SolidFieldNames::fragmentIDs, int(1));
  const auto pTypes = state.fields(SolidFieldNames::particleTypes, int(0));
  const auto corrections = state.fields(RKFieldNames::rkCorrections(order), RKCoefficients<Dimension>());
  const auto surfacePoint = state.fields(HydroFieldNames::surfacePoint, 0);
  const auto fClQ = state.fields(HydroFieldNames::ArtificialViscousClMultiplier, 0.0, true);
  const auto fCqQ = state.fields(HydroFieldNames::ArtificialViscousCqMultiplier, 0.0, true);
  const auto DvDxQ = state.fields(HydroFieldNames::ArtificialViscosityVelocityGradient, Tensor::zero(), true);
  CHECK(mass.size() == numNodeLists);
  CHECK(position.size() == numNodeLists);
  CHECK(velocity.size() == numNodeLists);
  CHECK(massDensity.size() == numNodeLists);
  CHECK(specificThermalEnergy.size() == numNodeLists);
  CHECK(H.size() == numNodeLists);
  CHECK(pressure.size() == numNodeLists);
  CHECK(soundSpeed.size() == numNodeLists);
  CHECK(S.size() == numNodeLists);
  CHECK(mu.size() == numNodeLists);
  CHECK(damage.size() == numNodeLists);
  CHECK(fragIDs.size() == numNodeLists);
  CHECK(pTypes.size() == numNodeLists);
  CHECK(corrections.size() == numNodeLists);
  CHECK(surfacePoint.size() == numNodeLists);
  CHECK(fClQ.size() == 0 or fClQ.size() == numNodeLists);
  CHECK(fCqQ.size() == 0 or fCqQ.size() == numNodeLists);
  CHECK(DvDxQ.size() == 0 or DvDxQ.size() == numNodeLists);

  // Derivative FieldLists.
  auto  DxDt = derivs.fields(IncrementState<Dimension, Vector>::prefix() + HydroFieldNames::position, Vector::zero());
  auto  DrhoDt = derivs.fields(IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, 0.0);
  auto  DvDt = derivs.fields(HydroFieldNames::hydroAcceleration, Vector::zero());
  auto  DepsDt = derivs.fields(IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
  auto  DvDx = derivs.fields(HydroFieldNames::velocityGradient, Tensor::zero());
  auto  localDvDx = derivs.fields(HydroFieldNames::internalVelocityGradient, Tensor::zero());
  auto  maxViscousPressure = derivs.fields(HydroFieldNames::maxViscousPressure, 0.0);
  auto  effViscousPressure = derivs.fields(HydroFieldNames::effectiveViscousPressure, 0.0);
  auto* pairAccelerationsPtr = derivs.template getPtr<PairAccelerationsType>(HydroFieldNames::pairAccelerations);
  auto  selfAccelerations = derivs.fields(HydroFieldNames::selfAccelerations, Vector::zero(), true);
  auto  XSPHDeltaV = derivs.fields(HydroFieldNames::XSPHDeltaV, Vector::zero());
  auto  DSDt = derivs.fields(IncrementState<Dimension, SymTensor>::prefix() + SolidFieldNames::deviatoricStress, SymTensor::zero());
  CHECK(DxDt.size() == numNodeLists);
  CHECK(DrhoDt.size() == numNodeLists);
  CHECK(DvDt.size() == numNodeLists);
  CHECK(DepsDt.size() == numNodeLists);
  CHECK(DvDx.size() == numNodeLists);
  CHECK(localDvDx.size() == numNodeLists);
  CHECK(maxViscousPressure.size() == numNodeLists);
  CHECK(effViscousPressure.size() == numNodeLists);
  CHECK(XSPHDeltaV.size() == numNodeLists);
  CHECK(DSDt.size() == numNodeLists);
  CHECK((compatibleEnergy and pairAccelerationsPtr->size() == npairs) or not compatibleEnergy);
  CHECK((compatibleEnergy     and selfAccelerations.size() == 0u) or
        (not compatibleEnergy and selfAccelerations.size() == numNodeLists));

  // Build the functor we use to compute the effective coupling between nodes.
  const NodeCoupling coupling;

  // Walk all the interacting pairs.
#pragma omp parallel
  {
    // Thread private scratch variables
    int i, j, nodeListi, nodeListj;
    Scalar Wi, Wj, Qi, Qj;
    QPiType QPiij, QPiji;
    Vector gradWi, gradWj;
    Vector deltagrad, forceij, forceji;
    Vector xij, vij, etai, etaj;
    SymTensor sigmai, sigmaj;

    typename SpheralThreads<Dimension>::FieldListStack threadStack;
    auto DvDt_thread = DvDt.threadCopy(threadStack);
    auto DepsDt_thread = DepsDt.threadCopy(threadStack);
    auto DvDx_thread = DvDx.threadCopy(threadStack);
    auto localDvDx_thread = localDvDx.threadCopy(threadStack);
    auto maxViscousPressure_thread = maxViscousPressure.threadCopy(threadStack, ThreadReduction::MAX);
    auto effViscousPressure_thread = effViscousPressure.threadCopy(threadStack);
    auto XSPHDeltaV_thread = XSPHDeltaV.threadCopy(threadStack);

#pragma omp for
    for (auto kk = 0u; kk < npairs; ++kk) {
      i = pairs[kk].i_node;
      j = pairs[kk].j_node;
      nodeListi = pairs[kk].i_list;
      nodeListj = pairs[kk].j_list;

      // Get the state for node i.
      const auto& posi = position(nodeListi, i);
      const auto  ri = abs(posi.y());
      const auto  circi = 2.0*M_PI*ri;
      const auto  mi = mass(nodeListi, i);
      const auto  mRZi = mi/circi;
      const auto& vi = velocity(nodeListi, i);
      const auto  rhoi = massDensity(nodeListi, i);
      //const auto  epsi = specificThermalEnergy(nodeListi, i);
      const auto  Pi = pressure(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto  ci = soundSpeed(nodeListi, i);
      const auto  Si = S(nodeListi, i);
      const auto  pTypei = pTypes(nodeListi, i);
      const auto& correctionsi = corrections(nodeListi, i);
      const auto  weighti = volume(nodeListi, i);  // Change CRKSPH weights here if need be!
      CHECK(mi > 0.0);
      CHECK(rhoi > 0.0);
      CHECK(weighti > 0.0);

      //auto& DrhoDti = DrhoDt(nodeListi, i);
      auto& DvDti = DvDt_thread(nodeListi, i);
      auto& DepsDti = DepsDt_thread(nodeListi, i);
      auto& DvDxi = DvDx_thread(nodeListi, i);
      auto& localDvDxi = localDvDx_thread(nodeListi, i);
      auto& maxViscousPressurei = maxViscousPressure_thread(nodeListi, i);
      auto& effViscousPressurei = effViscousPressure_thread(nodeListi, i);
      auto& XSPHDeltaVi = XSPHDeltaV_thread(nodeListi, i);

      // Get the state for node j
      const auto& posj = position(nodeListj, j);
      const auto  rj = abs(posj.y());
      const auto  circj = 2.0*M_PI*rj;
      const auto  mj = mass(nodeListj, j);
      const auto  mRZj = mj/circj;
      const auto& vj = velocity(nodeListj, j);
      const auto  rhoj = massDensity(nodeListj, j);
      //const auto  epsj = specificThermalEnergy(nodeListj, j);
      const auto  Pj = pressure(nodeListj, j);
      const auto& Hj = H(nodeListj, j);
      const auto  cj = soundSpeed(nodeListj, j);
      const auto  pTypej = pTypes(nodeListj, j);
      const auto& correctionsj = corrections(nodeListj, j);
      const auto& Sj = S(nodeListj, j);
      const auto  weightj = volume(nodeListj, j);     // Change CRKSPH weights here if need be!
      CHECK(mj > 0.0);
      CHECK(rhoj > 0.0);
      CHECK(weightj > 0.0);

      auto& DvDtj = DvDt_thread(nodeListj, j);
      auto& DepsDtj = DepsDt_thread(nodeListj, j);
      auto& DvDxj = DvDx_thread(nodeListj, j);
      auto& localDvDxj = localDvDx_thread(nodeListj, j);
      auto& maxViscousPressurej = maxViscousPressure_thread(nodeListj, j);
      auto& effViscousPressurej = effViscousPressure_thread(nodeListj, j);
      auto& XSPHDeltaVj = XSPHDeltaV_thread(nodeListj, j);

      // Node displacement.
      xij = posi - posj;
      vij = vi - vj;
      etai = Hi*xij;
      etaj = Hj*xij;

      // Flag if this is a contiguous material pair or not.
      const auto sameMatij = true; // nodeListi == nodeListj; // (nodeListi == nodeListj and fragIDi == fragIDj);

      // Flag if at least one particle is free (0).
      const auto freeParticle = (pTypei == 0 or pTypej == 0);

      // Symmetrized kernel weight and gradient.
      std::tie(Wj, gradWj) = WR.evaluateKernelAndGradient( xij, Hj, correctionsi);  // Hj because we compute RK using scatter formalism
      std::tie(Wi, gradWi) = WR.evaluateKernelAndGradient(-xij, Hi, correctionsj);
      deltagrad = gradWj - gradWi;

      // Find the damaged pair weighting scaling.
      const auto fDij = coupling(pairs[kk]);
      CHECK(fDij >= 0.0 and fDij <= 1.0);

      // Compute the artificial viscous pressure (Pi = P/rho^2 actually).
      Q.QPiij(QPiij, QPiji, Qi, Qj,
              nodeListi, i, nodeListj, j,
              posi, Hi, etai, vi, rhoi, ci,  
              posj, Hj, etaj, vj, rhoj, cj,
              fClQ, fCqQ, DvDxQ); 
      const auto Qaccij = (rhoi*rhoi*QPiij + rhoj*rhoj*QPiji)*deltagrad;
      const auto workQi = (rhoj*rhoj*QPiji*vij).dot(deltagrad);          // CRK
      const auto workQj = (rhoi*rhoi*QPiij*vij).dot(deltagrad);          // CRK
      maxViscousPressurei = max(maxViscousPressurei, 4.0*Qi);               // We need tighter timestep controls on the Q with CRK
      maxViscousPressurej = max(maxViscousPressurej, 4.0*Qj);
      effViscousPressurei += weightj * Qi * Wj;
      effViscousPressurej += weighti * Qj * Wi;

      // Velocity gradient.
      DvDxi -= fDij * weightj*vij.dyad(gradWj);
      DvDxj += fDij * weighti*vij.dyad(gradWi);
      if (sameMatij) {
        localDvDxi -= fDij * weightj*vij.dyad(gradWj);
        localDvDxj += fDij * weighti*vij.dyad(gradWi);
      }

      // Compute the stress tensors.
      if (sameMatij) {
        sigmai = fDij*Si - Pi * SymTensor::one();
        sigmaj = fDij*Sj - Pj * SymTensor::one();
      } else {
        sigmai = -Pi * SymTensor::one();
        sigmaj = -Pj * SymTensor::one();
      }

      // We decide between RK and CRK for the momentum and energy equations based on the surface condition.
      // Momentum
      forceij = (true ? // surfacePoint(nodeListi, i) <= 1 ? 
                 0.5*weighti*weightj*(-(sigmai + sigmaj)*deltagrad + Qaccij) :                                               // Type III CRK interpoint force.
                 mi*weightj*((-(sigmaj - sigmai)*gradWj)/rhoi + rhoi*QPiij*gradWj));                                         // RK
      forceji = (true ? // surfacePoint(nodeListj, j) <= 1 ?
                 0.5*weighti*weightj*(-(sigmai + sigmaj)*deltagrad + Qaccij) :                                               // Type III CRK interpoint force.
                 mj*weighti*((-(sigmaj - sigmai)*gradWi)/rhoj - rhoj*QPiji*gradWi));                                         // RK
      if (freeParticle) {
        DvDti -= forceij/mRZi;
        DvDtj += forceji/mRZj;
      }
      if (compatibleEnergy) {
        (*pairAccelerationsPtr)[kk][0] = -forceij/mRZi;
        (*pairAccelerationsPtr)[kk][1] =  forceji/mRZj;
      }

      // Energy
      DepsDti += (true ? // surfacePoint(nodeListi, i) <= 1 ?
                  0.5*weighti*weightj*(-sigmaj.dot(vij).dot(deltagrad) + workQi)/mRZi :                                      // CRK
                  (weightj*rhoi*QPiij*vij).dot(gradWj));                                                                     // RK, Q term only -- adiabatic portion added later
      DepsDtj += (true ? // surfacePoint(nodeListj, j) <= 1 ?
                  0.5*weighti*weightj*(-sigmai.dot(vij).dot(deltagrad) + workQj)/mRZj :                                      // CRK
                  (-weighti*rhoj*QPiji*vij).dot(gradWi));                                                                    // RK, Q term only -- adiabatic portion added later

      // Estimate of delta v (for XSPH).
      XSPHDeltaVi -= fDij*weightj*Wj*vij;
      XSPHDeltaVj += fDij*weighti*Wi*vij;
    }

    // Reduce the thread values to the master.
    threadReduceFieldLists<Dimension>(threadStack);

  }   // OMP parallel

  // Finish up the derivatives for each point.
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto& nodeList = mass[nodeListi]->nodeList();
    const auto ni = nodeList.numInternalNodes();

    // Check if we can identify a reference density.
    auto rho0 = 0.0;
    try {
      rho0 = dynamic_cast<const SolidEquationOfState<Dimension>&>(dynamic_cast<const SolidNodeList<Dimension>&>(nodeList).equationOfState()).referenceDensity();
      // cerr << "Setting reference density to " << rho0 << endl;
    } catch(...) {
      // cerr << "BLAGO!" << endl;
    }

#pragma omp parallel for
    for (auto i = 0u; i < ni; ++i) {

      // Get the state for node i.
      const auto& posi = position(nodeListi, i);
      const auto  ri = abs(posi.y());
      //const auto  circi = 2.0*M_PI*ri;
      const auto  mi = mass(nodeListi, i);
      //const auto  mRZi = mi/circi;
      const auto& vi = velocity(nodeListi, i);
      const auto  rhoi = massDensity(nodeListi, i);
      //const auto  epsi = specificThermalEnergy(nodeListi, i);
      const auto  Pi = pressure(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      //const auto  ci = soundSpeed(nodeListi, i);
      const auto  Si = S(nodeListi, i);
      const auto  STTi = -Si.Trace();
      const auto  mui = mu(nodeListi, i);
      const auto  zetai = abs((Hi*posi).y());
      const auto  hri = ri*safeInv(zetai);
      const auto  riInv = safeInv(ri, 0.25*hri);

      auto& DxDti = DxDt(nodeListi, i);
      auto& DrhoDti = DrhoDt(nodeListi, i);
      auto& DvDti = DvDt(nodeListi, i);
      auto& DepsDti = DepsDt(nodeListi, i);
      // auto& DvDxi = DvDx(nodeListi, i);
      auto& localDvDxi = localDvDx(nodeListi, i);
      auto& XSPHDeltaVi = XSPHDeltaV(nodeListi, i);
      auto& DSDti = DSDt(nodeListi, i);

      // Determine the position evolution, based on whether we're doing XSPH or not.
      if (XSPH) {
        DxDti = vi + XSPHDeltaVi;
      } else {
        DxDti = vi;
      }

      // Finish the acceleration.
      const Vector deltaDvDti(Si(1,0)/rhoi*riInv,
                              (Si(1,1) - STTi)/rhoi*riInv);
      DvDti += deltaDvDti;
      if (compatibleEnergy) selfAccelerations(nodeListi, i) = deltaDvDti;

      // Time evolution of the mass density.
      const auto vri = vi.y(); // + XSPHDeltaVi.y();
      DrhoDti = -rhoi*(localDvDxi.Trace() + vri*riInv);

      // For a surface point, add the RK thermal energy evolution.
      // if (surfacePoint(nodeListi, i) > 1) DepsDti += (Si - Pi*SymTensor::one()).doubledot(DvDxi)/rhoi;

      // Finish the specific thermal energy evolution.
      DepsDti += (STTi - Pi)/rhoi*vri*riInv;

      // If needed finish the total energy derivative.
      if (this->evolveTotalEnergy()) DepsDti = mi*(vi.dot(DvDti) + DepsDti);

      // Optionally use damage to ramp down stress on damaged material.
      const auto Di = (damageRelieveRubble ? 
                       max(0.0, min(1.0, damage(nodeListi, i).Trace() - 1.0)) :
                       0.0);
      // Hideali = (1.0 - Di)*Hideali + Di*Hfield0(nodeListi, i);
      // DHDti = (1.0 - Di)*DHDti + Di*(Hfield0(nodeListi, i) - Hi)*0.25/dt;

      // We also adjust the density evolution in the presence of damage.
      if (rho0 > 0.0) DrhoDti = (1.0 - Di)*DrhoDti - 0.25/dt*Di*(rhoi - rho0);

      // Determine the deviatoric stress evolution.
      const auto deformation = localDvDxi.Symmetric();
      const auto deformationTT = vi.y()*riInv;
      const auto spin = localDvDxi.SkewSymmetric();
      const auto deviatoricDeformation = deformation - ((deformation.Trace() + deformationTT)/3.0)*SymTensor::one();
      const auto spinCorrection = (spin*Si + Si*spin).Symmetric();
      DSDti = spinCorrection + 2.0*mui*deviatoricDeformation;

      // In the presence of damage, add a term to reduce the stress on this point.
      DSDti = (1.0 - Di)*DSDti - Di*Si*0.25/dt;
    }
  }
}

//------------------------------------------------------------------------------
// Apply the ghost boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
void
SolidCRKSPHRZ::
applyGhostBoundaries(State<Dimension>& state,
                     StateDerivatives<Dimension>& derivs) {

  // Convert the mass to mass/length before BCs are applied.
  FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  const FieldList<Dimension, Vector> pos = state.fields(HydroFieldNames::position, Vector::zero());
  const unsigned numNodeLists = mass.numFields();
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned n = mass[nodeListi]->numElements();
    for (unsigned i = 0; i != n; ++i) {
      const Scalar circi = 2.0*M_PI*abs(pos(nodeListi, i).y());
      CHECK(circi > 0.0);
#ifdef WIN32
      if (circi > 0.0) mass(nodeListi, i) /= circi;
#else
      mass(nodeListi, i) /= circi;
#endif
    }
  }

  // Apply ordinary BCs.
  SolidCRKSPH<Dimension>::applyGhostBoundaries(state, derivs);
  for (auto boundaryPtr: range(this->boundaryBegin(), this->boundaryEnd())) boundaryPtr->finalizeGhostBoundary();

  // Scale back to mass.
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned n = mass[nodeListi]->numElements();
    for (unsigned i = 0; i != n; ++i) {
      const Scalar circi = 2.0*M_PI*abs(pos(nodeListi, i).y());
      CHECK(circi > 0.0);
#ifdef WIN32
      if (circi > 0.0) mass(nodeListi, i) *= circi;
#else
      mass(nodeListi, i) *= circi;
#endif
    }
  }
}

//------------------------------------------------------------------------------
// Enforce the boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
void
SolidCRKSPHRZ::
enforceBoundaries(State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {

  // Convert the mass to mass/length before BCs are applied.
  FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  FieldList<Dimension, Vector> pos = state.fields(HydroFieldNames::position, Vector::zero());
  const unsigned numNodeLists = mass.numFields();
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned n = mass[nodeListi]->numInternalElements();
    for (unsigned i = 0; i != n; ++i) {
      const Scalar circi = 2.0*M_PI*abs(pos(nodeListi, i).y());
      CHECK(circi > 0.0);
#ifdef WIN32
      if (circi > 0.0) mass(nodeListi, i) /= circi;
#else
      mass(nodeListi, i) /= circi;
#endif
    }
  }

  // Apply ordinary BCs.
  SolidCRKSPH<Dimension>::enforceBoundaries(state, derivs);

  // Scale back to mass.
  // We also ensure no point approaches the z-axis too closely.
  FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero());
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned n = mass[nodeListi]->numInternalElements();
    //const Scalar nPerh = mass[nodeListi]->nodeList().nodesPerSmoothingScale();
    for (unsigned i = 0; i != n; ++i) {
      Vector& posi = pos(nodeListi, i);
      const Scalar circi = 2.0*M_PI*abs(posi.y());
#ifdef WIN32
      if (circi > 0.0) mass(nodeListi, i) *= circi;
#else
      mass(nodeListi, i) *= circi;
#endif
    }
  }
}

}
