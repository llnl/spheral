//---------------------------------Spheral++----------------------------------//
// SPH -- The classic SPH/ASPH hydrodynamic packages for Spheral++.
// 
// Created by JMO, Thu Nov 21 16:36:40 PST 2024
//----------------------------------------------------------------------------//
#include "config.hh"
#include "SPH/SPH.hh"
#include "FileIO/FileIO.hh"
#include "DataBase/State.hh"
#include "Physics/Physics.hh"
#include "Physics/GenericHydro.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Hydro/SpecificThermalEnergyPolicy.hh"
#include "Hydro/SpecificFromTotalThermalEnergyPolicy.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/ReplaceState.hh"
#include "DataBase/IncrementBoundedState.hh"
#include "DataBase/ReplaceBoundedState.hh"
#include "ArtificialViscosity/ArtificialViscosity.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Neighbor/PairwiseField.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/range.hh"
#include "Utilities/Timer.hh"
#include "Utilities/timingUtilities.hh"

#include <algorithm>
#include <fstream>
#include <map>
#include <vector>
using std::vector;
using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {

//------------------------------------------------------------------------------
// Construct with the given artificial viscosity and kernels.
//------------------------------------------------------------------------------
template<typename Dimension>
SPH<Dimension>::
SPH(DataBase<Dimension>& dataBase,
    ArtificialViscosity<Dimension>& Q,
    const TableKernel<Dimension>& W,
    const TableKernel<Dimension>& WPi,
    const double cfl,
    const bool useVelocityMagnitudeForDt,
    const bool compatibleEnergyEvolution,
    const bool evolveTotalEnergy,
    const bool gradhCorrection,
    const bool XSPH,
    const bool correctVelocityGradient,
    const bool sumMassDensityOverAllNodeLists,
    const MassDensityType densityUpdate,
    const double epsTensile,
    const double nTensile,
    const Vector& xmin,
    const Vector& xmax):
  SPHBase<Dimension>(dataBase,
                     Q,
                     W,
                     WPi,
                     cfl,
                     useVelocityMagnitudeForDt,
                     compatibleEnergyEvolution,
                     evolveTotalEnergy,
                     gradhCorrection,
                     XSPH,
                     correctVelocityGradient,
                     sumMassDensityOverAllNodeLists,
                     densityUpdate,
                     epsTensile,
                     nTensile,
                     xmin,
                     xmax),
  mPairAccelerationsPtr() {
}

//------------------------------------------------------------------------------
// Register the state we need/are going to evolve.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SPH<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {
  TIME_BEGIN("SPHregister");

  SPHBase<Dimension>::registerState(dataBase, state);

  // We have to choose either compatible or total energy evolution.
  const auto compatibleEnergy = this->compatibleEnergyEvolution();
  const auto evolveTotalEnergy = this->evolveTotalEnergy();
  VERIFY2(not (compatibleEnergy and evolveTotalEnergy),
          "SPH error : you cannot simultaneously use both compatibleEnergyEvolution and evolveTotalEnergy");

  // Register the specific thermal energy.
  auto specificThermalEnergy = dataBase.fluidSpecificThermalEnergy();
  if (compatibleEnergy) {
    state.enroll(specificThermalEnergy, make_policy<SpecificThermalEnergyPolicy<Dimension>>(dataBase));

  } else if (evolveTotalEnergy) {
    // If we're doing total energy, we register the specific energy to advance with the
    // total energy policy.
    state.enroll(specificThermalEnergy, make_policy<SpecificFromTotalThermalEnergyPolicy<Dimension>>());

  } else {
    // Otherwise we're just time-evolving the specific energy.
    state.enroll(specificThermalEnergy, make_policy<IncrementState<Dimension, Scalar>>());
  }

  TIME_END("SPHregister");
}

//------------------------------------------------------------------------------
// Register the state derivative fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SPH<Dimension>::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {
  TIME_BEGIN("SPHregisterDerivs");
  SPHBase<Dimension>::registerDerivatives(dataBase, derivs);
  const auto compatibleEnergy = this->compatibleEnergyEvolution();
  if (compatibleEnergy) {
    const auto& connectivityMap = dataBase.connectivityMap();
    mPairAccelerationsPtr = std::make_unique<PairAccelerationsType>(connectivityMap);
    derivs.enroll(HydroFieldNames::pairAccelerations, *mPairAccelerationsPtr);
  }
  TIME_END("SPHregisterDerivs");
}

//------------------------------------------------------------------------------
// Determine the principle derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SPH<Dimension>::
evaluateDerivatives(const typename Dimension::Scalar time,
                    const typename Dimension::Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivatives) const {

  // Depending on the type of the ArtificialViscosityView, dispatch the call to
  // the secondDerivativesLoop
  auto& Qhandle = this->artificialViscosity();
  if (Qhandle.QPiTypeIndex() == std::type_index(typeid(Scalar))) {
    chai::managed_ptr<ArtificialViscosityView<Dimension, Scalar>> Q = Qhandle.getScalarView();
    this->evaluateDerivativesImpl(time, dt, dataBase, state, derivatives, Q);
  } else {
    CHECK(Qhandle.QPiTypeIndex() == std::type_index(typeid(Tensor)));
    chai::managed_ptr<ArtificialViscosityView<Dimension, Tensor>> Q = Qhandle.getTensorView();
    this->evaluateDerivativesImpl(time, dt, dataBase, state, derivatives, Q);
  }
}
  
//------------------------------------------------------------------------------
// evaluateDerivatives
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename QType>
void
SPH<Dimension>::
evaluateDerivativesImpl(const typename Dimension::Scalar time,
                        const typename Dimension::Scalar dt,
                        const DataBase<Dimension>& dataBase,
                        const State<Dimension>& state,
                        StateDerivatives<Dimension>& derivs,
                        chai::managed_ptr<QType> Q) const {
  TIME_BEGIN("SPHevalDerivs");
  TIME_BEGIN("SPHevalDerivs_initial");

  using QPiType = typename QType::ReturnType;

  //static double totalLoopTime = 0.0;

  // The kernels and such.
  auto& W = this->kernel();
  auto& WQ = this->PiKernel();
  auto W_view = W.view();
  auto WQ_view = WQ.view();
  const auto  oneKernel = (W == WQ);

  // A few useful constants we'll use in the following loop.
  const auto tiny = 1.0e-30;
  const auto W0 = W(0.0, 1.0);
  const auto epsTensile = this->epsilonTensile();
  //const auto compatibleEnergy = this->compatibleEnergyEvolution();
  const auto evolveTotalEnergy = this->evolveTotalEnergy();
  const auto XSPH = this->XSPH();

  // The connectivity.
  auto& connectivityMap = dataBase.connectivityMap();
  const auto& nodeLists = connectivityMap.nodeLists();
  const auto numNodeLists = nodeLists.size();

  // The set of interacting node pairs.
  const auto& pairs_v = connectivityMap.nodePairList();
  const auto  pairs = pairs_v.view();
  const auto  npairs = pairs.size();

  // Get the state and derivative FieldLists.
  // State FieldLists.
  auto mass_v = state.fields(HydroFieldNames::mass, 0.0);
  auto position_v = state.fields(HydroFieldNames::position, Vector::zero());
  auto velocity_v = state.fields(HydroFieldNames::velocity, Vector::zero());
  auto massDensity_v = state.fields(HydroFieldNames::massDensity, 0.0);
  auto H_v = state.fields(HydroFieldNames::H, SymTensor::zero());
  auto pressure_v = state.fields(HydroFieldNames::pressure, 0.0);
  auto soundSpeed_v = state.fields(HydroFieldNames::soundSpeed, 0.0);
  auto omega_v = state.fields(HydroFieldNames::omegaGradh, 0.0);
  auto mass = mass_v.view();
  auto position = position_v.view();
  auto velocity = velocity_v.view();
  auto massDensity = massDensity_v.view();
  auto H = H_v.view();
  auto pressure = pressure_v.view();
  auto soundSpeed = soundSpeed_v.view();
  auto omega = omega_v.view();
  auto fClQ = state.fields(HydroFieldNames::ArtificialViscousClMultiplier, 0.0, true);
  auto fCqQ = state.fields(HydroFieldNames::ArtificialViscousCqMultiplier, 0.0, true);
  auto DvDxQ = state.fields(HydroFieldNames::ArtificialViscosityVelocityGradient, Tensor::zero(), true);
  auto DvDxQView = DvDxQ.view();
  auto fClQView = fClQ.view();
  auto fCqQView = fCqQ.view();
  CHECK(mass_v.size() == numNodeLists);
  CHECK(position_v.size() == numNodeLists);
  CHECK(velocity_v.size() == numNodeLists);
  CHECK(massDensity_v.size() == numNodeLists);
  CHECK(H_v.size() == numNodeLists);
  CHECK(pressure_v.size() == numNodeLists);
  CHECK(soundSpeed_v.size() == numNodeLists);
  CHECK(omega_v.size() == numNodeLists);
  CHECK(fClQ.size() == 0 or fClQ.size() == numNodeLists);
  CHECK(fCqQ.size() == 0 or fCqQ.size() == numNodeLists);
  CHECK(DvDxQ.size() == 0 or DvDxQ.size() == numNodeLists);

  // Derivative FieldLists.
  auto  rhoSum_v = derivs.fields(ReplaceState<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, 0.0);
  auto  normalization_v = derivs.fields(HydroFieldNames::normalization, 0.0);
  auto  DxDt_v = derivs.fields(IncrementState<Dimension, Vector>::prefix() + HydroFieldNames::position, Vector::zero());
  auto  DrhoDt_v = derivs.fields(IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, 0.0);
  auto  DvDt_v = derivs.fields(HydroFieldNames::hydroAcceleration, Vector::zero());
  auto  DepsDt_v = derivs.fields(IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
  auto  DvDx_v = derivs.fields(HydroFieldNames::velocityGradient, Tensor::zero());
  auto  localDvDx_v = derivs.fields(HydroFieldNames::internalVelocityGradient, Tensor::zero());
  auto  gradRho_v = derivs.fields(HydroFieldNames::massDensityGradient, Vector::zero());
  auto  M_v = derivs.fields(HydroFieldNames::M_SPHCorrection, Tensor::zero());
  auto  localM_v = derivs.fields("local " + HydroFieldNames::M_SPHCorrection, Tensor::zero());
  auto  maxViscousPressure_v = derivs.fields(HydroFieldNames::maxViscousPressure, 0.0);
  auto  effViscousPressure_v = derivs.fields(HydroFieldNames::effectiveViscousPressure, 0.0);
  //auto* pairAccelerationsPtr_v = derivs.template getPtr<PairAccelerationsType>(HydroFieldNames::pairAccelerations);
  auto  XSPHWeightSum_v = derivs.fields(HydroFieldNames::XSPHWeightSum, 0.0);
  auto  XSPHDeltaV_v = derivs.fields(HydroFieldNames::XSPHDeltaV, Vector::zero());
  auto rhoSum = rhoSum_v.view();
  auto normalization = normalization_v.view();
  auto DxDt = DxDt_v.view();
  auto DrhoDt = DrhoDt_v.view();
  auto DvDt = DvDt_v.view();
  auto DepsDt = DepsDt_v.view();
  auto DvDx = DvDx_v.view();
  auto localDvDx = localDvDx_v.view();
  auto gradRho = gradRho_v.view();
  auto M = M_v.view();
  auto localM = localM_v.view();
  auto maxViscousPressure = maxViscousPressure_v.view();
  auto effViscousPressure = effViscousPressure_v.view();
  auto XSPHWeightSum = XSPHWeightSum_v.view();
  auto XSPHDeltaV = XSPHDeltaV_v.view();
  CHECK(rhoSum_v.size() == numNodeLists);
  CHECK(normalization_v.size() == numNodeLists);
  CHECK(DxDt_v.size() == numNodeLists);
  CHECK(DrhoDt_v.size() == numNodeLists);
  CHECK(DvDt_v.size() == numNodeLists);
  CHECK(DepsDt_v.size() == numNodeLists);
  CHECK(DvDx_v.size() == numNodeLists);
  CHECK(localDvDx_v.size() == numNodeLists);
  CHECK(gradRho_v.size() == numNodeLists);
  CHECK(M_v.size() == numNodeLists);
  CHECK(localM_v.size() == numNodeLists);
  CHECK(maxViscousPressure_v.size() == numNodeLists);
  CHECK(effViscousPressure_v.size() == numNodeLists);
  CHECK(XSPHWeightSum_v.size() == numNodeLists);
  CHECK(XSPHDeltaV_v.size() == numNodeLists);
  //CHECK((compatibleEnergy and pairAccelerationsPtr->size() == npairs) or not compatibleEnergy);

  // The scale for the tensile correction.
  const auto& nodeList = mass_v[0]->nodeList();
  const auto  nPerh = nodeList.nodesPerSmoothingScale();
  const auto  WnPerh = W(1.0/nPerh, 1.0);
  bool CorrectVelocityGradient = this->mCorrectVelocityGradient;
  TIME_END("SPHevalDerivs_initial");

  // Walk all the interacting pairs.
  TIME_BEGIN("SPHevalDerivs_pairs");
  //RAJA::region<RAJA::seq_region>([=]()
  //#pragma omp parallel
  {
    // Thread private scratch variables
    // unsigned i, j, nodeListi, nodeListj;
    // Vector gradWi, gradWj, gradWQi, gradWQj;
    // Scalar Wi, gWi, WQi, gWQi, Wj, gWj, WQj, gWQj, Qi, Qj;
    // QPiType QPiij, QPiji;

    // typename SpheralThreads<Dimension>::FieldListStack threadStack;
    // auto rhoSum_thread = rhoSum.threadCopy(threadStack);
    // auto normalization_thread = normalization.threadCopy(threadStack);
    // auto DvDt_thread = DvDt.threadCopy(threadStack);
    // auto DepsDt_thread = DepsDt.threadCopy(threadStack);
    // auto DvDx_thread = DvDx.threadCopy(threadStack);
    // auto localDvDx_thread = localDvDx.threadCopy(threadStack);
    // auto gradRho_thread = gradRho.threadCopy(threadStack);
    // auto M_thread = M.threadCopy(threadStack);
    // auto localM_thread = localM.threadCopy(threadStack);
    // auto maxViscousPressure_thread = maxViscousPressure.threadCopy(threadStack, ThreadReduction::MAX);
    // auto effViscousPressure_thread = effViscousPressure.threadCopy(threadStack);
    // auto XSPHWeightSum_thread = XSPHWeightSum.threadCopy(threadStack);
    // auto XSPHDeltaV_thread = XSPHDeltaV.threadCopy(threadStack);

// #pragma omp for
//     for (auto kk = 0u; kk < npairs; ++kk) {
    RAJA::forall<EXEC_POLICY>(TRS_UINT(0u, npairs),
    [=] SPHERAL_HOST_DEVICE (size_t kk) {
      Vector gradWi, gradWj, gradWQi, gradWQj;
      Scalar Wi, gWi, WQi, gWQi, Wj, gWj, WQj, gWQj, Qi, Qj;
      QPiType QPiij, QPiji;
      size_t i = pairs[kk].i_node;
      size_t j = pairs[kk].j_node;
      size_t nodeListi = pairs[kk].i_list;
      size_t nodeListj = pairs[kk].j_list;

      // Get the state for node i.
      const auto& ri = position(nodeListi, i);
      const auto  mi = mass(nodeListi, i);
      const auto& vi = velocity(nodeListi, i);
      const auto  rhoi = massDensity(nodeListi, i);
      const auto  Pi = pressure(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto  ci = soundSpeed(nodeListi, i);
      const auto  omegai = omega(nodeListi, i);
      const auto  safeOmegai = safeInv(omegai, tiny);
      const auto  Hdeti = Hi.Determinant();
      CHECK(mi > 0.0);
      CHECK(rhoi > 0.0);
      CHECK(Hdeti > 0.0);

      auto& rhoSumi = rhoSum(nodeListi, i);
      auto& normi = normalization(nodeListi, i);
      auto& DvDti = DvDt(nodeListi, i);
      auto& DepsDti = DepsDt(nodeListi, i);
      auto& DvDxi = DvDx(nodeListi, i);
      auto& localDvDxi = localDvDx(nodeListi, i);
      auto& gradRhoi = gradRho(nodeListi, i);
      auto& Mi = M(nodeListi, i);
      auto& localMi = localM(nodeListi, i);
      auto& maxViscousPressurei = maxViscousPressure(nodeListi, i);
      auto& effViscousPressurei = effViscousPressure(nodeListi, i);
      auto& XSPHWeightSumi = XSPHWeightSum(nodeListi, i);
      auto& XSPHDeltaVi = XSPHDeltaV(nodeListi, i);

      // Get the state for node j
      const auto& rj = position(nodeListj, j);
      const auto  mj = mass(nodeListj, j);
      const auto& vj = velocity(nodeListj, j);
      const auto  rhoj = massDensity(nodeListj, j);
      const auto  Pj = pressure(nodeListj, j);
      const auto& Hj = H(nodeListj, j);
      const auto  cj = soundSpeed(nodeListj, j);
      const auto  omegaj = omega(nodeListj, j);
      const auto  safeOmegaj = safeInv(omegaj, tiny);
      const auto  Hdetj = Hj.Determinant();
      CHECK(mj > 0.0);
      CHECK(rhoj > 0.0);
      CHECK(Hdetj > 0.0);

      auto& rhoSumj = rhoSum(nodeListj, j);
      auto& normj = normalization(nodeListj, j);
      auto& DvDtj = DvDt(nodeListj, j);
      auto& DepsDtj = DepsDt(nodeListj, j);
      auto& DvDxj = DvDx(nodeListj, j);
      auto& localDvDxj = localDvDx(nodeListj, j);
      auto& gradRhoj = gradRho(nodeListj, j);
      auto& Mj = M(nodeListj, j);
      auto& localMj = localM(nodeListj, j);
      auto& maxViscousPressurej = maxViscousPressure(nodeListj, j);
      auto& effViscousPressurej = effViscousPressure(nodeListj, j);
      auto& XSPHWeightSumj = XSPHWeightSum(nodeListj, j);
      auto& XSPHDeltaVj = XSPHDeltaV(nodeListj, j);

      // Flag if this is a contiguous material pair or not.
      const bool sameMatij = true; // (nodeListi == nodeListj and fragIDi == fragIDj);

      // Node displacement.
      const auto rij = ri - rj;
      const auto etai = Hi*rij;
      const auto etaj = Hj*rij;
      const auto etaMagi = etai.magnitude();
      const auto etaMagj = etaj.magnitude();
      const auto etaiUnit = etai*safeInvVar(etaMagi);
      const auto etajUnit = etaj*safeInvVar(etaMagj);
      CHECK(etaMagi >= 0.0);
      CHECK(etaMagj >= 0.0);

      // Symmetrized kernel weight and gradient.
      W_view.kernelAndGradValue(etaMagi, Hdeti, Wi, gWi);
      W_view.kernelAndGradValue(etaMagj, Hdetj, Wj, gWj);
      gradWi = gWi*Hi*etaiUnit;
      gradWj = gWj*Hj*etajUnit;
      if (oneKernel) {
        WQi = Wi;
        WQj = Wj;
        gradWQi = gradWi;
        gradWQj = gradWj;
      } else {
        WQ_view.kernelAndGradValue(etaMagi, Hdeti, WQi, gWQi);
        WQ_view.kernelAndGradValue(etaMagj, Hdetj, WQj, gWQj);
        gradWQi = gWQi*Hi*etaiUnit;
        gradWQj = gWQj*Hj*etajUnit;
      }

      // Contribution to the sum density.
      if (nodeListi == nodeListj) {
        RAJA::atomicAdd<RAJA::auto_atomic>(&rhoSumi, mj*Wi);
        RAJA::atomicAdd<RAJA::auto_atomic>(&rhoSumj, mi*Wj);
        RAJA::atomicAdd<RAJA::auto_atomic>(&normi, mi/rhoi*Wi);
        RAJA::atomicAdd<RAJA::auto_atomic>(&normj, mj/rhoj*Wj);
      }

      // Compute the pair-wise artificial viscosity.
      const auto vij = vi - vj;
      Q->QPiij(QPiij, QPiji, Qi, Qj,
               nodeListi, i, nodeListj, j,
               ri, Hi, etai, vi, rhoi, ci,  
               rj, Hj, etaj, vj, rhoj, cj,
               fClQView, fCqQView, DvDxQView);

      const auto Qacci = 0.5*(QPiij*gradWQi);
      const auto Qaccj = 0.5*(QPiji*gradWQj);
      // const auto workQi = 0.5*(QPiij*vij).dot(gradWQi);
      // const auto workQj = 0.5*(QPiji*vij).dot(gradWQj);
      const auto workQi = vij.dot(Qacci);
      const auto workQj = vij.dot(Qaccj);
      RAJA::atomicMax<RAJA::auto_atomic>(&maxViscousPressurei, Qi);
      RAJA::atomicMax<RAJA::auto_atomic>(&maxViscousPressurej, Qj);
      RAJA::atomicAdd<RAJA::auto_atomic>(&effViscousPressurei, mj*Qi*WQi/rhoj);
      RAJA::atomicAdd<RAJA::auto_atomic>(&effViscousPressurej, mi*Qj*WQj/rhoi);

      // Determine an effective pressure including a term to fight the tensile instability.
      const auto Ri = epsTensile*FastMath::pow4(Wi/(Hdeti*WnPerh))*(Pi < 0.0 ? -Pi : 0.0);
      const auto Rj = epsTensile*FastMath::pow4(Wj/(Hdetj*WnPerh))*(Pj < 0.0 ? -Pj : 0.0);
      const auto Peffi = Pi + Ri;
      const auto Peffj = Pj + Rj;

      //Acceleration.
      CHECK(rhoi > 0.0);
      CHECK(rhoj > 0.0);
      const auto Prhoi = safeOmegai*Peffi/(rhoi*rhoi);
      const auto Prhoj = safeOmegaj*Peffj/(rhoj*rhoj);
      const auto deltaDvDt = Prhoi*gradWi + Prhoj*gradWj + Qacci + Qaccj;
      DvDti.atomicSub(mj*deltaDvDt);
      DvDtj.atomicAdd(mi*deltaDvDt);
      //if (compatibleEnergy) (*pairAccelerationsPtr)[kk] = -mj*deltaDvDt;  // Acceleration for i (j anti-symmetric)

      // Specific thermal energy evolution.
      // const Scalar workQij = 0.5*(mj*workQi + mi*workQj);
      RAJA::atomicAdd<RAJA::auto_atomic>(&DepsDti, mj*(Prhoi*vij.dot(gradWi) + workQi));
      RAJA::atomicAdd<RAJA::auto_atomic>(&DepsDtj, mi*(Prhoj*vij.dot(gradWj) + workQj));

      // Velocity gradient.
      const auto deltaDvDxi = mj*vij.dyad(gradWi);
      const auto deltaDvDxj = mi*vij.dyad(gradWj);
      DvDxi.atomicSub(deltaDvDxi);
      DvDxj.atomicSub(deltaDvDxj);
      if (sameMatij) {
        localDvDxi.atomicSub(deltaDvDxi);
        localDvDxj.atomicSub(deltaDvDxj);
      }

      // Estimate of delta v (for XSPH).
      if (XSPH and (sameMatij)) {
        const auto wXSPHij = 0.5*(mi/rhoi*Wi + mj/rhoj*Wj);
        RAJA::atomicAdd<RAJA::auto_atomic>(&XSPHWeightSumi, wXSPHij);
        RAJA::atomicAdd<RAJA::auto_atomic>(&XSPHWeightSumj, wXSPHij);
        XSPHDeltaVi.atomicSub(wXSPHij*vij);
        XSPHDeltaVj.atomicAdd(wXSPHij*vij);
      }

      // Mass density gradient
      if (sameMatij) {
        gradRhoi.atomicAdd(mj*(rhoj - rhoi)*gradWi);
        gradRhoj.atomicAdd(mi*(rhoj - rhoi)*gradWj);  // negatives cancel (rhoji and gradWj)
      }

      // Linear gradient correction term.
      Mi.atomicSub(mj*rij.dyad(gradWi));
      Mj.atomicSub(mi*rij.dyad(gradWj));
      if (sameMatij) {
        localMi.atomicSub(mj*rij.dyad(gradWi));
        localMj.atomicSub(mi*rij.dyad(gradWj));
      }

    }); // loop over pairs

    // Reduce the thread values to the master.
    //threadReduceFieldLists<Dimension>(threadStack);

  }   // OpenMP parallel region
  TIME_END("SPHevalDerivs_pairs");
  // Finish up the derivatives for each point.
  TIME_BEGIN("SPHevalDerivs_final");
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto& nodeList = mass_v[nodeListi]->nodeList();
    const auto ni = nodeList.numInternalNodes();
// #pragma omp parallel for
//     for (auto i = 0u; i < ni; ++i) {
    RAJA::forall<EXEC_POLICY>(TRS_UINT(0u, ni),
    [=] SPHERAL_HOST_DEVICE (size_t i) {

      // Get the state for node i.
      const auto& mi = mass(nodeListi, i);
      const auto& vi = velocity(nodeListi, i);
      const auto& rhoi = massDensity(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();
      //const auto  numNeighborsi = connectivityMap.numNeighborsForNode(nodeListi, i);
      CHECK(mi > 0.0);
      CHECK(rhoi > 0.0);
      CHECK(Hdeti > 0.0);

      auto& rhoSumi = rhoSum(nodeListi, i);
      auto& normi = normalization(nodeListi, i);
      auto& DxDti = DxDt(nodeListi, i);
      auto& DrhoDti = DrhoDt(nodeListi, i);
      auto& DvDti = DvDt(nodeListi, i);
      auto& DepsDti = DepsDt(nodeListi, i);
      auto& DvDxi = DvDx(nodeListi, i);
      auto& localDvDxi = localDvDx(nodeListi, i);
      auto& gradRhoi = gradRho(nodeListi, i);
      auto& Mi = M(nodeListi, i);
      auto& localMi = localM(nodeListi, i);
      auto& XSPHWeightSumi = XSPHWeightSum(nodeListi, i);
      auto& XSPHDeltaVi = XSPHDeltaV(nodeListi, i);

      // Add the self-contribution to density sum.
      rhoSumi += mi*W0*Hdeti;
      normi += mi/rhoi*W0*Hdeti;

      // Finish the gradient of the velocity.
      CHECK(rhoi > 0.0);
      if (CorrectVelocityGradient and
          std::abs(Mi.Determinant()) > 1.0e-10) {
          //and numNeighborsi > Dimension::pownu(2)) {
        Mi = Mi.Inverse();
        DvDxi = DvDxi*Mi;
      } else {
        DvDxi /= rhoi;
      }
      if (CorrectVelocityGradient and
          std::abs(localMi.Determinant()) > 1.0e-10) {
          //and numNeighborsi > Dimension::pownu(2)) {
        localMi = localMi.Inverse();
        localDvDxi = localDvDxi*localMi;
      } else {
        localDvDxi /= rhoi;
      }

      // Finish the mass density gradient
      gradRhoi /= rhoi;

      // Evaluate the continuity equation.
      DrhoDti = -rhoi*DvDxi.Trace();

      // If needed finish the total energy derivative.
      if (evolveTotalEnergy) DepsDti = mi*(vi.dot(DvDti) + DepsDti);

      // Determine the position evolution, based on whether we're doing XSPH or not.
      if (XSPH) {
        XSPHWeightSumi += Hdeti*mi/rhoi*W0;
        CHECK2(XSPHWeightSumi != 0.0, i << " " << XSPHWeightSumi);
        DxDti = vi + XSPHDeltaVi/max(tiny, XSPHWeightSumi);
      } else {
        DxDti = vi;
      }
    });
  }
  rhoSum.move(chai::CPU);
  normalization.move(chai::CPU);
  DxDt.move(chai::CPU);
  DrhoDt.move(chai::CPU);
  DvDt.move(chai::CPU);
  DepsDt.move(chai::CPU);
  DvDx.move(chai::CPU);
  localDvDx.move(chai::CPU);
  gradRho.move(chai::CPU);
  M.move(chai::CPU);
  localM.move(chai::CPU);
  maxViscousPressure.move(chai::CPU);
  effViscousPressure.move(chai::CPU);
  XSPHWeightSum.move(chai::CPU);
  XSPHDeltaV.move(chai::CPU);
  TIME_END("SPHevalDerivs_final");
  TIME_END("SPHevalDerivs");
}

}
