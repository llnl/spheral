//---------------------------------Spheral++----------------------------------//
// SPH -- The classic SPH/ASPH hydrodynamic packages for Spheral++.
// 
// Created by JMO, Thu Nov 21 16:36:40 PST 2024
//----------------------------------------------------------------------------//
#include "config.hh"
#include "SPH/SPH_RAJA.hh"
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
#include "Threading/ViewManager.hh"

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
SPH_RAJA<Dimension>::
SPH_RAJA(DataBase<Dimension>& dataBase,
         ArtificialViscosity<Dimension>& Q,
         const TableKernel<Dimension>& W,
         const TableKernel<Dimension>& WPi,
         const double cfl,
         const bool useVelocityMagnitudeForDt,
         const bool useNewAccelerationMagnitudeForDt,
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
  SPH<Dimension>(dataBase,
                 Q,
                 W,
                 WPi,
                 cfl,
                 useVelocityMagnitudeForDt,
                 useNewAccelerationMagnitudeForDt,
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
                 xmax) {
}

//------------------------------------------------------------------------------
// Determine the principle derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SPH_RAJA<Dimension>::
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
SPH_RAJA<Dimension>::
evaluateDerivativesImpl(const typename Dimension::Scalar time,
                        const typename Dimension::Scalar dt,
                        const DataBase<Dimension>& dataBase,
                        const State<Dimension>& state,
                        StateDerivatives<Dimension>& derivs,
                        chai::managed_ptr<QType>& Q) const {
  TIME_BEGIN("SPHevalDerivs_RAJA");
  TIME_BEGIN("SPHevalDerivs_initial_RAJA");

  using QPiType = typename QType::ReturnType;

  //static double totalLoopTime = 0.0;

  // The kernels and such.
  auto& W = this->kernel();
  auto& WQ = this->PiKernel();
  auto W_view = W.view();
  auto WQ_view = WQ.view();
  const auto oneKernel = (W == WQ);

  // A few useful constants we'll use in the following loop.
  const auto tiny = 1.0e-30;
  const auto W0 = W(0.0, 1.0);
  const auto epsTensile = this->epsilonTensile();
  const auto compatibleEnergy = this->compatibleEnergyEvolution();
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

  // Prepare the ViewManagers
  ViewManager<Dimension> stateMgr(state);
  ViewManager<Dimension> derivsMgr(derivs);

  // Get the state and derivative FieldLists.
  // State FieldLists.
  auto mass = stateMgr.fields(HydroFieldNames::mass, 0.0);
  auto position = stateMgr.fields(HydroFieldNames::position, Vector::zero());
  auto velocity = stateMgr.fields(HydroFieldNames::velocity, Vector::zero());
  auto massDensity = stateMgr.fields(HydroFieldNames::massDensity, 0.0);
  auto H = stateMgr.fields(HydroFieldNames::H, SymTensor::zero());
  auto pressure = stateMgr.fields(HydroFieldNames::pressure, 0.0);
  auto soundSpeed = stateMgr.fields(HydroFieldNames::soundSpeed, 0.0);
  auto omega = stateMgr.fields(HydroFieldNames::omegaGradh, 0.0);
  auto fClQ = stateMgr.fields(HydroFieldNames::ArtificialViscousClMultiplier, 0.0, true);
  auto fCqQ = stateMgr.fields(HydroFieldNames::ArtificialViscousCqMultiplier, 0.0, true);
  auto DvDxQ = stateMgr.fields(HydroFieldNames::ArtificialViscosityVelocityGradient, Tensor::zero(), true);
  CHECK(mass.size() == numNodeLists);
  CHECK(position.size() == numNodeLists);
  CHECK(velocity.size() == numNodeLists);
  CHECK(massDensity.size() == numNodeLists);
  CHECK(H.size() == numNodeLists);
  CHECK(pressure.size() == numNodeLists);
  CHECK(soundSpeed.size() == numNodeLists);
  CHECK(omega.size() == numNodeLists);
  CHECK(fClQ.size() == 0 or fClQ.size() == numNodeLists);
  CHECK(fCqQ.size() == 0 or fCqQ.size() == numNodeLists);
  CHECK(DvDxQ.size() == 0 or DvDxQ.size() == numNodeLists);

  // Derivative FieldLists.
  auto rhoSum = derivsMgr.fields(ReplaceState<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, 0.0);
  auto normalization = derivsMgr.fields(HydroFieldNames::normalization, 0.0);
  auto DxDt = derivsMgr.fields(IncrementState<Dimension, Vector>::prefix() + HydroFieldNames::position, Vector::zero());
  auto DrhoDt = derivsMgr.fields(IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, 0.0);
  auto DvDt = derivsMgr.fields(HydroFieldNames::hydroAcceleration, Vector::zero());
  auto DepsDt = derivsMgr.fields(IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
  auto DvDx = derivsMgr.fields(HydroFieldNames::velocityGradient, Tensor::zero());
  auto localDvDx = derivsMgr.fields(HydroFieldNames::internalVelocityGradient, Tensor::zero());
  auto gradRho = derivsMgr.fields(HydroFieldNames::massDensityGradient, Vector::zero());
  auto M = derivsMgr.fields(HydroFieldNames::M_SPHCorrection, Tensor::zero());
  auto localM = derivsMgr.fields("local " + HydroFieldNames::M_SPHCorrection, Tensor::zero());
  auto maxViscousPressure = derivsMgr.fields(HydroFieldNames::maxViscousPressure, 0.0);
  auto effViscousPressure = derivsMgr.fields(HydroFieldNames::effectiveViscousPressure, 0.0);
  auto XSPHWeightSum = derivsMgr.fields(HydroFieldNames::XSPHWeightSum, 0.0);
  auto XSPHDeltaV = derivsMgr.fields(HydroFieldNames::XSPHDeltaV, Vector::zero());
  auto pairAccelerations = derivsMgr.template get<PairAccelerationsType>(HydroFieldNames::pairAccelerations);
  CHECK(rhoSum.size() == numNodeLists);
  CHECK(normalization.size() == numNodeLists);
  CHECK(DxDt.size() == numNodeLists);
  CHECK(DrhoDt.size() == numNodeLists);
  CHECK(DvDt.size() == numNodeLists);
  CHECK(DepsDt.size() == numNodeLists);
  CHECK(DvDx.size() == numNodeLists);
  CHECK(localDvDx.size() == numNodeLists);
  CHECK(gradRho.size() == numNodeLists);
  CHECK(M.size() == numNodeLists);
  CHECK(localM.size() == numNodeLists);
  CHECK(maxViscousPressure.size() == numNodeLists);
  CHECK(effViscousPressure.size() == numNodeLists);
  CHECK(XSPHWeightSum.size() == numNodeLists);
  CHECK(XSPHDeltaV.size() == numNodeLists);
  CHECK((compatibleEnergy and pairAccelerations.size() == npairs) or not compatibleEnergy);

  // The scale for the tensile correction.
  const auto  nPerh = nodeLists[0]->nodesPerSmoothingScale();
  const auto  WnPerh = W(1.0/nPerh, 1.0);
  bool CorrectVelocityGradient = this->mCorrectVelocityGradient;
  TIME_END("SPHevalDerivs_initial_RAJA");

  // Walk all the interacting pairs.
  TIME_BEGIN("SPHevalDerivs_pairs_RAJA");
  {
    RAJA::forall<EXEC_POLICY>(TRS_UINT(0u, npairs),
    [=] SPHERAL_HOST_DEVICE (size_t kk) {
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
      Scalar Wi, Wj, gWi, gWj;
      W_view.kernelAndGradValue(etaMagi, Hdeti, Wi, gWi);
      W_view.kernelAndGradValue(etaMagj, Hdetj, Wj, gWj);
      Vector gradWi = gWi*Hi*etaiUnit;
      Vector gradWj = gWj*Hj*etajUnit;
      Scalar WQi, WQj;
      Vector gradWQi, gradWQj;
      if (oneKernel) {
        WQi = Wi;
        WQj = Wj;
        gradWQi = gradWi;
        gradWQj = gradWj;
      } else {
        Scalar gWQi, gWQj;
        WQ_view.kernelAndGradValue(etaMagi, Hdeti, WQi, gWQi);
        WQ_view.kernelAndGradValue(etaMagj, Hdetj, WQj, gWQj);
        gradWQi = gWQi*Hi*etaiUnit;
        gradWQj = gWQj*Hj*etajUnit;
      }

      // Compute the pair-wise artificial viscosity.
      const auto vij = vi - vj;
      QPiType QPiij(0.0);
      QPiType QPiji(0.0);
      Scalar Qi = 0.0;
      Scalar Qj = 0.0;
      Q->QPiij(QPiij, QPiji, Qi, Qj,
               nodeListi, i, nodeListj, j,
               ri, Hi, etai, vi, rhoi, ci,  
               rj, Hj, etaj, vj, rhoj, cj,
               fClQ, fCqQ, DvDxQ);

      // Contribution to the sum density.
      if (nodeListi == nodeListj) {
        GPUUtils::AtomicAddOp::apply(&rhoSumi, mj*Wi);
        GPUUtils::AtomicAddOp::apply(&rhoSumj, mi*Wj);
        GPUUtils::AtomicAddOp::apply(&normi, mi/rhoi*Wi);
        GPUUtils::AtomicAddOp::apply(&normj, mj/rhoj*Wj);
      }

      const auto Qacci = 0.5*(QPiij*gradWQi);
      const auto Qaccj = 0.5*(QPiji*gradWQj);
      const auto workQi = vij.dot(Qacci);
      const auto workQj = vij.dot(Qaccj);
      GPUUtils::AtomicMaxOp::apply(&maxViscousPressurei, Qi);
      GPUUtils::AtomicMaxOp::apply(&maxViscousPressurej, Qj);
      GPUUtils::AtomicAddOp::apply(&effViscousPressurei, mj*Qi*WQi/rhoj);
      GPUUtils::AtomicAddOp::apply(&effViscousPressurej, mi*Qj*WQj/rhoi);

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
      if (compatibleEnergy) pairAccelerations[kk] = -mj*deltaDvDt;  // Acceleration for i (j anti-symmetric)

      // Specific thermal energy evolution.
      GPUUtils::AtomicAddOp::apply(&DepsDti, mj*(Prhoi*vij.dot(gradWi) + workQi));
      GPUUtils::AtomicAddOp::apply(&DepsDtj, mi*(Prhoj*vij.dot(gradWj) + workQj));

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
        GPUUtils::AtomicAddOp::apply(&XSPHWeightSumi, wXSPHij);
        GPUUtils::AtomicAddOp::apply(&XSPHWeightSumj, wXSPHij);
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
    GPU_ERROR_CHECK
  }
  TIME_END("SPHevalDerivs_pairs_RAJA");

  // Finish up the derivatives for each point.
  TIME_BEGIN("SPHevalDerivs_final_RAJA");
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto ni = nodeLists[nodeListi]->numInternalNodes();
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
    GPU_ERROR_CHECK
  }
  derivsMgr.move(chai::CPU);
  TIME_END("SPHevalDerivs_final_RAJA");
  TIME_END("SPHevalDerivs_RAJA");
}

}
