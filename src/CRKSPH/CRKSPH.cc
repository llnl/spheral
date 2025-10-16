//---------------------------------Spheral++----------------------------------//
// CRKSPHBase -- Base class for the CRKSPH/ACRKSPH hydrodynamic packages
//
// Created by JMO, Mon Jul 19 21:52:29 PDT 2010
//----------------------------------------------------------------------------//
#include "CRKSPH/CRKSPH.hh"

#include "FileIO/FileIO.hh"
#include "RK/ReproducingKernel.hh"
#include "RK/RKFieldNames.hh"
#include "CRKSPH/computeCRKSPHSumMassDensity.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Hydro/entropyWeightingFunction.hh"
#include "Strength/SolidFieldNames.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/IncrementBoundedState.hh"
#include "DataBase/ReplaceState.hh"
#include "DataBase/ReplaceBoundedState.hh"
#include "DataBase/updateStateFields.hh"
#include "Hydro/SpecificThermalEnergyPolicy.hh"
#include "Hydro/SpecificFromTotalThermalEnergyPolicy.hh"
#include "Hydro/PressurePolicy.hh"
#include "Hydro/SoundSpeedPolicy.hh"
#include "Hydro/EntropyPolicy.hh"
#include "RK/ContinuityVolumePolicy.hh"
#include "ArtificialViscosity/ArtificialViscosity.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "Boundary/Boundary.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Neighbor/PairwiseField.hh"
#include "Utilities/Timer.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/range.hh"
#include "Utilities/newtonRaphson.hh"
#include "Utilities/SpheralFunctions.hh"
#include "Utilities/computeShepardsInterpolation.hh"
#include "Geometry/innerProduct.hh"
#include "Geometry/outerProduct.hh"

#include <algorithm>
#include <fstream>
#include <map>
#include <vector>

using std::vector;
using std::string;
using std::pair;
using std::make_pair;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {

//------------------------------------------------------------------------------
// Construct with the given artificial viscosity and kernels.
//------------------------------------------------------------------------------
template<typename Dimension>
CRKSPH<Dimension>::
CRKSPH(DataBase<Dimension>& dataBase,
       ArtificialViscosity<Dimension>& Q,
       const RKOrder order,
       const double cfl,
       const bool useVelocityMagnitudeForDt,
       const bool compatibleEnergyEvolution,
       const bool evolveTotalEnergy,
       const bool XSPH,
       const MassDensityType densityUpdate,
       const double epsTensile,
       const double nTensile):
  CRKSPHBase<Dimension>(dataBase,
                        Q,
                        order,
                        cfl,
                        useVelocityMagnitudeForDt,
                        compatibleEnergyEvolution,
                        evolveTotalEnergy,
                        XSPH,
                        densityUpdate,
                        epsTensile,
                        nTensile),
  mPairAccelerationsPtr() {
}

//------------------------------------------------------------------------------
// Register the state we need/are going to evolve.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPH<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {
  TIME_BEGIN("CRKregisterState");

  CRKSPHBase<Dimension>::registerState(dataBase, state);

  // We have to choose either compatible or total energy evolution.
  const auto compatibleEnergy = this->compatibleEnergyEvolution();
  const auto evolveTotalEnergy = this->evolveTotalEnergy();
  VERIFY2(not (compatibleEnergy and evolveTotalEnergy),
          "CRKSPH error : you cannot simultaneously use both compatibleEnergyEvolution and evolveTotalEnergy");

  // Register the specific thermal energy
  auto specificThermalEnergy = dataBase.fluidSpecificThermalEnergy();
  if (compatibleEnergy) {
    // The compatible energy update.
    state.enroll(specificThermalEnergy, make_policy<SpecificThermalEnergyPolicy<Dimension>>(dataBase));

  } else if (evolveTotalEnergy) {
    // If we're doing total energy, we register the specific energy to advance with the
    // total energy policy.
    state.enroll(specificThermalEnergy, make_policy<SpecificFromTotalThermalEnergyPolicy<Dimension>>());

  } else {
    // Otherwise we're just time-evolving the specific energy.
    state.enroll(specificThermalEnergy, make_policy<IncrementState<Dimension, Scalar>>());
  }
  TIME_END("CRKregisterState");
}

//------------------------------------------------------------------------------
// Register the state derivative fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPH<Dimension>::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {
  TIME_BEGIN("CRKregisterDerivatives");

  CRKSPHBase<Dimension>::registerDerivatives(dataBase, derivs);
  const auto compatibleEnergy = this->compatibleEnergyEvolution();
  if (compatibleEnergy) {
    const auto& connectivityMap = dataBase.connectivityMap();
    mPairAccelerationsPtr = std::make_unique<PairAccelerationsType>(connectivityMap);
    derivs.enroll(HydroFieldNames::pairAccelerations, *mPairAccelerationsPtr);
  }
  TIME_END("CRKregisterDerivatives");
}

//------------------------------------------------------------------------------
// Determine the principle derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPH<Dimension>::
evaluateDerivatives(const typename Dimension::Scalar time,
                    const typename Dimension::Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivatives) const {
  TIME_BEGIN("CRKevaluateDerivatives");

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
  TIME_END("CRKevaluateDerivatives");
}
  
//------------------------------------------------------------------------------
// Determine the principle derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename QType>
void
CRKSPH<Dimension>::
evaluateDerivativesImpl(const typename Dimension::Scalar /*time*/,
                        const typename Dimension::Scalar /*dt*/,
                        const DataBase<Dimension>& dataBase,
                        const State<Dimension>& state,
                        StateDerivatives<Dimension>& derivs,
                        chai::managed_ptr<QType> Q) const {
  TIME_BEGIN("CRKevaluateDerivativesImpl");

  using QPiType = typename QType::ReturnType;

  // The kernels and such.
  const auto order = this->correctionOrder();
  const auto& WR = state.template get<ReproducingKernel<Dimension>>(RKFieldNames::reproducingKernel(order));

  // A few useful constants we'll use in the following loop.
  //const double tiny = 1.0e-30;
  const auto compatibleEnergy = this->compatibleEnergyEvolution();
  const auto evolveTotalEnergy = this->evolveTotalEnergy();
  const auto XSPH = this->XSPH();

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
  const auto position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const auto massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  const auto H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const auto pressure = state.fields(HydroFieldNames::pressure, 0.0);
  const auto soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  const auto corrections = state.fields(RKFieldNames::rkCorrections(order), RKCoefficients<Dimension>());
  const auto surfacePoint = state.fields(HydroFieldNames::surfacePoint, 0);
  const auto fClQ = state.fields(HydroFieldNames::ArtificialViscousClMultiplier, 0.0, true);
  const auto fCqQ = state.fields(HydroFieldNames::ArtificialViscousCqMultiplier, 0.0, true);
  const auto DvDxQ = state.fields(HydroFieldNames::ArtificialViscosityVelocityGradient, Tensor::zero, true);
  CHECK(mass.size() == numNodeLists);
  CHECK(position.size() == numNodeLists);
  CHECK(velocity.size() == numNodeLists);
  CHECK(massDensity.size() == numNodeLists);
  CHECK(specificThermalEnergy.size() == numNodeLists);
  CHECK(H.size() == numNodeLists);
  CHECK(pressure.size() == numNodeLists);
  CHECK(soundSpeed.size() == numNodeLists);
  CHECK(corrections.size() == numNodeLists);
  CHECK(surfacePoint.size() == numNodeLists);
  CHECK(fClQ.size() == 0 or fClQ.size() == numNodeLists);
  CHECK(fCqQ.size() == 0 or fCqQ.size() == numNodeLists);
  CHECK(DvDxQ.size() == 0 or DvDxQ.size() == numNodeLists);

  // Derivative FieldLists.
  auto  DxDt = derivs.fields(IncrementState<Dimension, Vector>::prefix() + HydroFieldNames::position, Vector::zero);
  auto  DrhoDt = derivs.fields(IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, 0.0);
  auto  DvDt = derivs.fields(HydroFieldNames::hydroAcceleration, Vector::zero);
  auto  DepsDt = derivs.fields(IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
  auto  DvDx = derivs.fields(HydroFieldNames::velocityGradient, Tensor::zero);
  auto  localDvDx = derivs.fields(HydroFieldNames::internalVelocityGradient, Tensor::zero);
  auto  maxViscousPressure = derivs.fields(HydroFieldNames::maxViscousPressure, 0.0);
  auto  effViscousPressure = derivs.fields(HydroFieldNames::effectiveViscousPressure, 0.0);
  auto  XSPHDeltaV = derivs.fields(HydroFieldNames::XSPHDeltaV, Vector::zero);
  auto* pairAccelerationsPtr = derivs.template getPtr<PairAccelerationsType>(HydroFieldNames::pairAccelerations);
  CHECK(DxDt.size() == numNodeLists);
  CHECK(DrhoDt.size() == numNodeLists);
  CHECK(DvDt.size() == numNodeLists);
  CHECK(DepsDt.size() == numNodeLists);
  CHECK(DvDx.size() == numNodeLists);
  CHECK(localDvDx.size() == numNodeLists);
  CHECK(maxViscousPressure.size() == numNodeLists);
  CHECK(effViscousPressure.size() == numNodeLists);
  CHECK(XSPHDeltaV.size() == numNodeLists);
  CHECK((compatibleEnergy and pairAccelerationsPtr->size() == npairs) or not compatibleEnergy);

  // Walk all the interacting pairs.
#pragma omp parallel
  {
    // Thread private scratch variables
    int i, j, nodeListi, nodeListj;
    Scalar Wi, Wj, Qi, Qj;
    QPiType QPiij, QPiji;
    Vector gradWi, gradWj;
    Vector deltagrad, forceij, forceji;
    Vector rij, vij, etai, etaj;

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
      const auto& ri = position(nodeListi, i);
      const auto  mi = mass(nodeListi, i);
      const auto& vi = velocity(nodeListi, i);
      const auto  rhoi = massDensity(nodeListi, i);
      //const auto  epsi = specificThermalEnergy(nodeListi, i);
      const auto  Pi = pressure(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto  ci = soundSpeed(nodeListi, i);
      const auto& correctionsi = corrections(nodeListi, i);
      const auto  weighti = volume(nodeListi, i);  // Change CRKSPH weights here if need be!
      CHECK2(mi > 0.0, i << " " << mi);
      CHECK2(rhoi > 0.0, i << " " << rhoi);
      CHECK2(weighti > 0.0, i << " " << weighti);

      auto& DvDti = DvDt_thread(nodeListi, i);
      auto& DepsDti = DepsDt_thread(nodeListi, i);
      auto& DvDxi = DvDx_thread(nodeListi, i);
      auto& localDvDxi = localDvDx_thread(nodeListi, i);
      auto& maxViscousPressurei = maxViscousPressure_thread(nodeListi, i);
      auto& effViscousPressurei = effViscousPressure_thread(nodeListi, i);
      auto& XSPHDeltaVi = XSPHDeltaV_thread(nodeListi, i);

      // Get the state for node j
      const auto& rj = position(nodeListj, j);
      const auto  mj = mass(nodeListj, j);
      const auto& vj = velocity(nodeListj, j);
      const auto  rhoj = massDensity(nodeListj, j);
      //const auto  epsj = specificThermalEnergy(nodeListj, j);
      const auto  Pj = pressure(nodeListj, j);
      const auto& Hj = H(nodeListj, j);
      const auto  cj = soundSpeed(nodeListj, j);
      const auto& correctionsj = corrections(nodeListj, j);
      const auto  Hdetj = Hj.Determinant();
      const auto  weightj = volume(nodeListj, j);     // Change CRKSPH weights here if need be!
      CONTRACT_VAR(Hdetj);
      CHECK(mj > 0.0);
      CHECK(rhoj > 0.0);
      CHECK(Hdetj > 0.0);
      CHECK(weightj > 0.0);

      auto& DvDtj = DvDt_thread(nodeListj, j);
      auto& DepsDtj = DepsDt_thread(nodeListj, j);
      auto& DvDxj = DvDx_thread(nodeListj, j);
      auto& localDvDxj = localDvDx_thread(nodeListj, j);
      auto& maxViscousPressurej = maxViscousPressure_thread(nodeListj, j);
      auto& effViscousPressurej = effViscousPressure_thread(nodeListj, j);
      auto& XSPHDeltaVj = XSPHDeltaV_thread(nodeListj, j);

      // Node displacement.
      rij = ri - rj;
      vij = vi - vj;
      etai = Hi*rij;
      etaj = Hj*rij;

      // Symmetrized kernel weight and gradient.
      std::tie(Wj, gradWj) = WR.evaluateKernelAndGradient( rij, Hj, correctionsi);  // Hj because we compute RK using scatter formalism
      std::tie(Wi, gradWi) = WR.evaluateKernelAndGradient(-rij, Hi, correctionsj);
      deltagrad = gradWj - gradWi;

      // Compute the artificial viscous pressure (Pi = P/rho^2 actually).
      Q->QPiij(QPiij, QPiji, Qi, Qj,
               nodeListi, i, nodeListj, j,
               ri, Hi, etai, vi, rhoi, ci,  
               rj, Hj, etaj, vj, rhoj, cj,
               fClQ, fCqQ, DvDxQ);
      const auto Qaccij = (rhoi*rhoi*QPiij + rhoj*rhoj*QPiji)*deltagrad;
      // const auto workQij = 0.5*(vij.dot(Qaccij));
      const auto workQi = (rhoj*rhoj*QPiji*vij).dot(deltagrad);                   // CRK
      const auto workQj = (rhoi*rhoi*QPiij*vij).dot(deltagrad);                   // CRK
      // const auto workQVi =  vij.dot((rhoj*rhoj*QPiji).dot(gradWj));            //RK V and RK I Work
      // const auto workQVj =  vij.dot((rhoi*rhoi*QPiij).dot(gradWi));            //RK V and RK I Work
      maxViscousPressurei = max(maxViscousPressurei, 4.0*Qi);                     // We need tighter timestep controls on the Q with CRK
      maxViscousPressurej = max(maxViscousPressurej, 4.0*Qj);
      effViscousPressurei += weightj * Qi * Wj;
      effViscousPressurej += weighti * Qj * Wi;

      // Velocity gradient.
      DvDxi -= weightj*vij.dyad(gradWj);
      DvDxj += weighti*vij.dyad(gradWi);
      if (nodeListi == nodeListj) {
        localDvDxi -= weightj*vij.dyad(gradWj);
        localDvDxj += weighti*vij.dyad(gradWi);
      }

      // // Mass density gradient.
      // gradRhoi += weightj*(rhoj - rhoi)*gradWj;
      // gradRhoj += weighti*(rhoi - rhoj)*gradWi;

      // We decide between RK and CRK for the momentum and energy equations based on the surface condition.
      // Momentum
      forceij = (true ? // surfacePoint(nodeListi, i) <= 1 ? 
                 0.5*weighti*weightj*((Pi + Pj)*deltagrad + Qaccij) :                // Type III CRK interpoint force.
                 mi*weightj*((Pj - Pi)/rhoi*gradWj + rhoi*QPiij*gradWj));            // RK
      forceji = (true ? // surfacePoint(nodeListj, j) <= 1 ? 
                 0.5*weighti*weightj*((Pi + Pj)*deltagrad + Qaccij) :                // Type III CRK interpoint force.
                 mj*weighti*((Pj - Pi)/rhoj*gradWi - rhoj*QPiji*gradWi));            // RK
      DvDti -= forceij/mi;
      DvDtj += forceji/mj; 
      if (compatibleEnergy) (*pairAccelerationsPtr)[kk] = -forceij/mi;               // Acceleration for i (j anti-symmetric)

      // Energy
      DepsDti += (true ? // surfacePoint(nodeListi, i) <= 1 ? 
                  0.5*weighti*weightj*(Pj*vij.dot(deltagrad) + workQi)/mi :          // CRK
                  (weightj*rhoi*QPiij*vij).dot(gradWj));                             // RK
      DepsDtj += (true ? // surfacePoint(nodeListj, j) <= 1 ? 
                  0.5*weighti*weightj*(Pi*vij.dot(deltagrad) + workQj)/mj :          // CRK
                  (-weighti*rhoj*QPiji*vij).dot(gradWi));                            // RK

      // Estimate of delta v (for XSPH).
      if (XSPH and (nodeListi == nodeListj)) {
        XSPHDeltaVi -= weightj*Wj*vij;
        XSPHDeltaVj += weighti*Wi*vij;
      }
    }

    // Reduce the thread values to the master.
    threadReduceFieldLists<Dimension>(threadStack);

  }   // OMP parallel

  // Finish up the derivatives for each point.
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto& nodeList = mass[nodeListi]->nodeList();
    const auto ni = nodeList.numInternalNodes();
#pragma omp parallel for
    for (auto i = 0u; i < ni; ++i) {

      // Get the state for node i.
      const auto& mi = mass(nodeListi, i);
      const auto& vi = velocity(nodeListi, i);
      const auto& rhoi = massDensity(nodeListi, i);
      CHECK(mi > 0.0);
      CHECK(rhoi > 0.0);

      auto& DxDti = DxDt(nodeListi, i);
      auto& DrhoDti = DrhoDt(nodeListi, i);
      auto& DvDti = DvDt(nodeListi, i);
      auto& DepsDti = DepsDt(nodeListi, i);
      auto& DvDxi = DvDx(nodeListi, i);
      auto& XSPHDeltaVi = XSPHDeltaV(nodeListi, i);

      // Determine the position evolution, based on whether we're doing XSPH or not.
      if (XSPH) {
        DxDti = vi + XSPHDeltaVi;
      } else {
        DxDti = vi;
      }

      // Time evolution of the mass density.
      DrhoDti = -rhoi*DvDxi.Trace();

      // If needed finish the total energy derivative.
      if (evolveTotalEnergy) DepsDti = mi*(vi.dot(DvDti) + DepsDti);
    }
  }
  TIME_END("CRKevaluateDerivativesImpl");
}

}

