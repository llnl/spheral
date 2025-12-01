//---------------------------------Spheral++----------------------------------//
// SolidSphericalSPHHydroBase -- The SPH/ASPH solid material SPH hydrodynamic
//                               specialized for 1D Spherical (r) geometry.
//
// Based on the algorithm described in
// Omang, M., Børve, S., & Trulsen, J. (2006). SPH in spherical and cylindrical coordinates.
// Journal of Computational Physics, 213(1), 391412. https://doi.org/10.1016/j.jcp.2005.08.023
//
// Note this version is currently abusing our ordinary 1D geometric types,
// implicitly mapping x->r.
//
// Created by JMO, Tue Apr 26 16:28:55 PDT 2022
//----------------------------------------------------------------------------//
#include "SPH/SolidSphericalSPH.hh"
#include "FileIO/FileIO.hh"
#include "computeSPHSumMassDensity.hh"
#include "correctSPHSumMassDensity.hh"
#include "Utilities/NodeCoupling.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Hydro/SpecificThermalEnergyPolicy.hh"
#include "Hydro/SpecificFromTotalThermalEnergyPolicy.hh"
#include "Strength/SolidFieldNames.hh"
#include "NodeList/SolidNodeList.hh"
#include "Strength/DeviatoricStressPolicy.hh"
#include "Strength/BulkModulusPolicy.hh"
#include "Strength/PlasticStrainPolicy.hh"
#include "Strength/ShearModulusPolicy.hh"
#include "Strength/YieldStrengthPolicy.hh"
#include "Hydro/SphericalPositionPolicy.hh"
#include "Hydro/NonSymmetricSpecificThermalEnergyPolicy.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/IncrementBoundedState.hh"
#include "DataBase/ReplaceState.hh"
#include "DataBase/ReplaceBoundedState.hh"
#include "ArtificialViscosity/ArtificialViscosity.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "Boundary/Boundary.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Neighbor/PairwiseField.hh"
#include "Utilities/timingUtilities.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/range.hh"
#include "SolidMaterial/SolidEquationOfState.hh"
#include "Utilities/Timer.hh"

#include <algorithm>
#include <fstream>
#include <vector>
using std::vector;
using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {

namespace {

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
    return Dim<1>::SymTensor::zero();
  }
}

//------------------------------------------------------------------------------
// Compute a scalar sum for scalars or tensors
//------------------------------------------------------------------------------
Dim<1>::Scalar scalarSum(const Dim<1>::Scalar& x, const Dim<1>::Scalar& y) { return x + y; }
Dim<1>::Scalar scalarSum(const Dim<1>::Tensor& x, const Dim<1>::Tensor& y) { return x.xx() + y.xx(); }

}

//------------------------------------------------------------------------------
// Construct with the given artificial viscosity and kernels.
//------------------------------------------------------------------------------
SolidSphericalSPH::
SolidSphericalSPH(DataBase<Dimension>& dataBase,
                  ArtificialViscosityHandle<Dimension>& Q,
                  const SphericalKernel& W,
                  const SphericalKernel& WPi,
                  const SphericalKernel& WGrad,
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
                  const bool damageRelieveRubble,
                  const bool strengthInDamage,
                  const Vector& xmin,
                  const Vector& xmax):
  SolidSPH<Dim<1>>(dataBase,
                   Q,
                   W.baseKernel1d(),
                   WPi.baseKernel1d(),
                   WGrad.baseKernel1d(),
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
                   damageRelieveRubble,
                   strengthInDamage,
                   xmin,
                   xmax),
  mQself(2.0),
  mKernel(W),
  mPiKernel(WPi),
  mGradKernel(WGrad),
  mPairAccelerationsPtr(),
  mSelfAccelerations(FieldStorageType::CopyFields) {
}

//------------------------------------------------------------------------------
// Register the state we need/are going to evolve.
//------------------------------------------------------------------------------
void
SolidSphericalSPH::
registerState(DataBase<Dim<1>>& dataBase,
              State<Dim<1>>& state) {

  // The base class does most of it.
  SolidSPH<Dim<1>>::registerState(dataBase, state);

  // Re-regsiter the position update to prevent things going through the origin
  auto position = dataBase.fluidPosition();
  state.enroll(position, make_policy<SphericalPositionPolicy>());

  // We have to choose either compatible or total energy evolution.
  const auto compatibleEnergy = this->compatibleEnergyEvolution();
  const auto evolveTotalEnergy = this->evolveTotalEnergy();
  VERIFY2(not (compatibleEnergy and evolveTotalEnergy),
          "SPH error : you cannot simultaneously use both compatibleEnergyEvolution and evolveTotalEnergy");

  // Register the specific thermal energy.
  auto specificThermalEnergy = dataBase.fluidSpecificThermalEnergy();
  if (compatibleEnergy) {
    state.enroll(specificThermalEnergy, make_policy<NonSymmetricSpecificThermalEnergyPolicy<Dim<1>>>(dataBase));

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
SolidSphericalSPH::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {
  SolidSPH<Dimension>::registerDerivatives(dataBase, derivs);
  const auto compatibleEnergy = this->compatibleEnergyEvolution();
  if (compatibleEnergy) {
    const auto& connectivityMap = dataBase.connectivityMap();
    mPairAccelerationsPtr = std::make_unique<PairAccelerationsType>(connectivityMap);
    dataBase.resizeFluidFieldList(mSelfAccelerations, Vector::zero(), HydroFieldNames::selfAccelerations, false);
    derivs.enroll(HydroFieldNames::pairAccelerations, *mPairAccelerationsPtr);
    derivs.enroll(mSelfAccelerations);
  }
}

//------------------------------------------------------------------------------
// Stuff that occurs the beginning of a timestep
//------------------------------------------------------------------------------
void
SolidSphericalSPH::
preStepInitialize(const DataBase<Dimension>& dataBase, 
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {

  // We have to override the density summation to use the correct kernel
  switch(densityUpdate()) {

  case MassDensityType::IntegrateDensity:
    break;

  case MassDensityType::RigorousSumDensity:
  case MassDensityType::CorrectedSumDensity:
    {
      const auto& connectivityMap = dataBase.connectivityMap();
      const auto  position = state.fields(HydroFieldNames::position, Vector::zero());
      const auto  mass = state.fields(HydroFieldNames::mass, 0.0);
      const auto  H = state.fields(HydroFieldNames::H, SymTensor::zero());
      auto        massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
      computeSPHSumMassDensity(connectivityMap, this->kernel(), mSumMassDensityOverAllNodeLists, position, mass, H, massDensity);
      for (auto boundaryItr = this->boundaryBegin(); boundaryItr < this->boundaryEnd(); ++boundaryItr) (*boundaryItr)->applyFieldListGhostBoundary(massDensity);
      for (auto boundaryItr = this->boundaryBegin(); boundaryItr < this->boundaryEnd(); ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();
      if (densityUpdate() == MassDensityType::CorrectedSumDensity) {
        correctSPHSumMassDensity(connectivityMap, this->kernel(), mSumMassDensityOverAllNodeLists, position, mass, H, massDensity);
        for (auto boundaryItr = this->boundaryBegin(); boundaryItr < this->boundaryEnd(); ++boundaryItr) (*boundaryItr)->applyFieldListGhostBoundary(massDensity);
        for (auto boundaryItr = this->boundaryBegin(); boundaryItr < this->boundaryEnd(); ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();
      }
    }
    break;

  case MassDensityType::SumDensity:
    {
      auto       massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
      const auto massDensitySum = derivs.fields(ReplaceState<Dimension, Scalar>::prefix() + 
                                                HydroFieldNames::massDensity, 0.0);
      massDensity.assignFields(massDensitySum);
      for (auto boundaryItr = this->boundaryBegin(); boundaryItr < this->boundaryEnd(); ++boundaryItr) (*boundaryItr)->applyFieldListGhostBoundary(massDensity);
      for (auto boundaryItr = this->boundaryBegin(); boundaryItr < this->boundaryEnd(); ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();
    }
    break;

  default:
    VERIFY2(false, "Unsupported mass density definition for Spherical SPH");
    break;
  }
}

//------------------------------------------------------------------------------
// Determine the principle derivatives.
//------------------------------------------------------------------------------
void
SolidSphericalSPH::
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
SolidSphericalSPH::
evaluateDerivativesImpl(const Dim<1>::Scalar /*time*/,
                        const Dim<1>::Scalar dt,
                        const DataBase<Dim<1>>& dataBase,
                        const State<Dim<1>>& state,
                        StateDerivatives<Dim<1>>& derivs,
                        const QType& Q) const {

  TIME_BEGIN("SolidSphericalSPHevalDerivs");
  TIME_BEGIN("SolidSphericalSPHevalDerivs_initial");

  using QPiType = typename QType::ReturnType;

  // The kernels and such.
  const auto& W = this->kernel();
  const auto& WQ = this->PiKernel();
  const auto& WG = this->GradKernel();
  const auto& W1d = W.baseKernel1d();
  const auto  oneKernelQ = (W == WQ);
  const auto  oneKernelG = (W == WG);
  const auto  etaMax = W.etamax();
  const auto  W0 = W1d(0.0, 1.0);

  // A few useful constants we'll use in the following loop.
  const auto tiny = 1.0e-30;
  const auto epsTensile = this->epsilonTensile();
  const auto compatibleEnergy = this->compatibleEnergyEvolution();
  const auto evolveTotalEnergy = this->evolveTotalEnergy();
  const auto XSPH = this->XSPH();
  const auto correctVelocityGradient = this->correctVelocityGradient();

  // The connectivity.
  const auto& connectivityMap = dataBase.connectivityMap();
  const auto& nodeLists = connectivityMap.nodeLists();
  const auto numNodeLists = nodeLists.size();

  // The set of interacting node pairs.
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();

  // Get the state and derivative FieldLists.
  // State FieldLists.
  const auto mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto position = state.fields(HydroFieldNames::position, Vector::zero());
  const auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero());
  const auto massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  const auto H = state.fields(HydroFieldNames::H, SymTensor::zero());
  const auto pressure = state.fields(HydroFieldNames::pressure, 0.0);
  const auto soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  const auto omega = state.fields(HydroFieldNames::omegaGradh, 0.0);
  const auto S = state.fields(SolidFieldNames::deviatoricStress, SymTensor::zero());
  const auto mu = state.fields(SolidFieldNames::shearModulus, 0.0);
  const auto damage = state.fields(SolidFieldNames::tensorDamage, SymTensor::zero());
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
  CHECK(omega.size() == numNodeLists);
  CHECK(S.size() == numNodeLists);
  CHECK(mu.size() == numNodeLists);
  CHECK(damage.size() == numNodeLists);
  CHECK(fClQ.size() == 0 or fClQ.size() == numNodeLists);
  CHECK(fCqQ.size() == 0 or fCqQ.size() == numNodeLists);
  CHECK(DvDxQ.size() == 0 or DvDxQ.size() == numNodeLists);

  // Derivative FieldLists.
  auto  rhoSum = derivs.fields(ReplaceState<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, 0.0);
  auto  DxDt = derivs.fields(IncrementState<Dimension, Vector>::prefix() + HydroFieldNames::position, Vector::zero());
  auto  DrhoDt = derivs.fields(IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, 0.0);
  auto  DvDt = derivs.fields(HydroFieldNames::hydroAcceleration, Vector::zero());
  auto  DepsDt = derivs.fields(IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
  auto  DvDx = derivs.fields(HydroFieldNames::velocityGradient, Tensor::zero());
  auto  localDvDx = derivs.fields(HydroFieldNames::internalVelocityGradient, Tensor::zero());
  auto  M = derivs.fields(HydroFieldNames::M_SPHCorrection, Tensor::zero());
  auto  localM = derivs.fields("local " + HydroFieldNames::M_SPHCorrection, Tensor::zero());
  auto  maxViscousPressure = derivs.fields(HydroFieldNames::maxViscousPressure, 0.0);
  auto  effViscousPressure = derivs.fields(HydroFieldNames::effectiveViscousPressure, 0.0);
  auto  rhoSumCorrection = derivs.fields(HydroFieldNames::massDensityCorrection, 0.0);
  auto* pairAccelerationsPtr = derivs.template getPtr<PairAccelerationsType>(HydroFieldNames::pairAccelerations);
  auto  selfAccelerations = derivs.fields(HydroFieldNames::selfAccelerations, Vector::zero(), true);
  auto  XSPHWeightSum = derivs.fields(HydroFieldNames::XSPHWeightSum, 0.0);
  auto  XSPHDeltaV = derivs.fields(HydroFieldNames::XSPHDeltaV, Vector::zero());
  auto  DSDt = derivs.fields(IncrementState<Dimension, SymTensor>::prefix() + SolidFieldNames::deviatoricStress, SymTensor::zero());
  CHECK(rhoSum.size() == numNodeLists);
  CHECK(DxDt.size() == numNodeLists);
  CHECK(DrhoDt.size() == numNodeLists);
  CHECK(DvDt.size() == numNodeLists);
  CHECK(DepsDt.size() == numNodeLists);
  CHECK(DvDx.size() == numNodeLists);
  CHECK(localDvDx.size() == numNodeLists);
  CHECK(M.size() == numNodeLists);
  CHECK(localM.size() == numNodeLists);
  CHECK(maxViscousPressure.size() == numNodeLists);
  CHECK(effViscousPressure.size() == numNodeLists);
  CHECK(rhoSumCorrection.size() == numNodeLists);
  CHECK(XSPHWeightSum.size() == numNodeLists);
  CHECK(XSPHDeltaV.size() == numNodeLists);
  CHECK(DSDt.size() == numNodeLists);
  CHECK((compatibleEnergy and pairAccelerationsPtr->size() == npairs) or not compatibleEnergy);
  CHECK((compatibleEnergy     and selfAccelerations.size() == numNodeLists) or
        (not compatibleEnergy and selfAccelerations.size() == 0u));

  // The scale for the tensile correction.
  const auto& nodeList = mass[0]->nodeList();
  const auto  nPerh = nodeList.nodesPerSmoothingScale();
  const auto  WnPerh = W1d(1.0/nPerh, 1.0);
  TIME_END("SolidSphericalSPHevalDerivs_initial");

  // Walk all the interacting pairs.
  TIME_BEGIN("SolidSphericalSPHevalDerivs_pairs");
#pragma omp parallel
  {
    // Thread private  scratch variables.
    int i, j, nodeListi, nodeListj;
    Vector gradWii, gradWij, gradWQii, gradWQij, gradWGii, gradWGij;
    Vector gradWji, gradWjj, gradWQji, gradWQjj, gradWGji, gradWGjj;
    Scalar Wii, gWii, WQii, gWQii, Wji, gWji, WQji, gWQji, WGii, gWGii, WGji, gWGji;
    Scalar Wij, gWij, WQij, gWQij, Wjj, gWjj, WQjj, gWQjj, WGij, gWGij, WGjj, gWGjj;
    Scalar Qi, Qj;
    QPiType QPiij, QPiji;
    SymTensor sigmai, sigmaj, sigmarhoi, sigmarhoj;

    typename SpheralThreads<Dimension>::FieldListStack threadStack;
    auto rhoSum_thread = rhoSum.threadCopy(threadStack);
    auto DvDt_thread = DvDt.threadCopy(threadStack);
    auto DepsDt_thread = DepsDt.threadCopy(threadStack);
    auto DvDx_thread = DvDx.threadCopy(threadStack);
    auto localDvDx_thread = localDvDx.threadCopy(threadStack);
    auto M_thread = M.threadCopy(threadStack);
    auto localM_thread = localM.threadCopy(threadStack);
    auto maxViscousPressure_thread = maxViscousPressure.threadCopy(threadStack, ThreadReduction::MAX);
    auto effViscousPressure_thread = effViscousPressure.threadCopy(threadStack);
    auto rhoSumCorrection_thread = rhoSumCorrection.threadCopy(threadStack);
    auto XSPHWeightSum_thread = XSPHWeightSum.threadCopy(threadStack);
    auto XSPHDeltaV_thread = XSPHDeltaV.threadCopy(threadStack);
    auto DSDt_thread = DSDt.threadCopy(threadStack);

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
      const auto  Pi = pressure(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto  ci = soundSpeed(nodeListi, i);
      const auto& Si = S(nodeListi, i);
      const auto  safeOmegai = 1.0; // safeInv(omega(nodeListi, i), tiny);
      CHECK(mi > 0.0);
      CHECK(rhoi > 0.0);

      auto& rhoSumi = rhoSum_thread(nodeListi, i);
      auto& DvDti = DvDt_thread(nodeListi, i);
      auto& DepsDti = DepsDt_thread(nodeListi, i);
      auto& DvDxi = DvDx_thread(nodeListi, i);
      auto& localDvDxi = localDvDx_thread(nodeListi, i);
      auto& Mi = M_thread(nodeListi, i);
      auto& localMi = localM_thread(nodeListi, i);
      auto& maxViscousPressurei = maxViscousPressure_thread(nodeListi, i);
      auto& effViscousPressurei = effViscousPressure_thread(nodeListi, i);
      auto& rhoSumCorrectioni = rhoSumCorrection_thread(nodeListi, i);
      auto& XSPHWeightSumi = XSPHWeightSum_thread(nodeListi, i);
      auto& XSPHDeltaVi = XSPHDeltaV_thread(nodeListi, i);

      // Get the state for node j
      const auto& rj = position(nodeListj, j);
      const auto  mj = mass(nodeListj, j);
      const auto& vj = velocity(nodeListj, j);
      const auto  rhoj = massDensity(nodeListj, j);
      const auto  Pj = pressure(nodeListj, j);
      const auto& Hj = H(nodeListj, j);
      const auto  cj = soundSpeed(nodeListj, j);
      const auto& Sj = S(nodeListj, j);
      const auto  safeOmegaj = 1.0; // safeInv(omega(nodeListj, j), tiny);
      CHECK(mj > 0.0);
      CHECK(rhoj > 0.0);

      auto& rhoSumj = rhoSum_thread(nodeListj, j);
      auto& DvDtj = DvDt_thread(nodeListj, j);
      auto& DepsDtj = DepsDt_thread(nodeListj, j);
      auto& DvDxj = DvDx_thread(nodeListj, j);
      auto& localDvDxj = localDvDx_thread(nodeListj, j);
      auto& Mj = M_thread(nodeListj, j);
      auto& localMj = localM_thread(nodeListj, j);
      auto& maxViscousPressurej = maxViscousPressure_thread(nodeListj, j);
      auto& effViscousPressurej = effViscousPressure_thread(nodeListj, j);
      auto& rhoSumCorrectionj = rhoSumCorrection_thread(nodeListj, j);
      auto& XSPHWeightSumj = XSPHWeightSum_thread(nodeListj, j);
      auto& XSPHDeltaVj = XSPHDeltaV_thread(nodeListj, j);

      // Flag if this is a contiguous material pair or not.
      const auto sameMatij = true; // (nodeListi == nodeListj and fragIDi == fragIDj);
      const auto rij = ri - rj;

      // Determine how we're applying damage.
      const auto fDij = pairs[kk].f_couple;

      // Normalized node coordinates
      //   first subscript -> node position
      //  second subscript -> smoothing scale
      const auto etaii = Hi*ri;
      const auto etaji = Hi*rj;
      const auto etaij = Hj*ri;
      const auto etajj = Hj*rj;

      // Symmetrized kernel weight and gradient.
      // Note our subscript rules here:
      //   first subscript -> first node position in call: i -> (ri,rj) ; j -> (rj,ri)
      //  second subscript -> smoothing scale used
      W.kernelAndGrad(etaji, etaii, Hi, Wji, gradWji, gWji);
      W.kernelAndGrad(etajj, etaij, Hj, Wjj, gradWjj, gWjj);
      W.kernelAndGrad(etaii, etaji, Hi, Wii, gradWii, gWii);
      W.kernelAndGrad(etaij, etajj, Hj, Wij, gradWij, gWij);
      auto Wlookup = [](const bool copyVals, const SphericalKernel& W,
                        const Vector& etaj, const Vector& etai, const SymTensor& H,
                        const Scalar& Wj0, const Vector& gradWj0, const Scalar& gWj0,
                        Scalar& Wj, Vector& gradWj, Scalar& gWj) {
        if (copyVals) {
          Wj = Wj0;
          gWj = gWj0;
          gradWj = gradWj0;
        } else {
          W.kernelAndGrad(etaj, etai, H, Wj, gradWj, gWj);
        }
      };
      Wlookup(oneKernelQ, WQ,  etaji, etaii, Hi, Wji, gradWji, gWji, WQji, gradWQji, gWQji);
      Wlookup(oneKernelQ, WQ,  etajj, etaij, Hj, Wjj, gradWjj, gWjj, WQjj, gradWQjj, gWQjj);
      Wlookup(oneKernelQ, WQ,  etaii, etaji, Hi, Wii, gradWii, gWii, WQii, gradWQii, gWQii);
      Wlookup(oneKernelQ, WQ,  etaij, etajj, Hj, Wij, gradWij, gWij, WQij, gradWQij, gWQij);
      Wlookup(oneKernelG, WG,  etaji, etaii, Hi, Wji, gradWji, gWji, WGji, gradWGji, gWGji);
      Wlookup(oneKernelG, WG,  etajj, etaij, Hj, Wjj, gradWjj, gWjj, WGjj, gradWGjj, gWGjj);
      Wlookup(oneKernelG, WG,  etaii, etaji, Hi, Wii, gradWii, gWii, WGii, gradWGii, gWGii);
      Wlookup(oneKernelG, WG,  etaij, etajj, Hj, Wij, gradWij, gWij, WGij, gradWGij, gWGij);

      // Contribution to the sum density (only if the same material).
      if (nodeListi == nodeListj) {
        rhoSumi += mj*Wji;
        rhoSumj += mi*Wij;
      }

      // Contribution to the sum density correction
      rhoSumCorrectioni += mj * Wji / rhoj ;
      rhoSumCorrectionj += mi * Wij / rhoi ;

      // Compute the pair-wise artificial viscosity.
      const auto vij = vi - vj;
      Q.QPiij(QPiij, QPiji, Qi, Qj,
              nodeListi, i, nodeListj, j,
              ri, Hi, etaii - etaji, vi, rhoi, ci,  
              rj, Hj, etaij - etajj, vj, rhoj, cj,
              fClQ, fCqQ, DvDxQ);
      maxViscousPressurei = max(maxViscousPressurei, Qi);
      maxViscousPressurej = max(maxViscousPressurej, Qj);
      effViscousPressurei += mj*Qi*WQii/rhoj;
      effViscousPressurej += mi*Qj*WQij/rhoi;

      // Compute the stress tensors.
      if (sameMatij) {
        sigmai = fDij*Si - Pi * SymTensor::one();
        sigmaj = fDij*Sj - Pj * SymTensor::one();
      } else {
        sigmai = -Pi * SymTensor::one();
        sigmaj = -Pj * SymTensor::one();
      }

      // Compute the tensile correction to add to the stress as described in 
      // Gray, Monaghan, & Swift (Comput. Methods Appl. Mech. Eng., 190, 2001)
      if (epsTensile > 0.0) {
        const auto Wi = W1d((etaji - etaii).magnitude(), 1.0);
        const auto Wj = W1d((etajj - etaij).magnitude(), 1.0);
        const auto Ri = epsTensile*FastMath::pow4(Wi/(WnPerh))*tensileStressCorrection(sigmai);
        const auto Rj = epsTensile*FastMath::pow4(Wj/(WnPerh))*tensileStressCorrection(sigmaj);
        sigmai += Ri;
        sigmaj += Rj;
      }

      // Acceleration.
      CHECK(rhoi > 0.0);
      CHECK(rhoj > 0.0);
      sigmarhoi = safeOmegai*sigmai/(rhoi*rhoi);
      sigmarhoj = safeOmegaj*sigmaj/(rhoj*rhoj);
      const auto fQi = rj.x()*safeInv(ri.x() + rj.x());
      const auto deltaDvDti = mj*(sigmarhoi*gradWji + sigmarhoj*gradWjj - fQi        *(QPiij*gradWQji + QPiji*gradWQjj));
      const auto deltaDvDtj = mi*(sigmarhoj*gradWij + sigmarhoi*gradWii - (1.0 - fQi)*(QPiji*gradWQij + QPiij*gradWQii));
      DvDti += deltaDvDti;
      DvDtj += deltaDvDtj;
      if (compatibleEnergy) {
        (*pairAccelerationsPtr)[kk][0] = deltaDvDti;
        (*pairAccelerationsPtr)[kk][1] = deltaDvDtj;
      }

      // Specific thermal energy evolution.
      const auto QPisum = scalarSum(QPiij, QPiji);
      DepsDti -= mj*(sigmarhoi.xx() - 0.25*QPisum)*(vj.dot(gradWij) + vi.dot(gradWjj));
      DepsDtj -= mi*(sigmarhoj.xx() - 0.25*QPisum)*(vi.dot(gradWji) + vj.dot(gradWii));

      // Velocity gradient.
      const auto deltaDvDxi = -fDij * mj*vij.dyad(gradWGjj);
      const auto deltaDvDxj =  fDij * mi*vij.dyad(gradWGii);
      DvDxi += deltaDvDxi;
      DvDxj += deltaDvDxj;
      if (sameMatij) {
        localDvDxi += deltaDvDxi;
        localDvDxj += deltaDvDxj;
      }

      // Estimate of delta v (for XSPH).
      if (XSPH and sameMatij) {
        const auto wXSPHij = 0.5*(mi/rhoi*Wji + mj/rhoj*Wij);
        XSPHWeightSumi += wXSPHij;
        XSPHWeightSumj += wXSPHij;
        XSPHDeltaVi -= wXSPHij*vij;
        XSPHDeltaVj += wXSPHij*vij;
      }

      // Linear gradient correction term.
      const auto deltaMi = -mj*rij.dyad(gradWGjj);
      const auto deltaMj =  mi*rij.dyad(gradWGii);
      Mi += deltaMi;
      Mj += deltaMj;
      if (sameMatij) {
        localMi += deltaMi;
        localMj += deltaMj;
      }

    } // loop over pairs

    // Reduce the thread values to the master.
    threadReduceFieldLists<Dimension>(threadStack);

  }   // OpenMP parallel region
  TIME_END("SolidSphericalSPHevalDerivs_pairs");

  // Finish up the derivatives for each point.
  TIME_BEGIN("SolidSphericalSPHevalDerivs_final");
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto& nodeList = mass[nodeListi]->nodeList();
    const auto ni = nodeList.numInternalNodes();

    // Check if we can identify a reference density.
    auto rho0 = 0.0;
    try {
      rho0 = dynamic_cast<const SolidEquationOfState<Dimension>&>(dynamic_cast<const FluidNodeList<Dimension>&>(nodeList).equationOfState()).referenceDensity();
      // cerr << "Setting reference density to " << rho0 << endl;
    } catch(...) {
      // cerr << "BLAGO!" << endl;
    }

#pragma omp parallel for
    for (auto i = 0u; i < ni; ++i) {

      // Get the state for node i.
      const auto& ri = position(nodeListi, i);
      const auto& mi = mass(nodeListi, i);
      const auto& vi = velocity(nodeListi, i);
      const auto& rhoi = massDensity(nodeListi, i);
      const auto& Pi = pressure(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto& Si = S(nodeListi, i);
      const auto& mui = mu(nodeListi, i);
      const auto& ci = soundSpeed(nodeListi, i);
      const auto  safeOmegai = 1.0; // safeInv(omega(nodeListi, i), tiny);
      const auto  Hdeti = Hi.Determinant();
      const auto  hi = 1.0/Hdeti;
      const auto  numNeighborsi = connectivityMap.numNeighborsForNode(nodeListi, i);
      const auto  riInv = safeInvVar(ri.x(), 0.01*hi);
      CHECK(mi > 0.0);
      CHECK(rhoi > 0.0);
      CHECK(Hdeti > 0.0);

      auto& rhoSumi = rhoSum(nodeListi, i);
      auto& DxDti = DxDt(nodeListi, i);
      auto& DrhoDti = DrhoDt(nodeListi, i);
      auto& DvDti = DvDt(nodeListi, i);
      auto& DepsDti = DepsDt(nodeListi, i);
      auto& DvDxi = DvDx(nodeListi, i);
      auto& localDvDxi = localDvDx(nodeListi, i);
      auto& Mi = M(nodeListi, i);
      auto& localMi = localM(nodeListi, i);
      auto& maxViscousPressurei = maxViscousPressure(nodeListi, i);
      auto& effViscousPressurei = effViscousPressure(nodeListi, i);
      auto& rhoSumCorrectioni = rhoSumCorrection(nodeListi, i);
      auto& XSPHWeightSumi = XSPHWeightSum(nodeListi, i);
      auto& XSPHDeltaVi = XSPHDeltaV(nodeListi, i);
      auto& DSDti = DSDt(nodeListi, i);

      // Symmetrized kernel weight and gradient.
      const auto etaii = Hi*ri;
      const Vector etaQii(std::max(0.01, etaii[0]));
      double Wii, gWii;
      Vector gradWii, gradWQii, gradWGii;
      W.kernelAndGrad(etaii, etaii, Hi, Wii, gradWii, gWii);
      gradWQii = WQ.grad(etaQii, etaQii, Hi);
      gradWGii = oneKernelG ? gradWii : WG.grad(etaii, etaii, Hi);

      // Add the self-contribution to density sum.
      rhoSumi += mi*Wii;

      // Add the self-contribution to density sum correction.
      rhoSumCorrectioni += mi*Wii/rhoi ;

      // Correct the effective viscous pressure.
      effViscousPressurei /= rhoSumCorrectioni ;

      // Finish the gradient of the velocity.
      CHECK(rhoi > 0.0);
      if (correctVelocityGradient and
          std::abs(Mi.Determinant()) > 1.0e-10 and
          numNeighborsi > Dimension::pownu(2)) {
        Mi = Mi.Inverse();
        DvDxi = DvDxi*Mi;
      } else {
        DvDxi /= rhoi;
      }
      if (correctVelocityGradient and
          std::abs(localMi.Determinant()) > 1.0e-10 and
          numNeighborsi > Dimension::pownu(2)) {
        localMi = localMi.Inverse();
        localDvDxi = localDvDxi*localMi;
      } else {
        localDvDxi /= rhoi;
      }

      // Evaluate the continuity equation.
      const auto divv = DvDxi.xx() + 2.0*vi.x()*riInv;
      DrhoDti = -rhoi*divv;

      // If we're in range of the origin, compute an effective Q.
      Scalar Qi = 0.0;
      if (etaii.x() < etaMax and vi.x() < 0.0) {
        // const auto mu = -std::min(0.0, divv);  // * W1d(etaii.x(), 1.0)/W0;
        const auto mu = -std::min(0.0, vi.x()*riInv) * W1d(etaii.x(), 1.0)/W0;
        Qi = mQself*rhoi*hi*mu*(hi*mu + ci);
        maxViscousPressurei = std::max(maxViscousPressurei, Qi);
      }

      // Self-interaction for momentum (cause curvilinear coordinates are weird)
      const auto sigmai = Si - Pi * SymTensor::one();
      const auto deltaDvDti = mi*safeOmegai/(rhoi*rhoi)*(2.0*sigmai*gradWii - Qi*gradWQii) +  // self-interaction because kernel is not symmetric
                              Vector(3.0*Si.xx()/rhoi*riInv);                                 // hoop terms from theta & phi directions
      DvDti += deltaDvDti;
      if (compatibleEnergy) selfAccelerations(nodeListi, i) = deltaDvDti;

      // Specific thermal energy
      DepsDti -= 2.0*mi/(rhoi*rhoi)*(sigmai.xx() - 0.5*Qi)*vi.dot(gradWii) +
                 2.0*(0.5*Si.xx() + Pi)/rhoi*vi.x()*riInv;

      // If needed finish the total energy derivative.
      if (evolveTotalEnergy) DepsDti = mi*(vi.dot(DvDti) + DepsDti);

      // Determine the position evolution, based on whether we're doing XSPH or not.
      if (mXSPH) {
        XSPHWeightSumi += mi/rhoi*Wii;
        CHECK2(XSPHWeightSumi != 0.0, i << " " << XSPHWeightSumi);
        DxDti = vi + XSPHDeltaVi/max(tiny, XSPHWeightSumi);
      } else {
        DxDti = vi;
      }

      // Determine the deviatoric stress evolution.
      // Note the spin term is always zero in spherical coordinates.
      // const auto deviatoricDeformation = SymTensor(2.0/3.0*(localDvDxi.xx() + vi.x()*riInv));
      const auto deviatoricDeformation = localDvDxi.xx() - divv/3.0;
      DSDti.xx(2.0*mui*deviatoricDeformation);

      // Optionally use damage to ramp down stress on damaged material.
      const auto Di = max(0.0, min(1.0, damage(nodeListi, i).Trace() - 1.0));
      // Hideali = (1.0 - Di)*Hideali + Di*mHfield0(nodeListi, i);
      // DHDti = (1.0 - Di)*DHDti + Di*(mHfield0(nodeListi, i) - Hi)*0.25/dt;

      // We also adjust the density evolution in the presence of damage.
      if (rho0 > 0.0) DrhoDti = (1.0 - Di)*DrhoDti - 0.5/dt*Di*(rhoi - rho0);

      // // In the presence of damage, add a term to reduce the stress on this point.
      // DSDti = (1.0 - Di)*DSDti - 0.25/dt*Di*Si;
    }
  }
  TIME_END("SolidSphericalSPHevalDerivs_final");
  TIME_END("SolidSphericalSPHevalDerivs");
}

//------------------------------------------------------------------------------
// Apply the ghost boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
void
SolidSphericalSPH::
applyGhostBoundaries(State<Dim<1>>& state,
                     StateDerivatives<Dim<1>>& derivs) {
  // Convert the mass to mass/length^2 before BCs are applied.
  auto       mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto pos = state.fields(HydroFieldNames::position, Vector::zero());
  const auto numNodeLists = mass.numFields();
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto n = mass[nodeListi]->numElements();
    for (auto i = 0u; i < n; ++i) {
      const auto dA = FastMath::pow2(pos(nodeListi, i).x());
      CHECK(dA > 0.0);
      mass(nodeListi, i) /= dA;
    }
  }

  // Apply ordinary SPH BCs.
  SolidSPH<Dim<1>>::applyGhostBoundaries(state, derivs);
  for (auto boundaryPtr: range(this->boundaryBegin(), this->boundaryEnd())) boundaryPtr->finalizeGhostBoundary();

  // Scale back to mass.
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto n = mass[nodeListi]->numElements();
    for (auto i = 0u; i < n; ++i) {
      const auto dA = FastMath::pow2(pos(nodeListi, i).x());
      CHECK(dA > 0.0);
      mass(nodeListi, i) *= dA;
    }
  }
}

//------------------------------------------------------------------------------
// Enforce the boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
void
SolidSphericalSPH::
enforceBoundaries(State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {
  // Convert the mass to mass/length^2 before BCs are applied.
  auto       mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto pos = state.fields(HydroFieldNames::position, Vector::zero());
  const auto numNodeLists = mass.numFields();
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto n = mass[nodeListi]->numElements();
    for (auto i = 0u; i < n; ++i) {
      const auto dA = FastMath::pow2(pos(nodeListi, i).x());
      CHECK(dA > 0.0);
      mass(nodeListi, i) /= dA;
    }
  }

  // Apply ordinary SPH BCs.
  SolidSPH<Dim<1>>::applyGhostBoundaries(state, derivs);
  for (auto boundaryPtr: range(this->boundaryBegin(), this->boundaryEnd())) boundaryPtr->finalizeGhostBoundary();

  // Scale back to mass.
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto n = mass[nodeListi]->numElements();
    for (auto i = 0u; i < n; ++i) {
      const auto dA = FastMath::pow2(pos(nodeListi, i).x());
      CHECK(dA > 0.0);
      mass(nodeListi, i) *= dA;
    }
  }
}

}
