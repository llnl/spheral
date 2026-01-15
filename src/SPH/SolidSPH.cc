//---------------------------------Spheral++----------------------------------//
// SolidSPHHydroBase -- The SPH/ASPH solid material hydrodynamic package for Spheral++.
//
// Created by JMO, Fri Jul 30 11:07:33 PDT 2010
//----------------------------------------------------------------------------//
#include "SPH/SolidSPH.hh"
#include "FileIO/FileIO.hh"
#include "Utilities/NodeCoupling.hh"
#include "SPH/SPH.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Strength/SolidFieldNames.hh"
#include "NodeList/SolidNodeList.hh"
#include "Strength/DeviatoricStressPolicy.hh"
#include "Strength/BulkModulusPolicy.hh"
#include "Strength/PlasticStrainPolicy.hh"
#include "Strength/ShearModulusPolicy.hh"
#include "Strength/YieldStrengthPolicy.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/ReplaceState.hh"
#include "DataBase/ReplaceBoundedState.hh"
#include "DataBase/updateStateFields.hh"
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
#include <map>
#include <vector>
#include <optional>
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
SPHERAL_HOST_DEVICE
inline
Dim<1>::SymTensor
tensileStressCorrection(const Dim<1>::SymTensor& sigma) {
  if (sigma.xx() > 0.0) {
    return -sigma;
  } else {
    return Dim<1>::SymTensor::zero();
  }
}

SPHERAL_HOST_DEVICE
inline
Dim<2>::SymTensor
tensileStressCorrection(const Dim<2>::SymTensor& sigma) {
  const EigenStruct<2> eigen = sigma.eigenVectors();
  const double lambdax = eigen.eigenValues.x();
  const double lambday = eigen.eigenValues.y();
  Dim<2>::SymTensor result((lambdax > 0.0 ? -lambdax : 0.0), 0.0,
                           0.0,                              (lambday > 0.0 ? -lambday : 0.0));
  result.rotationalTransform(eigen.eigenVectors);
  return result;
}

SPHERAL_HOST_DEVICE
inline
Dim<3>::SymTensor
tensileStressCorrection(const Dim<3>::SymTensor& sigma) {
  const EigenStruct<3> eigen = sigma.eigenVectors();
  const double lambdax = eigen.eigenValues.x();
  const double lambday = eigen.eigenValues.y();
  const double lambdaz = eigen.eigenValues.z();
  Dim<3>::SymTensor result((lambdax > 0.0 ? -lambdax : 0.0), 0.0,                              0.0,
                           0.0,                              (lambday > 0.0 ? -lambday : 0.0), 0.0,
                           0.0,                              0.0,                              (lambdaz > 0.0 ? -lambdaz : 0.0));
  result.rotationalTransform(eigen.eigenVectors);
  return result;
}

//------------------------------------------------------------------------------
// Compute one minus the SymTensor in it's principle frame
//------------------------------------------------------------------------------
inline Dim<1>::SymTensor oneMinusEigenvalues(const Dim<1>::SymTensor& x) {
  return Dim<1>::SymTensor(1.0 - x[0]);
}

inline Dim<2>::SymTensor oneMinusEigenvalues(const Dim<2>::SymTensor& x) {
  const auto eigen = x.eigenVectors();
  Dim<2>::SymTensor result(1.0 - eigen.eigenValues[0], 0.0,
                           0.0, 1.0 - eigen.eigenValues[1]);
  result.rotationalTransform(eigen.eigenVectors);
  return result;
}

inline Dim<3>::SymTensor oneMinusEigenvalues(const Dim<3>::SymTensor& x) {
  const auto eigen = x.eigenVectors();
  Dim<3>::SymTensor result(1.0 - eigen.eigenValues[0], 0.0, 0.0,
                           0.0, 1.0 - eigen.eigenValues[1], 0.0,
                           0.0, 0.0, 1.0 - eigen.eigenValues[2]);
  result.rotationalTransform(eigen.eigenVectors);
  return result;
}

}

//------------------------------------------------------------------------------
// Construct with the given artificial viscosity and kernels.
//------------------------------------------------------------------------------
template<typename Dimension>
SolidSPH<Dimension>::
SolidSPH(DataBase<Dimension>& dataBase,
         ArtificialViscosity<Dimension>& Q,
         const TableKernel<Dimension>& W,
         const TableKernel<Dimension>& WPi,
         const TableKernel<Dimension>& WGrad,
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
  SPH<Dimension>(dataBase,
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
  mDamageRelieveRubble(damageRelieveRubble),
  mStrengthInDamage(strengthInDamage),
  mGradKernel(WGrad),
  mDdeviatoricStressDt(FieldStorageType::CopyFields),
  mBulkModulus(FieldStorageType::CopyFields),
  mShearModulus(FieldStorageType::CopyFields),
  mYieldStrength(FieldStorageType::CopyFields),
  mPlasticStrain0(FieldStorageType::CopyFields) {

  // Create storage for the state we're holding.
  mDdeviatoricStressDt = dataBase.newSolidFieldList(SymTensor::zero(), IncrementState<Dimension, Vector>::prefix() + SolidFieldNames::deviatoricStress);
  mBulkModulus = dataBase.newSolidFieldList(0.0, SolidFieldNames::bulkModulus);
  mShearModulus = dataBase.newSolidFieldList(0.0, SolidFieldNames::shearModulus);
  mYieldStrength = dataBase.newSolidFieldList(0.0, SolidFieldNames::yieldStrength);
  mPlasticStrain0 = dataBase.newSolidFieldList(0.0, SolidFieldNames::plasticStrain + "0");
}

//------------------------------------------------------------------------------
// On problem start up, we need to initialize our internal data.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidSPH<Dimension>::
initializeProblemStartupDependencies(DataBase<Dimension>& dataBase,
                                     State<Dimension>& state,
                                     StateDerivatives<Dimension>& derivs) {
  TIME_BEGIN("SolidSPHinitializeStartup");

  // Call the ancestor.
  SPHBase<Dimension>::initializeProblemStartupDependencies(dataBase, state, derivs);

  // Set the moduli.
  updateStateFields(SolidFieldNames::bulkModulus, state, derivs);
  updateStateFields(SolidFieldNames::shearModulus, state, derivs);
  updateStateFields(SolidFieldNames::yieldStrength, state, derivs);

  TIME_END("SolidSPHinitializeStartup");
}


//------------------------------------------------------------------------------
// Register the state we need/are going to evolve.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidSPH<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {
  TIME_BEGIN("SolidSPHregister");

  // Invoke SPHHydro's state.
  SPH<Dimension>::registerState(dataBase, state);

  // Create the local storage.
  dataBase.resizeFluidFieldList(mBulkModulus, 0.0, SolidFieldNames::bulkModulus, false);
  dataBase.resizeFluidFieldList(mShearModulus, 0.0, SolidFieldNames::shearModulus, false);
  dataBase.resizeFluidFieldList(mYieldStrength, 0.0, SolidFieldNames::yieldStrength, false);
  dataBase.resizeFluidFieldList(mPlasticStrain0, 0.0, SolidFieldNames::plasticStrain + "0", false);

  // Register the deviatoric stress and plastic strain to be evolved.
  auto S = dataBase.solidDeviatoricStress();
  auto ps = dataBase.solidPlasticStrain();
  state.enroll(S, make_policy<DeviatoricStressPolicy<Dimension>>());
  state.enroll(ps, make_policy<PlasticStrainPolicy<Dimension>>());

  // Register the bulk modulus, shear modulus, and yield strength.
  state.enroll(mBulkModulus, make_policy<BulkModulusPolicy<Dimension>>());
  state.enroll(mShearModulus, make_policy<ShearModulusPolicy<Dimension>>());
  state.enroll(mYieldStrength, make_policy<YieldStrengthPolicy<Dimension>>());

  // Register the damage with a default no-op update.
  // If there are any damage models running they can override this choice.
  auto D = dataBase.solidDamage();
  state.enroll(D);

  // Register the fragment IDs.
  auto fragIDs = dataBase.solidFragmentIDs();
  state.enroll(fragIDs);

  // Register the particle types.
  auto pTypes = dataBase.solidParticleTypes();
  state.enroll(pTypes);

  // And finally the intial plastic strain.
  mPlasticStrain0.assignFields(ps);
  state.enroll(mPlasticStrain0);
  TIME_END("SolidSPHregister");
}

//------------------------------------------------------------------------------
// Register the state derivative fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidSPH<Dimension>::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {
  TIME_BEGIN("SolidSPHregisterDerivs");

  // Call the ancestor method.
  SPH<Dimension>::registerDerivatives(dataBase, derivs);

  // Create the scratch fields.
  // Note we deliberately do not zero out the derivatives here!  This is because the previous step
  // info here may be used by other algorithms (like the CheapSynchronousRK2 integrator or
  // the ArtificialVisocisity::initialize step).
  const auto DSDtName = IncrementState<Dimension, Vector>::prefix() + SolidFieldNames::deviatoricStress;
  dataBase.resizeFluidFieldList(mDdeviatoricStressDt, SymTensor::zero(), DSDtName, false);

  derivs.enroll(mDdeviatoricStressDt);
  for (auto [nodeListi, solidNodeListPtr]: enumerate(dataBase.solidNodeListBegin(), dataBase.solidNodeListEnd())) {
    derivs.enroll(solidNodeListPtr->plasticStrainRate());
  }
  TIME_END("SolidSPHregisterDerivs");
}

//------------------------------------------------------------------------------
// Determine the principle derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidSPH<Dimension>::
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
// Determine the principle derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename QType>
void
SolidSPH<Dimension>::
evaluateDerivativesImpl(const typename Dimension::Scalar /*time*/,
                        const typename Dimension::Scalar dt,
                        const DataBase<Dimension>& dataBase,
                        const State<Dimension>& state,
                        StateDerivatives<Dimension>& derivs,
                        chai::managed_ptr<QType> Q) const {

  TIME_BEGIN("SolidSPHevalDerivs");
  TIME_BEGIN("SolidSPHevalDerivs_initial");

  using QPiType = typename QType::ReturnType;

  // The kernels and such.
  const auto& W = this->kernel();
  const auto& WQ = this->PiKernel();
  const auto& WG = this->GradKernel();
  auto W_view = W.view();
  auto WQ_view = WQ.view();
  auto WG_view = WG.view();
  const auto  oneKernelQ = (W == WQ);
  const auto  oneKernelG = (W == WG);

  // A few useful constants we'll use in the following loop.
  const auto tiny = 1.0e-30;
  const auto W0 = W(0.0, 1.0);
  const auto WQ0 = WQ(0.0, 1.0);
  const auto epsTensile = this->epsilonTensile();
  //const auto compatibleEnergy = this->compatibleEnergyEvolution();
  const auto evolveTotalEnergy = this->evolveTotalEnergy();
  const auto XSPH = this->XSPH();

  // The connectivity.
  const auto& connectivityMap = dataBase.connectivityMap();
  const auto& nodeLists = connectivityMap.nodeLists();
  const auto numNodeLists = nodeLists.size();
  const auto& pairs_v = connectivityMap.nodePairList();
  const auto pairs = pairs_v.view();
  const auto npairs = pairs.size();
  // const auto& coupling = connectivityMap.coupling();

  // Get the state and derivative FieldLists.
  // State FieldLists.
  auto mass_v = state.fields(HydroFieldNames::mass, 0.0);
  auto position_v = state.fields(HydroFieldNames::position, Vector::zero());
  auto velocity_v = state.fields(HydroFieldNames::velocity, Vector::zero());
  auto massDensity_v = state.fields(HydroFieldNames::massDensity, 0.0);
  auto specificThermalEnergy_v = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  auto H_v = state.fields(HydroFieldNames::H, SymTensor::zero());
  auto pressure_v = state.fields(HydroFieldNames::pressure, 0.0);
  auto soundSpeed_v = state.fields(HydroFieldNames::soundSpeed, 0.0);
  auto S_v = state.fields(SolidFieldNames::deviatoricStress, SymTensor::zero());
  auto mu_v = state.fields(SolidFieldNames::shearModulus, 0.0);
  auto damage_v = state.fields(SolidFieldNames::tensorDamage, SymTensor::zero());
  auto pTypes_v = state.fields(SolidFieldNames::particleTypes, int(0));
  auto fragID_v = state.fields(SolidFieldNames::fragmentIDs, int(0));
  auto omega_v = state.fields(HydroFieldNames::omegaGradh, 0.0);
  auto mass = mass_v.view();
  auto position = position_v.view();
  auto velocity = velocity_v.view();
  auto massDensity = massDensity_v.view();
  auto H = H_v.view();
  auto pressure = pressure_v.view();
  auto soundSpeed = soundSpeed_v.view();
  auto omega = omega_v.view();
  auto S = S_v.view();
  auto mu = mu_v.view();
  auto specificThermalEnergy = specificThermalEnergy_v.view();
  auto damage = damage_v.view();
  auto pTypes = pTypes_v.view();
  auto fragID = fragID_v.view();
  auto fClQ = state.fields(HydroFieldNames::ArtificialViscousClMultiplier, 0.0, true);
  auto fCqQ = state.fields(HydroFieldNames::ArtificialViscousCqMultiplier, 0.0, true);
  auto DvDxQ = state.fields(HydroFieldNames::ArtificialViscosityVelocityGradient, Tensor::zero(), true);
  auto DvDxQView = DvDxQ.view();
  auto fClQView = fClQ.view();
  auto fCqQView = fCqQ.view();
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
  CHECK(pTypes.size() == numNodeLists);
  CHECK(fragID.size() == numNodeLists);
  CHECK(fClQ.size() == 0 or fClQ.size() == numNodeLists);
  CHECK(fCqQ.size() == 0 or fCqQ.size() == numNodeLists);
  CHECK(DvDxQ.size() == 0 or DvDxQ.size() == numNodeLists);

  // Derivative FieldLists.
  auto  rhoSum_v = derivs.fields(ReplaceState<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, 0.0);
  auto  DxDt_v = derivs.fields(IncrementState<Dimension, Vector>::prefix() + HydroFieldNames::position, Vector::zero());
  auto  DrhoDt_v = derivs.fields(IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, 0.0);
  auto  DvDt_v = derivs.fields(HydroFieldNames::hydroAcceleration, Vector::zero());
  auto  DepsDt_v = derivs.fields(IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
  auto  DvDx_v = derivs.fields(HydroFieldNames::velocityGradient, Tensor::zero());
  auto  localDvDx_v = derivs.fields(HydroFieldNames::internalVelocityGradient, Tensor::zero());
  auto  M_v = derivs.fields(HydroFieldNames::M_SPHCorrection, Tensor::zero());
  auto  localM_v = derivs.fields("local " + HydroFieldNames::M_SPHCorrection, Tensor::zero());
  auto  rhoSumCorrection_v = derivs.fields(HydroFieldNames::massDensityCorrection, 0.0);
  auto  maxViscousPressure_v = derivs.fields(HydroFieldNames::maxViscousPressure, 0.0);
  auto  effViscousPressure_v = derivs.fields(HydroFieldNames::effectiveViscousPressure, 0.0);
  auto  XSPHWeightSum_v = derivs.fields(HydroFieldNames::XSPHWeightSum, 0.0);
  auto  XSPHDeltaV_v = derivs.fields(HydroFieldNames::XSPHDeltaV, Vector::zero());
  auto  DSDt_v = derivs.fields(IncrementState<Dimension, SymTensor>::prefix() + SolidFieldNames::deviatoricStress, SymTensor::zero());
  auto rhoSum = rhoSum_v.view();
  auto DxDt = DxDt_v.view();
  auto DrhoDt = DrhoDt_v.view();
  auto DvDt = DvDt_v.view();
  auto DepsDt = DepsDt_v.view();
  auto DvDx = DvDx_v.view();
  auto localDvDx = localDvDx_v.view();
  auto M = M_v.view();
  auto localM = localM_v.view();
  auto maxViscousPressure = maxViscousPressure_v.view();
  auto effViscousPressure = effViscousPressure_v.view();
  auto XSPHWeightSum = XSPHWeightSum_v.view();
  auto XSPHDeltaV = XSPHDeltaV_v.view();
  auto DSDt = DSDt_v.view();
  auto rhoSumCorrection = rhoSumCorrection_v.view();
  CHECK(rhoSum_v.size() == numNodeLists);
  CHECK(DxDt_v.size() == numNodeLists);
  CHECK(DrhoDt_v.size() == numNodeLists);
  CHECK(DvDt_v.size() == numNodeLists);
  CHECK(DepsDt_v.size() == numNodeLists);
  CHECK(DvDx_v.size() == numNodeLists);
  CHECK(localDvDx_v.size() == numNodeLists);
  CHECK(M_v.size() == numNodeLists);
  CHECK(localM_v.size() == numNodeLists);
  CHECK(maxViscousPressure_v.size() == numNodeLists);
  CHECK(effViscousPressure_v.size() == numNodeLists);
  CHECK(XSPHWeightSum_v.size() == numNodeLists);
  CHECK(XSPHDeltaV_v.size() == numNodeLists);
  //auto* pairAccelerationsPtr = derivs.template getPtr<PairAccelerationsType>(HydroFieldNames::pairAccelerations);

  // The scale for the tensile correction.
  const auto& nodeList = mass_v[0]->nodeList();
  const auto  nPerh = nodeList.nodesPerSmoothingScale();
  const auto  WnPerh = W(1.0/nPerh, 1.0);
  bool CorrectVelocityGradient = this->mCorrectVelocityGradient;
  TIME_END("SolidSPHevalDerivs_initial");

  // Walk all the interacting pairs.
  TIME_BEGIN("SolidSPHevalDerivs_pairs");
  //#pragma omp parallel
  {
    // Thread private  scratch variables.
    // size_t i, j, nodeListi, nodeListj;
    // Scalar Wi, gWi, WQi, gWQi, Wj, gWj, WQj, gWQj;
    // Vector gradWi, gradWj, gradWQi, gradWQj, gradWGi, gradWGj;
    // Scalar Qi, Qj;
    // QPiType QPiij, QPiji;
    // SymTensor sigmai, sigmaj, sigmarhoi, sigmarhoj;

// #pragma omp for
//     for (auto kk = 0u; kk < npairs; ++kk) {
    RAJA::forall<EXEC_POLICY>(TRS_UINT(0u, npairs),
    [=] SPHERAL_HOST_DEVICE (size_t kk) {
      Scalar Wi, gWi, WQi, gWQi, Wj, gWj, WQj, gWQj;
      Vector gradWQi, gradWQj, gradWGi, gradWGj;
      Scalar Qi, Qj;
      //QPiType QPiij, QPiji;
      SymTensor sigmai, sigmaj, sigmarhoi, sigmarhoj;
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
      const auto& Si = S(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();
      const auto  safeOmegai = safeInv(omegai, tiny);
      const auto  pTypei = pTypes(nodeListi, i);
      const auto  fragIDi = fragID(nodeListi, i);
      CHECK(mi > 0.0);
      CHECK(rhoi > 0.0);
      CHECK(Hdeti > 0.0);

      auto& rhoSumi = rhoSum(nodeListi, i);
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

      // Get the state for node j
      const auto& rj = position(nodeListj, j);
      const auto  mj = mass(nodeListj, j);
      const auto& vj = velocity(nodeListj, j);
      const auto  rhoj = massDensity(nodeListj, j);
      const auto  Pj = pressure(nodeListj, j);
      const auto& Hj = H(nodeListj, j);
      const auto  cj = soundSpeed(nodeListj, j);
      const auto  omegaj = omega(nodeListj, j);
      const auto& Sj = S(nodeListj, j);
      const auto  Hdetj = Hj.Determinant();
      const auto  safeOmegaj = safeInv(omegaj, tiny);
      const auto  pTypej = pTypes(nodeListj, j);
      const auto  fragIDj = fragID(nodeListj, j);
      CHECK(mj > 0.0);
      CHECK(rhoj > 0.0);
      CHECK(Hdetj > 0.0);

      auto& rhoSumj = rhoSum(nodeListj, j);
      auto& DvDtj = DvDt(nodeListj, j);
      auto& DepsDtj = DepsDt(nodeListj, j);
      auto& DvDxj = DvDx(nodeListj, j);
      auto& localDvDxj = localDvDx(nodeListj, j);
      auto& Mj = M(nodeListj, j);
      auto& localMj = localM(nodeListj, j);
      auto& maxViscousPressurej = maxViscousPressure(nodeListj, j);
      auto& effViscousPressurej = effViscousPressure(nodeListj, j);
      auto& rhoSumCorrectionj = rhoSumCorrection(nodeListj, j);
      auto& XSPHWeightSumj = XSPHWeightSum(nodeListj, j);
      auto& XSPHDeltaVj = XSPHDeltaV(nodeListj, j);

      // Flag if this is a contiguous material pair or not.
      const auto sameMatij = true; // (nodeListi == nodeListj and fragIDi == fragIDj);

      // Flag if at least one particle is free (0).
      const auto freeParticle = (pTypei == 0 or pTypej == 0);

      // Flag to turn off forces between different fragments,
      // only if the two particles are moving away from each other.
      auto fragdir = true;
      if (fragIDi != fragIDj) {
        const auto rdiff = ri-rj;
        const auto vdiff = vi-vj;
        const auto vdot = vdiff.dot(rdiff);
        if (vdot > 0.0) {
          fragdir = false;
        }
      }

      // Determine how we're applying damage.
      const auto fDij = pairs[kk].f_couple;

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
      Vector gradWi = gWi*Hi*etaiUnit;
      Vector gradWj = gWj*Hj*etajUnit;
      if (oneKernelQ) {
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
      if (oneKernelG) {
        gradWGi = gradWi;
        gradWGj = gradWj;
      } else {
        gradWGi = Hi*etaiUnit * WG_view.gradValue(etaMagi, Hdeti);
        gradWGj = Hj*etajUnit * WG_view.gradValue(etaMagj, Hdetj);
      }

      // Contribution to the sum density (only if the same material).
      if (nodeListi == nodeListj) {
        RAJA::atomicAdd<RAJA::auto_atomic>(&rhoSumi, mj*Wi);
        RAJA::atomicAdd<RAJA::auto_atomic>(&rhoSumj, mi*Wj);
        // rhoSumi += mj*Wi;
        // rhoSumj += mi*Wj;
      }

      // Contribution to the sum density correction
      // rhoSumCorrectioni += mj*WQi / rhoj;
      // rhoSumCorrectionj += mi*WQj / rhoi;
      RAJA::atomicAdd<RAJA::auto_atomic>(&rhoSumCorrectioni, mj*WQi / rhoj);
      RAJA::atomicAdd<RAJA::auto_atomic>(&rhoSumCorrectionj, mi*WQj / rhoi);

      // Compute the pair-wise artificial viscosity.
      const auto vij = vi - vj;
      QPiType QPiij(0.0);
      QPiType QPiji(0.0);
      Qi = 0.0;
      Qj = 0.0;
      Q->QPiij(QPiij, QPiji, Qi, Qj,
               nodeListi, i, nodeListj, j,
               ri, Hi, etai, vi, rhoi, ci,  
               rj, Hj, etaj, vj, rhoj, cj,
               fClQView, fCqQView, DvDxQView);
      const auto Qacci = 0.5*(QPiij*gradWQi);
      const auto Qaccj = 0.5*(QPiji*gradWQj);
      const auto workQi = vij.dot(Qacci);
      const auto workQj = vij.dot(Qaccj);
      // maxViscousPressurei = std::max(maxViscousPressurei, Qi);
      // maxViscousPressurej = std::max(maxViscousPressurej, Qj);
      // effViscousPressurei += mj*Qi*WQi/rhoj;
      // effViscousPressurej += mi*Qj*WQj/rhoi;
      RAJA::atomicMax<RAJA::auto_atomic>(&maxViscousPressurei, Qi);
      RAJA::atomicMax<RAJA::auto_atomic>(&maxViscousPressurej, Qj);
      RAJA::atomicAdd<RAJA::auto_atomic>(&effViscousPressurei, mj*Qi*WQi/rhoj);
      RAJA::atomicAdd<RAJA::auto_atomic>(&effViscousPressurej, mi*Qj*WQj/rhoi);

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
      const auto Ri = epsTensile*FastMath::pow4(Wi/(Hdeti*WnPerh))*tensileStressCorrection(sigmai);
      const auto Rj = epsTensile*FastMath::pow4(Wj/(Hdetj*WnPerh))*tensileStressCorrection(sigmaj);
      sigmai += Ri;
      sigmaj += Rj;

      // Acceleration.
      CHECK(rhoi > 0.0);
      CHECK(rhoj > 0.0);
      sigmarhoi = safeOmegai*sigmai/(rhoi*rhoi);
      sigmarhoj = safeOmegaj*sigmaj/(rhoj*rhoj);
      const auto deltaDvDt = sigmarhoi*gradWi + sigmarhoj*gradWj - Qacci - Qaccj;
      if (freeParticle && fragdir) {
        DvDti.atomicAdd(mj*deltaDvDt);
        DvDtj.atomicSub(mi*deltaDvDt);
        // DvDti += mj*deltaDvDt;
        // DvDtj -= mi*deltaDvDt;
      }
      //if (compatibleEnergy) (*pairAccelerationsPtr)[kk] = mj*deltaDvDt;  // Acceleration for i (j anti-symmetric)

      // Pair-wise portion of grad velocity.
      const auto deltaDvDxi = fDij * vij.dyad(gradWGi);
      const auto deltaDvDxj = fDij * vij.dyad(gradWGj);

      // Specific thermal energy evolution.
      RAJA::atomicSub<RAJA::auto_atomic>(&DepsDti,
                                         mj*(sigmarhoi.doubledot(deltaDvDxi.Symmetric()) - workQi));
      RAJA::atomicSub<RAJA::auto_atomic>(&DepsDtj,
                                         mi*(sigmarhoj.doubledot(deltaDvDxj.Symmetric()) - workQj));
      // DepsDti -= mj*(sigmarhoi.doubledot(deltaDvDxi.Symmetric()) - workQi);
      // DepsDtj -= mi*(sigmarhoj.doubledot(deltaDvDxj.Symmetric()) - workQj);

      // Velocity gradient.
      DvDxi.atomicSub(mj*deltaDvDxi);
      DvDxj.atomicSub(mi*deltaDvDxj);
      // DvDxi -= mj*deltaDvDxi;
      // DvDxj -= mi*deltaDvDxj;
      if (sameMatij) {
        // localDvDxi -= mj*deltaDvDxi;
        // localDvDxj -= mi*deltaDvDxj;
        localDvDxi.atomicSub(mj*deltaDvDxi);
        localDvDxj.atomicSub(mi*deltaDvDxj);
      }

      // Estimate of delta v (for XSPH).
      if (XSPH and sameMatij) {
        const auto wXSPHij = 0.5*(mi/rhoi*Wi + mj/rhoj*Wj);
        RAJA::atomicAdd<RAJA::auto_atomic>(&XSPHWeightSumi, wXSPHij);
        RAJA::atomicAdd<RAJA::auto_atomic>(&XSPHWeightSumj, wXSPHij);
        XSPHDeltaVi.atomicSub(wXSPHij*vij);
        XSPHDeltaVj.atomicAdd(wXSPHij*vij);
        // XSPHWeightSumi += wXSPHij;
        // XSPHWeightSumj += wXSPHij;
        // XSPHDeltaVi -= wXSPHij*vij;
        // XSPHDeltaVj += wXSPHij*vij;
      }

      // Linear gradient correction term.
      Mi.atomicSub(mj*rij.dyad(gradWGi));
      Mj.atomicSub(mi*rij.dyad(gradWGj));
      // Mi -= mj*rij.dyad(gradWGi);
      // Mj -= mi*rij.dyad(gradWGj);
      if (sameMatij) {
        localMi.atomicSub(mj*rij.dyad(gradWGi));
        localMj.atomicSub(mi*rij.dyad(gradWGj));
        // localMi -= mj*rij.dyad(gradWGi);
        // localMj -= mi*rij.dyad(gradWGj);
      }

    }); // loop over pairs

    // Reduce the thread values to the master.
    //threadReduceFieldLists<Dimension>(threadStack);

  }   // OpenMP parallel region
  TIME_END("SolidSPHevalDerivs_pairs");
  // Finish up the derivatives for each point.
  TIME_BEGIN("SolidSPHevalDerivs_final");
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto& nodeList = mass_v[nodeListi]->nodeList();
    const auto ni = nodeList.numInternalNodes();

    // // Check if we can identify a reference density.
    // auto rho0 = 0.0;
    // try {
    //   rho0 = dynamic_cast<const SolidEquationOfState<Dimension>&>(dynamic_cast<const FluidNodeList<Dimension>&>(nodeList).equationOfState()).referenceDensity();
    //   // cerr << "Setting reference density to " << rho0 << endl;
    // } catch(...) {
    //   // cerr << "BLAGO!" << endl;
    // }

// #pragma omp parallel for
//     for (auto i = 0u; i < ni; ++i) {
    RAJA::forall<EXEC_POLICY>(TRS_UINT(0u, ni),
    [=] SPHERAL_HOST_DEVICE (size_t i) {

      // Get the state for node i.
      const auto& mi = mass(nodeListi, i);
      const auto& vi = velocity(nodeListi, i);
      const auto& rhoi = massDensity(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto& Si = S(nodeListi, i);
      const auto& mui = mu(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();
      //const auto  numNeighborsi = connectivityMap.numNeighborsForNode(nodeListi, i);
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
      auto& effViscousPressurei = effViscousPressure(nodeListi, i);
      auto& rhoSumCorrectioni = rhoSumCorrection(nodeListi, i);
      auto& XSPHWeightSumi = XSPHWeightSum(nodeListi, i);
      auto& XSPHDeltaVi = XSPHDeltaV(nodeListi, i);
      auto& DSDti = DSDt(nodeListi, i);

      // Add the self-contribution to density sum.
      rhoSumi += mi*W0*Hdeti;

      // Add the self-contribution to density sum correction.
      rhoSumCorrectioni += mi*WQ0*Hdeti/rhoi ;

      // Correct the effective viscous pressure.
      effViscousPressurei /= rhoSumCorrectioni ;

      // Finish the gradient of the velocity.
      CHECK(rhoi > 0.0);
      if (CorrectVelocityGradient and
          std::abs(Mi.Determinant()) > 1.0e-10) { // and
          //numNeighborsi > Dimension::pownu(2)) {
        Mi = Mi.Inverse();
        DvDxi = DvDxi*Mi;
      } else {
        DvDxi /= rhoi;
      }
      if (CorrectVelocityGradient and
          std::abs(localMi.Determinant()) > 1.0e-10) {// and
        //numNeighborsi > Dimension::pownu(2)) {
        localMi = localMi.Inverse();
        localDvDxi = localDvDxi*localMi;
      } else {
        localDvDxi /= rhoi;
      }

      // Evaluate the continuity equation.
      DrhoDti = -rhoi*DvDxi.Trace();

      // If needed finish the total energy derivative.
      if (evolveTotalEnergy) DepsDti = mi*(vi.dot(DvDti) + DepsDti);

      // Determine the position evolution, based on whether we're doing XSPH or not.
      DxDti = vi;
      if (XSPH) {
        CHECK(XSPHWeightSumi >= 0.0);
        XSPHWeightSumi += Hdeti*mi/rhoi*W0 + 1.0e-30;
        DxDti += XSPHDeltaVi/XSPHWeightSumi;
      }

      // Determine the deviatoric stress evolution.
      const auto deformation = localDvDxi.Symmetric();
      const auto spin = localDvDxi.SkewSymmetric();
      const auto deviatoricDeformation = deformation - deformation.Trace()/3.0*SymTensor::one();
      const auto spinCorrection = (spin*Si + Si*spin).Symmetric();
      DSDti = spinCorrection + 2.0*mui*deviatoricDeformation;

      // Optionally use damage to ramp down stress on damaged material.
      // const auto Di = max(0.0, min(1.0, damage(nodeListi, i).Trace() - 1.0));
      // Hideali = (1.0 - Di)*Hideali + Di*mHfield0(nodeListi, i);
      // DHDti = (1.0 - Di)*DHDti + Di*(mHfield0(nodeListi, i) - Hi)*0.25/dt;

      // // We also adjust the density evolution in the presence of damage.
      // if (rho0 > 0.0) DrhoDti = (1.0 - Di)*DrhoDti - Di * 0.05*(rhoi - rho0)*ci*Hi.Trace()/Dimension::nDim;

      // // In the presence of damage, add a term to reduce the stress on this point.
      // DSDti = (1.0 - Di)*DSDti - 0.25/dt*Di*Si;
    });
  }
  rhoSum.move(chai::CPU);
  DxDt.move(chai::CPU);
  DrhoDt.move(chai::CPU);
  DvDt.move(chai::CPU);
  DepsDt.move(chai::CPU);
  DvDx.move(chai::CPU);
  localDvDx.move(chai::CPU);
  M.move(chai::CPU);
  localM.move(chai::CPU);
  maxViscousPressure.move(chai::CPU);
  effViscousPressure.move(chai::CPU);
  XSPHWeightSum.move(chai::CPU);
  XSPHDeltaV.move(chai::CPU);
  DSDt.move(chai::CPU);
  rhoSumCorrection.move(chai::CPU);
  TIME_END("SolidSPHevalDerivs_final");
  TIME_END("SolidSPHevalDerivs");
}

//------------------------------------------------------------------------------
// Apply the ghost boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidSPH<Dimension>::
applyGhostBoundaries(State<Dimension>& state,
                     StateDerivatives<Dimension>& derivs) {
  TIME_BEGIN("SolidSPHghostBounds");

  // Ancestor method.
  SPH<Dimension>::applyGhostBoundaries(state, derivs);

  // Apply boundary conditions to our extra strength variables.
  auto S = state.fields(SolidFieldNames::deviatoricStress, SymTensor::zero());
  auto K = state.fields(SolidFieldNames::bulkModulus, 0.0);
  auto mu = state.fields(SolidFieldNames::shearModulus, 0.0);
  auto Y = state.fields(SolidFieldNames::yieldStrength, 0.0);
  auto fragIDs = state.fields(SolidFieldNames::fragmentIDs, int(1));
  auto pTypes = state.fields(SolidFieldNames::particleTypes, int(0));

  for (auto boundaryPtr: range(this->boundaryBegin(), this->boundaryEnd())) {
    boundaryPtr->applyFieldListGhostBoundary(S);
    boundaryPtr->applyFieldListGhostBoundary(K);
    boundaryPtr->applyFieldListGhostBoundary(mu);
    boundaryPtr->applyFieldListGhostBoundary(Y);
    boundaryPtr->applyFieldListGhostBoundary(fragIDs);
    boundaryPtr->applyFieldListGhostBoundary(pTypes);
  }
  TIME_END("SolidSPHghostBounds");
}

//------------------------------------------------------------------------------
// Enforce the boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidSPH<Dimension>::
enforceBoundaries(State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {
  TIME_BEGIN("SolidSPHenforceBounds");

  // Ancestor method.
  SPH<Dimension>::enforceBoundaries(state, derivs);

  // Enforce boundary conditions on the extra strength variable.s
  auto S = state.fields(SolidFieldNames::deviatoricStress, SymTensor::zero());
  auto K = state.fields(SolidFieldNames::bulkModulus, 0.0);
  auto mu = state.fields(SolidFieldNames::shearModulus, 0.0);
  auto Y = state.fields(SolidFieldNames::yieldStrength, 0.0);
  auto fragIDs = state.fields(SolidFieldNames::fragmentIDs, int(1));
  auto pTypes = state.fields(SolidFieldNames::particleTypes, int(0));

  for (auto boundaryPtr: range(this->boundaryBegin(), this->boundaryEnd())) {
    boundaryPtr->enforceFieldListBoundary(S);
    boundaryPtr->enforceFieldListBoundary(K);
    boundaryPtr->enforceFieldListBoundary(mu);
    boundaryPtr->enforceFieldListBoundary(Y);
    boundaryPtr->enforceFieldListBoundary(fragIDs);
    boundaryPtr->enforceFieldListBoundary(pTypes);
  }
  TIME_END("SolidSPHenforceBounds");
}

//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidSPH<Dimension>::
dumpState(FileIO& file, const string& pathName) const {

  // Ancestor method.
  SPH<Dimension>::dumpState(file, pathName);

  file.write(mDdeviatoricStressDt, pathName + "/DdeviatoricStressDt");
  file.write(mBulkModulus, pathName + "/bulkModulus");
  file.write(mShearModulus, pathName + "/shearModulus");
  file.write(mYieldStrength, pathName + "/yieldStrength");
  file.write(mPlasticStrain0, pathName + "/plasticStrain0");
}

//------------------------------------------------------------------------------
// Restore the state from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidSPH<Dimension>::
restoreState(const FileIO& file, const string& pathName) {
 
  // Ancestor method.
  SPH<Dimension>::restoreState(file, pathName);

  file.read(mDdeviatoricStressDt, pathName + "/DdeviatoricStressDt");
  file.read(mBulkModulus, pathName + "/bulkModulus");
  file.read(mShearModulus, pathName + "/shearModulus");
  file.read(mYieldStrength, pathName + "/yieldStrength");
  file.read(mPlasticStrain0, pathName + "/plasticStrain0");
}

}
