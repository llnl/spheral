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
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

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
    return Dim<1>::SymTensor::zero;
  }
}

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
  mDdeviatoricStressDt = dataBase.newSolidFieldList(SymTensor::zero, IncrementState<Dimension, Vector>::prefix() + SolidFieldNames::deviatoricStress);
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
  dataBase.resizeFluidFieldList(mDdeviatoricStressDt, SymTensor::zero, DSDtName, false);

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
  const auto  oneKernelQ = (W == WQ);
  const auto  oneKernelG = (W == WG);

  // A few useful constants we'll use in the following loop.
  const auto tiny = 1.0e-30;
  const auto W0 = W(0.0, 1.0);
  const auto WQ0 = WQ(0.0, 1.0);
  const auto epsTensile = this->epsilonTensile();
  const auto compatibleEnergy = this->compatibleEnergyEvolution();
  const auto evolveTotalEnergy = this->evolveTotalEnergy();
  const auto XSPH = this->XSPH();

  // The connectivity.
  const auto& connectivityMap = dataBase.connectivityMap();
  const auto& nodeLists = connectivityMap.nodeLists();
  const auto numNodeLists = nodeLists.size();
  auto&       pairs = const_cast<NodePairList&>(connectivityMap.nodePairList());
  const auto  npairs = pairs.size();
  // const auto& coupling = connectivityMap.coupling();

  // Get the state and derivative FieldLists.
  // State FieldLists.
  const auto mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const auto massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  const auto H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const auto pressure = state.fields(HydroFieldNames::pressure, 0.0);
  const auto soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  const auto omega = state.fields(HydroFieldNames::omegaGradh, 0.0);
  const auto S = state.fields(SolidFieldNames::deviatoricStress, SymTensor::zero);
  const auto mu = state.fields(SolidFieldNames::shearModulus, 0.0);
  const auto damage = state.fields(SolidFieldNames::tensorDamage, SymTensor::zero);
  const auto pTypes = state.fields(SolidFieldNames::particleTypes, int(0));
  const auto fragID = state.fields(SolidFieldNames::fragmentIDs, int(0));
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
  auto  rhoSum = derivs.fields(ReplaceState<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, 0.0);
  auto  DxDt = derivs.fields(IncrementState<Dimension, Vector>::prefix() + HydroFieldNames::position, Vector::zero);
  auto  DrhoDt = derivs.fields(IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, 0.0);
  auto  DvDt = derivs.fields(HydroFieldNames::hydroAcceleration, Vector::zero);
  auto  DepsDt = derivs.fields(IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
  auto  DvDx = derivs.fields(HydroFieldNames::velocityGradient, Tensor::zero);
  auto  localDvDx = derivs.fields(HydroFieldNames::internalVelocityGradient, Tensor::zero);
  auto  M = derivs.fields(HydroFieldNames::M_SPHCorrection, Tensor::zero);
  auto  localM = derivs.fields("local " + HydroFieldNames::M_SPHCorrection, Tensor::zero);
  auto  maxViscousPressure = derivs.fields(HydroFieldNames::maxViscousPressure, 0.0);
  auto  effViscousPressure = derivs.fields(HydroFieldNames::effectiveViscousPressure, 0.0);
  auto  rhoSumCorrection = derivs.fields(HydroFieldNames::massDensityCorrection, 0.0);
  auto  XSPHWeightSum = derivs.fields(HydroFieldNames::XSPHWeightSum, 0.0);
  auto  XSPHDeltaV = derivs.fields(HydroFieldNames::XSPHDeltaV, Vector::zero);
  auto  DSDt = derivs.fields(IncrementState<Dimension, SymTensor>::prefix() + SolidFieldNames::deviatoricStress, SymTensor::zero);
  auto* pairAccelerationsPtr = derivs.template getPtr<PairAccelerationsType>(HydroFieldNames::pairAccelerations);
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

  // The scale for the tensile correction.
  const auto& nodeList = mass[0]->nodeList();
  const auto  nPerh = nodeList.nodesPerSmoothingScale();
  const auto  WnPerh = W(1.0/nPerh, 1.0);
  TIME_END("SolidSPHevalDerivs_initial");

  // Walk all the interacting pairs.
  TIME_BEGIN("SolidSPHevalDerivs_pairs");
#pragma omp parallel
  {
    // Thread private  scratch variables.
    int i, j, nodeListi, nodeListj;
    Scalar Wi, gWi, WQi, gWQi, Wj, gWj, WQj, gWQj;
    Vector gradWi, gradWj, gradWQi, gradWQj, gradWGi, gradWGj;
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
      //const auto  epsi = specificThermalEnergy(nodeListi, i);
      const auto  Pi = pressure(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto  ci = soundSpeed(nodeListi, i);
      const auto  omegai = omega(nodeListi, i);
      const auto& Si = S(nodeListi, i);
      //const auto  mui = mu(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();
      const auto  safeOmegai = safeInv(omegai, tiny);
      const auto  pTypei = pTypes(nodeListi, i);
      const auto  fragIDi = fragID(nodeListi, i);
      CHECK(mi > 0.0);
      CHECK(rhoi > 0.0);
      CHECK(Hdeti > 0.0);

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
      //const auto epsj = specificThermalEnergy(nodeListj, j);
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
      W.kernelAndGradValue(etaMagi, Hdeti, Wi, gWi);
      W.kernelAndGradValue(etaMagj, Hdetj, Wj, gWj);
      gradWi = gWi*Hi*etaiUnit;
      gradWj = gWj*Hj*etajUnit;
      if (oneKernelQ) {
        WQi = Wi;
        WQj = Wj;
        gradWQi = gradWi;
        gradWQj = gradWj;
      } else {
        WQ.kernelAndGradValue(etaMagi, Hdeti, WQi, gWQi);
        WQ.kernelAndGradValue(etaMagj, Hdetj, WQj, gWQj);
        gradWQi = gWQi*Hi*etaiUnit;
        gradWQj = gWQj*Hj*etajUnit;
      }
      if (oneKernelG) {
        gradWGi = gradWi;
        gradWGj = gradWj;
      } else {
        gradWGi = Hi*etaiUnit * WG.gradValue(etaMagi, Hdeti);
        gradWGj = Hj*etajUnit * WG.gradValue(etaMagj, Hdetj);
      }

      // Contribution to the sum density (only if the same material).
      if (nodeListi == nodeListj) {
        rhoSumi += mj*Wi;
        rhoSumj += mi*Wj;
      }

      // Contribution to the sum density correction
      rhoSumCorrectioni += mj * WQi / rhoj ;
      rhoSumCorrectionj += mi * WQj / rhoi ;

      // Compute the pair-wise artificial viscosity.
      const auto vij = vi - vj;
      Q->QPiij(QPiij, QPiji, Qi, Qj,
               nodeListi, i, nodeListj, j,
               ri, Hi, etai, vi, rhoi, ci,  
               rj, Hj, etaj, vj, rhoj, cj,
               fClQ, fCqQ, DvDxQ); 
      const auto Qacci = 0.5*(QPiij*gradWQi);
      const auto Qaccj = 0.5*(QPiji*gradWQj);
      const auto workQi = vij.dot(Qacci);
      const auto workQj = vij.dot(Qaccj);
      maxViscousPressurei = max(maxViscousPressurei, Qi);
      maxViscousPressurej = max(maxViscousPressurej, Qj);
      effViscousPressurei += mj*Qi*WQi/rhoj;
      effViscousPressurej += mi*Qj*WQj/rhoi;

      // Compute the stress tensors.
      if (sameMatij) {
        sigmai = fDij*Si - Pi * SymTensor::one;
        sigmaj = fDij*Sj - Pj * SymTensor::one;
      } else {
        sigmai = -Pi * SymTensor::one;
        sigmaj = -Pj * SymTensor::one;
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
        DvDti += mj*deltaDvDt;
        DvDtj -= mi*deltaDvDt;
      }
      if (compatibleEnergy) (*pairAccelerationsPtr)[kk] = mj*deltaDvDt;  // Acceleration for i (j anti-symmetric)

      // Pair-wise portion of grad velocity.
      const auto deltaDvDxi = fDij * vij.dyad(gradWGi);
      const auto deltaDvDxj = fDij * vij.dyad(gradWGj);

      // Specific thermal energy evolution.
      DepsDti -= mj*(sigmarhoi.doubledot(deltaDvDxi.Symmetric()) - workQi);
      DepsDtj -= mi*(sigmarhoj.doubledot(deltaDvDxj.Symmetric()) - workQj);

      // Velocity gradient.
      DvDxi -= mj*deltaDvDxi;
      DvDxj -= mi*deltaDvDxj;
      if (sameMatij) {
        localDvDxi -= mj*deltaDvDxi;
        localDvDxj -= mi*deltaDvDxj;
      }

      // Estimate of delta v (for XSPH).
      if (XSPH and sameMatij) {
        const auto wXSPHij = 0.5*(mi/rhoi*Wi + mj/rhoj*Wj);
        XSPHWeightSumi += wXSPHij;
        XSPHWeightSumj += wXSPHij;
        XSPHDeltaVi -= wXSPHij*vij;
        XSPHDeltaVj += wXSPHij*vij;
      }

      // Linear gradient correction term.
      Mi -= mj*rij.dyad(gradWGi);
      Mj -= mi*rij.dyad(gradWGj);
      if (sameMatij) {
        localMi -= mj*rij.dyad(gradWGi);
        localMj -= mi*rij.dyad(gradWGj);
      }

    } // loop over pairs

    // Reduce the thread values to the master.
    threadReduceFieldLists<Dimension>(threadStack);

  }   // OpenMP parallel region
  TIME_END("SolidSPHevalDerivs_pairs");

  // Finish up the derivatives for each point.
  TIME_BEGIN("SolidSPHevalDerivs_final");
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto& nodeList = mass[nodeListi]->nodeList();
    const auto ni = nodeList.numInternalNodes();

    // // Check if we can identify a reference density.
    // auto rho0 = 0.0;
    // try {
    //   rho0 = dynamic_cast<const SolidEquationOfState<Dimension>&>(dynamic_cast<const FluidNodeList<Dimension>&>(nodeList).equationOfState()).referenceDensity();
    //   // cerr << "Setting reference density to " << rho0 << endl;
    // } catch(...) {
    //   // cerr << "BLAGO!" << endl;
    // }

#pragma omp parallel for
    for (auto i = 0u; i < ni; ++i) {

      // Get the state for node i.
      const auto& mi = mass(nodeListi, i);
      const auto& vi = velocity(nodeListi, i);
      const auto& rhoi = massDensity(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto& Si = S(nodeListi, i);
      const auto& mui = mu(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();
      const auto  numNeighborsi = connectivityMap.numNeighborsForNode(nodeListi, i);
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
      if (this->mCorrectVelocityGradient and
          std::abs(Mi.Determinant()) > 1.0e-10 and
          numNeighborsi > Dimension::pownu(2)) {
        Mi = Mi.Inverse();
        DvDxi = DvDxi*Mi;
      } else {
        DvDxi /= rhoi;
      }
      if (this->mCorrectVelocityGradient and
          std::abs(localMi.Determinant()) > 1.0e-10 and
          numNeighborsi > Dimension::pownu(2)) {
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
      const auto deviatoricDeformation = deformation - deformation.Trace()/3.0*SymTensor::one;
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
    }
  }
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
  auto S = state.fields(SolidFieldNames::deviatoricStress, SymTensor::zero);
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
  auto S = state.fields(SolidFieldNames::deviatoricStress, SymTensor::zero);
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
