//---------------------------------Spheral++----------------------------------//
// PSPHHydroBase -- The PSPH/APSPH hydrodynamic package for Spheral++.
//
// Created by JMO, Wed Dec 16 20:52:02 PST 2015
//----------------------------------------------------------------------------//
#include "PSPH.hh"
#include "FileIO/FileIO.hh"
#include "computeSPHSumMassDensity.hh"
#include "computeSumVoronoiCellMassDensity.hh"
#include "computePSPHCorrections.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Physics/GenericHydro.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/ReplaceState.hh"
#include "DataBase/ReplaceBoundedState.hh"
#include "DataBase/applyPolicyToFieldList.hh"
#include "Hydro/GammaPolicy.hh"
#include "Hydro/PressurePolicy.hh"
#include "Hydro/SoundSpeedPolicy.hh"
#include "Mesh/generateMesh.hh"
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
#include "Utilities/globalBoundingVolumes.hh"
#include "Mesh/Mesh.hh"
#include "CRKSPH/volumeSpacing.hh"


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
PSPH<Dimension>::
PSPH(DataBase<Dimension>& dataBase,
     ArtificialViscosity<Dimension>& Q,
     const TableKernel<Dimension>& W,
     const TableKernel<Dimension>& WPi,
     const double cfl,
     const bool useVelocityMagnitudeForDt,
     const bool compatibleEnergyEvolution,
     const bool evolveTotalEnergy,
     const bool XSPH,
     const bool correctVelocityGradient,
     const bool HopkinsConductivity,
     const bool sumMassDensityOverAllNodeLists,
     const MassDensityType densityUpdate,
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
                 true,
                 XSPH,
                 correctVelocityGradient,
                 sumMassDensityOverAllNodeLists,
                 densityUpdate,
                 0.0,
                 1.0,
                 xmin,
                 xmax),
  mHopkinsConductivity(HopkinsConductivity),
  mGamma(FieldStorageType::CopyFields),
  mPSPHcorrection(FieldStorageType::CopyFields) {

  // Create storage for our internal state.
  mGamma = dataBase.newFluidFieldList(0.0, HydroFieldNames::gamma);
  mPSPHcorrection = dataBase.newFluidFieldList(0.0, HydroFieldNames::PSPHcorrection);
  dataBase.fluidGamma(mGamma);
}

//------------------------------------------------------------------------------
// On problem start up, we need to initialize our internal data.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PSPH<Dimension>::
initializeProblemStartupDependencies(DataBase<Dimension>& dataBase,
                                     State<Dimension>& state,
                                     StateDerivatives<Dimension>& derivs) {

  // SPH<Dimension>::initializeProblemStartupDependencies(dataBase, state, derivs);

  // The SPH class tries to update these using the policies, but since PSPH
  // handles them differently we have to override that behavior here and do
  // it differently.
  applyPolicyToFieldList(this->mPressure, make_policy<PressurePolicy<Dimension>>(), state, derivs);
  applyPolicyToFieldList(this->mSoundSpeed, make_policy<SoundSpeedPolicy<Dimension>>(), state, derivs);
}

//------------------------------------------------------------------------------
// Register the state we need/are going to evolve.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PSPH<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {

  // SPH does most of it.
  SPH<Dimension>::registerState(dataBase, state);

  // Make sure we're the right size.
  dataBase.resizeFluidFieldList(mGamma, 0.0, HydroFieldNames::gamma, false);
  dataBase.resizeFluidFieldList(mPSPHcorrection, 0.0, HydroFieldNames::PSPHcorrection, false);

  // We also require the fluid gamma.
  state.enroll(mGamma, make_policy<GammaPolicy<Dimension>>());

  // Override the default policies for pressure and sound speed.  We'll compute those
  // specially in the postStateUpdate.  Same goes for registering the PSPHcorrections.
  state.enroll(mPSPHcorrection);
  state.removePolicy(this->mPressure, true);
  state.removePolicy(this->mSoundSpeed, true);
}

//------------------------------------------------------------------------------
// Pre-step initializations.  Since the topology has just been changed we need
// to recompute the PSPH corrections.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PSPH<Dimension>::
preStepInitialize(const DataBase<Dimension>& dataBase,
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {

  // Base class
  SPH<Dimension>::preStepInitialize(dataBase, state, derivs);

  // Do the PSPH corrections.
  const TableKernel<Dimension>& W = this->kernel();
  const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
  const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
  const FieldList<Dimension, Scalar> specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  const FieldList<Dimension, Scalar> gamma = state.fields(HydroFieldNames::gamma, 0.0);
  const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);
  FieldList<Dimension, Scalar> rho = state.fields(HydroFieldNames::massDensity, 0.0);
  FieldList<Dimension, Scalar> P = state.fields(HydroFieldNames::pressure, 0.0);
  FieldList<Dimension, Scalar> cs = state.fields(HydroFieldNames::soundSpeed, 0.0);
  FieldList<Dimension, Scalar> PSPHcorrection = state.fields(HydroFieldNames::PSPHcorrection, 0.0);
  computePSPHCorrections(connectivityMap, W, mass, position, specificThermalEnergy, gamma, H, 
                         (this->mDensityUpdate != MassDensityType::IntegrateDensity),
                         rho, P, cs, PSPHcorrection);
  for (auto boundaryPtr: range(this->boundaryBegin(), this->boundaryEnd())) {
    boundaryPtr->applyFieldListGhostBoundary(rho);
    boundaryPtr->applyFieldListGhostBoundary(P);
    boundaryPtr->applyFieldListGhostBoundary(cs);
    boundaryPtr->applyFieldListGhostBoundary(PSPHcorrection);
  }
  for (auto boundaryPtr: range(this->boundaryBegin(), this->boundaryEnd())) boundaryPtr->finalizeGhostBoundary();
}

//------------------------------------------------------------------------------
// Post-state update.  This is where we can recompute the PSPH pressure and
// corrections.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
PSPH<Dimension>::
postStateUpdate(const Scalar t, 
                const Scalar dt,
                const DataBase<Dimension>& dataBase, 
                State<Dimension>& state,
                StateDerivatives<Dimension>& derivs) {

  // Check if the base class needs anything
  SPH<Dimension>::postStateUpdate(t, dt, dataBase, state, derivs);

  // Do the PSPH corrections.
  const TableKernel<Dimension>& W = this->kernel();
  const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
  const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
  const FieldList<Dimension, Scalar> specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  const FieldList<Dimension, Scalar> gamma = state.fields(HydroFieldNames::gamma, 0.0);
  const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);
  FieldList<Dimension, Scalar> rho = state.fields(HydroFieldNames::massDensity, 0.0);
  FieldList<Dimension, Scalar> P = state.fields(HydroFieldNames::pressure, 0.0);
  FieldList<Dimension, Scalar> cs = state.fields(HydroFieldNames::soundSpeed, 0.0);
  FieldList<Dimension, Scalar> PSPHcorrection = state.fields(HydroFieldNames::PSPHcorrection, 0.0);
  computePSPHCorrections(connectivityMap, W, mass, position, specificThermalEnergy, gamma, H,
                         (this->mDensityUpdate != MassDensityType::IntegrateDensity),
                         rho, P, cs, PSPHcorrection);

  // We depend on the caller knowing to finalize the ghost boundaries!
  return true;
}

//------------------------------------------------------------------------------
// Determine the principle derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PSPH<Dimension>::
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
// The actual implementation
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename QType>
void
PSPH<Dimension>::
evaluateDerivativesImpl(const typename Dimension::Scalar /*time*/,
                        const typename Dimension::Scalar /*dt*/,
                        const DataBase<Dimension>& dataBase,
                        const State<Dimension>& state,
                        StateDerivatives<Dimension>& derivs,
                        chai::managed_ptr<QType> Q) const {

  using QPiType = typename QType::ReturnType;

  // The kernels and such.
  const auto& W = this->kernel();
  const auto& WQ = this->PiKernel();

  // A few useful constants we'll use in the following loop.
  const double tiny = 1.0e-10;
  const Scalar W0 = W(0.0, 1.0);
  const auto compatibleEnergy = this->compatibleEnergyEvolution();
  const auto XSPH = this->XSPH();

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
  const auto position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const auto massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  const auto H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const auto pressure = state.fields(HydroFieldNames::pressure, 0.0);
  const auto soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  const auto gamma = state.fields(HydroFieldNames::gamma, 0.0);
  const auto PSPHcorrection = state.fields(HydroFieldNames::PSPHcorrection, 0.0);
  const auto reducingViscosityMultiplierQ = state.fields(HydroFieldNames::ArtificialViscousCqMultiplier, 0.0, true);
  const auto reducingViscosityMultiplierL = state.fields(HydroFieldNames::ArtificialViscousClMultiplier, 0.0, true);
  const auto DvDxQ = state.fields(HydroFieldNames::ArtificialViscosityVelocityGradient, Tensor::zero);

  CHECK(mass.size() == numNodeLists);
  CHECK(position.size() == numNodeLists);
  CHECK(velocity.size() == numNodeLists);
  CHECK(massDensity.size() == numNodeLists);
  CHECK(specificThermalEnergy.size() == numNodeLists);
  CHECK(H.size() == numNodeLists);
  CHECK(pressure.size() == numNodeLists);
  CHECK(soundSpeed.size() == numNodeLists);
  CHECK(gamma.size() == numNodeLists);
  CHECK(PSPHcorrection.size() == numNodeLists);
  CHECK((not mHopkinsConductivity) or (reducingViscosityMultiplierQ.size() == numNodeLists));
  CHECK((not mHopkinsConductivity) or (reducingViscosityMultiplierL.size() == numNodeLists));

  // Derivative FieldLists.
  auto  rhoSum = derivs.fields(ReplaceState<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, 0.0);
  auto  normalization = derivs.fields(HydroFieldNames::normalization, 0.0);
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
  auto* pairAccelerationsPtr = derivs.template getPtr<PairAccelerationsType>(HydroFieldNames::pairAccelerations);
  auto  XSPHWeightSum = derivs.fields(HydroFieldNames::XSPHWeightSum, 0.0);
  auto  XSPHDeltaV = derivs.fields(HydroFieldNames::XSPHDeltaV, Vector::zero);
  CHECK(rhoSum.size() == numNodeLists);
  CHECK(normalization.size() == numNodeLists);
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
  CHECK(XSPHWeightSum.size() == numNodeLists);
  CHECK(XSPHDeltaV.size() == numNodeLists);
  CHECK((compatibleEnergy and pairAccelerationsPtr->size() == npairs) or not compatibleEnergy);

  // Walk all the interacting pairs.
#pragma omp parallel
  {
    // Thread private scratch variables
    unsigned i, j, nodeListi, nodeListj;
    Scalar Wi, gWi, WQi, gWQi, Wj, gWj, WQj, gWQj, Qi, Qj;
    QPiType QPiij, QPiji;

    typename SpheralThreads<Dimension>::FieldListStack threadStack;
    auto rhoSum_thread = rhoSum.threadCopy(threadStack);
    auto normalization_thread = normalization.threadCopy(threadStack);
    auto DvDt_thread = DvDt.threadCopy(threadStack);
    auto DepsDt_thread = DepsDt.threadCopy(threadStack);
    auto DvDx_thread = DvDx.threadCopy(threadStack);
    auto localDvDx_thread = localDvDx.threadCopy(threadStack);
    auto M_thread = M.threadCopy(threadStack);
    auto localM_thread = localM.threadCopy(threadStack);
    auto maxViscousPressure_thread = maxViscousPressure.threadCopy(threadStack, ThreadReduction::MAX);
    auto effViscousPressure_thread = effViscousPressure.threadCopy(threadStack);
    auto XSPHWeightSum_thread = XSPHWeightSum.threadCopy(threadStack);
    auto XSPHDeltaV_thread = XSPHDeltaV.threadCopy(threadStack);

#pragma omp for
    for (auto kk = 0u; kk < npairs; ++kk) {
      i = pairs[kk].i_node;
      j = pairs[kk].j_node;
      nodeListi = pairs[kk].i_list;
      nodeListj = pairs[kk].j_list;

      // Get the state for node i.
      const auto& ri = position(nodeListi, i);
      const auto& mi = mass(nodeListi, i);
      const auto& vi = velocity(nodeListi, i);
      const auto& rhoi = massDensity(nodeListi, i);
      const auto& epsi = specificThermalEnergy(nodeListi, i);
      const auto& Pi = pressure(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto& ci = soundSpeed(nodeListi, i);
      const auto& gammai = gamma(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();
      CHECK(mi > 0.0);
      CHECK(rhoi > 0.0);
      CHECK(Hdeti > 0.0);

      auto& rhoSumi = rhoSum_thread(nodeListi, i);
      auto& normi = normalization_thread(nodeListi, i);
      auto& DvDti = DvDt_thread(nodeListi, i);
      auto& DepsDti = DepsDt_thread(nodeListi, i);
      auto& DvDxi = DvDx_thread(nodeListi, i);
      auto& localDvDxi = localDvDx_thread(nodeListi, i);
      auto& Mi = M_thread(nodeListi, i);
      auto& localMi = localM_thread(nodeListi, i);
      auto& maxViscousPressurei = maxViscousPressure_thread(nodeListi, i);
      auto& effViscousPressurei = effViscousPressure_thread(nodeListi, i);
      auto& XSPHWeightSumi = XSPHWeightSum_thread(nodeListi, i);
      auto& XSPHDeltaVi = XSPHDeltaV_thread(nodeListi, i);

      // Get the state for node j
      const auto& rj = position(nodeListj, j);
      const auto& mj = mass(nodeListj, j);
      const auto& vj = velocity(nodeListj, j);
      const auto& rhoj = massDensity(nodeListj, j);
      const auto& epsj = specificThermalEnergy(nodeListj, j);
      const auto& Pj = pressure(nodeListj, j);
      const auto& Hj = H(nodeListj, j);
      const auto& cj = soundSpeed(nodeListj, j);
      const auto& gammaj = gamma(nodeListj, j);
      const auto  Hdetj = Hj.Determinant();
      CHECK(mj > 0.0);
      CHECK(rhoj > 0.0);
      CHECK(Hdetj > 0.0);

      auto& rhoSumj = rhoSum_thread(nodeListj, j);
      auto& normj = normalization_thread(nodeListj, j);
      auto& DvDtj = DvDt_thread(nodeListj, j);
      auto& DepsDtj = DepsDt_thread(nodeListj, j);
      auto& DvDxj = DvDx_thread(nodeListj, j);
      auto& localDvDxj = localDvDx_thread(nodeListj, j);
      auto& Mj = M_thread(nodeListj, j);
      auto& localMj = localM_thread(nodeListj, j);
      auto& maxViscousPressurej = maxViscousPressure_thread(nodeListj, j);
      auto& effViscousPressurej = effViscousPressure_thread(nodeListj, j);
      auto& XSPHWeightSumj = XSPHWeightSum_thread(nodeListj, j);
      auto& XSPHDeltaVj = XSPHDeltaV_thread(nodeListj, j);

      // Flag if this is a contiguous material pair or not.
      const bool sameMatij = true; // (nodeListi == nodeListj and fragIDi == fragIDj);

      // Node displacement.
      const auto rij = ri - rj;
      const auto etai = Hi*rij;
      const auto etaj = Hj*rij;
      const auto etaMagi = etai.magnitude();
      const auto etaMagj = etaj.magnitude();
      CHECK(etaMagi >= 0.0);
      CHECK(etaMagj >= 0.0);

      // Symmetrized kernel weight and gradient.
      W.kernelAndGradValue(etaMagi, Hdeti, Wi, gWi);
      WQ.kernelAndGradValue(etaMagi, Hdeti, WQi, gWQi);
      const auto Hetai = Hi*etai.unitVector();
      const auto gradWi = gWi*Hetai;
      const auto gradWQi = gWQi*Hetai;

      W.kernelAndGradValue(etaMagj, Hdetj, Wj, gWj);
      WQ.kernelAndGradValue(etaMagj, Hdetj, WQj, gWQj);
      const auto Hetaj = Hj*etaj.unitVector();
      const auto gradWj = gWj*Hetaj;
      const auto gradWQj = gWQj*Hetaj;

      // Contribution to the sum density.
      if (nodeListi == nodeListj) {
        rhoSumi += mj*Wj;
        rhoSumj += mi*Wi;
        normi += mi/rhoi*Wj;
        normj += mj/rhoj*Wi;
      }

      // Compute the pair-wise artificial viscosity.
      const auto vij = vi - vj;
      Q->QPiij(QPiij, QPiji, Qi, Qj,
               nodeListi, i, nodeListj, j,
               ri, Hi, etai, vi, rhoi, ci,  
               rj, Hj, etaj, vj, rhoj, cj,
               reducingViscosityMultiplierL, reducingViscosityMultiplierQ, DvDxQ); 
      const auto Qacci = 0.5*(QPiij*gradWQi);
      const auto Qaccj = 0.5*(QPiji*gradWQj);
      // const auto workQi = 0.5*(QPiij*vij).dot(gradWQi);
      // const auto workQj = 0.5*(QPiji*vij).dot(gradWQj);
      const auto workQi = vij.dot(Qacci);
      const auto workQj = vij.dot(Qaccj);
      maxViscousPressurei = max(maxViscousPressurei, Qi);
      maxViscousPressurej = max(maxViscousPressurej, Qj);
      effViscousPressurei += mj*Qi*WQi/rhoj;
      effViscousPressurej += mi*Qj*WQj/rhoi;

      // Acceleration.
      CHECK(rhoi > 0.0);
      CHECK(rhoj > 0.0);
      const auto& Fcorri=PSPHcorrection(nodeListi, i);
      const auto& Fcorrj=PSPHcorrection(nodeListj, j);
      const auto  Fij=1.0-Fcorri*safeInv(mj*epsj, tiny);
      const auto  Fji=1.0-Fcorrj*safeInv(mi*epsi, tiny);
      const auto  engCoef=(gammai-1)*(gammaj-1)*epsi*epsj;
      const auto  deltaDvDt = engCoef*(gradWi*Fij*safeInv(Pi, tiny) + gradWj*Fji*safeInv(Pj, tiny)) + Qacci + Qaccj;
      DvDti -= mj*deltaDvDt;
      DvDtj += mi*deltaDvDt;
      if (compatibleEnergy) (*pairAccelerationsPtr)[kk] = -mj*deltaDvDt;  // Acceleration for i (j anti-symmetric)

      // Specific thermal energy evolution.
      DepsDti += mj*(engCoef*Fij*safeInv(Pi, tiny)*vij.dot(gradWi) + workQi);
      DepsDtj += mi*(engCoef*Fji*safeInv(Pj, tiny)*vij.dot(gradWj) + workQj);

      //ADD ARITIFICIAL CONDUCTIVITY IN HOPKINS 2014A
      if (mHopkinsConductivity) {
        const auto  alph_c = 0.25;//Parameter = 0.25 in Hopkins 2014
        const auto  Vs = ci+cj-3.0*vij.dot(rij.unitVector());
        const auto& Qalpha_i = reducingViscosityMultiplierL(nodeListi, i); //Both L and Q corrections are the same for Cullen Viscosity
        const auto& Qalpha_j = reducingViscosityMultiplierL(nodeListj, j); //Both L and Q corrections are the same for Cullen Viscosity
        //DepsDti += (Vs > 0.0)*alph_c*mi*mj*(Qalpha_i+Qalpha_j)*0.5*Vs*(epsi-epsj)*abs(Pi-Pj)*((gradWi+gradWj).dot(rij.unitVector()))*safeInv((Pi+Pj+1e-30)*(rhoi+rhoj+1e-30));
        //DepsDtj += (Vs > 0.0)*alph_c*mi*mj*(Qalpha_i+Qalpha_j)*0.5*Vs*(epsi-epsj)*abs(Pi-Pj)*((gradWi+gradWj).dot(rij.unitVector()))*safeInv((Pi+Pj+1e-30)*(rhoi+rhoj+1e-30));
        //DepsDti += (Vs > 0.0) ? alph_c*mj*(Qalpha_i+Qalpha_j)*0.5*Vs*(epsi-epsj)*abs(Pi-Pj)*((gradWi).dot(rij.unitVector()))/max((Pi+Pj)*0.5*(rhoi+rhoj),tiny) : 0.0;
        //DepsDtj += (Vs > 0.0) ? alph_c*mi*(Qalpha_i+Qalpha_j)*0.5*Vs*(epsj-epsi)*abs(Pi-Pj)*((gradWj).dot(rij.unitVector()))/max((Pi+Pj)*0.5*(rhoi+rhoj),tiny) : 0.0;
        DepsDti += (Vs > 0.0) ? alph_c*mj*(Qalpha_i+Qalpha_j)*0.5*Vs*(epsi-epsj)*abs(Pi-Pj)*((gradWi+gradWj).dot(rij.unitVector()))*safeInv((Pi+Pj)*(rhoi+rhoj),tiny) : 0.0;
        DepsDtj += (Vs > 0.0) ? alph_c*mi*(Qalpha_i+Qalpha_j)*0.5*Vs*(epsj-epsi)*abs(Pi-Pj)*((gradWi+gradWj).dot(rij.unitVector()))*safeInv((Pi+Pj)*(rhoi+rhoj),tiny) : 0.0;
        //const Scalar tmpi = (Vs > 0.0) ? alph_c*mj*(Qalpha_i+Qalpha_j)*0.5*Vs*(epsi-epsj)*abs(Pi-Pj)*((gradWi+gradWj).dot(rij.unitVector()))*safeInv((Pi+Pj)*(rhoi+rhoj)) : 0.0;
        //const Scalar tmpj = (Vs > 0.0) ? alph_c*mi*(Qalpha_i+Qalpha_j)*0.5*Vs*(epsj-epsi)*abs(Pi-Pj)*((gradWi+gradWj).dot(rij.unitVector()))*safeInv((Pi+Pj)*(rhoi+rhoj)) : 0.0;
        //DepsDti += tmpi;
        //DepsDtj += tmpj;
        //DepsDti += (Vs > 0.0) ? alph_c*mj*(Qalpha_i+Qalpha_j)*0.5*Vs*(epsi-epsj)*abs(Pi-Pj)*((gradWi+gradWj).dot(rij.unitVector()))*safeInv((Pi+Pj)*(rhoi+rhoj)) : 0.0;
        //DepsDtj += (Vs > 0.0) ? alph_c*mi*(Qalpha_i+Qalpha_j)*0.5*Vs*(epsj-epsi)*abs(Pi-Pj)*((gradWi+gradWj).dot(rij.unitVector()))*safeInv((Pi+Pj)*(rhoi+rhoj)) : 0.0;
        //DepsDti += (Vs > 0.0 && Pi > 1e-4 && Pj > 1e-4) ? alph_c*mj*(Qalpha_i+Qalpha_j)*0.5*Vs*(epsi-epsj)*abs(Pi-Pj)*((gradWi+gradWj).dot(rij.unitVector()))*safeInv((Pi+Pj)*(rhoi+rhoj),tiny) : 0.0;
        //DepsDtj += (Vs > 0.0 && Pi > 1e-4 && Pj > 1e-4) ? alph_c*mi*(Qalpha_i+Qalpha_j)*0.5*Vs*(epsj-epsi)*abs(Pi-Pj)*((gradWi+gradWj).dot(rij.unitVector()))*safeInv((Pi+Pj)*(rhoi+rhoj),tiny) : 0.0;
      }

      // Velocity gradient.
      const auto deltaDvDxi = mj*vij.dyad(gradWi);
      const auto deltaDvDxj = mi*vij.dyad(gradWj);
      DvDxi -= deltaDvDxi; 
      DvDxj -= deltaDvDxj;
      if (sameMatij) {
        localDvDxi -= deltaDvDxi; 
        localDvDxj -= deltaDvDxj;
      }

      // Estimate of delta v (for XSPH).
      if (XSPH and (sameMatij)) {
        const auto wXSPHij = 0.5*(mi/rhoi*Wi + mj/rhoj*Wj);
        XSPHWeightSumi += wXSPHij;
        XSPHWeightSumj += wXSPHij;
        XSPHDeltaVi -= wXSPHij*vij;
        XSPHDeltaVj += wXSPHij*vij;
      }

      // Linear gradient correction term.
      Mi -= mj*rij.dyad(gradWi);
      Mj -= mi*rij.dyad(gradWj);
      if (sameMatij) {
        localMi -= mj*rij.dyad(gradWi);
        localMj -= mi*rij.dyad(gradWj);
      }

    } // loop over pairs

    // Reduce the thread values to the master.
    threadReduceFieldLists<Dimension>(threadStack);

  }   // OpenMP parallel region

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
      const auto& Hi = H(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();
      const auto numNeighborsi = connectivityMap.numNeighborsForNode(nodeListi, i);
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
      auto& Mi = M(nodeListi, i);
      auto& localMi = localM(nodeListi, i);
      auto& XSPHWeightSumi = XSPHWeightSum(nodeListi, i);
      auto& XSPHDeltaVi = XSPHDeltaV(nodeListi, i);

      // Add the self-contribution to density sum.
      rhoSumi += mi*W0*Hdeti;
      normi += mi/rhoi*W0*Hdeti;

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
      if (this->mEvolveTotalEnergy) DepsDti = mi*(vi.dot(DvDti) + DepsDti);

      // Determine the position evolution, based on whether we're doing XSPH or not.
      if (this->XSPH()) {
        XSPHWeightSumi += Hdeti*mi/rhoi*W0;
        CHECK2(XSPHWeightSumi != 0.0, i << " " << XSPHWeightSumi);
        DxDti = vi + XSPHDeltaVi*safeInvVar(XSPHWeightSumi, tiny);
      } else {
        DxDti = vi;
      }
    }
  }
}

//------------------------------------------------------------------------------
// Apply the ghost boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PSPH<Dimension>::
applyGhostBoundaries(State<Dimension>& state,
                     StateDerivatives<Dimension>& derivs) {

  // SPH does most of it.
  SPH<Dimension>::applyGhostBoundaries(state, derivs);

  // Apply boundary conditions to the PSPH corrections.
  FieldList<Dimension, Scalar> gamma = state.fields(HydroFieldNames::gamma, 0.0);
  FieldList<Dimension, Scalar> PSPHcorrection = state.fields(HydroFieldNames::PSPHcorrection, 0.0);

  for (auto boundaryPtr: range(this->boundaryBegin(), this->boundaryEnd())) {
    boundaryPtr->applyFieldListGhostBoundary(gamma);
    boundaryPtr->applyFieldListGhostBoundary(PSPHcorrection);
  }}

//------------------------------------------------------------------------------
// Enforce the boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PSPH<Dimension>::
enforceBoundaries(State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {

  // SPH does most of it.
  SPH<Dimension>::enforceBoundaries(state, derivs);

  // Enforce boundary conditions on the PSPH corrections.
  FieldList<Dimension, Scalar> gamma = state.fields(HydroFieldNames::gamma, 0.0);
  FieldList<Dimension, Scalar> PSPHcorrection = state.fields(HydroFieldNames::PSPHcorrection, 0.0);

  for (auto boundaryPtr: range(this->boundaryBegin(), this->boundaryEnd())) {
    boundaryPtr->enforceFieldListBoundary(gamma);
    boundaryPtr->enforceFieldListBoundary(PSPHcorrection);
  }
}

//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PSPH<Dimension>::
dumpState(FileIO& file, const string& pathName) const {

  // SPH does most of it.
  SPH<Dimension>::dumpState(file, pathName);

  file.write(mGamma, pathName + "/gamma");
  file.write(mPSPHcorrection, pathName + "/PSPHcorrection");
}

//------------------------------------------------------------------------------
// Restore the state from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PSPH<Dimension>::
restoreState(const FileIO& file, const string& pathName) {
 
  // SPH does most of it.
  SPH<Dimension>::restoreState(file, pathName);

  file.read(mGamma, pathName + "/gamma");
  file.read(mPSPHcorrection, pathName + "/PSPHcorrection");
}

}
