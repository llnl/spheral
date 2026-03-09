//---------------------------------Spheral++----------------------------------//
// SPHRZ -- An SPH/ASPH hydrodynamic package for Spheral++,
//          specialized for 2D RZ (cylindrical) geometry.
//
// This RZ version is a naive area-weighting implementation, nothing as
// highfalutin as the Garcia-Senz approach.
//
// Note this version is currently abusing our ordinary 2D geometric types,
// implicitly mapping x->z, y->r.
//
// Created by JMO, Fri May  6 16:18:36 PDT 2016
//----------------------------------------------------------------------------//
#include "SPHRZ.hh"
#include "FileIO/FileIO.hh"
#include "computeSPHSumMassDensity.hh"
#include "correctSPHSumMassDensity.hh"
#include "computeSumVoronoiCellMassDensity.hh"
#include "computeSPHOmegaGradhCorrection.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Hydro/SpecificThermalEnergyPolicy.hh"
#include "Hydro/SpecificFromTotalThermalEnergyPolicy.hh"
#include "Physics/GenericHydro.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/ReplaceState.hh"
#include "DataBase/ReplaceBoundedState.hh"
#include "Hydro/RZNonSymmetricSpecificThermalEnergyPolicy.hh"
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
#include "Utilities/Timer.hh"
#include "Mesh/Mesh.hh"
#include "CRKSPH/volumeSpacing.hh"
#include "Geometry/GeometryRegistrar.hh"
#include "Geometry/toroidalVolume.hh"

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
SPHRZ::
SPHRZ(DataBase<Dimension>& dataBase,
               ArtificialViscosityHandle<Dim<2>>& Q,
               const TableKernel<Dim<2>>& W,
               const TableKernel<Dim<2>>& WPi,
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
  SPHBase<Dim<2>>(dataBase,
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
  mPairAccelerationsPtr(std::make_unique<PairAccelerationsType>()),
  mPairWorkPtr(std::make_unique<PairWorkType>()),
  mSelfAccelerations(FieldStorageType::CopyFields),
  mMassRZ(FieldStorageType::CopyFields),
  mMassDensityRZ(FieldStorageType::CopyFields),
  mDmassDensityDtRZ(FieldStorageType::CopyFields) {
}

//------------------------------------------------------------------------------
// On problem start up, we need to initialize our internal data.
//------------------------------------------------------------------------------
void
SPHRZ::
initializeProblemStartup(DataBase<Dimension>& dataBase) {
  TIME_BEGIN("SPHRZInitializeStartup");
  SPHBase<Dimension>::initializeProblemStartup(dataBase);
  mMassRZ = dataBase.newFluidFieldList(0.0, HydroFieldNames::massRZ);
  mMassDensityRZ = dataBase.newFluidFieldList(0.0, HydroFieldNames::massDensityRZ);
  mDmassDensityDtRZ = dataBase.newFluidFieldList(0.0, IncrementBoundedState<Dimension, Scalar>::prefix() + HydroFieldNames::massDensityRZ);
  TIME_END("SPHRZInitializeStartup");
}

//------------------------------------------------------------------------------
// On problem start up, we need to initialize our internal data.
//------------------------------------------------------------------------------
void
SPHRZ::
initializeProblemStartupDependencies(DataBase<Dimension>& dataBase,
                                     State<Dimension>& state,
                                     StateDerivatives<Dimension>& derivs) {
  TIME_BEGIN("SPHRZInitializeStartupDependencies");
  SPHBase<Dimension>::initializeProblemStartupDependencies(dataBase, state, derivs);
  dataBase.resizeFluidFieldList(mMassRZ, 0.0, HydroFieldNames::massRZ);
  dataBase.resizeFluidFieldList(mMassDensityRZ, 0.0, HydroFieldNames::massDensityRZ);
  dataBase.resizeFluidFieldList(mDmassDensityDtRZ, 0.0, IncrementBoundedState<Dimension, Scalar>::prefix() + HydroFieldNames::massDensityRZ);

  // When we come in the initial conditions for mass and density are 2D areal
  // values, so we need to set up our areal and real 3D values appropriately.
  const auto pos = state.fields(HydroFieldNames::position, Vector::zero());
  auto       mass = state.fields(HydroFieldNames::mass, 0.0);
  auto       rho = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto nfields = mass.numFields();
  for (auto k = 0u; k < nfields; ++k) {
    const auto n = mass[k]->numInternalElements();
    for (auto i = 0u; i < n; ++i) {
      CHECK(rho(k,i) > 0.0);
      const auto ri = abs(pos(k,i).y());
      // const auto Ai = mass(k,i)/rho(k,i);
      // const auto Vi = 2.0*M_PI*ri*Ai;
      // const auto di = std::sqrt(Ai);
      // const auto Vi = cylindricalToroidalVolume(di, ri);
      // const auto Ri = std::sqrt(Ai/M_PI);
      // const auto Vi = circularToroidalVolume(Ri, ri);
      mMassRZ(k,i) = mass(k,i);
      mMassDensityRZ(k,i) = rho(k,i);
      mass(k,i) *= 2.0*M_PI*ri;
    }
  }
  TIME_END("SPHRZInitializeStartupDependencies");
}

//------------------------------------------------------------------------------
// Register the state we need/are going to evolve.
//------------------------------------------------------------------------------
void
SPHRZ::
registerState(DataBase<Dim<2>>& dataBase,
              State<Dim<2>>& state) {

  // The base class does most of it.
  SPHBase<Dimension>::registerState(dataBase, state);

  // RZ mass
  state.enroll(mMassRZ);

  // RZ mass density
  for (auto [nodeListi, fluidNodeListPtr]: enumerate(dataBase.fluidNodeListBegin(), dataBase.fluidNodeListEnd())) {
    state.enroll(*mMassDensityRZ[nodeListi], make_policy<IncrementBoundedState<Dimension, Scalar>>(fluidNodeListPtr->rhoMin(),
                                                                                                   fluidNodeListPtr->rhoMax()));
  }

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
SPHRZ::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {
  SPHBase<Dimension>::registerDerivatives(dataBase, derivs);
  derivs.enroll(mDmassDensityDtRZ);
  const auto compatibleEnergy = this->compatibleEnergyEvolution();
  if (compatibleEnergy) {
    const auto& connectivityMap = dataBase.connectivityMap();
    mPairAccelerationsPtr = std::make_unique<PairAccelerationsType>(connectivityMap);
    mPairWorkPtr = std::make_unique<PairWorkType>(connectivityMap);
    // dataBase.resizeFluidFieldList(mSelfAccelerations, Vector::zero(), HydroFieldNames::selfAccelerations, false);
  }
  derivs.enroll(HydroFieldNames::pairAccelerations, *mPairAccelerationsPtr);
  derivs.enroll(HydroFieldNames::pairWork, *mPairWorkPtr);
  // derivs.enroll(mSelfAccelerations);
}

//------------------------------------------------------------------------------
// Stuff that occurs the beginning of a timestep
//------------------------------------------------------------------------------
void
SPHRZ::
preStepInitialize(const DataBase<Dimension>& dataBase, 
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {

  TIME_BEGIN("SPHRZPreStepInitialize");

  if (densityUpdate() == MassDensityType::IntegrateDensity) return;

  // this->applyGhostBoundaries(state, derivs);
  // for (auto* boundPtr: this->boundaryConditions()) boundPtr->finalizeGhostBoundary();

  // We're going to do something expensive to update the mass density, so get ready.
  const auto& connectivityMap = state.connectivityMap();
  const auto  position = state.fields(HydroFieldNames::position, Vector::zero());
  const auto  mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto  massRZ = state.fields(HydroFieldNames::massRZ, 0.0);
  const auto  H = state.fields(HydroFieldNames::H, SymTensor::zero());
  auto        massDensityRZ = state.fields(HydroFieldNames::massDensityRZ, 0.0);


  switch(densityUpdate()) {

  case MassDensityType::RigorousSumDensity:
  case MassDensityType::CorrectedSumDensity:
    {
      computeSPHSumMassDensity(connectivityMap, this->kernel(), mSumMassDensityOverAllNodeLists, position, massRZ, H, massDensityRZ);
      if (densityUpdate() == MassDensityType::CorrectedSumDensity) {
        for (auto* boundPtr: this->boundaryConditions()) boundPtr->applyFieldListGhostBoundary(massDensityRZ);
        for (auto* boundPtr: this->boundaryConditions()) boundPtr->finalizeGhostBoundary();
        correctSPHSumMassDensity(connectivityMap, this->kernel(), mSumMassDensityOverAllNodeLists, position, massRZ, H, massDensityRZ);
      }
    }
    break;

  case MassDensityType::VoronoiCellDensity:
    {
      this->updateVolume(state, false);
      const auto volume = state.fields(HydroFieldNames::volume, 0.0);
      massDensityRZ = massRZ / volume;
    }
    break;

  case MassDensityType::SumVoronoiCellDensity:
    {
      this->updateVolume(state, true);
      const auto volume = state.fields(HydroFieldNames::volume, 0.0);
      computeSumVoronoiCellMassDensity(connectivityMap, this->kernel(), position, massRZ, volume, H, massDensityRZ);
    }
    break;

  default:
    VERIFY2(false, "SPHRZ::preStepInitialize did not handle a density update choice : " << static_cast<int>(densityUpdate()));
    break;
  }

  // Update the real mass density based on the areal (RZ) density
  auto massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  for (auto k = 0u; k < massDensity.numFields(); ++k) {
    const auto n = massDensity[k]->numInternalElements();
    for (auto i = 0u; i < n; ++i) {
      CHECK(massDensityRZ(k,i) > 0.0);
      const auto ri = abs(position(k,i).y());
      CHECK2(ri > 0.0, "Bad position for node " << i << " : " << position(k,i));
      const auto Ai = massRZ(k,i)/massDensityRZ(k,i);
      const auto Vi = 2.0*M_PI*ri*Ai;
      // const auto di = std::sqrt(massRZ(k,i)/massDensityRZ(k,i));
      // const auto Vi = cylindricalToroidalVolume(di, ri);
      // // const auto Ri = std::sqrt(massRZ(k,i)/(M_PI*massDensityRZ(k,i)));
      // // const auto Vi = circularToroidalVolume(Ri, ri);
      massDensity(k,i) = mass(k,i)/Vi;
    }
  }
  for (auto* boundPtr: this->boundaryConditions()) {
    boundPtr->applyFieldListGhostBoundary(massDensityRZ);
    boundPtr->applyFieldListGhostBoundary(massDensity);
  }
  for (auto* boundPtr: this->boundaryConditions()) boundPtr->finalizeGhostBoundary();

  TIME_END("SPHRZPreStepInitialize");
}

//------------------------------------------------------------------------------
// Determine the principle derivatives.
//------------------------------------------------------------------------------
void
SPHRZ::
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
SPHRZ::
evaluateDerivativesImpl(const Dim<2>::Scalar time,
                        const Dim<2>::Scalar dt,
                        const DataBase<Dim<2>>& dataBase,
                        const State<Dim<2>>& state,
                        StateDerivatives<Dim<2>>& derivs,
                        const QType& Q) const {

  using QPiType = typename QType::ReturnType;

  // The kernels and such.
  const auto& W = this->kernel();
  const auto& WQ = this->PiKernel();
  const auto  oneKernel = (W == WQ);

  // A few useful constants we'll use in the following loop.
  // const auto tiny = 1.0e-30;
  const Scalar W0 = W(0.0, 1.0);
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
  const auto massRZ = state.fields(HydroFieldNames::massRZ, 0.0);
  const auto position = state.fields(HydroFieldNames::position, Vector::zero());
  const auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero());
  const auto massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto massDensityRZ = state.fields(HydroFieldNames::massDensityRZ, 0.0);
  const auto H = state.fields(HydroFieldNames::H, SymTensor::zero());
  const auto pressure = state.fields(HydroFieldNames::pressure, 0.0);
  const auto soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  const auto omega = state.fields(HydroFieldNames::omegaGradh, 0.0);
  const auto fClQ = state.fields(HydroFieldNames::ArtificialViscousClMultiplier, 0.0, true);
  const auto fCqQ = state.fields(HydroFieldNames::ArtificialViscousCqMultiplier, 0.0, true);
  const auto DvDxQ = state.fields(HydroFieldNames::ArtificialViscosityVelocityGradient, Tensor::zero(), true);
  CHECK(mass.size() == numNodeLists);
  CHECK(massRZ.size() == numNodeLists);
  CHECK(position.size() == numNodeLists);
  CHECK(velocity.size() == numNodeLists);
  CHECK(massDensity.size() == numNodeLists);
  CHECK(massDensityRZ.size() == numNodeLists);
  CHECK(H.size() == numNodeLists);
  CHECK(pressure.size() == numNodeLists);
  CHECK(soundSpeed.size() == numNodeLists);
  CHECK(omega.size() == numNodeLists);
  CHECK(fClQ.size() == 0 or fClQ.size() == numNodeLists);
  CHECK(fCqQ.size() == 0 or fCqQ.size() == numNodeLists);
  CHECK(DvDxQ.size() == 0 or DvDxQ.size() == numNodeLists);

  // Derivative FieldLists.
  auto  rhoSum = derivs.fields(ReplaceState<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, 0.0);
  auto  normalization = derivs.fields(HydroFieldNames::normalization, 0.0);
  auto  DxDt = derivs.fields(IncrementState<Dimension, Vector>::prefix() + HydroFieldNames::position, Vector::zero());
  auto  DrhoDt = derivs.fields(IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, 0.0);
  auto  DrhoDtRZ = derivs.fields(IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::massDensityRZ, 0.0);
  auto  DvDt = derivs.fields(HydroFieldNames::hydroAcceleration, Vector::zero());
  auto  DepsDt = derivs.fields(IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
  auto  DvDx = derivs.fields(HydroFieldNames::velocityGradient, Tensor::zero());
  auto  localDvDx = derivs.fields(HydroFieldNames::internalVelocityGradient, Tensor::zero());
  auto  M = derivs.fields(HydroFieldNames::M_SPHCorrection, Tensor::zero());
  auto  localM = derivs.fields("local " + HydroFieldNames::M_SPHCorrection, Tensor::zero());
  auto  maxViscousPressure = derivs.fields(HydroFieldNames::maxViscousPressure, 0.0);
  auto  effViscousPressure = derivs.fields(HydroFieldNames::effectiveViscousPressure, 0.0);
  auto& pairAccelerations = derivs.template get<PairAccelerationsType>(HydroFieldNames::pairAccelerations);
  auto& pairWork = derivs.template get<PairWorkType>(HydroFieldNames::pairWork);
  // auto  selfAccelerations = derivs.fields(HydroFieldNames::selfAccelerations, Vector::zero(), true);
  auto  XSPHWeightSum = derivs.fields(HydroFieldNames::XSPHWeightSum, 0.0);
  auto  XSPHDeltaV = derivs.fields(HydroFieldNames::XSPHDeltaV, Vector::zero());
  CHECK(rhoSum.size() == numNodeLists);
  CHECK(normalization.size() == numNodeLists);
  CHECK(DxDt.size() == numNodeLists);
  CHECK(DrhoDt.size() == numNodeLists);
  CHECK(DrhoDtRZ.size() == numNodeLists);
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
  CHECK((compatibleEnergy and pairAccelerations.size() == npairs) or not compatibleEnergy);
  CHECK((compatibleEnergy and pairWork.size() == npairs) or not compatibleEnergy);
  // CHECK((compatibleEnergy     and selfAccelerations.size() == numNodeLists) or
  //       (not compatibleEnergy and selfAccelerations.size() == 0u));

  // Walk all the interacting pairs.
#pragma omp parallel
  {
    // Thread private scratch variables
    int i, j, nodeListi, nodeListj;
    Scalar Wi, gWi, WQi, gWQi, Wj, gWj, WQj, gWQj;
    Vector gradWi, gradWj, gradWQi, gradWQj;
    Scalar Qi, Qj;
    QPiType QPiij, QPiji;

    typename SpheralThreads<Dim<2>>::FieldListStack threadStack;
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
      const auto& posi = position(nodeListi, i);
      const auto  ri = std::abs(posi.y());
      const auto  mi = mass(nodeListi, i);
      const auto  mRZi = massRZ(nodeListi, i);
      const auto& vi = velocity(nodeListi, i);
      const auto  vri = vi.y();
      const auto  rhoi = massDensity(nodeListi, i);
      const auto  rhoRZi = massDensityRZ(nodeListi, i);
      const auto  Pi = pressure(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto  ci = soundSpeed(nodeListi, i);
      // const auto& omegai = omega(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();
      const auto  zetai = abs((Hi*posi).y());
      const auto  hri = ri*safeInv(zetai);
      const auto  riInv = safeInvVar(ri, 0.1*hri);
      // const auto  safeOmegai = safeInv(omegai, tiny);
      // const auto  Ai = mRZi/rhoRZi;
      // const auto  zetai = abs((Hi*posi).y());
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
      const auto& posj = position(nodeListj, j);
      const auto  rj = std::abs(posj.y());
      const auto  mj = mass(nodeListj, j);
      const auto  mRZj = massRZ(nodeListj, j);
      const auto& vj = velocity(nodeListj, j);
      const auto  vrj = vj.y();
      const auto  rhoj = massDensity(nodeListj, j);
      const auto  rhoRZj = massDensityRZ(nodeListj, j);
      const auto  Pj = pressure(nodeListj, j);
      const auto& Hj = H(nodeListj, j);
      const auto  cj = soundSpeed(nodeListj, j);
      // const auto& omegaj = omega(nodeListj, j);
      const auto  Hdetj = Hj.Determinant();
      const auto  zetaj = abs((Hj*posj).y());
      const auto  hrj = rj*safeInv(zetaj);
      const auto  rjInv = safeInvVar(rj, 0.1*hrj);
      // const auto  safeOmegaj = safeInv(omegaj, tiny);
      // const auto  Aj = mRZj/rhoRZj;
      // const auto  zetaj = abs((Hj*posj).y());
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
      const auto xij = posi - posj;
      const auto etai = Hi*xij;
      const auto etaj = Hj*xij;
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
      if (oneKernel) {
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

      // Contribution to the sum density.
      if (nodeListi == nodeListj) {
        rhoSumi += mRZj*Wi;
        rhoSumj += mRZi*Wj;
        normi += mRZi/rhoi*Wi;
        normj += mRZj/rhoj*Wj;
      }

      // Compute the pair-wise artificial viscosity.
      const auto vij = vi - vj;
      Q.QPiij(QPiij, QPiji, Qi, Qj,
              nodeListi, i, nodeListj, j,
              posi, Hi, etai, vi, rhoRZi, ci,  
              posj, Hj, etaj, vj, rhoRZj, cj,
              fClQ, fCqQ, DvDxQ); 
      const auto Qacci = 0.5*(QPiij*gradWQi);
      const auto Qaccj = 0.5*(QPiji*gradWQj);
      const auto workQi = 0.5*(QPiij*vij).dot(gradWQi);// - Qi*rhoRZi/rhoi*vri*riInv;
      const auto workQj = 0.5*(QPiji*vij).dot(gradWQj);// - Qj*rhoRZj/rhoj*vrj*rjInv;
      // const auto workQi = vij.dot(Qacci);
      // const auto workQj = vij.dot(Qaccj);
      maxViscousPressurei = max(maxViscousPressurei, Qi);
      maxViscousPressurej = max(maxViscousPressurej, Qj);
      effViscousPressurei += mRZj*Qi*WQi/rhoj;
      effViscousPressurej += mRZi*Qj*WQj/rhoi;

      // Acceleration.
      CHECK(rhoi > 0.0);
      CHECK(rhoj > 0.0);
      const auto Prhoi = Pi/(rhoi*rhoRZi);
      const auto Prhoj = Pj/(rhoj*rhoRZj);
      // const auto deltaDvDti = mRZj*(0.5*(Pj + Qj)*(gradWi + gradWj)/rhoRZj);// + Qacci + Qaccj);
      // const auto deltaDvDtj = mRZi*(0.5*(Pi + Qi)*(gradWi + gradWj)/rhoRZi);// + Qacci + Qaccj);
      // const auto deltaDvDti = mRZj*(rhoRZi/rhoi * (Pi/(rhoRZi*rhoRZi)*gradWi + Pj/(rhoRZj*rhoRZj)*gradWj) + Qacci + Qaccj);
      // const auto deltaDvDtj = mRZi*(rhoRZj/rhoj * (Pi/(rhoRZi*rhoRZi)*gradWi + Pj/(rhoRZj*rhoRZj)*gradWj) + Qacci + Qaccj);
      // const auto deltaDvDti = mRZj*(Prhoj*gradWj + Pi*rhoj/(rhoi*rhoi*rhoRZj)*gradWi + Qacci + Qaccj);
      // const auto deltaDvDtj = mRZi*(Prhoi*gradWi + Pj*rhoi/(rhoj*rhoj*rhoRZi)*gradWj + Qacci + Qaccj);
      const auto deltaDvDt = Prhoi*gradWi + Prhoj*gradWj + Qacci + Qaccj;
      const auto deltaDvDti = mRZj*deltaDvDt;
      const auto deltaDvDtj = mRZi*deltaDvDt;
      DvDti -= deltaDvDti;
      DvDtj += deltaDvDtj;

      // Specific thermal energy evolution.
      // const auto worki = mRZj*(Pi/(rhoi*rhoi)*vij.dot(gradWi) + workQi);
      // const auto workj = mRZi*(Pj/(rhoj*rhoj)*vij.dot(gradWj) + workQj);
      // const auto worki = mRZj*(Pi/(rhoi*rhoRZj)*vij.dot(gradWi));// + workQi);
      // const auto workj = mRZi*(Pj/(rhoj*rhoRZi)*vij.dot(gradWj));// + workQj);
      const auto worki = mRZj*(Prhoi*vij.dot(gradWi) + workQi);
      const auto workj = mRZi*(Prhoj*vij.dot(gradWj) + workQj);
      DepsDti += worki;
      DepsDtj += workj;

      // Update the history for compatible energy update
      if (compatibleEnergy) {
        pairAccelerations[kk][0] = -deltaDvDti;
        pairAccelerations[kk][1] =  deltaDvDtj;
        pairWork[kk][0] = worki;
        pairWork[kk][1] = workj;
      }

      // Velocity gradient.
      const auto deltaDvDxi = mRZj*vij.dyad(gradWi);
      const auto deltaDvDxj = mRZi*vij.dyad(gradWj);
      DvDxi -= deltaDvDxi;
      DvDxj -= deltaDvDxj;
      if (sameMatij) {
        localDvDxi -= deltaDvDxi;
        localDvDxj -= deltaDvDxj;
      }

      // Estimate of delta v (for XSPH).
      if (sameMatij) {// or min(zetai, zetaj) < 1.0) {
        const auto wXSPHij = 0.5*(mRZi/rhoi*Wi + mRZj/rhoj*Wj);
        XSPHWeightSumi += wXSPHij;
        XSPHWeightSumj += wXSPHij;
        XSPHDeltaVi -= wXSPHij*vij;
        XSPHDeltaVj += wXSPHij*vij;
      }

      // Linear gradient correction term.
      Mi -= mRZj*xij.dyad(gradWi);
      Mj -= mRZi*xij.dyad(gradWj);
      if (sameMatij) {
        localMi -= mRZj*xij.dyad(gradWi);
        localMj -= mRZi*xij.dyad(gradWj);
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
      const auto& posi = position(nodeListi, i);
      const auto  ri = posi.y();                  // Can be negative for ghost points!
      const auto  mi = mass(nodeListi, i);
      const auto  mRZi = massRZ(nodeListi, i);
      const auto& vi = velocity(nodeListi, i);
      const auto  rhoi = massDensity(nodeListi, i);
      const auto  rhoRZi = massDensityRZ(nodeListi, i);
      const auto  Pi = pressure(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();
      const auto  zetai = (Hi*posi).y();
      const auto  hri = ri*safeInv(zetai);
      CHECK(hri >= 0.0);
      const auto  riInv = safeInvVar(ri, 0.1*hri);
      // const auto  riInv = safeInv(ri, tiny);
      const auto  numNeighborsi = connectivityMap.numNeighborsForNode(nodeListi, i);
      // const auto  Ai = mRZi/rhoRZi;
      // const auto  Vi = mi/rhoi;
      // const auto  riInv = 2.0*M_PI*Ai/Vi;
      // const auto  riInv = safeInv(ri, 0.1*std::sqrt(Ai));
      CHECK(rhoi > 0.0);
      CHECK(rhoRZi > 0.0);
      CHECK(Hdeti > 0.0);

      // auto& rhoSumi = rhoSum(nodeListi, i);
      auto& normi = normalization(nodeListi, i);
      auto& DxDti = DxDt(nodeListi, i);
      auto& DrhoDti = DrhoDt(nodeListi, i);
      auto& DrhoDtRZi = DrhoDtRZ(nodeListi, i);
      auto& DvDti = DvDt(nodeListi, i);
      auto& DepsDti = DepsDt(nodeListi, i);
      auto& DvDxi = DvDx(nodeListi, i);
      auto& localDvDxi = localDvDx(nodeListi, i);
      auto& Mi = M(nodeListi, i);
      auto& localMi = localM(nodeListi, i);
      auto& XSPHWeightSumi = XSPHWeightSum(nodeListi, i);
      auto& XSPHDeltaVi = XSPHDeltaV(nodeListi, i);

      // // Add the self-contribution to density sum.
      // rhoSumi += mRZi*W0*Hdeti;
      // rhoSumi /= circi;
      normi += mRZi*W0*Hdeti;

      // // Finish the acceleration -- self hoop strain.
      // const Vector deltaDvDti(0.0, Pi/rhoRZi*riInv);
      // DvDti += deltaDvDti;
      // if (compatibleEnergy) selfAccelerations(nodeListi, i) = deltaDvDti;

      // Finish the gradient of the velocity.
      CHECK(rhoi > 0.0);
      if (correctVelocityGradient and
          std::abs(Mi.Determinant()) > 1.0e-10 and
          numNeighborsi > Dimension::pownu(2)) {
        Mi = Mi.Inverse();
        DvDxi = DvDxi*Mi;
      } else {
        DvDxi /= rhoRZi;
      }
      if (correctVelocityGradient and
          std::abs(localMi.Determinant()) > 1.0e-10 and
          numNeighborsi > Dimension::pownu(2)) {
        localMi = localMi.Inverse();
        localDvDxi = localDvDxi*localMi;
      } else {
        localDvDxi /= rhoRZi;
      }

      // Finish the continuity equation.
      XSPHWeightSumi += Hdeti*mRZi*W0;
      CHECK2(XSPHWeightSumi != 0.0, i << " " << XSPHWeightSumi);
      XSPHDeltaVi /= XSPHWeightSumi;
      const auto  vri = vi.y(); // + 0.5*dt*DvDti.y();
      DrhoDtRZi = -rhoRZi*DvDxi.Trace();
      DrhoDti = -rhoi*(DvDxi.Trace() + vri*riInv);

      // Finish the specific thermal energy evolution.
      DepsDti -= Pi/rhoi*vri*riInv;

      // If needed finish the total energy derivative.
      if (evolveTotalEnergy) DepsDti = mi*(vi.dot(DvDti) + DepsDti);

      // Determine the position evolution, based on whether we're doing XSPH or not.
      if (XSPH) {
        DxDti = vi + XSPHDeltaVi;
      } else {
        DxDti = vi;
      }
    }
  }
}

//------------------------------------------------------------------------------
// Finalize derivatives
//------------------------------------------------------------------------------
void
SPHRZ::
finalizeDerivatives(const Scalar time,
                    const Scalar dt,
                    const DataBase<Dim<2>>& dataBase,
                    const State<Dim<2>>& state,
                    StateDerivatives<Dim<2>>& derivs) const {

  // If we're using compatible energy mode we need to apply BCs to DepsDt
  const auto compatibleEnergy = this->compatibleEnergyEvolution();
  if (compatibleEnergy) {
    auto DvDt = derivs.fields(HydroFieldNames::hydroAcceleration, Vector::zero());
    auto DepsDt = derivs.fields(IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
    for (auto* bptr: this->boundaryConditions()) {
      bptr->applyFieldListGhostBoundary(DvDt);
      bptr->applyFieldListGhostBoundary(DepsDt);
    }
    for (auto* bptr: this->boundaryConditions()) bptr->finalizeGhostBoundary();
  }

}

//------------------------------------------------------------------------------
// Apply the ghost boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
void
SPHRZ::
applyGhostBoundaries(State<Dim<2>>& state,
                     StateDerivatives<Dim<2>>& derivs) {

  // Our state
  auto massRZ = state.fields(HydroFieldNames::massRZ, 0.0);
  auto rhoRZ = state.fields(HydroFieldNames::massDensityRZ, 0.0);
  for (auto boundaryPtr: this->boundaryConditions()) {
    boundaryPtr->applyFieldListGhostBoundary(massRZ);
    boundaryPtr->applyFieldListGhostBoundary(rhoRZ);
  }

  // // Convert the mass to mass/length before BCs are applied.
  // // const auto pos = state.fields(HydroFieldNames::position, Vector::zero());
  // auto mass = state.fields(HydroFieldNames::mass, 0.0);
  // mass /= massRZ;
  // // const auto numNodeLists = mass.numFields();
  // // for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
  // //   const auto n = mass[nodeListi]->numElements();
  // //   for (auto i = 0u; i < n; ++i) {
  // //     const Scalar circi = 2.0*M_PI*abs(pos(nodeListi, i).y());
  // //     CHECK(circi > 0.0);
  // //     mass(nodeListi, i) /= circi;
  // //   }
  // // }

  // Apply ordinary SPH BCs.
  SPHBase<Dim<2>>::applyGhostBoundaries(state, derivs);
  // for (auto boundaryPtr: this->boundaryConditions()) boundaryPtr->finalizeGhostBoundary();

  // // Scale back to mass.
  // mass *= massRZ;
  // // for (unsigned nodeListi = 0; nodeListi < numNodeLists; ++nodeListi) {
  // //   const unsigned n = mass[nodeListi]->numElements();
  // //   for (unsigned i = 0; i < n; ++i) {
  // //     const Scalar circi = 2.0*M_PI*abs(pos(nodeListi, i).y());
  // //     CHECK(circi > 0.0);
  // //     mass(nodeListi, i) *= circi;
  // //   }
  // // }

}

//------------------------------------------------------------------------------
// Enforce the boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
void
SPHRZ::
enforceBoundaries(State<Dim<2>>& state,
                  StateDerivatives<Dim<2>>& derivs) {

  // Our state
  auto massRZ = state.fields(HydroFieldNames::massRZ, 0.0);
  auto rhoRZ = state.fields(HydroFieldNames::massDensityRZ, 0.0);
  for (auto boundaryPtr: this->boundaryConditions()) {
    boundaryPtr->enforceFieldListBoundary(massRZ);
    boundaryPtr->enforceFieldListBoundary(rhoRZ);
  }

  // // Convert the mass to mass/length before BCs are applied.
  // FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  // FieldList<Dimension, Vector> pos = state.fields(HydroFieldNames::position, Vector::zero());
  // const unsigned numNodeLists = mass.numFields();
  // for (unsigned nodeListi = 0; nodeListi < numNodeLists; ++nodeListi) {
  //   const unsigned n = mass[nodeListi]->numElements();
  //   for (unsigned i = 0; i < n; ++i) {
  //     const Scalar circi = 2.0*M_PI*abs(pos(nodeListi, i).y());
  //     CHECK(circi > 0.0);
  //     mass(nodeListi, i) /= circi;
  //   }
  // }

  // Apply ordinary SPH BCs.
  SPHBase<Dim<2>>::enforceBoundaries(state, derivs);

  // // Scale back to mass.
  // // We also ensure no point approaches the z-axis too closely.
  // FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero());
  // for (unsigned nodeListi = 0; nodeListi < numNodeLists; ++nodeListi) {
  //   const unsigned n = mass[nodeListi]->numElements();
  //   for (unsigned i = 0; i < n; ++i) {
  //     Vector& posi = pos(nodeListi, i);
  //     const Scalar circi = 2.0*M_PI*abs(posi.y());
  //     mass(nodeListi, i) *= circi;
  //   }
  // }
}

}
