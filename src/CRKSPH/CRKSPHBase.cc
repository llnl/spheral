//---------------------------------Spheral++----------------------------------//
// CRKSPHBase -- Base class for the CRKSPH/ACRKSPH hydrodynamic packages
//
// Created by JMO, Mon Jul 19 21:52:29 PDT 2010
//----------------------------------------------------------------------------//
#include "CRKSPH/CRKSPHBase.hh"

#include "FileIO/FileIO.hh"
#include "RK/ReproducingKernel.hh"
#include "RK/RKFieldNames.hh"
#include "CRKSPH/computeCRKSPHSumMassDensity.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Hydro/entropyWeightingFunction.hh"
#include "Strength/SolidFieldNames.hh"
#include "Physics/GenericHydro.hh"
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

namespace Spheral {

//------------------------------------------------------------------------------
// Construct with the given artificial viscosity and kernels.
//------------------------------------------------------------------------------
template<typename Dimension>
CRKSPHBase<Dimension>::
CRKSPHBase(DataBase<Dimension>& dataBase,
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
  GenericHydro<Dimension>(Q, cfl, useVelocityMagnitudeForDt),
  mOrder(order),
  mDensityUpdate(densityUpdate),
  mCompatibleEnergyEvolution(compatibleEnergyEvolution),
  mEvolveTotalEnergy(evolveTotalEnergy),
  mXSPH(XSPH),
  mEpsTensile(epsTensile),
  mnTensile(nTensile),
  mTimeStepMask(FieldStorageType::CopyFields),
  mPressure(FieldStorageType::CopyFields),
  mSoundSpeed(FieldStorageType::CopyFields),
  mEntropy(FieldStorageType::CopyFields),
  mXSPHDeltaV(FieldStorageType::CopyFields),
  mDxDt(FieldStorageType::CopyFields),
  mDvDt(FieldStorageType::CopyFields),
  mDmassDensityDt(FieldStorageType::CopyFields),
  mDspecificThermalEnergyDt(FieldStorageType::CopyFields),
  mDvDx(FieldStorageType::CopyFields),
  mInternalDvDx(FieldStorageType::CopyFields),
  mRestart(registerWithRestart(*this)) {

  // Create storage for our internal state.
  mTimeStepMask = dataBase.newFluidFieldList(int(0), HydroFieldNames::timeStepMask);
  mPressure = dataBase.newFluidFieldList(0.0, HydroFieldNames::pressure);
  mSoundSpeed = dataBase.newFluidFieldList(0.0, HydroFieldNames::soundSpeed);
  mEntropy = dataBase.newFluidFieldList(0.0, HydroFieldNames::entropy);
  mXSPHDeltaV = dataBase.newFluidFieldList(Vector::zero(), HydroFieldNames::XSPHDeltaV);
  mDxDt = dataBase.newFluidFieldList(Vector::zero(), IncrementState<Dimension, Vector>::prefix() + HydroFieldNames::position);
  mDvDt = dataBase.newFluidFieldList(Vector::zero(), HydroFieldNames::hydroAcceleration);
  mDmassDensityDt = dataBase.newFluidFieldList(0.0, IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity);
  mDspecificThermalEnergyDt = dataBase.newFluidFieldList(0.0, IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy);
  mDvDx = dataBase.newFluidFieldList(Tensor::zero(), HydroFieldNames::velocityGradient);
  mInternalDvDx = dataBase.newFluidFieldList(Tensor::zero(), HydroFieldNames::internalVelocityGradient);
}

//------------------------------------------------------------------------------
// On problem start up, we need to initialize our internal data.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHBase<Dimension>::
initializeProblemStartupDependencies(DataBase<Dimension>& dataBase,
                                     State<Dimension>& state,
                                     StateDerivatives<Dimension>& derivs) {
  TIME_BEGIN("CRKBaseInitializeProblemStartupDependencies");

  // Initialize the pressure, sound speed, and entropy.
  updateStateFields<Dimension>(HydroFieldNames::pressure, state, derivs);
  updateStateFields<Dimension>(HydroFieldNames::soundSpeed, state, derivs);
  updateStateFields<Dimension>(HydroFieldNames::entropy, state, derivs);
  TIME_END("CRKBaseInitializeProblemStartupDependencies");
}

//------------------------------------------------------------------------------
// Register the state we need/are going to evolve.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHBase<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {
  TIME_BEGIN("CRKBaseRegisterState");

  // Create the local storage for time step mask, pressure, sound speed, and correction fields.
  dataBase.resizeFluidFieldList(mTimeStepMask, 1, HydroFieldNames::timeStepMask);
  dataBase.resizeFluidFieldList(mEntropy,    0.0,                   HydroFieldNames::entropy, false);
  dataBase.resizeFluidFieldList(mPressure,   0.0,                   HydroFieldNames::pressure, false);
  dataBase.resizeFluidFieldList(mSoundSpeed, 0.0,                   HydroFieldNames::soundSpeed, false);

  // We have to choose either compatible or total energy evolution.
  VERIFY2(not (mCompatibleEnergyEvolution and mEvolveTotalEnergy),
          "CRKSPH error : you cannot simultaneously use both compatibleEnergyEvolution and evolveTotalEnergy");

  // Now register away.
  // Mass.
  auto mass = dataBase.fluidMass();
  state.enroll(mass);

  // Volume.
  // Note: we depend on RKCorrections having already registered the volume, but it was registered
  // without an update policy.
  auto vol = state.fields(HydroFieldNames::volume, 0.0);
  CHECK(vol.size() == dataBase.numFluidNodeLists());
  state.enroll(vol, make_policy<ContinuityVolumePolicy<Dimension>>());

  // Mass density (with NodeList dependent limits)
  auto massDensity = dataBase.fluidMassDensity();
  for (auto [nodeListi, fluidNodeListPtr]: enumerate(dataBase.fluidNodeListBegin(), dataBase.fluidNodeListEnd())) {
    state.enroll(*massDensity[nodeListi], make_policy<IncrementBoundedState<Dimension, Scalar>>(fluidNodeListPtr->rhoMin(),
                                                                                                fluidNodeListPtr->rhoMax()));
  }

  // Position
  // We make this dependent on the thermal energy in case we're using the compatible energy update in RZ
  auto position = dataBase.fluidPosition();
  state.enroll(position, make_policy<IncrementState<Dimension, Vector>>({HydroFieldNames::specificThermalEnergy}));

  // Register the velocity
  // We make this dependent on the thermal energy in case we're using the compatible energy update
  auto velocity = dataBase.fluidVelocity();
  state.enroll(velocity, make_policy<IncrementState<Dimension, Vector>>({HydroFieldNames::position,
                                                                         HydroFieldNames::specificThermalEnergy},
                                                                        true));  // Use all DvDt sources (wildcard)

  // Register the entropy.
  state.enroll(mEntropy, make_policy<EntropyPolicy<Dimension>>());

  // Register the time step mask, initialized to 1 so that everything defaults to being
  // checked.
  state.enroll(mTimeStepMask);

  // Register the pressure and sound speed.
  state.enroll(mPressure, make_policy<PressurePolicy<Dimension>>());
  state.enroll(mSoundSpeed, make_policy<SoundSpeedPolicy<Dimension>>());
  TIME_END("CRKBaseRegisterState");
}

//------------------------------------------------------------------------------
// Register the state derivative fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHBase<Dimension>::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {
  TIME_BEGIN("CRKBaseRegisterDerivatives");

  const auto DxDtName = IncrementState<Dimension, Vector>::prefix() + HydroFieldNames::position;
  const auto DvDtName = HydroFieldNames::hydroAcceleration;

  // Create the scratch fields.
  // Note we deliberately do not zero out the derivatives here!  This is because the previous step
  // info here may be used by other algorithms (like the CheapSynchronousRK2 integrator or
  // the ArtificialViscosity::initialize step).
  dataBase.resizeFluidFieldList(mXSPHDeltaV, Vector::zero(), HydroFieldNames::XSPHDeltaV, false);
  dataBase.resizeFluidFieldList(mDxDt, Vector::zero(), IncrementState<Dimension, Vector>::prefix() + HydroFieldNames::position, false);
  dataBase.resizeFluidFieldList(mDvDt, Vector::zero(), HydroFieldNames::hydroAcceleration, false);
  dataBase.resizeFluidFieldList(mDmassDensityDt, 0.0, IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, false);
  dataBase.resizeFluidFieldList(mDspecificThermalEnergyDt, 0.0, IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy, false);
  dataBase.resizeFluidFieldList(mDvDx, Tensor::zero(), HydroFieldNames::velocityGradient, false);
  dataBase.resizeFluidFieldList(mInternalDvDx, Tensor::zero(), HydroFieldNames::internalVelocityGradient, false);

  derivs.enroll(mXSPHDeltaV);

  // These two (the position and velocity updates) may be registered
  // by other physics packages as well, so we need to be careful
  // not to duplicate if so.
  if (not derivs.registered(mDxDt)) derivs.enroll(mDxDt);
  if (not derivs.registered(mDvDt)) derivs.enroll(mDvDt);

  derivs.enroll(mDmassDensityDt);
  derivs.enroll(mDspecificThermalEnergyDt);
  derivs.enroll(mDvDx);
  derivs.enroll(mInternalDvDx);
  TIME_END("CRKBaseRegisterDerivatives");
}

//------------------------------------------------------------------------------
// This method is called once at the beginning of a timestep, after all state registration.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHBase<Dimension>::
preStepInitialize(const DataBase<Dimension>& dataBase, 
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& /*derivs*/) {
  TIME_BEGIN("CRKBasePreStepInitialize");

  // Depending on the mass density advancement selected, we may want to replace the 
  // mass density.
  // Note: we depend on RKCorrections having already updated the volume!
  if (mDensityUpdate == MassDensityType::RigorousSumDensity or
      mDensityUpdate == MassDensityType::VoronoiCellDensity) {
    auto        massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
    const auto& WR = state.template get<ReproducingKernel<Dimension>>(RKFieldNames::reproducingKernel(mOrder));
    const auto& W = WR.kernel();
    const auto& connectivityMap = dataBase.connectivityMap();
    const auto  mass = state.fields(HydroFieldNames::mass, 0.0);
    const auto  H = state.fields(HydroFieldNames::H, SymTensor::zero());
    const auto  position = state.fields(HydroFieldNames::position, Vector::zero());
    const auto  vol = state.fields(HydroFieldNames::volume, 0.0);
    if (densityUpdate() == MassDensityType::RigorousSumDensity) {
      computeCRKSPHSumMassDensity(connectivityMap, W, position, mass, vol, H, massDensity);
    } else {
      massDensity.assignFields(mass/vol);
    }
    for (auto boundaryPtr: range(this->boundaryBegin(), this->boundaryEnd())) boundaryPtr->applyFieldListGhostBoundary(massDensity);
    for (auto boundaryPtr: range(this->boundaryBegin(), this->boundaryEnd())) boundaryPtr->finalizeGhostBoundary();
  }
  TIME_END("CRKBasePreStepInitialize");
}

//------------------------------------------------------------------------------
// Finalize the derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHBase<Dimension>::
finalizeDerivatives(const typename Dimension::Scalar /*time*/,
                    const typename Dimension::Scalar /*dt*/,
                    const DataBase<Dimension>& /*dataBase*/,
                    const State<Dimension>& /*state*/,
                    StateDerivatives<Dimension>& derivs) const {
  TIME_BEGIN("CRKBaseFinalizeDerivatives");

  // If we're using the compatible energy discretization, we need to enforce
  // boundary conditions on the accelerations.
  if (compatibleEnergyEvolution()) {
    auto accelerations = derivs.fields(HydroFieldNames::hydroAcceleration, Vector::zero());
    auto DepsDt = derivs.fields(IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
    for (auto boundaryPtr: range(this->boundaryBegin(), this->boundaryEnd())) {
      boundaryPtr->applyFieldListGhostBoundary(accelerations);
      boundaryPtr->applyFieldListGhostBoundary(DepsDt);
    }
    for (auto boundaryPtr: range(this->boundaryBegin(), this->boundaryEnd())) boundaryPtr->finalizeGhostBoundary();
  }
  TIME_END("CRKBaseFinalizeDerivatives");
}

//------------------------------------------------------------------------------
// Apply the ghost boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHBase<Dimension>::
applyGhostBoundaries(State<Dimension>& state,
                     StateDerivatives<Dimension>& /*derivs*/) {
  TIME_BEGIN("CRKBaseApplyGhostBoundaries");

  // Apply boundary conditions to the basic fluid state Fields.
  // volume, mass, and massDensity handled by RKCorrections
  auto specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero());
  auto pressure = state.fields(HydroFieldNames::pressure, 0.0);
  auto soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  auto entropy = state.fields(HydroFieldNames::entropy, 0.0);

  for (auto boundaryPtr: range(this->boundaryBegin(), this->boundaryEnd())) {
    boundaryPtr->applyFieldListGhostBoundary(specificThermalEnergy);
    boundaryPtr->applyFieldListGhostBoundary(velocity);
    boundaryPtr->applyFieldListGhostBoundary(pressure);
    boundaryPtr->applyFieldListGhostBoundary(soundSpeed);
    boundaryPtr->applyFieldListGhostBoundary(entropy);
  }
  TIME_END("CRKBaseApplyGhostBoundaries");
}

//------------------------------------------------------------------------------
// Enforce the boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHBase<Dimension>::
enforceBoundaries(State<Dimension>& state,
                  StateDerivatives<Dimension>& /*derivs*/) {
  TIME_BEGIN("CRKBaseEnforceBoundaries");

  // Enforce boundary conditions on the fluid state Fields.
  // volume, mass, and massDensity handled by RKCorrections
  auto specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero());
  auto pressure = state.fields(HydroFieldNames::pressure, 0.0);
  auto soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  auto entropy = state.fields(HydroFieldNames::entropy, 0.0);

  for (auto boundaryPtr: range(this->boundaryBegin(), this->boundaryEnd())) {
    boundaryPtr->enforceFieldListBoundary(specificThermalEnergy);
    boundaryPtr->enforceFieldListBoundary(velocity);
    boundaryPtr->enforceFieldListBoundary(pressure);
    boundaryPtr->enforceFieldListBoundary(soundSpeed);
    boundaryPtr->enforceFieldListBoundary(entropy);
  }
  TIME_END("CRKBaseEnforceBoundaries");
}

//------------------------------------------------------------------------------
// Return the RK orders we want to use
//------------------------------------------------------------------------------
template<typename Dimension>
std::set<RKOrder>
CRKSPHBase<Dimension>::
requireReproducingKernels() const {
  return std::set<RKOrder>({RKOrder::ZerothOrder, mOrder});
}

//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHBase<Dimension>::
dumpState(FileIO& file, const string& pathName) const {
  file.write(mTimeStepMask, pathName + "/timeStepMask");
  file.write(mPressure, pathName + "/pressure");
  file.write(mSoundSpeed, pathName + "/soundSpeed");
  file.write(mEntropy, pathName + "/entropy");
  file.write(mXSPHDeltaV, pathName + "/XSPHDeltaV");

  file.write(mDxDt, pathName + "/DxDt");
  file.write(mDvDt, pathName + "/DvDt");
  file.write(mDmassDensityDt, pathName + "/DmassDensityDt");
  file.write(mDspecificThermalEnergyDt, pathName + "/DspecificThermalEnergyDt");
  file.write(mDvDx, pathName + "/DvDx");
  file.write(mInternalDvDx, pathName + "/internalDvDx");
}

//------------------------------------------------------------------------------
// Restore the state from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHBase<Dimension>::
restoreState(const FileIO& file, const string& pathName) {
  file.read(mTimeStepMask, pathName + "/timeStepMask");
  file.read(mPressure, pathName + "/pressure");
  file.read(mSoundSpeed, pathName + "/soundSpeed");
  file.read(mEntropy, pathName + "/entropy");
  file.read(mXSPHDeltaV, pathName + "/XSPHDeltaV");

  file.read(mDxDt, pathName + "/DxDt");
  file.read(mDvDt, pathName + "/DvDt");
  file.read(mDmassDensityDt, pathName + "/DmassDensityDt");
  file.read(mDspecificThermalEnergyDt, pathName + "/DspecificThermalEnergyDt");
  file.read(mDvDx, pathName + "/DvDx");
  file.read(mInternalDvDx, pathName + "/internalDvDx");
}

}

