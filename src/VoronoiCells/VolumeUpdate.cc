//---------------------------------Spheral++----------------------------------//
// VolumeUpdate
//
// Base class for volume computation packages. Handles non-Voronoi volume
// types (MassOverDensity, SumVolume, HullVolume, HVolume). VoronoiCells
// inherits from this for full Voronoi geometry computation.
//----------------------------------------------------------------------------//
#include "VoronoiCells/VolumeUpdate.hh"
#include "VoronoiCells/GeometryScaling.hh"
#include "Boundary/Boundary.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "FileIO/FileIO.hh"
#include "Geometry/Dimension.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Kernel/TableKernel.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "RK/computeRKSumVolume.hh"
#include "RK/computeHullVolumes.hh"
#include "RK/computeHVolumes.hh"
#include "Utilities/DBC.hh"

#include <limits>

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename Dimension>
VolumeUpdate<Dimension>::
VolumeUpdate(const VolumeType volumeType,
             const TableKernel<Dimension>& W,
             const bool updateInStep,
             const bool updateInFinalize):
  mVolume(FieldStorageType::CopyFields),
  mVolume3d(FieldStorageType::CopyFields),
  mVolumeType(volumeType),
  mUpdateInStep(updateInStep),
  mUpdateInFinalize(updateInFinalize),
  mKernel(W),
  mRestart(registerWithRestart(*this)) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
VolumeUpdate<Dimension>::
~VolumeUpdate() {
}

//------------------------------------------------------------------------------
// Compute volumes based on volume type
//------------------------------------------------------------------------------
template<typename Dimension>
void
VolumeUpdate<Dimension>::
computeVolume(const DataBase<Dimension>& dataBase,
              State<Dimension>& state) {
  auto vol = state.fields(HydroFieldNames::volume, 0.0);
  auto vol3d = state.fields(HydroFieldNames::volume3d, 0.0);
  const auto mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto H = state.fields(HydroFieldNames::H, SymTensor::zero());
  const auto position = state.fields(HydroFieldNames::position, Vector::zero());

  switch (mVolumeType) {
  case VolumeType::MassOverDensity:
    CHECK(mass.size() == massDensity.size());
    // mass/rho gives V_phys in non-Cartesian; unscale to get V_coord
    vol3d.assignFields(mass / massDensity);
    unscaleFromGeometry(position, vol3d, vol);
    break;

  case VolumeType::SumVolume:
    {
      const auto& connectivityMap = dataBase.connectivityMap();
      computeRKSumVolume(connectivityMap, mKernel, position, mass, H, vol);
      scaleForGeometry(position, vol, vol3d);
    }
    break;

  case VolumeType::HullVolume:
    {
      const auto& connectivityMap = dataBase.connectivityMap();
      computeHullVolumes(connectivityMap, mKernel.kernelExtent(), position, H, vol);
      scaleForGeometry(position, vol, vol3d);
    }
    break;

  case VolumeType::HVolume:
    {
      const auto nPerh = vol.nodeListPtrs()[0]->nodesPerSmoothingScale();
      computeHVolumes(nPerh, H, vol);
      scaleForGeometry(position, vol, vol3d);
    }
    break;

  case VolumeType::VoronoiVolume:
    VERIFY2(false, "VoronoiVolume requires the VoronoiCells package, not VolumeUpdate");
    break;
  }
}

//------------------------------------------------------------------------------
// Size up our FieldLists on problem startup
//------------------------------------------------------------------------------
template<typename Dimension>
void
VolumeUpdate<Dimension>::
initializeProblemStartup(DataBase<Dimension>& dataBase) {
  mVolume = dataBase.newFluidFieldList(0.0, HydroFieldNames::volume);
  mVolume3d = dataBase.newFluidFieldList(0.0, HydroFieldNames::volume3d);
}

//------------------------------------------------------------------------------
// Compute initial volumes
//------------------------------------------------------------------------------
template<typename Dimension>
void
VolumeUpdate<Dimension>::
initializeProblemStartupDependencies(DataBase<Dimension>& dataBase,
                                     State<Dimension>& state,
                                     StateDerivatives<Dimension>& derivs) {
  dataBase.resizeFluidFieldList(mVolume, 0.0, HydroFieldNames::volume, false);
  dataBase.resizeFluidFieldList(mVolume3d, 0.0, HydroFieldNames::volume3d, false);
  computeVolume(dataBase, state);
  for (auto* boundaryPtr: this->boundaryConditions()) {
    boundaryPtr->initializeProblemStartup(false);
  }

  // Apply boundaries
  this->applyGhostBoundaries(state, derivs);
  for (auto* boundaryPtr: this->boundaryConditions()) {
    boundaryPtr->finalizeGhostBoundary();
  }
}

//------------------------------------------------------------------------------
// Register the state
//------------------------------------------------------------------------------
template<typename Dimension>
void
VolumeUpdate<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {
  state.enroll(mVolume);
  state.enroll(mVolume3d);
  
  // Stuff VolumeUpdate needs that might have been enrolled elsewhere
  auto position = dataBase.fluidPosition();
  auto mass = dataBase.fluidMass();
  auto massDensity = dataBase.fluidMassDensity();
  auto H = dataBase.fluidHfield();
  if (not state.registered(position)) state.enroll(position);
  if (not state.registered(mass)) state.enroll(mass);
  if (not state.registered(massDensity)) state.enroll(massDensity);
  if (not state.registered(H)) state.enroll(H);
}

//------------------------------------------------------------------------------
// No derivatives to register
//------------------------------------------------------------------------------
template<typename Dimension>
void
VolumeUpdate<Dimension>::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {
}

//------------------------------------------------------------------------------
// Apply the ghost boundary conditions.
// Volume is stored as coordinate-plane values and applied directly.
// Mass carries a geometric factor; unscale before BCs, rescale afterward.
//------------------------------------------------------------------------------
template<typename Dimension>
void
VolumeUpdate<Dimension>::
applyGhostBoundaries(State<Dimension>& state,
                     StateDerivatives<Dimension>& derivs) {
  auto vol = state.fields(HydroFieldNames::volume, 0.0);
  auto vol3d = state.fields(HydroFieldNames::volume3d, 0.0);
  auto mass = state.fields(HydroFieldNames::mass, 0.0);
  auto rho = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto pos = state.fields(HydroFieldNames::position, Vector::zero());

  {
    auto guard = unscaledRegion(pos, mass, vol3d);
    for (auto* boundaryPtr: this->boundaryConditions()) {
      boundaryPtr->applyFieldListGhostBoundary(vol);
      boundaryPtr->applyFieldListGhostBoundary(vol3d);
      boundaryPtr->applyFieldListGhostBoundary(mass);
      boundaryPtr->applyFieldListGhostBoundary(rho);
    }
    if (guard.applied()) {
      for (auto* boundaryPtr: this->boundaryConditions()) {
        boundaryPtr->finalizeGhostBoundary();
      }
    }
  }
}

//------------------------------------------------------------------------------
// Enforce the boundary conditions (internal nodes only).
//------------------------------------------------------------------------------
template<typename Dimension>
void
VolumeUpdate<Dimension>::
enforceBoundaries(State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {
  auto vol = state.fields(HydroFieldNames::volume, 0.0);
  auto vol3d = state.fields(HydroFieldNames::volume3d, 0.0);
  auto mass = state.fields(HydroFieldNames::mass, 0.0);
  auto rho = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto pos = state.fields(HydroFieldNames::position, Vector::zero());

  {
    auto guard = unscaledInternalRegion(pos, mass, vol3d);
    for (auto* boundaryPtr: this->boundaryConditions()) {
      boundaryPtr->enforceFieldListBoundary(vol);
      boundaryPtr->enforceFieldListBoundary(vol3d);
      boundaryPtr->enforceFieldListBoundary(mass);
      boundaryPtr->enforceFieldListBoundary(rho);
    }
  }
}

//------------------------------------------------------------------------------
// No time step vote
//------------------------------------------------------------------------------
template<typename Dimension>
typename VolumeUpdate<Dimension>::TimeStepType
VolumeUpdate<Dimension>::
dt(const DataBase<Dimension>& /*dataBase*/,
   const State<Dimension>& /*state*/,
   const StateDerivatives<Dimension>& /*derivs*/,
   const Scalar /*currentTime*/) const {
  return std::make_pair(std::numeric_limits<double>::max(), std::string("VolumeUpdate: no vote"));
}

//------------------------------------------------------------------------------
// No derivatives to evaluate
//------------------------------------------------------------------------------
template<typename Dimension>
void
VolumeUpdate<Dimension>::
evaluateDerivatives(const Scalar /*time*/,
                    const Scalar /*dt*/,
                    const DataBase<Dimension>& /*dataBase*/,
                    const State<Dimension>& /*state*/,
                    StateDerivatives<Dimension>& /*derivatives*/) const {
}

//------------------------------------------------------------------------------
// Recompute volumes at start of step
//------------------------------------------------------------------------------
template<typename Dimension>
void
VolumeUpdate<Dimension>::
preStepInitialize(const DataBase<Dimension>& dataBase,
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {
  if (mUpdateInStep) {
    computeVolume(dataBase, state);
    this->applyGhostBoundaries(state, derivs);
    for (auto* boundaryPtr: this->boundaryConditions()) {
      boundaryPtr->finalizeGhostBoundary();
    }
  }
}

//------------------------------------------------------------------------------
// Post state update — no-op (volumes recomputed at start of step)
//------------------------------------------------------------------------------
template<typename Dimension>
bool
VolumeUpdate<Dimension>::
postStateUpdate(const Scalar /*time*/,
                const Scalar /*dt*/,
                const DataBase<Dimension>& dataBase,
                State<Dimension>& state,
                StateDerivatives<Dimension>& derivs) {
  return true;
}

//------------------------------------------------------------------------------
// Implicit finalize — recompute volumes for implicit consumers
//------------------------------------------------------------------------------
template<typename Dimension>
void
VolumeUpdate<Dimension>::
finalize(const Scalar /*time*/,
                 const Scalar /*dt*/,
                 DataBase<Dimension>& dataBase,
                 State<Dimension>& state,
                 StateDerivatives<Dimension>& derivs) {
  if (mUpdateInFinalize) {
    computeVolume(dataBase, state);
    this->applyGhostBoundaries(state, derivs);
    for (auto* boundaryPtr: this->boundaryConditions()) {
      boundaryPtr->finalizeGhostBoundary();
    }
  }
}

//------------------------------------------------------------------------------
// Dump the current state to the given file
//------------------------------------------------------------------------------
template<typename Dimension>
void
VolumeUpdate<Dimension>::
dumpState(FileIO& file, const std::string& pathName) const {
  file.write(mVolume, pathName + "/Volume");
  file.write(mVolume3d, pathName + "/PhysicalVolume");
}

//------------------------------------------------------------------------------
// Restore the state from the given file
//------------------------------------------------------------------------------
template<typename Dimension>
void
VolumeUpdate<Dimension>::
restoreState(const FileIO& file, const std::string& pathName) {
  file.read(mVolume, pathName + "/Volume");
  file.read(mVolume3d, pathName + "/PhysicalVolume");
}

} // end namespace Spheral
