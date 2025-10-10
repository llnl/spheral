//---------------------------------Spheral++----------------------------------//
// RKCorrections
//
// Computes RK corrections for other physics packages
//----------------------------------------------------------------------------//
#include "RK/RKCorrections.hh"
#include "RK/RKUtilities.hh"
#include "RK/RKFieldNames.hh"
#include "RK/computeRKVolumes.hh"
#include "Boundary/Boundary.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "FileIO/FileIO.hh"
#include "Geometry/Dimension.hh"
#include "Kernel/TableKernel.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Strength/SolidFieldNames.hh"

#include <limits>

namespace Spheral {

using std::vector;

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename Dimension>
RKCorrections<Dimension>::
RKCorrections(const std::set<RKOrder> orders,
              const DataBase<Dimension>& dataBase,
              const TableKernel<Dimension>& W,
              const RKVolumeType volumeType,
              const bool needHessian,
              const bool updateInFinalize):
  mOrders(orders),
  mDataBase(dataBase),
  mVolumeType(volumeType),
  mNeedHessian(needHessian),
  mUpdateInFinalize(updateInFinalize),
  mWR(),
  mVolume(FieldStorageType::CopyFields),
  mSurfaceArea(FieldStorageType::CopyFields),
  mNormal(FieldStorageType::CopyFields),
  mCorrections(),
  mSurfacePoint(FieldStorageType::CopyFields),
  mEtaVoidPoints(FieldStorageType::CopyFields),
  mCells(FieldStorageType::CopyFields),
  mCellFaceFlags(FieldStorageType::CopyFields),
  mRestart(registerWithRestart(*this)) {

  mOrders.insert(RKOrder::ZerothOrder);     // We always at least want ZerothOrder
  for (auto order: mOrders) {
    mWR.emplace(std::make_pair(order, ReproducingKernel<Dimension>(W, order)));   // cause ReproducingKernel is not default constructible
    mCorrections.emplace(std::make_pair(order, FieldList<Dimension, RKCoefficients<Dimension>>(FieldStorageType::CopyFields)));
  }

  mVolume = dataBase.newFluidFieldList(0.0, HydroFieldNames::volume);
  mSurfaceArea = dataBase.newFluidFieldList(0.0, HydroFieldNames::surfaceArea);
  mNormal = dataBase.newFluidFieldList(Vector::zero, HydroFieldNames::normal);
  mSurfacePoint = dataBase.newFluidFieldList(0, HydroFieldNames::surfacePoint);
  mEtaVoidPoints = dataBase.newFluidFieldList(std::vector<Vector>(), HydroFieldNames::etaVoidPoints);
  if (mVolumeType == RKVolumeType::RKVoronoiVolume) {
    mCells = dataBase.newFluidFieldList(FacetedVolume(), HydroFieldNames::cells);
    mCellFaceFlags = dataBase.newFluidFieldList(std::vector<CellFaceFlag>(), HydroFieldNames::cellFaceFlags);
  }
  mDeltaCentroid = dataBase.newFluidFieldList(Vector::zero, "delta centroid");
  
  ENSURE(mWR.size() == mOrders.size());
  ENSURE(mCorrections.size() == mOrders.size());
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
RKCorrections<Dimension>::
~RKCorrections() {
}

//------------------------------------------------------------------------------
// Optional hook to initialize once when the problem is starting up
//------------------------------------------------------------------------------
template<typename Dimension>
void
RKCorrections<Dimension>::
initializeProblemStartup(DataBase<Dimension>& dataBase) {
  for (auto order: mOrders) {
    mCorrections[order] = dataBase.newFluidFieldList(RKCoefficients<Dimension>(), RKFieldNames::rkCorrections(order));
  }
}

//------------------------------------------------------------------------------
// Optional hook to initialize once when the problem is starting up
//------------------------------------------------------------------------------
template<typename Dimension>
void
RKCorrections<Dimension>::
initializeProblemStartupDependencies(DataBase<Dimension>& dataBase,
                                     State<Dimension>& state,
                                     StateDerivatives<Dimension>& derivs) {
  // Get some more data
  const auto& connectivityMap = dataBase.connectivityMap();
  const auto& W = mWR.begin()->second.kernel();
  const auto  mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto  H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const auto  position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto  massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto  damage = state.fields(SolidFieldNames::tensorDamage, SymTensor::zero, true);
  
  // Compute the volumes
  computeRKVolumes(connectivityMap, W,
                   position, mass, massDensity, H, damage,
                   mFacetedBoundaries, mFacetedHoles, this->boundaryConditions(), mVolumeType,
                   mSurfacePoint, mDeltaCentroid, mEtaVoidPoints, mCells, mCellFaceFlags,
                   mVolume);
  
  // Propagate volume to constant ghost nodes
  for (auto boundItr = this->boundaryBegin(); boundItr < this->boundaryEnd(); ++boundItr) (*boundItr)->initializeProblemStartup(false);

  // Apply boundaries to newly computed terms
  for (auto boundItr = this->boundaryBegin(); boundItr < this->boundaryEnd(); ++boundItr) {
    (*boundItr)->applyFieldListGhostBoundary(mVolume);
    if (mVolumeType == RKVolumeType::RKVoronoiVolume) {
      (*boundItr)->applyFieldListGhostBoundary(mSurfacePoint);
      (*boundItr)->applyFieldListGhostBoundary(mEtaVoidPoints);
    }
  }
  for (auto boundItr = this->boundaryBegin(); boundItr < this->boundaryEnd(); ++boundItr) {
    (*boundItr)->finalizeGhostBoundary();
  }
  
  // // Allocate correction fields
  // for (auto order: mOrders) mCorrections[order] = dataBase.newFluidFieldList(RKCoefficients<Dimension>(), RKFieldNames::rkCorrections(order));
  
  // Compute corrections
  for (auto order: mOrders) {
    if (mOrders.size() == 1 or order != RKOrder::ZerothOrder) {  // Zeroth order is always computed anyway
      mWR[order].computeCorrections(connectivityMap, mVolume, position, H,
                                    mNeedHessian, mCorrections[RKOrder::ZerothOrder], mCorrections[order]);
    }
  }

  // Boundaries may need to be reinitialized.
  for (auto boundItr = this->boundaryBegin(); boundItr < this->boundaryEnd(); ++boundItr) (*boundItr)->initializeProblemStartup(false);

  // Apply boundaries to corrections before computing normal
  for (auto boundItr = this->boundaryBegin(); boundItr < this->boundaryEnd(); ++boundItr) {
    for (auto order: mOrders) {
      (*boundItr)->applyFieldListGhostBoundary(mCorrections[order]);
    }
  }
  for (auto boundItr = this->boundaryBegin(); boundItr < this->boundaryEnd(); ++boundItr) {
    (*boundItr)->finalizeGhostBoundary();
  }

  // Compute normal direction
  mWR[RKOrder::ZerothOrder].computeNormal(connectivityMap, mVolume, position, H,
                                          mCorrections[RKOrder::ZerothOrder], mSurfaceArea, mNormal);
}

//------------------------------------------------------------------------------
// Register the state
//------------------------------------------------------------------------------
template<typename Dimension>
void
RKCorrections<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {
  // Stuff RKCorrections owns
  state.enroll(RKFieldNames::rkOrders, mOrders);
  for (auto order: mOrders) {
    state.enroll(RKFieldNames::reproducingKernel(order), mWR[order]);
    state.enroll(mCorrections[order]);
  }
  state.enroll(mVolume);
  state.enroll(mSurfaceArea);
  state.enroll(mNormal);
  
  state.enroll(mSurfacePoint);
  state.enroll(mEtaVoidPoints);
  if (mVolumeType == RKVolumeType::RKVoronoiVolume) {
    state.enroll(mCells);
    state.enroll(mCellFaceFlags);
  }

  // Stuff RKCorrections needs that might have been enrolled elsewhere
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
RKCorrections<Dimension>::
registerDerivatives(DataBase<Dimension>& /*dataBase*/,
                    StateDerivatives<Dimension>& /*derivs*/) {
}

//------------------------------------------------------------------------------
// Apply the ghost boundary conditions
//------------------------------------------------------------------------------
template<typename Dimension>
void
RKCorrections<Dimension>::
applyGhostBoundaries(State<Dimension>& state,
                     StateDerivatives<Dimension>& derivs) {
  CONTRACT_VAR(derivs);
  // Get state variables
  auto vol = state.fields(HydroFieldNames::volume, 0.0);
  auto mass = state.fields(HydroFieldNames::mass, 0.0);
  auto massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  auto surfaceArea = state.fields(HydroFieldNames::surfaceArea, 0.0);
  auto normal = state.fields(HydroFieldNames::normal, Vector::zero);
  auto surfacePoint = state.fields(HydroFieldNames::surfacePoint, 0);
  auto etaVoidPoints = state.fields(HydroFieldNames::etaVoidPoints, std::vector<Vector>());
  auto zerothCorrections = state.fields(RKFieldNames::rkCorrections(RKOrder::ZerothOrder), RKCoefficients<Dimension>());

  // Apply ghost boundary conditions
  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->applyFieldListGhostBoundary(vol);
    (*boundaryItr)->applyFieldListGhostBoundary(mass);
    (*boundaryItr)->applyFieldListGhostBoundary(massDensity);
    (*boundaryItr)->applyFieldListGhostBoundary(surfaceArea);
    (*boundaryItr)->applyFieldListGhostBoundary(normal);
    (*boundaryItr)->applyFieldListGhostBoundary(surfacePoint);
    (*boundaryItr)->applyFieldListGhostBoundary(etaVoidPoints);
    for (auto order: mOrders) {
      auto corrections = state.fields(RKFieldNames::rkCorrections(order), RKCoefficients<Dimension>());
      (*boundaryItr)->applyFieldListGhostBoundary(corrections);
    }
  }
}

//------------------------------------------------------------------------------
// Enforce the boundary conditions
//------------------------------------------------------------------------------
template<typename Dimension>
void
RKCorrections<Dimension>::
enforceBoundaries(State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {
  CONTRACT_VAR(derivs);
  // Get state variables
  auto vol = state.fields(HydroFieldNames::volume, 0.0);
  auto mass = state.fields(HydroFieldNames::mass, 0.0);
  auto massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  auto surfaceArea = state.fields(HydroFieldNames::surfaceArea, 0.0);
  auto normal = state.fields(HydroFieldNames::normal, Vector::zero);
  auto surfacePoint = state.fields(HydroFieldNames::surfacePoint, 0);
  auto etaVoidPoints = state.fields(HydroFieldNames::etaVoidPoints, std::vector<Vector>());

  // Enforce boundary conditions
  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->enforceFieldListBoundary(vol);
    (*boundaryItr)->enforceFieldListBoundary(mass);
    (*boundaryItr)->enforceFieldListBoundary(massDensity);
    (*boundaryItr)->enforceFieldListBoundary(surfaceArea);
    (*boundaryItr)->enforceFieldListBoundary(normal);
    (*boundaryItr)->enforceFieldListBoundary(surfacePoint);
  } 
}

//------------------------------------------------------------------------------
// No time step vote
//------------------------------------------------------------------------------
template<typename Dimension>
typename RKCorrections<Dimension>::TimeStepType
RKCorrections<Dimension>::
dt(const DataBase<Dimension>& /*dataBase*/, 
   const State<Dimension>& /*state*/,
   const StateDerivatives<Dimension>& /*derivs*/,
   const Scalar /*currentTime*/) const {
  return std::make_pair(std::numeric_limits<double>::max(), std::string("RKCorrections: no vote"));
}

//------------------------------------------------------------------------------
// Compute new volumes
//------------------------------------------------------------------------------
template<typename Dimension>
void
RKCorrections<Dimension>::
preStepInitialize(const DataBase<Dimension>& dataBase, 
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {
  CONTRACT_VAR(derivs);
  // Get data
  const auto& W = mWR.begin()->second.kernel();
  const auto& connectivityMap = dataBase.connectivityMap();
  const auto  mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto  H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const auto  position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto  damage = state.fields(SolidFieldNames::tensorDamage, SymTensor::zero, true);
  const auto  massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  auto        volume = state.fields(HydroFieldNames::volume, 0.0);
  auto        surfacePoint = state.fields(HydroFieldNames::surfacePoint, 0);
  FieldList<Dimension, FacetedVolume> cells;
  FieldList<Dimension, std::vector<CellFaceFlag>> cellFaceFlags;
  if (mVolumeType == RKVolumeType::RKVoronoiVolume) {
    cells = state.fields(HydroFieldNames::cells, FacetedVolume());
    cellFaceFlags = state.fields(HydroFieldNames::cellFaceFlags, std::vector<CellFaceFlag>());
  }
  
  // Compute volumes
  computeRKVolumes(connectivityMap, W,
                   position, mass, massDensity, H, damage,
                   mFacetedBoundaries, mFacetedHoles, this->boundaryConditions(), mVolumeType,
                   surfacePoint, mDeltaCentroid, mEtaVoidPoints, cells, cellFaceFlags,
                   volume);
  
  // Apply ghost boundaries to Voronoi stuff
  for (auto boundItr = this->boundaryBegin(); boundItr < this->boundaryEnd(); ++boundItr) {
    (*boundItr)->applyFieldListGhostBoundary(volume);
    if (mVolumeType == RKVolumeType::RKVoronoiVolume) {
      (*boundItr)->applyFieldListGhostBoundary(cells);
      (*boundItr)->applyFieldListGhostBoundary(cellFaceFlags);
      (*boundItr)->applyFieldListGhostBoundary(surfacePoint);
      (*boundItr)->applyFieldListGhostBoundary(mEtaVoidPoints);
    }
  }
  for (auto boundItr = this->boundaryBegin(); boundItr < this->boundaryEnd(); ++boundItr) {
    (*boundItr)->finalizeGhostBoundary();
  }
}

//------------------------------------------------------------------------------
// Compute new RK corrections
//------------------------------------------------------------------------------
template<typename Dimension>
bool
RKCorrections<Dimension>::
initialize(const typename Dimension::Scalar time,
           const typename Dimension::Scalar dt,
           const DataBase<Dimension>& dataBase,
           State<Dimension>& state,
           StateDerivatives<Dimension>& derivs) {
  CONTRACT_VAR(time);
  CONTRACT_VAR(dt);
  CONTRACT_VAR(derivs);
  // Get data
  const auto& connectivityMap = dataBase.connectivityMap();
  const auto  H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const auto  position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto  volume = state.fields(HydroFieldNames::volume, 0.0);
  auto        zerothCorrections = state.fields(RKFieldNames::rkCorrections(RKOrder::ZerothOrder), RKCoefficients<Dimension>());
  auto        surfaceArea = state.fields(HydroFieldNames::surfaceArea, 0.0);
  auto        normal = state.fields(HydroFieldNames::normal, Vector::zero);
  
  // Compute corrections
  for (auto order: mOrders) {
    if (mOrders.size() == 1 or order != RKOrder::ZerothOrder) {  // Zeroth order is always computed anyway
      auto corrections = state.fields(RKFieldNames::rkCorrections(order), RKCoefficients<Dimension>());
      mWR[order].computeCorrections(connectivityMap, volume, position, H,
                                    mNeedHessian, zerothCorrections, corrections);
    }
  }
  
  // Apply ghost boundaries to corrections
  return true;

  // for (auto boundaryItr = this->boundaryBegin(); boundaryItr < this->boundaryEnd(); ++boundaryItr) {
  //   (*boundaryItr)->applyFieldListGhostBoundary(zerothCorrections);
  //   (*boundaryItr)->applyFieldListGhostBoundary(surfaceArea);
  //   (*boundaryItr)->applyFieldListGhostBoundary(normal);
  //   for (auto order: mOrders) {
  //     if (order != RKOrder::ZerothOrder) {
  //       auto corrections = state.fields(RKFieldNames::rkCorrections(order), RKCoefficients<Dimension>());
  //       (*boundaryItr)->applyFieldListGhostBoundary(corrections);
  //     }
  //   }
  // }
  // for (auto boundItr = this->boundaryBegin(); boundItr < this->boundaryEnd(); ++boundItr) {
  //   (*boundItr)->finalizeGhostBoundary();
  // }
}

//------------------------------------------------------------------------------
// No derivatives to evaluate
//------------------------------------------------------------------------------
template<typename Dimension>
void
RKCorrections<Dimension>::
evaluateDerivatives(const Scalar /*time*/,
                    const Scalar /*dt*/,
                    const DataBase<Dimension>& /*dataBase*/,
                    const State<Dimension>& /*state*/,
                    StateDerivatives<Dimension>& /*derivatives*/) const {
}

//------------------------------------------------------------------------------
// Nothing to finalize
//------------------------------------------------------------------------------
template<typename Dimension>
void
RKCorrections<Dimension>::
finalize(const Scalar time, 
         const Scalar dt,
         DataBase<Dimension>& dataBase, 
         State<Dimension>& state,
         StateDerivatives<Dimension>& derivs) {
  if (mUpdateInFinalize) {
    // Calculate new volumes
    preStepInitialize(dataBase, state, derivs);
    
    // Calculate new corrections
    initialize(time, dt, dataBase, state, derivs);
  }
}

//------------------------------------------------------------------------------
// Add a faceted boundary
//------------------------------------------------------------------------------
template<typename Dimension>
void
RKCorrections<Dimension>::
addFacetedBoundary(const FacetedVolume& cell,
                   const std::vector<FacetedVolume>& holes) {
  const auto numExisting = mFacetedBoundaries.size();
  for (auto i = 0u; i < numExisting; ++i) {
    if (cell == mFacetedBoundaries[i] and holes == mFacetedHoles[i]) {
      std::cerr << "tried to add same faceted boundary twice" << std::endl;
      return;
    }
  }
  mFacetedBoundaries.push_back(cell);
  mFacetedHoles.push_back(holes);
}

//------------------------------------------------------------------------------
// Dump the current state to the given file
//------------------------------------------------------------------------------
template<typename Dimension>
void
RKCorrections<Dimension>::
dumpState(FileIO& file, const std::string& pathName) const {
  file.write(mVolume, pathName + "/Volume");
  file.write(mSurfacePoint, pathName + "/surfacePoint");
  for (const auto& corr: mCorrections) {
    file.write(corr.second, pathName + "/" + RKFieldNames::rkCorrections(corr.first));
  }
}

//------------------------------------------------------------------------------
// Restore the state from the given file
//------------------------------------------------------------------------------
template<typename Dimension>
void
RKCorrections<Dimension>::
restoreState(const FileIO& file, const std::string& pathName) {
  file.read(mVolume, pathName + "/Volume");
  file.read(mSurfacePoint, pathName + "/surfacePoint");
  for (auto order: mOrders) {
    file.read(mCorrections[order], pathName + "/" + RKFieldNames::rkCorrections(order));
  }
}

//------------------------------------------------------------------------------
// WR
//------------------------------------------------------------------------------
template<typename Dimension>
const ReproducingKernel<Dimension>&
RKCorrections<Dimension>::
WR(const RKOrder order) const {
  const auto itr = mWR.find(order);
  VERIFY2(itr != mWR.end(),
          "RKCorrections::WR error: attempt to access for unknown correction");
  return itr->second;
}

//------------------------------------------------------------------------------
// corrections
//------------------------------------------------------------------------------
template<typename Dimension>
const FieldList<Dimension, RKCoefficients<Dimension>>&
RKCorrections<Dimension>::
corrections(const RKOrder order) const {
  const auto itr = mCorrections.find(order);
  VERIFY2(itr != mCorrections.end(),
          "RKCorrections::corrections error: attempt to access for unknown correction");
  return itr->second;
}

} // end namespace Spheral
