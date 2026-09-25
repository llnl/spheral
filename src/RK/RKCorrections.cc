//---------------------------------Spheral++----------------------------------//
// RKCorrections
//
// Computes RK corrections for other physics packages.
//----------------------------------------------------------------------------//
#include "RK/RKCorrections.hh"
#include "RK/RKUtilities.hh"
#include "RK/RKFieldNames.hh"
#include "Boundary/Boundary.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "FileIO/FileIO.hh"
#include "Geometry/Dimension.hh"
#include "Kernel/TableKernel.hh"
#include "Hydro/HydroFieldNames.hh"

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
              const bool needHessian,
              const bool updateInStep,
              const bool updateInFinalize):
  mOrders(orders),
  mDataBase(dataBase),
  mNeedHessian(needHessian),
  mUpdateInStep(updateInStep),
  mUpdateInFinalize(updateInFinalize),
  mWR(),
  mSurfaceArea(FieldStorageType::CopyFields),
  mNormal(FieldStorageType::CopyFields),
  mCorrections(),
  mRestart(registerWithRestart(*this)) {

  mOrders.insert(RKOrder::ZerothOrder);     // We always at least want ZerothOrder
  for (auto order: mOrders) {
    mWR.emplace(std::make_pair(order, ReproducingKernel<Dimension>(W, order)));
    mCorrections.emplace(std::make_pair(order, FieldList<Dimension, RKCoefficients<Dimension>>(FieldStorageType::CopyFields)));
  }

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
// Allocate correction fields on problem startup
//------------------------------------------------------------------------------
template<typename Dimension>
void
RKCorrections<Dimension>::
initializeProblemStartup(DataBase<Dimension>& dataBase) {
  mSurfaceArea = dataBase.newFluidFieldList(0.0, HydroFieldNames::surfaceArea);
  mNormal = dataBase.newFluidFieldList(Vector::zero(), HydroFieldNames::normal);
  for (auto order: mOrders) {
    mCorrections[order] = dataBase.newFluidFieldList(RKCoefficients<Dimension>(), RKFieldNames::rkCorrections(order));
  }
  ENSURE(mWR.size() == mOrders.size());
  ENSURE(mCorrections.size() == mOrders.size());
}

//------------------------------------------------------------------------------
// Compute initial corrections
//------------------------------------------------------------------------------
template<typename Dimension>
void
RKCorrections<Dimension>::
initializeProblemStartupDependencies(DataBase<Dimension>& dataBase,
                                     State<Dimension>& state,
                                     StateDerivatives<Dimension>& derivs) {
  const auto& connectivityMap = dataBase.connectivityMap();
  const auto  H = state.fields(HydroFieldNames::H, SymTensor::zero());
  const auto  position = state.fields(HydroFieldNames::position, Vector::zero());
  const auto  volume = state.fields(HydroFieldNames::volume, 0.0);

  // Compute corrections
  for (auto order: mOrders) {
    if (mOrders.size() == 1 or order != RKOrder::ZerothOrder) {
      mWR[order].computeCorrections(connectivityMap, volume, position, H,
                                    mNeedHessian, mCorrections[RKOrder::ZerothOrder], mCorrections[order]);
    }
  }

  // Boundaries may need to be reinitialized.
  for (auto boundItr = this->boundaryBegin(); boundItr < this->boundaryEnd(); ++boundItr) (*boundItr)->initializeProblemStartup(false);

  // Apply boundaries to corrections before computing normal
  this->applyGhostBoundaries(state, derivs);
  for (auto boundItr = this->boundaryBegin(); boundItr < this->boundaryEnd(); ++boundItr) {
    (*boundItr)->finalizeGhostBoundary();
  }

  // Compute normal direction
  // mWR[RKOrder::ZerothOrder].computeNormal(connectivityMap, volume, position, H,
  //                                         mCorrections[RKOrder::ZerothOrder], mSurfaceArea, mNormal);
}

//------------------------------------------------------------------------------
// Register the state
//------------------------------------------------------------------------------
template<typename Dimension>
void
RKCorrections<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {
  // RK-owned state
  state.enroll(RKFieldNames::rkOrders, mOrders);
  for (auto order: mOrders) {
    state.enroll(RKFieldNames::reproducingKernel(order), mWR[order]);
    state.enroll(mCorrections[order]);
  }
  state.enroll(mSurfaceArea);
  state.enroll(mNormal);
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
  auto surfaceArea = state.fields(HydroFieldNames::surfaceArea, 0.0);
  auto normal = state.fields(HydroFieldNames::normal, Vector::zero());

  for (auto* boundaryPtr: this->boundaryConditions()) {
    boundaryPtr->applyFieldListGhostBoundary(surfaceArea);
    boundaryPtr->applyFieldListGhostBoundary(normal);
    for (auto order: mOrders) {
      auto corrections = state.fields(RKFieldNames::rkCorrections(order), RKCoefficients<Dimension>());
      boundaryPtr->applyFieldListGhostBoundary(corrections);
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
  auto surfaceArea = state.fields(HydroFieldNames::surfaceArea, 0.0);
  auto normal = state.fields(HydroFieldNames::normal, Vector::zero());

  for (auto* boundaryPtr: this->boundaryConditions()) {
    boundaryPtr->enforceFieldListBoundary(surfaceArea);
    boundaryPtr->enforceFieldListBoundary(normal);
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
// Recompute RK corrections from current state
//------------------------------------------------------------------------------
template<typename Dimension>
void
RKCorrections<Dimension>::
updateCorrections(const DataBase<Dimension>& dataBase,
                  State<Dimension>& state) {
  const auto& connectivityMap = dataBase.connectivityMap();
  const auto  H = state.fields(HydroFieldNames::H, SymTensor::zero());
  const auto  position = state.fields(HydroFieldNames::position, Vector::zero());
  const auto  volume = state.fields(HydroFieldNames::volume, 0.0);
  auto        zerothCorrections = state.fields(RKFieldNames::rkCorrections(RKOrder::ZerothOrder), RKCoefficients<Dimension>());

  for (auto order: mOrders) {
    if (mOrders.size() == 1 or order != RKOrder::ZerothOrder) {
      auto corrections = state.fields(RKFieldNames::rkCorrections(order), RKCoefficients<Dimension>());
      mWR[order].computeCorrections(connectivityMap, volume, position, H,
                                    mNeedHessian, zerothCorrections, corrections);
    }
  }
}

//------------------------------------------------------------------------------
// Compute new RK corrections
//------------------------------------------------------------------------------
template<typename Dimension>
bool
RKCorrections<Dimension>::
initialize(const typename Dimension::Scalar /*time*/,
           const typename Dimension::Scalar /*dt*/,
           const DataBase<Dimension>& dataBase,
           State<Dimension>& state,
           StateDerivatives<Dimension>& derivs) {
  if (mUpdateInStep) {
    updateCorrections(dataBase, state);
  }
  return true;
}

//------------------------------------------------------------------------------
// Post state update — no-op (corrections recomputed in initialize)
//------------------------------------------------------------------------------
template<typename Dimension>
bool
RKCorrections<Dimension>::
postStateUpdate(const Scalar /*time*/,
                const Scalar /*dt*/,
                const DataBase<Dimension>& dataBase,
                State<Dimension>& state,
                StateDerivatives<Dimension>& derivs) {
  return true;
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
// Implicit finalize — recompute corrections for implicit consumers
//------------------------------------------------------------------------------
template<typename Dimension>
void
RKCorrections<Dimension>::
finalize(const Scalar /*time*/, 
                 const Scalar /*dt*/,
                 DataBase<Dimension>& dataBase, 
                 State<Dimension>& state,
                 StateDerivatives<Dimension>& derivs) {
  if (mUpdateInFinalize) {
    updateCorrections(dataBase, state);
    this->applyGhostBoundaries(state, derivs);
    for (auto boundItr = this->boundaryBegin(); boundItr < this->boundaryEnd(); ++boundItr) {
      (*boundItr)->finalizeGhostBoundary();
    }
  }
}

//------------------------------------------------------------------------------
// Dump the current state to the given file
//------------------------------------------------------------------------------
template<typename Dimension>
void
RKCorrections<Dimension>::
dumpState(FileIO& file, const std::string& pathName) const {
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
