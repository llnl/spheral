//---------------------------------Spheral++----------------------------------//
// VoronoiCells
//
// Computes polytopes for each point similar to the Voronoi tessellation
//----------------------------------------------------------------------------//
#include "VoronoiCells/VoronoiCells.hh"
#include "VoronoiCells/computeVoronoiVolume.hh"
#include "VoronoiCells/UpdateVoronoiCells.hh"
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
VoronoiCells<Dimension>::
VoronoiCells(const Scalar kernelExtent,
             const vector<FacetedVolume>& facetedBoundaries,
             const vector<vector<FacetedVolume>>& facetedHoles):
  mEtaMax(kernelExtent),
  mVolume(FieldStorageType::CopyFields),
  mWeight(FieldStorageType::CopyFields),
  mSurfacePoint(FieldStorageType::CopyFields),
  mEtaVoidPoints(FieldStorageType::CopyFields),
  mCells(FieldStorageType::CopyFields),
  mCellFaceFlags(FieldStorageType::CopyFields),
  mDeltaCentroid(FieldStorageType::CopyFields),
  mFacetedBoundaries(facetedBoundaries),
  mFacetedHoles(facetedHoles),
  mRestart(registerWithRestart(*this)) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
VoronoiCells<Dimension>::
~VoronoiCells() {
}

//------------------------------------------------------------------------------
// Size up our FieldLists on problem startup
//------------------------------------------------------------------------------
template<typename Dimension>
void
VoronoiCells<Dimension>::
initializeProblemStartup(DataBase<Dimension>& dataBase) {
  mVolume = dataBase.newFluidFieldList(0.0, HydroFieldNames::volume);
  // mWeight = dataBase.newFluidFieldList(0.0, "Voronoi weight");
  mSurfacePoint = dataBase.newFluidFieldList(0, HydroFieldNames::surfacePoint);
  mEtaVoidPoints = dataBase.newFluidFieldList(std::vector<Vector>(), HydroFieldNames::etaVoidPoints);
  mCells = dataBase.newFluidFieldList(FacetedVolume(), HydroFieldNames::cells);
  mCellFaceFlags = dataBase.newFluidFieldList(std::vector<CellFaceFlag>(), HydroFieldNames::cellFaceFlags);
  mDeltaCentroid = dataBase.newFluidFieldList(Vector::zero, "delta centroid");
}

//------------------------------------------------------------------------------
// On problem initialization we need to compute the cells.  Onace a calculation
// is going we rely on the cells being updated at the end of the prior step.
//------------------------------------------------------------------------------
template<typename Dimension>
void
VoronoiCells<Dimension>::
initializeProblemStartupDependencies(DataBase<Dimension>& dataBase,
                                     State<Dimension>& state,
                                     StateDerivatives<Dimension>& derivs) {

  // Ensure our state is sized correctly
  dataBase.resizeFluidFieldList(mVolume, 0.0, HydroFieldNames::volume, false);
  // dataBase.resizeFluidFieldList(mWeight, 0.0, "Voronoi weight", false);
  dataBase.resizeFluidFieldList(mSurfacePoint, 0, HydroFieldNames::surfacePoint, false);
  dataBase.resizeFluidFieldList(mEtaVoidPoints, vector<Vector>(), HydroFieldNames::etaVoidPoints, false);
  dataBase.resizeFluidFieldList(mCells, FacetedVolume(), HydroFieldNames::cells, false);
  dataBase.resizeFluidFieldList(mCellFaceFlags, vector<CellFaceFlag>(), HydroFieldNames::cellFaceFlags, false);
  dataBase.resizeFluidFieldList(mDeltaCentroid, Vector::zero, "delta centroid", false);

  // Use our finalize method to compute the cell geometry
  this->preStepInitialize(dataBase, state, derivs);

  // Propagate our state to constant any ghost nodes
  for (auto* boundaryPtr: this->boundaryConditions()) boundaryPtr->initializeProblemStartup(false);
}

//------------------------------------------------------------------------------
// Register the state
//------------------------------------------------------------------------------
template<typename Dimension>
void
VoronoiCells<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {
  state.enroll(mVolume);
  state.enroll(mSurfacePoint);
  state.enroll(mCellFaceFlags);
  state.enroll(mCells);
  // state.enroll(mCells, make_policy<UpdateVoronoiCells<Dimension>>(mVolume,
  //                                                                 mWeight,
  //                                                                 mDeltaCentroid,
  //                                                                 mEtaVoidPoints,
  //                                                                 this->boundaryConditions(),
  //                                                                 mFacetedBoundaries,
  //                                                                 mFacetedHoles));
}

//------------------------------------------------------------------------------
// No derivatives to register
//------------------------------------------------------------------------------
template<typename Dimension>
void
VoronoiCells<Dimension>::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {
}

//------------------------------------------------------------------------------
// Apply the ghost boundary conditions
//------------------------------------------------------------------------------
template<typename Dimension>
void
VoronoiCells<Dimension>::
applyGhostBoundaries(State<Dimension>& state,
                     StateDerivatives<Dimension>& derivs) {
  auto cells = state.template fields<FacetedVolume>(HydroFieldNames::cells);
  auto vol = state.fields(HydroFieldNames::volume, 0.0);
  auto surfacePoint = state.fields(HydroFieldNames::surfacePoint, 0);
  for (auto* boundaryPtr: this->boundaryConditions()) {
    boundaryPtr->applyFieldListGhostBoundary(cells);
    boundaryPtr->applyFieldListGhostBoundary(vol);
    boundaryPtr->applyFieldListGhostBoundary(surfacePoint);
  }
}

//------------------------------------------------------------------------------
// Enforce the boundary conditions
//------------------------------------------------------------------------------
template<typename Dimension>
void
VoronoiCells<Dimension>::
enforceBoundaries(State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {
  auto cells = state.template fields<FacetedVolume>(HydroFieldNames::cells);
  auto vol = state.fields(HydroFieldNames::volume, 0.0);
  auto surfacePoint = state.fields(HydroFieldNames::surfacePoint, 0);
  for (auto* boundaryPtr: this->boundaryConditions()) {
    boundaryPtr->enforceFieldListBoundary(cells);
    boundaryPtr->enforceFieldListBoundary(vol);
    boundaryPtr->enforceFieldListBoundary(surfacePoint);
  } 
}

//------------------------------------------------------------------------------
// No time step vote
//------------------------------------------------------------------------------
template<typename Dimension>
typename VoronoiCells<Dimension>::TimeStepType
VoronoiCells<Dimension>::
dt(const DataBase<Dimension>& /*dataBase*/, 
   const State<Dimension>& /*state*/,
   const StateDerivatives<Dimension>& /*derivs*/,
   const Scalar /*currentTime*/) const {
  return std::make_pair(std::numeric_limits<double>::max(), std::string("VoronoiCells: no vote"));
}

//------------------------------------------------------------------------------
// No derivatives to evaluate
//------------------------------------------------------------------------------
template<typename Dimension>
void
VoronoiCells<Dimension>::
evaluateDerivatives(const Scalar /*time*/,
                    const Scalar /*dt*/,
                    const DataBase<Dimension>& /*dataBase*/,
                    const State<Dimension>& /*state*/,
                    StateDerivatives<Dimension>& /*derivatives*/) const {
}

//------------------------------------------------------------------------------
// Initialize at the start of a physics cycle.
// This is when we do the expensive operation of computing the Voronoi cell
// geometry from scratch.
//------------------------------------------------------------------------------
template<typename Dimension>
void
VoronoiCells<Dimension>::
preStepInitialize(const DataBase<Dimension>& dataBase, 
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {

  // State we need to compute the Voronoi cells
  const auto& cm = state.connectivityMap();
  const auto  pos = state.fields(HydroFieldNames::position, Vector::zero);
  const auto  H = state.fields(HydroFieldNames::H, SymTensor::zero);
  // const auto  mass = state.fields(HydroFieldNames::mass, 0.0);
  // const auto  rho = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto  D = state.fields(SolidFieldNames::tensorDamage, SymTensor::zero, true);
  auto& boundaries = this->boundaryConditions();

//   // Use m/rho to estimate our weighting to roughly match cell volumes
//   const auto numNodeLists = dataBase.numFluidNodeLists();
//   for (auto k = 0u; k < numNodeLists; ++k) {
//     const auto n = mass[k]->numInternalElements();
// #pragma omp parallel for
//     for (auto i = 0u; i < n; ++i) {
//       CHECK(rho(k,i) > 0.0);
//       mVolume(k,i) = mass(k,i)/rho(k,i);
//     }
//   }
  
//   // Enforce boundaries on the volume
//   for (auto* bcPtr: boundaries) bcPtr->applyFieldListGhostBoundary(mVolume);
//   for (auto* bcPtr: boundaries) bcPtr->finalizeGhostBoundary();

//   // We can now compute the weights from our volumes (including ghosts)
//   for (auto k = 0u; k < numNodeLists; ++k) {
//     const auto n = mass[k]->numElements();    // ghosts as well!
// #pragma omp parallel for
//     for (auto i = 0u; i < n; ++i) {
//       CHECK(mVolume(k,i) > 0.0);
//       mWeight(k,i) = 1.0/Dimension::rootnu(mVolume(k,i));
//     }
//   }

  // Compute the cell data.  Note we are using the fact the state versions of the things
  // we're updating (mSurfacePoint, mCells, etc.) are just pointing at our internal fields.
  computeVoronoiVolume(pos, H, cm, D, mFacetedBoundaries, mFacetedHoles, boundaries, mWeight,
                       mSurfacePoint, mVolume, mDeltaCentroid, mEtaVoidPoints, mCells, mCellFaceFlags);
}

//------------------------------------------------------------------------------
// Provide a hook to be called after the state has been updated and 
// boundary conditions have been enforced.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
VoronoiCells<Dimension>::
postStateUpdate(const Scalar time,
                const Scalar dt,
                const DataBase<Dimension>& dataBase, 
                State<Dimension>& state,
                StateDerivatives<Dimension>& derivs) {
  this->preStepInitialize(dataBase, state, derivs);
  return true;
}

//------------------------------------------------------------------------------
// Add a faceted boundary
//------------------------------------------------------------------------------
template<typename Dimension>
void
VoronoiCells<Dimension>::
addFacetedBoundary(const FacetedVolume& bound,
                   const std::vector<FacetedVolume>& holes) {
  const auto numExisting = mFacetedBoundaries.size();
  for (auto i = 0u; i < numExisting; ++i) {
    if (bound == mFacetedBoundaries[i] and holes == mFacetedHoles[i]) {
      std::cerr << "tried to add same faceted boundary twice" << std::endl;
      return;
    }
  }
  mFacetedBoundaries.push_back(bound);
  mFacetedHoles.push_back(holes);
}

//------------------------------------------------------------------------------
// Dump the current state to the given file
//------------------------------------------------------------------------------
template<typename Dimension>
void
VoronoiCells<Dimension>::
dumpState(FileIO& file, const std::string& pathName) const {
  // file.write(mVolume, pathName + "/Voronoi_volume");
  // file.write(mWeight, pathName + "/weight");
  // file.write(mSurfacePoint, pathName + "/surfacePoint");
  // file.write(mCellFaceFlags, pathName + "/cellFaceFlags");
}

//------------------------------------------------------------------------------
// Restore the state from the given file
//------------------------------------------------------------------------------
template<typename Dimension>
void
VoronoiCells<Dimension>::
restoreState(const FileIO& file, const std::string& pathName) {
}

} // end namespace Spheral
