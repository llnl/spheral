//---------------------------------Spheral++----------------------------------//
// VoronoiCells
//
// Computes polytopes for each point similar to the Voronoi tessellation.
// Inherits from VolumeUpdate; overrides computeVolume for Voronoi geometry.
//----------------------------------------------------------------------------//
#include "VoronoiCells/VoronoiCells.hh"
#include "VoronoiCells/GeometryScaling.hh"
#include "VoronoiCells/computeVoronoiVolume.hh"
#include "VoronoiCells/UpdateVoronoiCells.hh"
#include "Boundary/Boundary.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "FileIO/FileIO.hh"
#include "Geometry/Dimension.hh"
#include "Geometry/GeometryRegistrar.hh"
#include "Kernel/TableKernel.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Strength/SolidFieldNames.hh"
#include "Utilities/SpheralMessage.hh"

#include <limits>
#include <algorithm>

namespace Spheral {

using std::vector;

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename Dimension>
VoronoiCells<Dimension>::
VoronoiCells(const VolumeType volumeType,
             const TableKernel<Dimension>& W,
             const vector<FacetedVolume>& facetedBoundaries,
             const vector<vector<FacetedVolume>>& facetedHoles,
             const bool updateInStep,
             const bool updateInFinalize):
  VolumeUpdate<Dimension>(volumeType,
                          W,
                          updateInStep,
                          updateInFinalize),
  mEtaMax(W.kernelExtent()),
  mWeight(FieldStorageType::CopyFields),
  mSurfacePoint(FieldStorageType::CopyFields),
  mEtaVoidPoints(FieldStorageType::CopyFields),
  mCells(FieldStorageType::CopyFields),
  mCellFaceFlags(FieldStorageType::CopyFields),
  mDeltaCentroid(FieldStorageType::CopyFields),
  mFacetedBoundaries(facetedBoundaries),
  mFacetedHoles(facetedHoles) {
  if (facetedHoles.empty()) mFacetedHoles.resize(mFacetedBoundaries.size());
  ENSURE(mFacetedBoundaries.size() == mFacetedHoles.size());
}

//------------------------------------------------------------------------------
// Compute the Voronoi cell geometry (overrides VolumeUpdate::computeVolume).
// Always computes the full Voronoi tessellation (cells, surfacePoint, etc.)
// but overwrites the volume field if the user chose a non-Voronoi volume type.
//------------------------------------------------------------------------------
template<typename Dimension>
void
VoronoiCells<Dimension>::
computeVolume(const DataBase<Dimension>& dataBase,
              State<Dimension>& state) {
  const auto& cm = state.connectivityMap();
  const auto  pos = state.fields(HydroFieldNames::position, Vector::zero());
  const auto  H = state.fields(HydroFieldNames::H, SymTensor::zero());
  const auto  D = state.fields(SolidFieldNames::tensorDamage, SymTensor::zero(), true);
  const auto  mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto  rho = state.fields(HydroFieldNames::massDensity, 0.0);
  auto vol = state.fields(HydroFieldNames::volume, 0.0);
  auto vol3d = state.fields(HydroFieldNames::volume3d, 0.0);
  auto surfacePoint = state.fields(HydroFieldNames::surfacePoint, 0);
  auto cells = state.fields(HydroFieldNames::cells, FacetedVolume());
  auto cellFaceFlags = state.fields(HydroFieldNames::cellFaceFlags, std::vector<CellFaceFlag>());
  auto etaVoidPoints = state.fields(HydroFieldNames::etaVoidPoints, std::vector<Vector>());
  
  // Pre-seed volumes with mass/density, then unscale the geometric factor
  // so the seed is in raw coordinate-plane units.  computeVoronoiVolume1d
  // reads these as fallback extents for boundary nodes.
  vol.assignFields(mass / rho);
  unscaleFromGeometry(pos, vol);

  auto& boundaries = this->boundaryConditions();

  // Always compute full Voronoi tessellation for cell geometry
  computeVoronoiVolume(pos, H, cm, D, mFacetedBoundaries, mFacetedHoles, boundaries, mWeight,
                       surfacePoint, vol, mDeltaCentroid, etaVoidPoints, cells, cellFaceFlags);

  // If the user chose a non-Voronoi volume type, overwrite both volume fields
  // (base class computeVolume handles both volume and volume3d)
  if (this->volumeType() != VolumeType::VoronoiVolume) {
    VolumeUpdate<Dimension>::computeVolume(dataBase, state);
  }
  else {
    // Voronoi volume is coordinate-plane; scale into volume3d
    scaleForGeometry(pos, vol, vol3d);
  }
}

//------------------------------------------------------------------------------
// Size up our FieldLists on problem startup
//------------------------------------------------------------------------------
template<typename Dimension>
void
VoronoiCells<Dimension>::
initializeProblemStartup(DataBase<Dimension>& dataBase) {
  VolumeUpdate<Dimension>::initializeProblemStartup(dataBase);
  // mWeight = dataBase.newFluidFieldList(0.0, "Voronoi weight");
  mSurfacePoint = dataBase.newFluidFieldList(0, HydroFieldNames::surfacePoint);
  mEtaVoidPoints = dataBase.newFluidFieldList(std::vector<Vector>(), HydroFieldNames::etaVoidPoints);
  mCells = dataBase.newFluidFieldList(FacetedVolume(), HydroFieldNames::cells);
  mCellFaceFlags = dataBase.newFluidFieldList(std::vector<CellFaceFlag>(), HydroFieldNames::cellFaceFlags);
  mDeltaCentroid = dataBase.newFluidFieldList(Vector::zero(), "delta centroid");
}

//------------------------------------------------------------------------------
// Compute initial cells
//------------------------------------------------------------------------------
template<typename Dimension>
void
VoronoiCells<Dimension>::
initializeProblemStartupDependencies(DataBase<Dimension>& dataBase,
                                     State<Dimension>& state,
                                     StateDerivatives<Dimension>& derivs) {
  // Resize Voronoi-specific fields
  // dataBase.resizeFluidFieldList(mWeight, 0.0, "Voronoi weight", false);
  dataBase.resizeFluidFieldList(mSurfacePoint, 0, HydroFieldNames::surfacePoint, false);
  dataBase.resizeFluidFieldList(mEtaVoidPoints, vector<Vector>(), HydroFieldNames::etaVoidPoints, false);
  dataBase.resizeFluidFieldList(mCells, FacetedVolume(), HydroFieldNames::cells, false);
  dataBase.resizeFluidFieldList(mCellFaceFlags, vector<CellFaceFlag>(), HydroFieldNames::cellFaceFlags, false);
  dataBase.resizeFluidFieldList(mDeltaCentroid, Vector::zero(), "delta centroid", false);

  // Delegate to base for volume resize + computeVolume + ghost BCs
  VolumeUpdate<Dimension>::initializeProblemStartupDependencies(dataBase, state, derivs);

  // Apply boundaries to newly computed terms
  auto surfacePoint = state.fields(HydroFieldNames::surfacePoint, 0);
  auto etaVoidPoints = state.fields(HydroFieldNames::etaVoidPoints, std::vector<Vector>());
  for (auto boundItr = this->boundaryBegin(); boundItr < this->boundaryEnd(); ++boundItr) {
    (*boundItr)->applyFieldListGhostBoundary(surfacePoint);
    (*boundItr)->applyFieldListGhostBoundary(etaVoidPoints);
  }
  for (auto boundItr = this->boundaryBegin(); boundItr < this->boundaryEnd(); ++boundItr) {
    (*boundItr)->finalizeGhostBoundary();
  }
}

//------------------------------------------------------------------------------
// Register the state — volume from base, plus Voronoi-specific fields
//------------------------------------------------------------------------------
template<typename Dimension>
void
VoronoiCells<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {
  VolumeUpdate<Dimension>::registerState(dataBase, state);
  state.enroll(mSurfacePoint);
  state.enroll(mEtaVoidPoints);
  state.enroll(mCellFaceFlags);
  state.enroll(mCells);
}

//------------------------------------------------------------------------------
// Apply the ghost boundary conditions
//------------------------------------------------------------------------------
template<typename Dimension>
void
VoronoiCells<Dimension>::
applyGhostBoundaries(State<Dimension>& state,
                     StateDerivatives<Dimension>& derivs) {
  VolumeUpdate<Dimension>::applyGhostBoundaries(state, derivs);
  auto cells = state.template fields<FacetedVolume>(HydroFieldNames::cells);
  auto surfacePoint = state.fields(HydroFieldNames::surfacePoint, 0);
  auto etaVoidPoints = state.fields(HydroFieldNames::etaVoidPoints, std::vector<Vector>());
  for (auto* bcPtr: this->boundaryConditions()) {
    bcPtr->applyFieldListGhostBoundary(etaVoidPoints);
    bcPtr->applyFieldListGhostBoundary(cells);
    bcPtr->applyFieldListGhostBoundary(surfacePoint);
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
  VolumeUpdate<Dimension>::enforceBoundaries(state, derivs);
  auto cells = state.template fields<FacetedVolume>(HydroFieldNames::cells);
  auto surfacePoint = state.fields(HydroFieldNames::surfacePoint, 0);
  auto etaVoidPoints = state.fields(HydroFieldNames::etaVoidPoints, std::vector<Vector>());
  for (auto* bcPtr: this->boundaryConditions()) {
    bcPtr->enforceFieldListBoundary(etaVoidPoints);
    bcPtr->enforceFieldListBoundary(cells);
    bcPtr->enforceFieldListBoundary(surfacePoint);
  }
}

//------------------------------------------------------------------------------
// Add a faceted boundary
//------------------------------------------------------------------------------
template<typename Dimension>
void
VoronoiCells<Dimension>::
addFacetedBoundary(const FacetedVolume& bound,
                   const std::vector<FacetedVolume>& holes) {
  if (std::ranges::find(mFacetedBoundaries, bound) == mFacetedBoundaries.end()) {
    mFacetedBoundaries.push_back(bound);
  } else {
    SpheralWarning << "tried to add same faceted boundary twice" << std::endl;
  }
  if (std::ranges::find(mFacetedHoles, holes) == mFacetedHoles.end()) {
    mFacetedHoles.push_back(holes);
  } else {
    SpheralWarning << "tried to add same faceted holes twice" << std::endl;
  }
  ENSURE(mFacetedBoundaries.size() == mFacetedHoles.size());
}

//------------------------------------------------------------------------------
// Dump the current state to the given file
//------------------------------------------------------------------------------
template<typename Dimension>
void
VoronoiCells<Dimension>::
dumpState(FileIO& file, const std::string& pathName) const {
  VolumeUpdate<Dimension>::dumpState(file, pathName);
  file.write(mSurfacePoint, pathName + "/surfacePoint");
}

//------------------------------------------------------------------------------
// Restore the state from the given file
//------------------------------------------------------------------------------
template<typename Dimension>
void
VoronoiCells<Dimension>::
restoreState(const FileIO& file, const std::string& pathName) {
  VolumeUpdate<Dimension>::restoreState(file, pathName);
  file.read(mSurfacePoint, pathName + "/surfacePoint");
}

} // end namespace Spheral
