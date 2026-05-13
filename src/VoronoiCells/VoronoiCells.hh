//---------------------------------Spheral++----------------------------------//
// VoronoiCells
//
// Computes polytopes for each point similar to the Voronoi tessellation.
// Inherits volume ownership from VolumeUpdate, overrides computeVolume
// to perform full Voronoi geometry computation.
//----------------------------------------------------------------------------//
#ifndef __Spheral_VoronoiCells__
#define __Spheral_VoronoiCells__

#include "VoronoiCells/VolumeUpdate.hh"
#include "Geometry/CellFaceFlag.hh"

namespace Spheral {

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class DataBase;
template<typename Dimension> class Boundary;

template<typename Dimension>
class VoronoiCells : public VolumeUpdate<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using FacetedVolume = typename Dimension::FacetedVolume;
  
  using BoundaryIterator = typename std::vector<Boundary<Dimension>*>::iterator;
  using ConstBoundaryIterator = typename std::vector<Boundary<Dimension>*>::const_iterator;
  using TimeStepType = typename std::pair<double, std::string>;

  // Constructor
  VoronoiCells(const VolumeType volumeType,
               const TableKernel<Dimension>& W,
               const std::vector<FacetedVolume>& facetedBoundaries = std::vector<FacetedVolume>(),
               const std::vector<std::vector<FacetedVolume>>& facetedHoles = std::vector<std::vector<FacetedVolume>>(),
               const bool updateInStep = true,
               const bool updateInFinalize = false);

  // Destructor.
  virtual ~VoronoiCells();

  // Override computeVolume for Voronoi geometry
  virtual void computeVolume(const DataBase<Dimension>& dataBase,
                             State<Dimension>& state) override;

  // Size up our FieldLists on problem startup
  virtual void initializeProblemStartup(DataBase<Dimension>& dataBase) override;

  // Compute initial cells
  virtual void initializeProblemStartupDependencies(DataBase<Dimension>& dataBase,
                                                     State<Dimension>& state,
                                                     StateDerivatives<Dimension>& derivs) override;

  // Register additional Voronoi state (surfacePoint, cells, cellFaceFlags)
  virtual void registerState(DataBase<Dimension>& dataBase,
                             State<Dimension>& state) override;

  // Apply boundary conditions to ghost points
  virtual void applyGhostBoundaries(State<Dimension>& state,
                                    StateDerivatives<Dimension>& derivs) override;
  
  // Enforce boundary conditions for internal points
  virtual void enforceBoundaries(State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivs) override;

  // Add a faceted boundary
  virtual void addFacetedBoundary(const FacetedVolume& bound,
                                  const std::vector<FacetedVolume>& holes = std::vector<FacetedVolume>());
  
  // Methods required for restarting.
  virtual std::string label() const override { return "VoronoiCells"; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIO& file, const std::string& pathName);

  // Parameters
  Scalar kernelExtent() const { return mEtaMax; }

  // The Voronoi-specific state field lists
  const FieldList<Dimension, Scalar>&                    weight()            const { return mWeight; }       
  const FieldList<Dimension, int>&                       surfacePoint()      const { return mSurfacePoint; } 
  const FieldList<Dimension, std::vector<Vector>>&       etaVoidPoints()     const { return mEtaVoidPoints; }
  const FieldList<Dimension, FacetedVolume>&             cells()             const { return mCells; }        
  const FieldList<Dimension, std::vector<CellFaceFlag>>& cellFaceFlags()     const { return mCellFaceFlags; }
  const FieldList<Dimension, Vector>&                    deltaCentroid()     const { return mDeltaCentroid; }
  const std::vector<FacetedVolume>&                      facetedBoundaries() const { return mFacetedBoundaries; }
  const std::vector<std::vector<FacetedVolume>>&         facetedHoles()      const { return mFacetedHoles; }

  // No default constructor, copying, or assignment.
  VoronoiCells() = delete;
  VoronoiCells(const VoronoiCells&) = delete;
  VoronoiCells& operator=(const VoronoiCells&) = delete;

private:
  //--------------------------- Private Interface ---------------------------//
  Scalar mEtaMax;
  FieldList<Dimension, Scalar> mWeight;
  FieldList<Dimension, int> mSurfacePoint;
  FieldList<Dimension, std::vector<Vector>> mEtaVoidPoints;
  FieldList<Dimension, FacetedVolume> mCells;
  FieldList<Dimension, std::vector<CellFaceFlag>> mCellFaceFlags;
  FieldList<Dimension, Vector> mDeltaCentroid;
  std::vector<FacetedVolume> mFacetedBoundaries;
  std::vector<std::vector<FacetedVolume>> mFacetedHoles;
};

}

#endif
