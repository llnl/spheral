//---------------------------------Spheral++----------------------------------//
// VoronoiCells
//
// Computes polytopes for each point similar to the Voronoi tessellation
//----------------------------------------------------------------------------//
#ifndef __Spheral_VoronoiCells__
#define __Spheral_VoronoiCells__

#include "DataOutput/registerWithRestart.hh"
#include "Field/FieldList.hh"
#include "Geometry/CellFaceFlag.hh"
#include "Physics/Physics.hh"

namespace Spheral {

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class DataBase;
template<typename Dimension> class Boundary;

template<typename Dimension>
class VoronoiCells : public Physics<Dimension> {
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
  VoronoiCells(const Scalar kernelExtent,
               const std::vector<FacetedVolume>& facetedBoundaries = std::vector<FacetedVolume>(),
               const std::vector<std::vector<FacetedVolume>>& facetedHoles = std::vector<std::vector<FacetedVolume>>());

  // Destructor.
  virtual ~VoronoiCells();

  //******************************************************************************//
  // An optional hook to initialize once when the problem is starting up.
  // This is called after the materials and NodeLists are created. This method
  // should set the sizes of all arrays owned by the physics package and initialize
  // independent variables.
  // It is assumed after this method has been called it is safe to call
  // Physics::registerState to create full populated State objects.
  virtual void initializeProblemStartup(DataBase<Dimension>& dataBase) override;

  // A second optional method to be called on startup, after Physics::initializeProblemStartup
  // has been called.
  // This method is called after independent variables have been initialized and put into
  // the state and derivatives. During this method, the dependent state, such as
  // temperature and pressure, is initialized so that all the fields in the initial
  // state and derivatives objects are valid.
  virtual void initializeProblemStartupDependencies(DataBase<Dimension>& dataBase,
                                                    State<Dimension>& state,
                                                    StateDerivatives<Dimension>& derivs) override;

  // Evaluate derivatives
  virtual void evaluateDerivatives(const Scalar time,
                                   const Scalar dt,
                                   const DataBase<Dimension>& dataBase,
                                   const State<Dimension>& state,
                                   StateDerivatives<Dimension>& derivatives) const override;
  
  // Optional hook to be called at the beginning of a time step.
  virtual void preStepInitialize(const DataBase<Dimension>& dataBase, 
                                 State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivs) override;

  // Vote on a time step.
  virtual TimeStepType dt(const DataBase<Dimension>& dataBase, 
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs,
                          const Scalar currentTime) const override;

  // Register the state
  virtual void registerState(DataBase<Dimension>& dataBase,
                             State<Dimension>& state) override;

  // Register the state derivatives
  virtual void registerDerivatives(DataBase<Dimension>& dataBase,
                                   StateDerivatives<Dimension>& derivs) override;

  // Apply boundary conditions to ghost points
  virtual void applyGhostBoundaries(State<Dimension>& state,
                                    StateDerivatives<Dimension>& derivs) override;
  
  // Enforce boundary conditions for internal points
  virtual void enforceBoundaries(State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivs) override;
  
  // Provide a hook to be called after the state has been updated and 
  // boundary conditions have been enforced.
  virtual bool postStateUpdate(const Scalar time, 
                               const Scalar dt,
                               const DataBase<Dimension>& dataBase, 
                               State<Dimension>& state,
                               StateDerivatives<Dimension>& derivatives) override;

 // Add a faceted boundary
  virtual void addFacetedBoundary(const FacetedVolume& bound,
                                  const std::vector<FacetedVolume>& holes);
  
  // We do require the connecitivity
  virtual bool requireConnectivity() const override { return true; }
  
  // Methods required for restarting.
  virtual std::string label() const override { return "VoronoiCells"; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIO& file, const std::string& pathName);

  // Parameters
  Scalar kernelExtent() const { return mEtaMax; }

  // The state field lists we're maintaining.
  const FieldList<Dimension, Scalar>&                    volume()            const { return mVolume; }       
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
  FieldList<Dimension, Scalar> mVolume, mWeight;
  FieldList<Dimension, int> mSurfacePoint;
  FieldList<Dimension, std::vector<Vector>> mEtaVoidPoints;
  FieldList<Dimension, FacetedVolume> mCells;
  FieldList<Dimension, std::vector<CellFaceFlag>> mCellFaceFlags;
  FieldList<Dimension, Vector> mDeltaCentroid;
  std::vector<FacetedVolume> mFacetedBoundaries;
  std::vector<std::vector<FacetedVolume>> mFacetedHoles;
  
  // The restart registration.
  RestartRegistrationType mRestart;
};

}

#endif
