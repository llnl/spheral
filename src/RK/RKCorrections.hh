//---------------------------------Spheral++----------------------------------//
// RKCorrections
//
// Computes RK corrections for other physics packages
//----------------------------------------------------------------------------//
#ifndef __Spheral_RKCorrections__
#define __Spheral_RKCorrections__

#include "RK/RKCorrectionParams.hh"
#include "RK/ReproducingKernel.hh"
#include "DataOutput/registerWithRestart.hh"
#include "Field/FieldList.hh"
#include "Geometry/CellFaceFlag.hh"
#include "Physics/Physics.hh"

#include <unordered_map>
#include <set>

namespace Spheral {

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class DataBase;
template<typename Dimension> class Boundary;

template<typename Dimension>
class RKCorrections : public Physics<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::FacetedVolume FacetedVolume;
  
  typedef typename std::vector<Boundary<Dimension>*>::iterator BoundaryIterator;
  typedef typename std::vector<Boundary<Dimension>*>::const_iterator ConstBoundaryIterator;
  typedef typename std::pair<double, std::string> TimeStepType;

  // Constructor
  RKCorrections(const std::set<RKOrder> orders,
                const DataBase<Dimension>& dataBase,
                const TableKernel<Dimension>& W,
                const RKVolumeType volumeType,
                const bool needHessian,
                const bool updateInFinalize);

  // Destructor.
  virtual ~RKCorrections();

  // Evaluate derivatives
  virtual void evaluateDerivatives(const Scalar time,
                                   const Scalar dt,
                                   const DataBase<Dimension>& dataBase,
                                   const State<Dimension>& state,
                                   StateDerivatives<Dimension>& derivatives) const override;
  
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

  // Label for the package
  virtual std::string label() const override { return "RKCorrections"; }

  // Apply boundary conditions to ghost points
  virtual void applyGhostBoundaries(State<Dimension>& state,
                                    StateDerivatives<Dimension>& derivs) override;
  
  // Enforce boundary conditions for internal points
  virtual void enforceBoundaries(State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivs) override;
  
  // Initialize field lists and calculate initial RK corrections
  virtual void initializeProblemStartup(DataBase<Dimension>& dataBase) override;
  
  // Initialize field lists and calculate initial RK corrections
  virtual void initializeProblemStartupDependencies(DataBase<Dimension>& dataBase,
                                                    State<Dimension>& state,
                                                    StateDerivatives<Dimension>& derivs) override;
  
  // Compute the volumes
  virtual void preStepInitialize(const DataBase<Dimension>& dataBase, 
                                 State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivs) override;
  
  // Compute RK corrections
  virtual bool initialize(const Scalar time, 
                          const Scalar dt,
                          const DataBase<Dimension>& dataBase, 
                          State<Dimension>& state,
                          StateDerivatives<Dimension>& derivs) override;

  // Finalize
  virtual void finalize(const Scalar time, 
                        const Scalar dt,
                        DataBase<Dimension>& dataBase, 
                        State<Dimension>& state,
                        StateDerivatives<Dimension>& derivs) override;

  // Add a faceted boundary
  virtual void addFacetedBoundary(const FacetedVolume& cell,
                                  const std::vector<FacetedVolume>& holes);
  
  // We do require the connecitivity
  virtual bool requireConnectivity() const override { return true; }
  
  // Methods required for restarting.
  virtual void dumpState(FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIO& file, const std::string& pathName);

  // Parameters
  std::set<RKOrder> correctionOrders() const { return mOrders; }
  RKVolumeType      volumeType()       const { return mVolumeType; }
  bool              needHessian()      const { return mNeedHessian; }

  // The state field lists we're maintaining.
  const FieldList<Dimension, Scalar>&                    volume()        const { return mVolume; }
  const FieldList<Dimension, Scalar>&                    surfaceArea()   const { return mSurfaceArea; }
  const FieldList<Dimension, Vector>&                    normal()        const { return mNormal; }
  const FieldList<Dimension, int>&                       surfacePoint()  const { return mSurfacePoint; }
  const FieldList<Dimension, std::vector<Vector>>&       etaVoidPoints() const { return mEtaVoidPoints; }
  const FieldList<Dimension, FacetedVolume>&             cells()         const { return mCells; }        
  const FieldList<Dimension, std::vector<CellFaceFlag>>& cellFaceFlags() const { return mCellFaceFlags; }
  const FieldList<Dimension, Vector>&                    deltaCentroid() const { return mDeltaCentroid; }

  // RKOrder dependent state
  const ReproducingKernel<Dimension>&                    WR(const RKOrder order)          const;
  const FieldList<Dimension, RKCoefficients<Dimension>>& corrections(const RKOrder order) const;

private:
  //--------------------------- Private Interface ---------------------------//

  // Data
  std::set<RKOrder> mOrders;
  const DataBase<Dimension>& mDataBase;
  const RKVolumeType mVolumeType;
  const bool mNeedHessian;
  const bool mUpdateInFinalize;
  std::unordered_map<RKOrder, ReproducingKernel<Dimension>> mWR;

  // State
  FieldList<Dimension, Scalar> mVolume;
  
  // Corrections
  FieldList<Dimension, Scalar> mSurfaceArea;
  FieldList<Dimension, Vector> mNormal;
  std::unordered_map<RKOrder, FieldList<Dimension, RKCoefficients<Dimension>>> mCorrections;
  
  // Voronoi stuff
  FieldList<Dimension, int> mSurfacePoint;
  FieldList<Dimension, std::vector<Vector>> mEtaVoidPoints;
  FieldList<Dimension, FacetedVolume> mCells;
  FieldList<Dimension, std::vector<CellFaceFlag>> mCellFaceFlags;
  FieldList<Dimension, Vector> mDeltaCentroid;
  std::vector<FacetedVolume> mFacetedBoundaries;
  std::vector<std::vector<FacetedVolume>> mFacetedHoles;
  
  // The restart registration.
  RestartRegistrationType mRestart;

  // No default constructor, copying, or assignment.
  RKCorrections();
  RKCorrections(const RKCorrections&);
  RKCorrections& operator=(const RKCorrections&);
}; // end RKCorrections

} // end namespace Spheral

#endif
