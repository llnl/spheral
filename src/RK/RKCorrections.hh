//---------------------------------Spheral++----------------------------------//
// RKCorrections
//
// Computes RK corrections for other physics packages.
//----------------------------------------------------------------------------//
#ifndef __Spheral_RKCorrections__
#define __Spheral_RKCorrections__

#include "RK/RKCorrectionParams.hh"
#include "RK/ReproducingKernel.hh"
#include "DataOutput/registerWithRestart.hh"
#include "Field/FieldList.hh"
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
  using VolumeRequirements = typename Physics<Dimension>::VolumeRequirements;

  // Constructor
  RKCorrections(const std::set<RKOrder> orders,
                const DataBase<Dimension>& dataBase,
                const TableKernel<Dimension>& W,
                const bool needHessian,
                const bool updateInStep,
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
  
  // Initialize field lists
  virtual void initializeProblemStartup(DataBase<Dimension>& dataBase) override;
  
  // Compute initial corrections
  virtual void initializeProblemStartupDependencies(DataBase<Dimension>& dataBase,
                                                     State<Dimension>& state,
                                                     StateDerivatives<Dimension>& derivs) override;
  
  // Compute RK corrections
  virtual bool initialize(const Scalar time, 
                          const Scalar dt,
                          const DataBase<Dimension>& dataBase, 
                          State<Dimension>& state,
                          StateDerivatives<Dimension>& derivs) override;

  // Recompute corrections after state update
  virtual bool postStateUpdate(const Scalar time,
                               const Scalar dt,
                               const DataBase<Dimension>& dataBase,
                               State<Dimension>& state,
                               StateDerivatives<Dimension>& derivs) override;

  // Finalize — recompute corrections at end of step
  virtual void finalize(const Scalar time, 
                                const Scalar dt,
                                DataBase<Dimension>& dataBase, 
                                State<Dimension>& state,
                                StateDerivatives<Dimension>& derivs) override;

  // RK needs volumes (not Voronoi) to compute corrections
  virtual VolumeRequirements requireVolumes() const override {
    return {mUpdateInStep, mUpdateInFinalize, false};
  }
  
  // Methods required for restarting.
  virtual void dumpState(FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIO& file, const std::string& pathName);

  // Parameters
  std::set<RKOrder> correctionOrders() const { return mOrders; }
  bool              needHessian()      const { return mNeedHessian; }
  bool              updateInStep()     const { return mUpdateInStep; }
  bool              updateInFinalize() const { return mUpdateInFinalize; }

  // RK-specific state
  const FieldList<Dimension, Scalar>&                    surfaceArea()   const { return mSurfaceArea; }
  const FieldList<Dimension, Vector>&                    normal()        const { return mNormal; }

  // RKOrder dependent state
  const ReproducingKernel<Dimension>&                    WR(const RKOrder order)          const;
  const FieldList<Dimension, RKCoefficients<Dimension>>& corrections(const RKOrder order) const;

private:
  //--------------------------- Private Interface ---------------------------//

  // Data
  std::set<RKOrder> mOrders;
  const DataBase<Dimension>& mDataBase;
  const bool mNeedHessian;
  const bool mUpdateInStep;
  const bool mUpdateInFinalize;
  std::unordered_map<RKOrder, ReproducingKernel<Dimension>> mWR;

  // Corrections
  FieldList<Dimension, Scalar> mSurfaceArea;
  FieldList<Dimension, Vector> mNormal;
  std::unordered_map<RKOrder, FieldList<Dimension, RKCoefficients<Dimension>>> mCorrections;
  
  // Private helper for correction computation
  void updateCorrections(const DataBase<Dimension>& dataBase, State<Dimension>& state);

  // The restart registration.
  RestartRegistrationType mRestart;

  // No default constructor, copying, or assignment.
  RKCorrections();
  RKCorrections(const RKCorrections&);
  RKCorrections& operator=(const RKCorrections&);
}; // end RKCorrections

} // end namespace Spheral

#endif
