//---------------------------------Spheral++----------------------------------//
// Integrator -- The topmost abstract base class for all integrator classes
// Spheral++.  Integrator classes take a list of Physics packages and advance
// them in time.
//
// Created by JMO, Wed May 31 21:58:08 PDT 2000
//----------------------------------------------------------------------------//
#ifndef __Spheral_Integrator__
#define __Spheral_Integrator__

#include "DataOutput/registerWithRestart.hh"
#include "NodeList/NodeListRegistrar.hh"
#include "Physics/Physics.hh"
#include "Utilities/SpheralMessage.hh"
#include "Utilities/DBC.hh"

#include <string>
#include <vector>

namespace Spheral {

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class DataBase;
template<typename Dimension, typename DataType> class FieldList;
template<typename Dimension> class Boundary;
class FileIO;

template<typename Dimension>
class Integrator {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  using PackageIterator = typename std::vector<Physics<Dimension>*>::iterator;
  using ConstPackageIterator = typename std::vector<Physics<Dimension>*>::const_iterator;

  using BoundaryIterator = typename std::vector<Boundary<Dimension>*>::iterator;
  using ConstBoundaryIterator = typename std::vector<Boundary<Dimension>*>::const_iterator;

  using TimeStepType = typename Physics<Dimension>::TimeStepType;

  // Constructors.
  Integrator(DataBase<Dimension>& dataBase,
             const std::vector<Physics<Dimension>*>& physicsPackages = std::vector<Physics<Dimension>*>());
  Integrator& operator=(const Integrator& rhs) = default;
  virtual ~Integrator() = default;

  // All Integrator classes must define the dt and step methods.
  virtual bool step(Scalar maxTime,
                    State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) = 0;
  virtual bool step(Scalar maxTime);

  // Provide a method of looping over the physics packages and picking a
  // time step.
  virtual Scalar selectDt(const Scalar dtMin, 
                          const Scalar dtMax,
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs) const;

  // Perform generic initializations at the beginning of a timestep.
  // To be called once per advance cycle.
  virtual void preStepInitialize(State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivs) const;

  // Prepare all physics packages for calls to evaluateDerivatives.
  // To be called before any call to Physics::evaluateDerivatives, therefore potentially
  // several times during a time step.
  virtual void initializeDerivatives(const double t,
                                     const double dt,
                                     State<Dimension>& state,
                                     StateDerivatives<Dimension>& derivs) const;

  // Iterate over all physics packages and call evaluateDerivatives.
  virtual void evaluateDerivatives(const Scalar t,
                                   const Scalar dt,
                                   const DataBase<Dimension>& dataBase,
                                   const State<Dimension>& state,
                                   StateDerivatives<Dimension>& derivs) const;

  // Iterate over all physics packages and call finalizeDerivatives.
  virtual void finalizeDerivatives(const Scalar t,
                                   const Scalar dt,
                                   const DataBase<Dimension>& dataBase,
                                   const State<Dimension>& state,
                                   StateDerivatives<Dimension>& derivs) const;

  // Iterate over all physics packages and call postStateUpdate
  virtual void postStateUpdate(const Scalar t,
                               const Scalar dt,
                               const DataBase<Dimension>& dataBase,
                               State<Dimension>& state,
                               StateDerivatives<Dimension>& derivs) const;

  // Finalize at the end a timestep, therefore called once at the end of a timestep.
  virtual void postStepFinalize(const double t,
                                const double dt,
                                State<Dimension>& state,
                                StateDerivatives<Dimension>& derivs) const;

  // Add a Physics package.
  void appendPhysicsPackage(Physics<Dimension>& package);

  // Reset the Physics packages to a new set
  void resetPhysicsPackages(std::vector<Physics<Dimension>*>& packages);

  // Test if the given Physics package is listed in the integrator.
  bool havePhysicsPackage(const Physics<Dimension>& package) const;

  // Get the unique set of boundary conditions across all physics packages.
  std::vector<Boundary<Dimension>*> uniqueBoundaryConditions() const;

  // Set the ghost nodes for all node lists according to the boundary 
  // conditions.
  void setGhostNodes() const;

  // Set the ghost node values on the Fields of the nodes lists in the
  // data base.
  void applyGhostBoundaries(State<Dimension>& state,
                            StateDerivatives<Dimension>& derivs) const;

  // Finalize the ghost node boundary conditions.
  void finalizeGhostBoundaries() const;

  // Find the nodes in violation of the boundary conditions.
  void setViolationNodes() const;

  // Reset any internal nodes in violation of boundary conditions to be brought 
  // into compliance.
  void enforceBoundaries(State<Dimension>& state,
                         StateDerivatives<Dimension>& derivs) const;

  // Copy the ghost positions and H's from one state to another.
  void copyGhostState(const State<Dimension>& state0,
                      State<Dimension>& state1) const;

  // Access the current time.
  Scalar currentTime() const                                                        { return mCurrentTime; }
  void currentTime(const Scalar x)                                                  { mCurrentTime = x; }

  // Access the current cycle.
  int currentCycle() const                                                          { return mCurrentCycle; }
  void currentCycle(const int x)                                                    { mCurrentCycle = x; }

  // Access the minimum allowed timestep.
  Scalar dtMin() const                                                              { return mDtMin; }
  void dtMin(const Scalar x)                                                        { mDtMin = x; }

  // Access the maximum allowed timestep.
  Scalar dtMax() const                                                              { return mDtMax; }
  void dtMax(const Scalar x)                                                        { mDtMax = x; }

  // Access the last timestep.
  Scalar lastDt() const                                                             { return mLastDt; }
  void lastDt(const Scalar x)                                                       { mLastDt = x; }

  // Access the timestep growth factor.
  Scalar dtGrowth() const                                                           { return mDtGrowth; }
  void dtGrowth(const Scalar x)                                                     { mDtGrowth = x; }

  // The fraction of the timestep we consider when checking for stable behavior.
  Scalar dtCheckFrac() const                                                        { return mDtCheckFrac; }
  void dtCheckFrac(const Scalar x)                                                  { mDtCheckFrac = x; }

  // Public const access to the DataBase.
  const DataBase<Dimension>& dataBase() const                                       { return mDataBase.get(); }

  // Access the list of physics packages.
  const std::vector<Physics<Dimension>*>& physicsPackages() const                   { return mPhysicsPackages; }

  // Provide standard iterator methods over the physics package list.
  PackageIterator physicsPackagesBegin()                                            { return mPhysicsPackages.begin(); }
  PackageIterator physicsPackagesEnd()                                              { return mPhysicsPackages.end(); }

  ConstPackageIterator physicsPackagesBegin() const                                 { return mPhysicsPackages.begin(); }
  ConstPackageIterator physicsPackagesEnd() const                                   { return mPhysicsPackages.end(); }

  // Flag to determine whether or not to be rigorous about about boundaries.
  bool rigorousBoundaries() const                                                   { DeprecationWarning("Integrator::rigorousBoundaries"); return false; }
  void rigorousBoundaries(const bool x)                                             { DeprecationWarning("Integrator::rigorousBoundaries"); }

  // If we're not being rigorous about boundary conditions, how frequently
  // do we update them?
  int updateBoundaryFrequency() const                                               { return mUpdateBoundaryFrequency; }
  void updateBoundaryFrequency(const int x)                                         { mUpdateBoundaryFrequency = x; }

  // Select whether the integrator is verbose or not during a cycle.
  bool verbose() const                                                              { return mVerbose; }
  void verbose(const bool x)                                                        { mVerbose = x; }

  // Should the integrator check interim timestep votes and abort steps?
  bool allowDtCheck() const                                                         { return mAllowDtCheck; }
  void allowDtCheck(const bool x)                                                   { mAllowDtCheck = x; }

  // Select whether we should run in a mode the ensures domain decomposition independence.
  // Possibly some performance impact.
  bool domainDecompositionIndependent() const                                       { return NodeListRegistrar<Dimension>::instance().domainDecompositionIndependent(); }
  void domainDecompositionIndependent(const bool x)                                 { NodeListRegistrar<Dimension>::instance().domainDecompositionIndependent(x); }

  // Select whether we're going to enforce culling of ghost nodes or not.
  bool cullGhostNodes() const                                                       { return mCullGhostNodes; }
  void cullGhostNodes(const bool x)                                                 { mCullGhostNodes = x; }

  // The timestep multiplier
  Scalar dtMultiplier()                                                       const { return mDtMultiplier; }

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const { return "Integrator"; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIO& file, const std::string& pathName);
  //****************************************************************************

  // Forbidden methods
  Integrator() = delete;

protected:
  //-------------------------- Protected Interface --------------------------//
  // Allow write access to the DataBase for descendent classes.
  DataBase<Dimension>& accessDataBase()                                     const  { return mDataBase.get(); }

  // How should we query a physics package for the time step?
  virtual TimeStepType dt(const Physics<Dimension>* pkg,
                          const DataBase<Dimension>& dataBase,
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs,
                          const Scalar currentTime) const;

  mutable Scalar mDtMultiplier;

private:
  //--------------------------- Private Interface ---------------------------//
  Scalar mDtMin, mDtMax, mDtGrowth, mLastDt, mDtCheckFrac, mCurrentTime;
  int mCurrentCycle, mUpdateBoundaryFrequency;
  bool mVerbose, mAllowDtCheck, mRequireConnectivity, mRequireGhostConnectivity, mRequireOverlapConnectivity, mRequireIntersectionConnectivity;
  std::reference_wrapper<DataBase<Dimension>> mDataBase;
  std::vector<Physics<Dimension>*> mPhysicsPackages;
  bool mCullGhostNodes;

  // The restart registration.
  RestartRegistrationType mRestart;
};

}

#endif
