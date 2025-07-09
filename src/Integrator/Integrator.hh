//---------------------------------Spheral++----------------------------------//
// Integrator -- The topmost abstract base class for all integrator classes
// Spheral++.  Integrator classes take a list of Physics packages and advance
// them in time.
//
// Created by JMO, Wed May 31 21:58:08 PDT 2000
//----------------------------------------------------------------------------//
#ifndef Integrator_HH
#define Integrator_HH

#include "DataOutput/registerWithRestart.hh"

#ifdef USE_MPI
#include <mpi.h>
#endif

#include <string>
#include <vector>

namespace Spheral {

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class DataBase;
template<typename Dimension> class Physics;
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

  // Constructors.
  Integrator(DataBase<Dimension>& dataBase);
  Integrator(DataBase<Dimension>& dataBase,
             const std::vector<Physics<Dimension>*>& physicsPackages);

  // Destructor.
  virtual ~Integrator() = default;

  // Assignment.
  Integrator& operator=(const Integrator& rhs) = default;

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
  Scalar currentTime() const;
  void currentTime(Scalar time);

  // Access the current cycle.
  int currentCycle() const;
  void currentCycle(int cycle);

  // Access the minimum allowed timestep.
  Scalar dtMin() const;
  void dtMin(Scalar dt);

  // Access the maximum allowed timestep.
  Scalar dtMax() const;
  void dtMax(Scalar dt);

  // Access the last timestep.
  Scalar lastDt() const;
  void lastDt(Scalar dt);

  // Access the timestep growth factor.
  Scalar dtGrowth() const;
  void dtGrowth(Scalar fraction);

  // The fraction of the timestep we consider when checking for stable behavior.
  Scalar dtCheckFrac() const;
  void dtCheckFrac(Scalar fraction);

  // Public const access to the DataBase.
  const DataBase<Dimension>& dataBase() const;

  // Access the list of physics packages.
  const std::vector<Physics<Dimension>*>& physicsPackages() const;

  // Provide standard iterator methods over the physics package list.
  PackageIterator physicsPackagesBegin();
  PackageIterator physicsPackagesEnd();

  ConstPackageIterator physicsPackagesBegin() const;
  ConstPackageIterator physicsPackagesEnd() const;

  // Flag to determine whether or not to be rigorous about about boundaries.
  bool rigorousBoundaries() const;
  void rigorousBoundaries(bool value);

  // If we're not being rigorous about boundary conditions, how frequently
  // do we update them?
  int updateBoundaryFrequency() const;
  void updateBoundaryFrequency(int value);

  // Select whether the integrator is verbose or not during a cycle.
  bool verbose() const;
  void verbose(bool x);

  // Should the integrator check interim timestep votes and abort steps?
  bool allowDtCheck() const;
  void allowDtCheck(bool x);

  // Select whether we should run in a mode the ensures domain decomposition independence.
  // Possibly some performance impact.
  bool domainDecompositionIndependent() const;
  void domainDecompositionIndependent(bool x);

  // Select whether we're going to enforce culling of ghost nodes or not.
  bool cullGhostNodes() const;
  void cullGhostNodes(bool x);

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const { return "Integrator"; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIO& file, const std::string& pathName);
  //****************************************************************************

  // Forbiddent methods
  Integrator() = delete;

protected:
  //-------------------------- Protected Interface --------------------------//
  // Allow write access to the DataBase for descendent classes.
  DataBase<Dimension>& accessDataBase() const;

private:
  //--------------------------- Private Interface ---------------------------//
  Scalar mDtMin, mDtMax, mDtGrowth, mLastDt, mDtMultiplier, mDtCheckFrac, mCurrentTime;
  int mCurrentCycle, mUpdateBoundaryFrequency;
  bool mVerbose, mAllowDtCheck, mRequireConnectivity, mRequireGhostConnectivity, mRequireOverlapConnectivity, mRequireIntersectionConnectivity;
  std::reference_wrapper<DataBase<Dimension>> mDataBase;
  std::vector<Physics<Dimension>*> mPhysicsPackages;
  bool mCullGhostNodes;

  // The restart registration.
  RestartRegistrationType mRestart;
};

}

#include "IntegratorInline.hh"

#endif
