//---------------------------------Spheral++----------------------------------//
// Integrator -- The topmost abstract base class for all integrator classes
// Spheral++.  Integrator classes take a list of Physics packages and advance
// them in time.
//
// Created by JMO, Wed May 31 21:58:08 PDT 2000
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
#include "DataOutput/Restart.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/FieldList.hh"
#include "Physics/Physics.hh"
#include "Boundary/Boundary.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Utilities/range.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Distributed/allReduce.hh"
#include "Distributed/Communicator.hh"
#include "Utilities/DBC.hh"
#include "Integrator.hh"

#include <limits.h>
#include <float.h>
#include <algorithm>

using std::vector;
using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {

//------------------------------------------------------------------------------
// Construct with the given DataBase and Physics packages.
//------------------------------------------------------------------------------
template<typename Dimension>
Integrator<Dimension>::
Integrator(DataBase<Dimension>& dataBase,
           const vector<Physics<Dimension>*>& physicsPackages):
  mDtMultiplier(1.0),
  mDtMin(0.0),
  mDtMax(std::numeric_limits<Scalar>::max()),
  mDtGrowth(2.0),
  mLastDt(1e-5),
  mDtCheckFrac(0.5),
  mCurrentTime(0.0),
  mCurrentCycle(0),
  mUpdateBoundaryFrequency(1),
  mVerboseStep(1),
  mVerbose(false),
  mAllowDtCheck(false),
  mRequireConnectivity(true),
  mRequireGhostConnectivity(false),
  mRequireOverlapConnectivity(false),
  mRequireIntersectionConnectivity(false),
  mDataBase(dataBase),
  mPhysicsPackages(),
  mCullGhostNodes(true),
  mRestart(registerWithRestart(*this)) {
  for (auto* pkg: physicsPackages) this->appendPhysicsPackage(*pkg);
}

//------------------------------------------------------------------------------
// step
//------------------------------------------------------------------------------
template<typename Dimension>
bool
Integrator<Dimension>::
step(const typename Dimension::Scalar maxTime) {
  DataBase<Dimension>& db = this->accessDataBase();

  // Check if we need to construct connectivity.
  mRequireConnectivity = false;
  mRequireGhostConnectivity = false;
  mRequireOverlapConnectivity = false;
  mRequireIntersectionConnectivity = false;
  for (auto* physicsPtr: mPhysicsPackages) {
    const auto [conn, ghost, overlap, intersect] = physicsPtr->requireConnectivity();
    mRequireConnectivity = (mRequireConnectivity or conn);
    mRequireGhostConnectivity = (mRequireGhostConnectivity or ghost);
    mRequireOverlapConnectivity = (mRequireOverlapConnectivity or overlap);
    mRequireIntersectionConnectivity = (mRequireIntersectionConnectivity or intersect);
  }

  // Set the ghost nodes (this updates the ConnectivityMap as well in the DataBase)
  if (mCurrentCycle % mUpdateBoundaryFrequency == 0) setGhostNodes();

  // Build the state and derivatives
  State<Dimension> state(db, this->physicsPackagesBegin(), this->physicsPackagesEnd());
  StateDerivatives<Dimension> derivs(db, this->physicsPackagesBegin(), this->physicsPackagesEnd());

  // Set boundary properties
  applyGhostBoundaries(state, derivs);

  // Try to advance using the derived class step method
  auto success = false;
  auto count = 0u;
  auto maxIterations = 10u;
  while (not success and count < maxIterations) {
    ++count;
    if (count == maxIterations) mAllowDtCheck = false;
    success = this->step(maxTime, state, derivs);
    if (count == maxIterations) mAllowDtCheck = true;
    if (not success) {
      mDtMultiplier *= 0.5;
      if (Process::getRank() == 0) cerr << "Integrator::step reported unstable timestep -- cutting dt and trying again: " << count << "/10" << endl;
    }
  }
  mDtMultiplier = 1.0;
  return success;
}

//------------------------------------------------------------------------------
// Loop over the stored physics packages and pick the minimum timestep.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
Integrator<Dimension>::
selectDt(const typename Dimension::Scalar dtMin,
         const typename Dimension::Scalar dtMax,
         const State<Dimension>& state,
         const StateDerivatives<Dimension>& derivs) const {

  REQUIRE(dtMin >= 0 and dtMax > 0);
  REQUIRE(dtMin <= dtMax);

  // Get the current time and data base.
  auto        t = currentTime();
  const auto& db = dataBase();

  // Loop over each package, and pick their timesteps.
  TimeStepType dt(dtMax, "");
  for (auto* pkg: mPhysicsPackages) {
    auto dtVote = this->dt(pkg, db, state, derivs, t);
    if (dtVote.first > 0.0 and dtVote.first < dt.first) dt = dtVote;
  }

  // Apply any dt scaling due to iteration
  dt.first *= mDtMultiplier;

  // We also require that the timestep is not allowed to grow faster than a
  // prescribed fraction.
  dt.first = std::min(dt.first, dtGrowth()*lastDt());

  // Enforce the timestep boundaries.
  dt.first = std::min(dtMax, std::max(dtMin, dt.first));

  CHECK(dt.first >= 0.0 and
        dt.first >= dtMin and dt.first <= dtMax);

  // In the parallel case we need to find the minimum timestep across all processors.
#ifdef SPHERAL_ENABLE_GLOBALDT_REDUCTION
  const auto [globalDt, rank] = allReduceLoc(dt.first, SPHERAL_OP_MINLOC);
#else
  const auto globalDt = dt.first;
  const auto rank = Process::getRank();
#endif

  // Are we verbose?
  if (rank == Process::getRank() and (verbose() and currentCycle() % verboseStep() == 0)) {
    std::cout << "Selected timestep of "
              << dt.first << " on rank " << rank << std::endl
              << dt.second << std::endl;
  }
  std::cout.flush();
  dt.first = globalDt;

  return dt.first;
}

//------------------------------------------------------------------------------
// Initializations for integrators.  To be called once at the beginning of an
// integrator step.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Integrator<Dimension>::preStepInitialize(State<Dimension>& state,
                                         StateDerivatives<Dimension>& derivs) const {
  auto& db = mDataBase.get();
  for (auto* physicsPtr: mPhysicsPackages) {
    physicsPtr->preStepInitialize(db, state, derivs);
  }
}

//------------------------------------------------------------------------------
// Initialize all physics packages before evaluating derivatives.
// Called before each call to Physics::evaluateDerivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Integrator<Dimension>::initializeDerivatives(const double t,
                                             const double dt,
                                             State<Dimension>& state,
                                             StateDerivatives<Dimension>& derivs) const {

  // Initialize the work fields.
  auto& db = mDataBase.get();
  auto work = db.globalWork();
  work = 0.0;

  // Loop over the physics packages and perform any necessary initializations.
  auto updateBoundaries = false;
  for (auto* physicsPtr: mPhysicsPackages) {
    updateBoundaries |= physicsPtr->initialize(t, dt, db, state, derivs);
  }

  // Apply boundaries if requested.
  if (updateBoundaries) {
    this->applyGhostBoundaries(state, derivs);
    this->finalizeGhostBoundaries();
  }
}

//------------------------------------------------------------------------------
// Iterate over all physics packages and call evaluateDerivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Integrator<Dimension>::evaluateDerivatives(const Scalar t,
                                           const Scalar dt,
                                           const DataBase<Dimension>& dataBase,
                                           const State<Dimension>& state,
                                           StateDerivatives<Dimension>& derivs) const {

  // Loop over the physics packages and have them evaluate their derivatives.
  for (auto* physicsPtr: mPhysicsPackages) {
    physicsPtr->evaluateDerivatives(t, dt, dataBase, state, derivs);
  }
}

//------------------------------------------------------------------------------
// Iterate over all physics packages and call finalizeDerivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Integrator<Dimension>::finalizeDerivatives(const Scalar t,
                                           const Scalar dt,
                                           const DataBase<Dimension>& dataBase,
                                           const State<Dimension>& state,
                                           StateDerivatives<Dimension>& derivs) const {

  // Loop over the physics packages and have them finalize their derivatives.
  for (auto* physicsPtr: mPhysicsPackages) {
    physicsPtr->finalizeDerivatives(t, dt, dataBase, state, derivs);
  }
}

//------------------------------------------------------------------------------
// Iterate over all physics packages and let them do any post state update 
// stuff.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Integrator<Dimension>::postStateUpdate(const Scalar t,
                                       const Scalar dt,
                                       const DataBase<Dimension>& dataBase, 
                                       State<Dimension>& state,
                                       StateDerivatives<Dimension>& derivs) const {

  // Loop over the physics packages.
  auto updateBoundaries = false;
  for (auto* physicsPtr: mPhysicsPackages) {
    updateBoundaries |= physicsPtr->postStateUpdate(t, dt, dataBase, state, derivs);
  }

  // Apply boundaries if requested.
  if (updateBoundaries) {
    this->applyGhostBoundaries(state, derivs);
    this->finalizeGhostBoundaries();
  }
}

//------------------------------------------------------------------------------
// Finalize at the completion of a time step.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Integrator<Dimension>::postStepFinalize(const double t,
                                        const double dt,
                                        State<Dimension>& state,
                                        StateDerivatives<Dimension>& derivs) const {

  // Loop over the physics packages and perform any necessary finalizations.
  auto& db = mDataBase.get();
  for (auto* physicsPtr: mPhysicsPackages) {
    physicsPtr->finalize(t, dt, db, state, derivs);
  }
}

//------------------------------------------------------------------------------
// Add a physics package.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Integrator<Dimension>::
appendPhysicsPackage(Physics<Dimension>& package) {
  if (!havePhysicsPackage(package)) {
    for (auto* packagePtr: package.preSubPackages()) this->appendPhysicsPackage(*packagePtr);
    mPhysicsPackages.push_back(&package);
    for (auto* packagePtr: package.postSubPackages()) this->appendPhysicsPackage(*packagePtr);
  } else {
    cerr << "Warning: attempt to append Physics package " << &package
         << "to Integrator " << this << " which already has it." << endl;
  }
}

//------------------------------------------------------------------------------
// Reset the Physics packages to a new set
//------------------------------------------------------------------------------
template<typename Dimension>
void
Integrator<Dimension>::
resetPhysicsPackages(std::vector<Physics<Dimension>*>& packages) {
  mPhysicsPackages = packages;
}

//------------------------------------------------------------------------------
// Test if the given physics package is listed in the integrator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
Integrator<Dimension>::
havePhysicsPackage(const Physics<Dimension>& package) const {
  return count(mPhysicsPackages.begin(), mPhysicsPackages.end(), &package) > 0;
}

//------------------------------------------------------------------------------
// Get the unique set of boundary conditions across all physics packages.
//------------------------------------------------------------------------------
template<typename Dimension>
vector<Boundary<Dimension>*>
Integrator<Dimension>::
uniqueBoundaryConditions() const {

  // The result we're going to build up.
  vector<Boundary<Dimension>*> result;

  // Iterate over each physics package.
  for (auto* physicsPtr: mPhysicsPackages) {

    // Iterate over the boundary conditions associated with this package.
    for (auto* boundaryPtr: physicsPtr->boundaryConditions()) {
      if (find(result.begin(), result.end(), boundaryPtr) == result.end())
        result.push_back(boundaryPtr);
    }

  }

  // Sort by priority (stable to preserve insertion order among equal priorities)
  std::stable_sort(result.begin(), result.end(),
                   [](const Boundary<Dimension>* a, const Boundary<Dimension>* b) {
                     return a->priority() < b->priority(); });

  BEGIN_CONTRACT_SCOPE
  // Ensure that all boundary conditions are included in the result
  for (auto* physicsPtr: mPhysicsPackages) {
    for (auto* boundaryPtr: physicsPtr->boundaryConditions()) {
      CONTRACT_VAR(boundaryPtr);
      ENSURE(find(result.begin(), result.end(), boundaryPtr) != result.end());
    }
  }

  // Also ensure that there are no duplicates.
  for (auto* boundaryPtr: result) {
    CONTRACT_VAR(boundaryPtr);
    ENSURE(count(result.begin(), result.end(), boundaryPtr) == 1);
  }

  // Ensure boundaries are sorted by priority.
  for (auto i = 1u; i < result.size(); ++i) {
    ENSURE2(result[i-1]->priority() <= result[i]->priority(),
            "Boundary priority ordering violated: " << result[i-1]->priority() << " > " << result[i]->priority());
  }
  END_CONTRACT_SCOPE

  // That's it.
  return result;
}

//------------------------------------------------------------------------------
// Set up the ghost nodes for boundary conditions.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Integrator<Dimension>::setGhostNodes() const {

  // Get that DataBase.
  auto& db = mDataBase.get();

  // Get the complete set of unique boundary conditions.
  const auto boundaries = uniqueBoundaryConditions();

  // Remove any old ghost node information from the NodeLists.
  for (auto* nodeListPtr: db.fluidNodeListPtrs()) nodeListPtr->numGhostNodes(0);
  for (auto* nodeListPtr: db.DEMNodeListPtrs()) nodeListPtr->numGhostNodes(0);

  // If we're need overlap connectivity, we need to double the kernel extent before setting ghost nodes.
  if (mRequireOverlapConnectivity) {
    auto maxeta = 2.0*db.kernelExtent();
    db.kernelExtent(maxeta); // Implies reinitializeNeighbors
  } else {
    db.reinitializeNeighbors();
  }

  // Iterate over the boundaries and set their ghost node info.
  for (auto* boundaryPtr: boundaries) {
    boundaryPtr->setAllGhostNodes(db);
    boundaryPtr->finalizeGhostBoundary();
    db.reinitializeNeighbors();
  }

  // If we doubled the kernel extents for overlap connectivity, put 'em back.
  if (mRequireOverlapConnectivity) {
    auto maxeta = 0.5*db.kernelExtent();
    db.kernelExtent(maxeta);  // Implies reinitializeNeighbors
  }

  if (mRequireConnectivity) {

    // Update the connectivity.
    db.updateConnectivityMap(mRequireGhostConnectivity, mRequireOverlapConnectivity, mRequireIntersectionConnectivity);

    // If we're culling ghost nodes, do it now.
    if (mCullGhostNodes and 
        (not this->domainDecompositionIndependent()) and
        (not mRequireGhostConnectivity) and
        (not mRequireOverlapConnectivity)) {
      const auto numNodeLists = db.numNodeLists();
      const auto& cm = db.connectivityMap();
      const auto& pairs = cm.nodePairList();

      // First build the set of flags indicating which nodes are used.
      FieldList<Dimension, size_t> flags = db.newGlobalFieldList(size_t(0u), "active nodes");
      for (auto p: pairs) {
        flags(p.i_list, p.i_node) = 1;
        flags(p.j_list, p.j_node) = 1;
      }
      
      // Look for points needed for boundary conditions
      for (auto [nodeListi, nodeListPtr]: enumerate(db.nodeListPtrs())) {
        const auto& nodeList = *nodeListPtr;

        // Ghost nodes that are control nodes for other ghost nodes we're keeping must
        // be kept as well.
        const auto firstGhostNode = nodeList.firstGhostNode();
        for (auto* boundaryPtr: boundaries) {
          const auto& controlNodes = boundaryPtr->controlNodes(nodeList);
          const auto& ghostNodes = boundaryPtr->ghostNodes(nodeList);
          // CHECK(controlNodes.size() == ghostNodes.size());  // Not true if this is a DistributedBoundary!
          for (auto i: controlNodes) {
            if (i >= firstGhostNode) flags(nodeListi, i) = 1;
          }

          // Boundary conditions are allowed to opt out of culling entirely.
          if (not boundaryPtr->allowGhostCulling()) {
            for (auto i: ghostNodes) flags(nodeListi, i) = 1;
          }
        }
      }

      // Create the index mapping from old to new node orderings.
      FieldList<Dimension, size_t> old2newIndexMap = db.newGlobalFieldList(size_t(0u), "index map");
      for (auto [nodeListi, nodeListPtr]: enumerate(db.nodeListPtrs())) {
        const auto numNodes = nodeListPtr->numNodes();
        for (auto i = 0u; i < numNodes; ++i) old2newIndexMap(nodeListi, i) = i;
      }

      // Now use these flags to cull the boundary conditions.
      vector<size_t> numNodesRemoved(numNodeLists, 0);
      for (auto* boundaryPtr: boundaries) {
        boundaryPtr->cullGhostNodes(flags, old2newIndexMap, numNodesRemoved);
      }

      // Patch up the connectivity map.
      db.patchConnectivityMap(flags, old2newIndexMap);

      // Now the boundary conditions have been updated, so we can go ahead and remove
      // the ghost nodes themselves from the NodeLists.
      for (auto [nodeListi, nodeListPtr]: enumerate(db.nodeListPtrs())) {
        vector<size_t> nodesToRemove;
        for (auto i = nodeListPtr->firstGhostNode(); i < nodeListPtr->numNodes(); ++i) {
          if (flags(nodeListi, i) == 0) nodesToRemove.push_back(i);
        }
        nodeListPtr->deleteNodes(nodesToRemove);
      }
      db.updateNeighbors();

      // All nodes should now be labeled as keepers.
      BEGIN_CONTRACT_SCOPE
      {
        for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
          ENSURE(flags[nodeListi]->numElements() == 0 or
                 *std::min_element(flags[nodeListi]->begin(), flags[nodeListi]->end()) == 1);
        }
      }
      END_CONTRACT_SCOPE

      // The ConnectivityMap should be valid too.
      ENSURE(db.connectivityMap().valid());
    }

  // } else {

  //   // We're not connectivity and don't need ghost nodes, so make sure all 
  //   // boundaries are empty.
  //   const vector<Boundary<Dimension>*> boundaries = uniqueBoundaryConditions();
  //   for (ConstBoundaryIterator boundaryItr = boundaries.begin(); 
  //        boundaryItr != boundaries.end();
  //        ++boundaryItr) (*boundaryItr)->reset(db);
  }
}

//------------------------------------------------------------------------------
// Apply the ghost boundary conditions.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Integrator<Dimension>::applyGhostBoundaries(State<Dimension>& state,
                                            StateDerivatives<Dimension>& derivs) const {

//   // Start our work timer.
//   typedef Timing::Time Time;
//   const Time start = Timing::currentTime();

  // Get that DataBase.
  auto& db = mDataBase.get();

  // If we're being rigorous about boundaries, we have to reset the ghost nodes.
  const auto boundaries = uniqueBoundaryConditions();

  // If we didn't call setGhostNodes, then make each boundary update it's 
  // ghost node info (position and H).
  for (auto* boundaryPtr: boundaries) {
    for (auto* nodeListPtr: db.nodeListPtrs()) boundaryPtr->updateGhostNodes(*nodeListPtr);
    boundaryPtr->finalizeGhostBoundary();
  }
  db.updateNeighbors();

  // Iterate over the physics packages, and have them apply ghost boundaries
  // for their state.
  for (auto* physicsPtr: mPhysicsPackages) {
    physicsPtr->applyGhostBoundaries(state, derivs);
  }

  //   // Update the work per node fields.
  //   double deltaLocal = Timing::difference(start, Timing::currentTime());
  //   int numLocalNodes = 0;
  //   for (typename DataBase<Dimension>::NodeListIterator nodeListItr = db.nodeListBegin();
  //        nodeListItr != db.nodeListEnd(); 
  //        ++nodeListItr) numLocalNodes += (*nodeListItr)->numInternalNodes();
  //   const double workPerNode = deltaLocal/(numLocalNodes + 1.0e-30);
  //   FieldList<Dimension, Scalar> work = db.globalWork();
  //   work += workPerNode;

  // Finalize the boundaries.
  this->finalizeGhostBoundaries();
}

//------------------------------------------------------------------------------
// Finalize the ghost boundary conditions.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Integrator<Dimension>::finalizeGhostBoundaries() const {

//   // Start our work timer.
//   typedef Timing::Time Time;
//   const Time start = Timing::currentTime();

  // Get that DataBase.
  //DataBase<Dimension>& db = accessDataBase();

  // If we're being rigorous about boundaries, we have to reset the ghost nodes.
  const vector<Boundary<Dimension>*> boundaries = uniqueBoundaryConditions();
  for (auto* boundaryPtr: boundaries) {
    boundaryPtr->finalizeGhostBoundary();
  }

  //   // Update the work per node fields.
  //   double deltaLocal = Timing::difference(start, Timing::currentTime());
  //   int numLocalNodes = 0;
  //   for (typename DataBase<Dimension>::NodeListIterator nodeListItr = db.nodeListBegin();
  //        nodeListItr != db.nodeListEnd(); 
  //        ++nodeListItr) numLocalNodes += (*nodeListItr)->numInternalNodes();
  //   const double workPerNode = deltaLocal/(numLocalNodes + 1.0e-30);
  //   FieldList<Dimension, Scalar> work = db.globalWork();
  //   work += workPerNode;
}

//------------------------------------------------------------------------------
// Identify the internal nodes in violation of the boundary conditions.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Integrator<Dimension>::setViolationNodes() const {

  // Get that DataBase.
  auto& db = mDataBase.get();

  // Get the complete set of unique boundary conditions.
  const vector<Boundary<Dimension>*> boundaries = uniqueBoundaryConditions();

  // Have each boundary identify the set of nodes that violate it.
  for (auto* boundaryPtr: boundaries) boundaryPtr->setAllViolationNodes(db);

  // Fix neighbor information.
  db.updateNeighbors();
}

//------------------------------------------------------------------------------
// Enforce the boundary conditions.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Integrator<Dimension>::enforceBoundaries(State<Dimension>& state,
                                         StateDerivatives<Dimension>& derivs) const {

  // Have each boundary identify the set of nodes in violation.  This also resets
  // the positions and H's of the nodes to be in compliance.
  setViolationNodes();

  // Iterate over the physics packages, and have them apply ghost boundaries
  // for their state.
  for (auto* physicsPtr: mPhysicsPackages) {
    physicsPtr->enforceBoundaries(state, derivs);
  }
}

//------------------------------------------------------------------------------
// Copy the ghost positions and H tensor from one state to another.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Integrator<Dimension>::copyGhostState(const State<Dimension>& state0,
                                      State<Dimension>& state1) const {
  const FieldList<Dimension, Vector> x0 = state0.fields(HydroFieldNames::position, Vector::zero());
  const FieldList<Dimension, SymTensor> H0 = state0.fields(HydroFieldNames::H, SymTensor::zero());
  FieldList<Dimension, Vector> x1 = state1.fields(HydroFieldNames::position, Vector::zero());
  FieldList<Dimension, SymTensor> H1 = state1.fields(HydroFieldNames::H, SymTensor::zero());
  const auto numNodeLists = x0.numFields();
  CHECK(H0.numFields() == numNodeLists);
  CHECK(x1.numFields() == numNodeLists);
  CHECK(H1.numFields() == numNodeLists);
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto firstGhostNode = x0[nodeListi]->numInternalElements();
    const auto numNodes = x0[nodeListi]->numElements();
    CHECK(H0[nodeListi]->nodeListPtr() == x0[nodeListi]->nodeListPtr());
    CHECK(x1[nodeListi]->nodeListPtr() == x0[nodeListi]->nodeListPtr());
    CHECK(H1[nodeListi]->nodeListPtr() == x0[nodeListi]->nodeListPtr());
    CHECK(H0[nodeListi]->numInternalElements() == firstGhostNode);
    CHECK(x1[nodeListi]->numInternalElements() == firstGhostNode);
    CHECK(H1[nodeListi]->numInternalElements() == firstGhostNode);
    CHECK(H0[nodeListi]->numElements() == numNodes);
    CHECK(x1[nodeListi]->numElements() == numNodes);
    CHECK(H1[nodeListi]->numElements() == numNodes);
    for (auto nodei = firstGhostNode; nodei < numNodes; ++nodei) {
      x1(nodeListi, nodei) = x0(nodeListi, nodei);
      H1(nodeListi, nodei) = H0(nodeListi, nodei);
    }
  }
}

//------------------------------------------------------------------------------
// Dump the current state of the Integrator to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Integrator<Dimension>::
dumpState(FileIO& file, const string& pathName) const {

  // file.write(dtMin(), pathName + "/dtMin");
  // file.write(dtMax(), pathName + "/dtMax");
  // file.write(dtGrowth(), pathName + "/dtGrowth");
  file.write(mLastDt, pathName + "/lastDt");
  file.write(mCurrentTime, pathName + "/currentTime");
  file.write(mCurrentCycle, pathName + "/currentCycle");
  file.write(mDtMultiplier, pathName + "/dtMultiplier");
}  

//------------------------------------------------------------------------------
// Restore the state of the NodeList from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Integrator<Dimension>::
restoreState(const FileIO& file, const string& pathName) {

  // file.read(mDtMin, pathName + "/dtMin");
  // file.read(mDtMax, pathName + "/dtMax");
  // file.read(mDtGrowth, pathName + "/dtGrowth");
  file.read(mLastDt, pathName + "/lastDt");
  file.read(mCurrentTime, pathName + "/currentTime");
  file.read(mCurrentCycle, pathName + "/currentCycle");
  file.read(mDtMultiplier, pathName + "/dtMultiplier");
}

//------------------------------------------------------------------------------
// How should we query a physics package for the time step?
//------------------------------------------------------------------------------
template<typename Dimension>
typename Integrator<Dimension>::TimeStepType
Integrator<Dimension>::
dt(const Physics<Dimension>* pkg,
   const DataBase<Dimension>& dataBase,
   const State<Dimension>& state,
   const StateDerivatives<Dimension>& derivs,
   const Scalar currentTime) const {
  return pkg->dt(dataBase, state, derivs, currentTime);
}

}
