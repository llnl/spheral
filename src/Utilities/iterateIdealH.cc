//------------------------------------------------------------------------------
// Iterate the ideal H algorithm to converge on a new H field.
// This routine replaces the H field in place.
//------------------------------------------------------------------------------
#include "iterateIdealH.hh"
#include "Field/FieldList.hh"
#include "Distributed/allReduce.hh"
#include "Utilities/range.hh"
#include "Distributed/Communicator.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/IncrementBoundedState.hh"
#include "DataBase/ReplaceBoundedState.hh"
#include "Geometry/GeometryRegistrar.hh"
#include "SmoothingScale/ASPHSmoothingScale.hh"
#include "Utilities/Timer.hh"
#include "Utilities/SpheralMessage.hh"

#include <ctime>
using std::vector;
using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {

template<typename Dimension>
void
iterateIdealH(DataBase<Dimension>& dataBase,
              vector<Physics<Dimension>*>& packages,     // Should include smoothing scale algorithm
              const vector<Boundary<Dimension>*>& boundaries,
              const int maxIterations,
              const double tolerance,
              const double nPerhForIteration,
              const bool sphericalStart,
              const bool fixDeterminant) {
  TIME_FUNCTION;

  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  // Start the timing.
  const auto t0 = clock();

  // Extract the state we care about.
  auto pos = dataBase.fluidPosition();
  auto m = dataBase.fluidMass();
  auto rho = dataBase.fluidMassDensity();
  auto H = dataBase.fluidHfield();

  // If we're using the fixDeterminant, take a snapshot of the input H determinants.
  auto Hdet0 = dataBase.newFluidFieldList(double());
  if (fixDeterminant) {
    for (auto itr = dataBase.internalNodeBegin();
         itr != dataBase.internalNodeEnd();
         ++itr) Hdet0(itr) = H(itr).Determinant();
  }

  // Store the input nperh for each NodeList.
  // If we're rescaling the nodes per h for our work, make a cut at it.
  vector<double> nperh0;
  // Pulled divide by nPerhForIteration out of loop to improve optimization
  if (distinctlyGreaterThan(nPerhForIteration, 0.0)) {
    for (auto* nodeListPtr: range(dataBase.fluidNodeListBegin(), dataBase.fluidNodeListEnd())) {
      const auto nperh = nodeListPtr->nodesPerSmoothingScale();
      nperh0.push_back(nperh);
      auto& Hfield = **(H.fieldForNodeList(*nodeListPtr));
      Hfield *= Dimension::rootnu(nperh / nPerhForIteration);
      nodeListPtr->nodesPerSmoothingScale(nPerhForIteration);
    }
  }
  else {
    for (auto* nodeListPtr: range(dataBase.fluidNodeListBegin(), dataBase.fluidNodeListEnd())) {
      const auto nperh = nodeListPtr->nodesPerSmoothingScale();
      nperh0.push_back(nperh);
    }
  }
  CHECK(nperh0.size() == dataBase.numFluidNodeLists());

  // If we are both fixing the H determinant and rescaling the nodes per h,
  // we need a snapshot of the rescaled H determinant as well.
  auto Hdet1 = dataBase.newFluidFieldList(double());
  if (fixDeterminant) {
    for (auto itr = dataBase.internalNodeBegin();
         itr != dataBase.internalNodeEnd();
         ++itr) Hdet1(itr) = H(itr).Determinant();
  }

  // If requested, start by making all the H's round (volume preserving).
  if (sphericalStart) {
    for (auto itr = dataBase.internalNodeBegin();
         itr != dataBase.internalNodeEnd();
         ++itr) {
      const auto Hdeti = H(itr).Determinant();
      const auto Hi = Dimension::rootnu(Hdeti) * SymTensor::one();
      H(itr) = Hi;
    }
  }

  // // Check if we're using ASPH and radialOnly.  If so we'll switch to fixShape for the iteration.
  // auto radialOnly = false;
  // ASPHSmoothingScale<Dimension>* asphPkg = nullptr;
  // for (auto* pkg: packages) {
  //   asphPkg = dynamic_cast<ASPHSmoothingScale<Dimension>*>(pkg);
  //   if (asphPkg != nullptr and asphPkg->radialOnly()) {
  //     radialOnly = true;
  //     asphPkg->radialOnly(false);
  //     asphPkg->fixShape(true);
  //     break;
  //   }
  // }

  // Build a list of flags to indicate which nodes have been completed.
  auto flagNodeDone = dataBase.newFluidFieldList(0, "node completed");

  // Since we don't have a hydro there are a few other fields we need registered.
  auto zerothMoment = dataBase.newFluidFieldList(0.0, HydroFieldNames::massZerothMoment);
  auto firstMoment = dataBase.newFluidFieldList(Vector::zero(), HydroFieldNames::massFirstMoment);
  auto DvDx = dataBase.newFluidFieldList(Tensor::zero(), HydroFieldNames::velocityGradient);
  auto DHDt = dataBase.newFluidFieldList(SymTensor::zero(), IncrementBoundedState<Dimension, SymTensor>::prefix() + HydroFieldNames::H);
  auto H1 = dataBase.newFluidFieldList(SymTensor::zero(), ReplaceBoundedState<Dimension, SymTensor>::prefix() + HydroFieldNames::H);

  // Iterate until we either hit the max iterations or the H's achieve convergence.
  auto maxDeltaH = 2.0*tolerance;
  auto itr = 0;
  while (itr < maxIterations and maxDeltaH > tolerance) {
    ++itr;
    maxDeltaH = 0.0;
    // flagNodeDone = 0;

    // Remove any old ghost node information from the NodeLists.
    for (auto* nodeListPtr: range(dataBase.fluidNodeListBegin(), dataBase.fluidNodeListEnd())) {
      nodeListPtr->numGhostNodes(0);
      nodeListPtr->neighbor().updateNodes();
    }

    // Enforce boundary conditions.
    for (auto* boundaryPtr: boundaries) {
      boundaryPtr->setAllGhostNodes(dataBase);
      boundaryPtr->applyFieldListGhostBoundary(m);
      boundaryPtr->applyFieldListGhostBoundary(rho);
      boundaryPtr->finalizeGhostBoundary();
      for (auto* nodeListPtr: range(dataBase.fluidNodeListBegin(), dataBase.fluidNodeListEnd())) nodeListPtr->neighbor().updateNodes();
    }

    // Update connectivity
    dataBase.updateConnectivityMap(false, false, false);

    // Prepare the state and derivatives
    State<Dimension> state(dataBase, packages);
    state.enrollConnectivityMap(dataBase.connectivityMapPtr(false, false, false));
    StateDerivatives<Dimension> derivs(dataBase, packages);
    state.enroll(pos);
    state.enroll(m);
    state.enroll(rho);
    derivs.enroll(zerothMoment);
    derivs.enroll(firstMoment);
    derivs.enroll(DvDx);
    derivs.enroll(DHDt);
    derivs.enroll(H1);

    // Some methods (ASPH) update both Hideal and H in the finalize, so we make a copy of the state
    // to give the methods
    auto state1 = state;
    state1.copyState();

    // Call the smoothing scale package to get a new vote on the ideal H
    for (auto* pkg: packages) pkg->preStepInitialize(dataBase, state1, derivs);
    for (auto* pkg: packages) pkg->initialize(0.0, 1.0, dataBase, state1, derivs);
    derivs.Zero();
    for (auto* pkg: packages) pkg->evaluateDerivatives(0.0, 1.0, dataBase, state1, derivs);
    for (auto* pkg: packages) pkg->finalizeDerivatives(0.0, 1.0, dataBase, state1, derivs);
    for (auto* pkg: packages) pkg->finalize(0.0, 1.0, dataBase, state1, derivs);
    
    // Set the new H and measure how much it changed
    for (auto [nodeListi_, nodeListPtr_] : enumerate(dataBase.fluidNodeListPtrs())) {     // __clang__
      const auto nodeListi = nodeListi_;                                                  // __clang__
      const auto nodeListPtr = nodeListPtr_;                                              // __clang__
      const auto ni = nodeListPtr->numInternalNodes();

#pragma omp parallel for
      for (auto i = 0u; i < ni; ++i) {
        if (flagNodeDone(nodeListi, i) == 0) {

          // If we are preserving the determinant, do it.
          if (fixDeterminant) {
            H1(nodeListi, i) *= Dimension::rootnu(Hdet1(nodeListi, i)/H1(nodeListi, i).Determinant());
            CHECK(fuzzyEqual(H1(nodeListi, i).Determinant(), Hdet1(nodeListi, i)));
          }

          // Check how much this H has changed.
          const auto H1sqrt = H1(nodeListi, i).sqrt();
          const auto phi = (H1sqrt*H(nodeListi, i).Inverse()*H1sqrt).Symmetric().eigenValues();
          const auto phimin = phi.minElement();
          const auto phimax = phi.maxElement();
          const auto deltaHi = max(abs(phimin - 1.0), abs(phimax - 1.0));
          if (deltaHi <= tolerance) flagNodeDone(nodeListi, i) = 1;
#pragma omp critical
          {
            maxDeltaH = max(maxDeltaH, deltaHi);
          }

// Assign the new H
          H(nodeListi, i) = H1(nodeListi, i);
        }
      }
    }

    // Globally reduce the max H change.
    maxDeltaH = allReduce(maxDeltaH, SPHERAL_OP_MAX);

    // Output the statitics.
    if (maxIterations > 1) SpheralMessage("iterateIdealH: (iteration, deltaH) = (" << itr << ", " << maxDeltaH << ")");
  }

  // If we have rescaled the nodes per h, now we have to iterate the H determinant
  // to convergence with the proper nperh.
  if (distinctlyGreaterThan(nPerhForIteration, 0.0)) {

    // Reset the nperh.
    for (auto [k, nodeListPtr]: enumerate(dataBase.fluidNodeListBegin(), dataBase.fluidNodeListEnd())) {
      CHECK(k < nperh0.size());
      //const double nperh = nperh0[k];
//       Field<Dimension, SymTensor>& Hfield = **(H.fieldForNodeList(**nodeListItr));
//       Hfield *= Dimension::rootnu(nPerhForIteration/nperh);
      nodeListPtr->nodesPerSmoothingScale(nperh0[k]);
    }
  }

  // If we're fixing the determinant, restore them.
  if (fixDeterminant) {
    for (auto itr = dataBase.internalNodeBegin();
         itr != dataBase.internalNodeEnd();
         ++itr) {
      H(itr) *= Dimension::rootnu(Hdet0(itr)/H(itr).Determinant());
      ENSURE(fuzzyEqual(H(itr).Determinant(), Hdet0(itr), 1.0e-10));
    }
  }

  // Leave the boundary conditions properly enforced.
  for (auto* nodeListPtr: range(dataBase.fluidNodeListBegin(), dataBase.fluidNodeListEnd())) {
    nodeListPtr->numGhostNodes(0);
    nodeListPtr->neighbor().updateNodes();
  }
  for (auto* boundaryPtr: range(boundaries.begin(), boundaries.end())) {
    boundaryPtr->setAllGhostNodes(dataBase);
    boundaryPtr->finalizeGhostBoundary();
    for (auto* nodeListPtr: range(dataBase.fluidNodeListBegin(), dataBase.fluidNodeListEnd())) {
      nodeListPtr->neighbor().updateNodes();
    }
  }

  for (auto* boundaryPtr: range(boundaries.begin(), boundaries.end())) boundaryPtr->applyFieldListGhostBoundary(m);
  for (auto* boundaryPtr: range(boundaries.begin(), boundaries.end())) boundaryPtr->finalizeGhostBoundary();

  // // Restore ASPH radialOnly choice if necessary
  // if (radialOnly) {
  //   CHECK(asphPkg != nullptr);
  //   asphPkg->radialOnly(true);
  //   asphPkg->fixShape(false);
  // }

  // Report the final timing.
  const auto t1 = clock();
  if (maxIterations > 1) SpheralMessage("iterateIdealH: required a total of " << ((t1 - t0)/CLOCKS_PER_SEC) << " seconds.");
}

}

