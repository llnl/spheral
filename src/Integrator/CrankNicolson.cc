//---------------------------------Spheral++----------------------------------//
// CrankNicolson -- Advance the set of Physics packages in time using first
// order Backward Euler method.  All packages are advanced at one
// timestep simultaneously each step, i.e., synchronously.
//
// Created by JMO, Mon Oct 21 14:32:05 PDT 2024
//----------------------------------------------------------------------------//
#include "Integrator/CrankNicolson.hh"
#include "DataOutput/Restart.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "Physics/Physics.hh"
#include "Utilities/DBC.hh"

using std::vector;
using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename Dimension>
CrankNicolson<Dimension>::
CrankNicolson(DataBase<Dimension>& dataBase,
              const vector<Physics<Dimension>*> physicsPackages,
              const Scalar beta,
              const Scalar tol,
              const size_t maxIterations):
  ImplicitIntegrator<Dimension>(dataBase, physicsPackages, tol),
  mBeta(beta),
  mMaxIterations(maxIterations) {
  this->allowDtCheck(true);
}

//------------------------------------------------------------------------------
// Take a step.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
CrankNicolson<Dimension>::
step(typename Dimension::Scalar maxTime,
     State<Dimension>& state,
     StateDerivatives<Dimension>& derivs) {

  // Get the current time and data base.
  const auto t = this->currentTime();
  DataBase<Dimension>& db = this->accessDataBase();
  const auto tol = this->convergenceTolerance();

  // Initalize the integrator.
  this->preStepInitialize(state, derivs);

  // Determine the minimum timestep across all packages.
  const Scalar dt = this->selectDt(min(this->dtMin(), maxTime - t),
                                   min(this->dtMax(), maxTime - t),
                                   state,
                                   derivs);
  const auto hdt = 0.5*dt;
  const auto hdt1 = hdt*(1.0 - mBeta);
  const auto hdt2 = hdt*mBeta;
  CHECK(hdt1 + hdt2 == hdt);

  // We now assume we can reuse the last steps derivative estimate as the
  // beginning of step for this iteration.  Similar assumption used in
  // CheapSynchronousRK2.  If you want to be more rigorous uncomment the
  // following block to explicitly determine the beginning of step derivatives.

  // Evaluate the beginning of step derivatives.
  this->initializeDerivatives(t, dt, state, derivs);
  derivs.Zero();
  this->evaluateDerivatives(t, dt, db, state, derivs);
  this->finalizeDerivatives(t, dt, db, state, derivs);

  // Make a copy of the initial state and derivatives
  State<Dimension> state0(state), state1(state);
  StateDerivatives<Dimension> derivs0(derivs), derivs1(derivs);
  state0.copyState();
  state1.copyState();
  derivs0.copyState();
  derivs1.copyState();
  CHECK(state0 == state);

  // If we have not yet accrued enough previous step information to make a
  // prediction about the next state, just advance explicitly
  if (mDtMultiplier < 1.0) {
    ++mNumExplicitSteps;

    //..........................................................................
    // Do a standard RK2 step
    // Trial advance the state to the mid timestep point.
    state.update(derivs, hdt, t, hdt);
    this->currentTime(t + hdt);
    this->applyGhostBoundaries(state, derivs);
    this->finalizeGhostBoundaries();
    this->postStateUpdate(t + hdt, hdt, db, state, derivs);

    // Evaluate the derivatives at the trial midpoint conditions.
    this->initializeDerivatives(t + hdt, hdt, state, derivs);
    derivs.Zero();
    this->evaluateDerivatives(t + hdt, hdt, db, state, derivs);
    this->finalizeDerivatives(t + hdt, hdt, db, state, derivs);

    // Advance the state from the beginning of the cycle using the midpoint 
    // derivatives.
    state.assign(state0);
    state.update(derivs, dt, t, dt);
    this->currentTime(t + dt);
    this->applyGhostBoundaries(state, derivs);
    this->finalizeGhostBoundaries();
    this->postStateUpdate(t + dt, dt, db, state, derivs);

  } else {

    // Initial Forward Euler prediction for the end of timestep state
    state.update(derivs, dt, t, dt, false);
    this->applyGhostBoundaries(state, derivs);
    this->finalizeGhostBoundaries();
    this->postStateUpdate(t + dt, dt, db, state, derivs);

    // Initial independent variable vector
    vector<double> solution1, solution;
    // state0.serializeIndependentData(solution0);
    // const auto n = solution0.size();

    // Iterate!
    auto done = false;
    size_t iterations = 0u;
    while (iterations++ < mMaxIterations and not done) {

      // Last pass on new state
      state1.assign(state);
      if (mBeta < 1.0) derivs1.assign(derivs);
      state1.serializeIndependentData(solution1);
      // const auto n = solution1.size();

      // Derivatives at the last prediction
      this->initializeDerivatives(t + dt, dt, state, derivs);
      derivs.Zero();
      this->evaluateDerivatives(t + dt, dt, db, state, derivs);
      this->finalizeDerivatives(t + dt, dt, db, state, derivs);

      // Estimate new n+1 solution
      // Are we blending old and new solutions?
      state.assign(state0);
      if (mBeta < 1.0) {
        // state.serializeIndependentData(solution);
        // CHECK(solution.size() == n);
        // for (auto i = 0u; i < n; ++i) solution[i] = (1.0 - mBeta)*solution1[i] + mBeta*solution[i];
        // state.deserializeIndependentData(solution);
        state.update(derivs0, hdt,    t,              hdt);
        state.update(derivs1, hdt1,   t + hdt,        hdt1);
        state.update(derivs,  hdt2,   t + hdt + hdt1, hdt2);
      } else {
        state.update(derivs0, hdt, t,       hdt);
        state.update(derivs,  hdt, t + hdt, hdt);
      }

      // Finish state update
      this->applyGhostBoundaries(state, derivs);
      this->finalizeGhostBoundaries();
      this->postStateUpdate(t + dt, dt, db, state, derivs);

      // Compare for convergence
      const auto maxResidual = this->computeResiduals(state, state1);
      done = maxResidual < tol;
      SpheralMessage("=============> CrankNicolson: " << iterations << "/" << mMaxIterations << " : " << maxResidual << "/" << tol);
    }

    // Did we succeed?
    if (not done) {
      if (not this->verbose()) this->computeResiduals(state, state1, true);
      state.assign(state0);
      return false;
    }
    ++mNumImplicitSteps;
  }
  
  // Apply any physics specific finalizations.
  this->postStepFinalize(t + dt, dt, state, derivs);

  // Enforce boundaries.
  this->enforceBoundaries(state, derivs);

  // Set the new current time and last time step.
  this->currentTime(t + dt);
  this->currentCycle(this->currentCycle() + 1);
  this->lastDt(dt);

  return true;
}

}
