//---------------------------------Spheral++----------------------------------//
// BackwardEuler -- Advance the set of Physics packages in time using first
// order Backward Euler method.  All packages are advanced at one
// timestep simultaneously each step, i.e., synchronously.
//
// Created by JMO, Mon Oct 21 14:32:05 PDT 2024
//----------------------------------------------------------------------------//
#include "Integrator/BackwardEuler.hh"
#include "Integrator/ImplicitIntegrationVectorOperator.hh"
#include "Solvers/KINSOL.hh"

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
BackwardEuler<Dimension>::
BackwardEuler(DataBase<Dimension>& dataBase,
              const vector<Physics<Dimension>*> physicsPackages,
              const Scalar beta,
              const Scalar ftol,
              const Scalar tol,
              const size_t maxIterations):
  ImplicitIntegrator<Dimension>(dataBase, physicsPackages, tol),
  mBeta(beta),
  mftol(ftol),
  mtM2(-1.0),
  mtM1(-1.0),
  mMaxIters(maxIterations),
  mSolutionM2(),
  mSolutionM1() {
  this->allowDtCheck(true);
}

//------------------------------------------------------------------------------
// Take a step.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
BackwardEuler<Dimension>::
step(typename Dimension::Scalar maxTime,
     State<Dimension>& state,
     StateDerivatives<Dimension>& derivs) {

  // Get the current time and data base.
  Scalar t = this->currentTime();
  DataBase<Dimension>& db = this->accessDataBase();

  // Initalize the integrator.
  this->preStepInitialize(state, derivs);

  // Determine the minimum timestep across all packages.
  const Scalar dt = this->selectDt(min(this->dtMin(), maxTime - t),
                                   min(this->dtMax(), maxTime - t),
                                   state,
                                   derivs);

  // Evaluate the beginning of step derivatives.
  this->initializeDerivatives(t, dt, state, derivs);
  derivs.Zero();
  this->evaluateDerivatives(t, dt, db, state, derivs);
  this->finalizeDerivatives(t, dt, db, state, derivs);

  // Make a copy of the initial state and derivatives
  State<Dimension> state0(state);
  StateDerivatives<Dimension> derivs0(derivs);
  state0.copyState();
  derivs0.copyState();
  CHECK(state0 == state);

  // If we have not yet accrued enough previous step information to make a
  // prediction about the next state, just advance explicitly
  if (mDtMultiplier < 1.0 or this->currentCycle() < 4) {

    //..........................................................................
    // Do a standard RK2 step
    ++mNumExplicitSteps;
    const auto hdt = 0.5*dt;
    
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

    //..........................................................................
    // Make a parabolic fitted prediction for the initial guess.  We use the Lagrange
    // interpolation formula for a second-order polynomial.
    vector<double> solution0;
    // state.serializeIndependentData(solution0);  // t_n
    // const auto n = solution0.size();
    // CHECK(mSolutionM2.size() == n);             // t_{n-2}
    // CHECK(mSolutionM1.size() == n);             // t_{n-1}
    // CHECK(mtM2 < mtM1 and mtM1 < t);
    // const auto x1 = mtM2, x2 = mtM1, x3 = t, x = t + dt;
    // // cerr << "Solution t_{n-2}:";
    // // for (auto i = 0u; i < n; ++i) cerr << " " << mSolutionM2[i];
    // // cerr << endl
    // //      << "Solution t_{n-1}:";
    // // for (auto i = 0u; i < n; ++i) cerr << " " << mSolutionM1[i];
    // // cerr << endl
    // //      << "Solution   t_{n}:";
    // // for (auto i = 0u; i < n; ++i) cerr << " " << solution0[i];
    // for (auto i = 0u; i < n; ++i) {
    //   solution0[i] = (mSolutionM2[i]*(x - x2)*(x - x3)/((x1 - x2)*(x1 - x3)) +
    //                   mSolutionM1[i]*(x - x1)*(x - x3)/((x2 - x1)*(x2 - x3)) +
    //                   solution0[i]  *(x - x1)*(x - x2)/((x3 - x1)*(x3 - x2)));
    // }
    // // cerr << endl
    // //      << "Solution t_{n+1}:";
    // // for (auto i = 0u; i < n; ++i) cerr << " " << solution0[i];
    // // cerr << endl;

    // Initial Forward Euler prediction for the end of timestep state
    // state.deserializeIndependentData(solution0);
    state.update(derivs, dt, t, dt, false);
    this->applyGhostBoundaries(state, derivs);
    this->finalizeGhostBoundaries();
    this->postStateUpdate(t + dt, dt, db, state, derivs);

    // Derivatives at the ForwardEuler prediction
    this->initializeDerivatives(t + dt, dt, state, derivs);
    derivs.Zero();
    this->evaluateDerivatives(t + dt, dt, db, state, derivs);
    this->finalizeDerivatives(t + dt, dt, db, state, derivs);

    // Now do a Crank-Nicolson prediction for the solution at t+dt
    state.assign(state0);
    state.update(derivs0, 0.5*dt, t,          0.5*dt);
    state.update(derivs,  0.5*dt, t + 0.5*dt, 0.5*dt);
    this->applyGhostBoundaries(state, derivs);
    this->finalizeGhostBoundaries();
    this->postStateUpdate(t + dt, dt, db, state, derivs);
    state.serializeIndependentData(solution0);

    // Build a solver
    KINSOL solver;
    solver.fnormtol(mftol);
    solver.scsteptol(this->convergenceTolerance());
    solver.numMaxIters(mMaxIters);

    // Build the VectorOperator
    ImplicitIntegrationVectorOperator op(t, dt, mBeta, state0, derivs0, state, derivs, *this);

    // Iterate on the new state
    // std::vector<double> solution;
    // state.serializeIndependentData(solution);
    auto numIters = solver.solve(op, solution0);
    state.assign(state0);
    if (numIters == size_t(mMaxIters)) {
      return false;
    }

    // Unpack the solution into the final state.
    ++mNumImplicitSteps;
    state.update(derivs, dt, t, dt, false);
    this->applyGhostBoundaries(state, derivs);
    this->finalizeGhostBoundaries();
    this->postStateUpdate(t + dt, dt, db, state, derivs);

  }
  
  // Apply any physics specific finalizations.
  this->postStepFinalize(t + dt, dt, state, derivs);

  // Enforce boundaries.
  this->enforceBoundaries(state, derivs);

  // Copy the final state for the next step
  // cerr << "TIMES: " << mtM2 << " " << mtM1 << " " << t << " " << t + dt << endl;
  mtM2 = mtM1;
  mtM1 = t;
  mSolutionM2 = mSolutionM1;
  state0.serializeIndependentData(mSolutionM1);

  // Set the new current time and last time step.
  this->currentTime(t + dt);
  this->currentCycle(this->currentCycle() + 1);
  this->lastDt(dt);

  return true;
}

}
