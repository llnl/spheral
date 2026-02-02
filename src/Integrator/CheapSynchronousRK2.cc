//---------------------------------Spheral++----------------------------------//
// CheapSynchronousRK2 -- Advance the set of Physics packages in time using 
// my standard cheat on the second order Runge-Kutta.  We just skip the end of
// step derivative evaluation, and assume that we can use the previous cycles
// mid-step derivatives to determine this steps mid-step state.
//
// Created by JMO, Mon Jun 12 17:52:27 PDT 2000
//----------------------------------------------------------------------------//
#include "DataOutput/Restart.hh"
#include "CheapSynchronousRK2.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "Physics/Physics.hh"
#include "Utilities/DBC.hh"
#include "Utilities/Timer.hh"

using std::vector;
using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {

//------------------------------------------------------------------------------
// Construct with the given DataBase and Physics packages.
//------------------------------------------------------------------------------
template<typename Dimension>
CheapSynchronousRK2<Dimension>::
CheapSynchronousRK2(DataBase<Dimension>& dataBase,
               const vector<Physics<Dimension>*>& physicsPackages):
  Integrator<Dimension>(dataBase, physicsPackages) {
}

//------------------------------------------------------------------------------
// Take a step.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
CheapSynchronousRK2<Dimension>::
step(typename Dimension::Scalar maxTime,
     State<Dimension>& state,
     StateDerivatives<Dimension>& derivs) {

  TIME_BEGIN("CheapRK2");

  // Get the current time and data base.
  Scalar t = this->currentTime();
  DataBase<Dimension>& db = this->accessDataBase();

  // Initalize the integrator.
  TIME_BEGIN("CheapRK2PreInit");
  this->preStepInitialize(state, derivs);
  this->initializeDerivatives(t, 0.0, state, derivs);
  TIME_END("CheapRK2PreInit");

  // Determine the minimum timestep across all packages.
  TIME_BEGIN("CheapRK2Dt");
  const Scalar dt = this->selectDt(min(this->dtMin(), maxTime - t),
                                   min(this->dtMax(), maxTime - t),
                                   state,
                                   derivs);
  const double hdt = 0.5*dt;
  TIME_END("CheapRK2Dt");

  // Copy the beginning of step state.
  TIME_BEGIN("CheapRK2CopyState");
  State<Dimension> state0(state);
  state0.copyState();
  TIME_END("CheapRK2CopyState");

  // Trial advance the state to the mid timestep point.
  TIME_BEGIN("CheapRK2MidStep");
  state.timeAdvanceOnly(true);
  state.update(derivs, hdt, t, hdt);
  this->currentTime(t + hdt);
  this->applyGhostBoundaries(state, derivs);
  this->finalizeGhostBoundaries();
  this->postStateUpdate(t + hdt, hdt, db, state, derivs);
  TIME_END("CheapRK2MidStep");

  // Evaluate the derivatives at the midpoint.
  TIME_BEGIN("CheapRK2EvalDerivs");
  this->initializeDerivatives(t + hdt, hdt, state, derivs);
  derivs.Zero();
  this->evaluateDerivatives(t + hdt, hdt, db, state, derivs);
  this->finalizeDerivatives(t + hdt, hdt, db, state, derivs);
  TIME_END("CheapRK2EvalDerivs");

  // Check if the timestep is still a good idea...
  if (this->allowDtCheck()) {
    const auto dtnew = this->selectDt(min(this->dtMin(), maxTime - t),
                                      min(this->dtMax(), maxTime - t),
                                      state,
                                      derivs);
    if (dtnew < 0.5*dt) {
      state.assign(state0);
      return false;
    }
  }

  // Advance the state from the beginning of the cycle using the midpoint derivatives.
  TIME_BEGIN("CheapRK2EndStep");
  state.timeAdvanceOnly(false);
  //this->copyGhostState(state, state0);
  state.assign(state0);
  state.update(derivs, dt, t, dt);
  this->currentTime(t + dt);
  this->applyGhostBoundaries(state, derivs);
  this->finalizeGhostBoundaries();
  this->postStateUpdate(t + dt, dt, db, state, derivs);
  TIME_END("CheapRK2EndStep");

  // Apply any physics specific finalizations.
  TIME_BEGIN("CheapRK2Finalize");
  this->postStepFinalize(t + dt, dt, state, derivs);
  TIME_END("CheapRK2Finalize");

  // Enforce boundaries.
  TIME_BEGIN("CheapRK2EnforceBound");
  this->enforceBoundaries(state, derivs);
  TIME_END("CheapRK2EnforceBound");

  // Set the new current time and last time step.
  this->currentCycle(this->currentCycle() + 1);
  this->lastDt(dt);
  TIME_END("CheapRK2");

  return true;
}

}
