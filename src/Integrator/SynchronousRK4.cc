//---------------------------------Spheral++----------------------------------//
// SynchronousRK4 -- Advance the set of Physics packages in time using fourth
// order Runge-Kutta.  All packages are advanced at one timestep simultaneously
// each step, i.e., synchronously.
//
// Created by JMO, Thu Jun 14 10:22:47 PDT 2012
//----------------------------------------------------------------------------//
#include "DataOutput/Restart.hh"
#include "SynchronousRK4.hh"
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
// Construct with the given DataBase and Physics packages.
//------------------------------------------------------------------------------
template<typename Dimension>
SynchronousRK4<Dimension>::
SynchronousRK4(DataBase<Dimension>& dataBase,
               const vector<Physics<Dimension>*>& physicsPackages):
  Integrator<Dimension>(dataBase, physicsPackages) {
}

//------------------------------------------------------------------------------
// Take a step.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
SynchronousRK4<Dimension>::
step(typename Dimension::Scalar maxTime,
     State<Dimension>& state,
     StateDerivatives<Dimension>& derivs1) {

  // Get the current time and data base.
  Scalar t = this->currentTime();
  DataBase<Dimension>& db = this->accessDataBase();

  // Initalize the integrator.
  this->preStepInitialize(state, derivs1);

  // Determine the minimum timestep across all packages.
  const Scalar dt = this->selectDt(min(this->dtMin(), maxTime - t),
                                   min(this->dtMax(), maxTime - t),
                                   state,
                                   derivs1);

  // Zero out the derivatives, and make some independent copies
  StateDerivatives<Dimension> derivs2(derivs1), derivs3(derivs1), derivs4(derivs1);
  derivs2.copyState();
  derivs3.copyState();
  derivs4.copyState();

  // Make a copy of the state we'll use for our intermediate estimates.
  State<Dimension> tmpstate(state);
  tmpstate.copyState();

  // Stage 1:
  // Get derivs1(t_n, state(t_n))
  this->initializeDerivatives(t, dt, state, derivs1);
  derivs1.Zero();
  this->evaluateDerivatives(t, dt, db, state, derivs1);
  this->finalizeDerivatives(t, dt, db, state, derivs1);

  // Stage 2: 
  // Get derivs2(t_n + 0.5*dt, state(t_n + 0.5*dt*derivs1))
  tmpstate.update(derivs1, 0.5*dt, t, 0.5*dt);
  this->applyGhostBoundaries(tmpstate, derivs1);
  this->finalizeGhostBoundaries();
  this->postStateUpdate(t + 0.5*dt, 0.5*dt, db, tmpstate, derivs1);
  this->initializeDerivatives(t + 0.5*dt, 0.5*dt, tmpstate, derivs2);
  derivs2.Zero();
  this->evaluateDerivatives(t + 0.5*dt, 0.5*dt, db, tmpstate, derivs2);
  this->finalizeDerivatives(t + 0.5*dt, 0.5*dt, db, tmpstate, derivs2);

  // Check if the timestep is still a good idea...
  if (this->allowDtCheck()) {
    const auto dtnew = this->selectDt(min(this->dtMin(), maxTime - t),
                                      min(this->dtMax(), maxTime - t),
                                      tmpstate,
                                      derivs2);
    if (dtnew < this->dtCheckFrac()*dt) {
      this->currentTime(t);
      return false;
    }
  }

  // Stage 3: 
  // Get derivs3(t_n + 0.5*dt, state(t_n + 0.5*dt*derivs2))
  tmpstate = state;
  tmpstate.copyState();
  tmpstate.update(derivs2, 0.5*dt, t, 0.5*dt);
  this->applyGhostBoundaries(tmpstate, derivs2);
  this->finalizeGhostBoundaries();
  this->postStateUpdate(t + 0.5*dt, 0.5*dt, db, tmpstate, derivs2);
  this->initializeDerivatives(t + 0.5*dt, 0.5*dt, tmpstate, derivs3);
  derivs3.Zero();
  this->evaluateDerivatives(t + 0.5*dt, 0.5*dt, db, tmpstate, derivs3);
  this->finalizeDerivatives(t + 0.5*dt, 0.5*dt, db, tmpstate, derivs3);

  // Check if the timestep is still a good idea...
  if (this->allowDtCheck()) {
    const auto dtnew = this->selectDt(min(this->dtMin(), maxTime - t),
                                      min(this->dtMax(), maxTime - t),
                                      tmpstate,
                                      derivs3);
    if (dtnew < this->dtCheckFrac()*dt) {
      this->currentTime(t);
      return false;
    }
  }

  // Stage 4: 
  // Get derivs3(t_n + dt, state(t_n + dt*derivs3))
  tmpstate = state;
  tmpstate.copyState();
  tmpstate.update(derivs3, dt, t, dt);
  this->applyGhostBoundaries(tmpstate, derivs3);
  this->finalizeGhostBoundaries();
  this->postStateUpdate(t + dt, dt, db, tmpstate, derivs3);
  this->initializeDerivatives(t + dt, dt, tmpstate, derivs4);
  derivs4.Zero();
  this->evaluateDerivatives(t + dt, dt, db, tmpstate, derivs4);
  this->finalizeDerivatives(t + dt, dt, db, tmpstate, derivs4);

  // Check if the timestep is still a good idea...
  if (this->allowDtCheck()) {
    const auto dtnew = this->selectDt(min(this->dtMin(), maxTime - t),
                                      min(this->dtMax(), maxTime - t),
                                      tmpstate,
                                      derivs4);
    if (dtnew < this->dtCheckFrac()*dt) {
      this->currentTime(t);
      return false;
    }
  }

  // Advance.
  // Now we can apply the RK4 algorithm to advance the actual state over the full time step.
  // Conceptually we are doing:
  //   state(t_n + dt) = state(t_n) + dt/6*(derivs1 + 2*derivs2 + 2*derivs3 + derivs4)
  state.update(derivs1, dt/6.0, t, dt);
  state.update(derivs2, dt/3.0, t, dt);
  state.update(derivs3, dt/3.0, t, dt);
  state.update(derivs4, dt/6.0, t, dt);
  this->applyGhostBoundaries(state, derivs4);
  this->finalizeGhostBoundaries();
  this->postStateUpdate(t + dt, dt, db, state, derivs4);

  // Apply any physics specific finalizations.
  this->postStepFinalize(t + dt, dt, state, derivs4);

  // Enforce boundaries.
  this->enforceBoundaries(state, derivs4);

  // Set the new current time and last time step.
  this->currentCycle(this->currentCycle() + 1);
  this->currentTime(t + dt);
  this->lastDt(dt);
  return true;
}

}
