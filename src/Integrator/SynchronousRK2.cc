//---------------------------------Spheral++----------------------------------//
// SynchronousRK2 -- Advance the set of Physics packages in time using second
// order Runge-Kutta.  All packages are advanced at one timestep simultaneously
// each step, i.e., synchronously.
//
// Created by JMO, Thu Jun  1 22:26:45 PDT 2000
//----------------------------------------------------------------------------//
#include "DataOutput/Restart.hh"
#include "SynchronousRK2.hh"
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
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {

//------------------------------------------------------------------------------
// Construct with the given DataBase and Physics packages.
//------------------------------------------------------------------------------
template<typename Dimension>
SynchronousRK2<Dimension>::
SynchronousRK2(DataBase<Dimension>& dataBase,
               const vector<Physics<Dimension>*>& physicsPackages):
  Integrator<Dimension>(dataBase, physicsPackages) {
}

//------------------------------------------------------------------------------
// Take a step.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
SynchronousRK2<Dimension>::
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
  const double hdt = 0.5*dt;

  // Evaluate the beginning of step derivatives.
  this->initializeDerivatives(t, hdt, state, derivs);
  derivs.Zero();
  this->evaluateDerivatives(t, hdt, db, state, derivs);
  this->finalizeDerivatives(t, hdt, db, state, derivs);

  // Copy the beginning of step state.
  State<Dimension> state0(state);
  state0.copyState();

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

  // Check if the timestep is still a good idea...
  if (this->allowDtCheck()) {
    const auto dtnew = this->selectDt(min(this->dtMin(), maxTime - t),
                                      min(this->dtMax(), maxTime - t),
                                      state,
                                      derivs);
    if (dtnew < this->dtCheckFrac()*dt) {
      this->currentTime(t);
      state.assign(state0);
      return false;
    }
  }

  // Advance the state from the beginning of the cycle using the midpoint 
  // derivatives.
  // this->copyGhostState(state, state0);
  state.assign(state0);
  state.update(derivs, dt, t, dt);
  this->currentTime(t + dt);
  this->applyGhostBoundaries(state, derivs);
  this->finalizeGhostBoundaries();
  this->postStateUpdate(t + dt, dt, db, state, derivs);

  // Apply any physics specific finalizations.
  this->postStepFinalize(t + dt, dt, state, derivs);

  // Enforce boundaries.
  this->enforceBoundaries(state, derivs);

  // Set the new current time and last time step.
  this->currentCycle(this->currentCycle() + 1);
  this->lastDt(dt);
  return true;
}

}
