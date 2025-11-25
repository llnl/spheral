//---------------------------------Spheral++----------------------------------//
// Second order predictor corrector time integrator.
// Based on the predictor corrector scheme described in Monaghans tensile
// instability paper (Monaghan 2000, JCP 159, 290.)
//
// Created by JMO, Wed Dec  4 21:39:36 PST 2002
//----------------------------------------------------------------------------//
#include "DataOutput/Restart.hh"
#include "PredictorCorrector.hh"
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
PredictorCorrector<Dimension>::
PredictorCorrector(DataBase<Dimension>& dataBase,
               const vector<Physics<Dimension>*>& physicsPackages):
  Integrator<Dimension>(dataBase, physicsPackages) {
}

//------------------------------------------------------------------------------
// Take a step.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
PredictorCorrector<Dimension>::
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

  // Zero out the derivatives.
  derivs.Zero();

  // Evaluate the beginning of step derivatives.
  this->initializeDerivatives(t, dt, state, derivs);
  this->evaluateDerivatives(t, dt, db, state, derivs);
  this->finalizeDerivatives(t, dt, db, state, derivs);

  // Copy the beginning of step state and derivatives.
  State<Dimension> state0(state);
  StateDerivatives<Dimension> derivs0(derivs);
  state0.copyState();
  derivs0.copyState();

  // Predictor step -- trial advance to the end of step.
  state.update(derivs, dt, t, dt);

  // Enforce Boundary conditions.
  this->enforceBoundaries(state, derivs);
  this->applyGhostBoundaries(state, derivs);
  this->finalizeGhostBoundaries();
                               
  // Do any physics specific stuff relating to the fact the state was just updated.
  this->postStateUpdate(t + dt, dt, db, state, derivs);

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

  // Loop over the physics packages and perform any necessary initializations.
  this->initializeDerivatives(t + dt, dt, state, derivs);

  // Zero out the stored derivatives.
  derivs.Zero();

  // Now loop over the packages and evaluate the derivatives at the
  // predicted end of step state.
  this->evaluateDerivatives(t + dt, dt, db, state, derivs);
  this->finalizeDerivatives(t + dt, dt, db, state, derivs);

  // Apply the corrector stage.
  const Scalar hdt = 0.5*dt;
  this->copyGhostState(state, state0);
  state.assign(state0);
  state.update(derivs0, hdt, t, hdt);
  this->applyGhostBoundaries(state, derivs);
  state.update(derivs, hdt, t + hdt, hdt);
  this->currentTime(t + dt);

  // Enforce boundaries.
  this->enforceBoundaries(state, derivs);
  this->applyGhostBoundaries(state, derivs);
  this->finalizeGhostBoundaries();

  // Do any physics specific stuff relating to the fact the state was just updated.
  this->postStateUpdate(t + dt, dt, db, state, derivs);

  // Apply any physics specific finalizations.
  this->postStepFinalize(t + dt, dt, state, derivs);

  // Set the new current time and last time step.
  this->currentCycle(this->currentCycle() + 1);
  this->lastDt(dt);
  return true;
}

}
