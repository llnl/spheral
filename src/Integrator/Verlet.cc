//---------------------------------Spheral++----------------------------------//
// Verlet -- Advance the set of Physics packages in time using the second
// order Verlet algorithm.  This method is symplectic in the absence of 
// dissipation.
//
// Based on the description in 
// Monaghan JJ. Smoothed particle hydrodynamics. Reports on progress in physics. 2005
//
// Created by JMO, Sat Aug 23 10:33:33 PDT 2014
//----------------------------------------------------------------------------//
#include "DataOutput/Restart.hh"
#include "Verlet.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "Physics/Physics.hh"
#include "Hydro/HydroFieldNames.hh"
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
Verlet<Dimension>::
Verlet(DataBase<Dimension>& dataBase,
               const vector<Physics<Dimension>*>& physicsPackages):
  Integrator<Dimension>(dataBase, physicsPackages) {
}

//------------------------------------------------------------------------------
// Take a step.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
Verlet<Dimension>::
step(typename Dimension::Scalar maxTime,
     State<Dimension>& state,
     StateDerivatives<Dimension>& derivs) {

  TIME_BEGIN("Verlet");

  // Get the current time and data base.
  auto  t = this->currentTime();
  auto& db = this->accessDataBase();

  // Initalize the integrator.
  TIME_BEGIN("VerletPreInit");
  this->preStepInitialize(state, derivs);
  TIME_END("VerletPreInit");

  // Copy the beginning of step positions.
  TIME_BEGIN("VerletCopyPos0");
  auto pos0 = state.fields(HydroFieldNames::position, Vector::zero());
  pos0.copyFields();
  TIME_END("VerletCopyPos0");

  // Determine the minimum timestep across all packages.
  TIME_BEGIN("VerletDt");
  const Scalar dtMin = min(this->dtMin(), maxTime - t);
  const Scalar dtMax = min(this->dtMax(), maxTime - t);
  const Scalar dt0 = this->selectDt(dtMin, dtMax, state, derivs);
  const Scalar hdt0 = 0.5*dt0;
  const auto dtcheck = this->allowDtCheck();
  const auto dtcheckFrac = this->dtCheckFrac();
  TIME_END("VerletDt");

  // If we're doing dt checking, we need to copy the initial state.
  std::unique_ptr<State<Dimension>> state0;
  if (dtcheck) {
    TIME_BEGIN("VerletCopyState0");
    state0 = std::make_unique<State<Dimension>>(state);
    state0->copyState();
    TIME_END("VerletCopyState0");
  }

  // Evaluate the beginning of step derivatives.
  TIME_BEGIN("VerletEvalDerivs1");
  this->initializeDerivatives(t, dt0, state, derivs);
  derivs.Zero();
  this->evaluateDerivatives(t, dt0, db, state, derivs);
  this->finalizeDerivatives(t, dt0, db, state, derivs);
  TIME_END("VerletEvalDerivs1");

  // Predict state at the mid-point.
  // state.timeAdvanceOnly(true);
  TIME_BEGIN("VerletPredict1");
  state.update(derivs, hdt0, t, dt0);
  this->enforceBoundaries(state, derivs);
  this->applyGhostBoundaries(state, derivs);
  this->finalizeGhostBoundaries();
  this->postStateUpdate(t + hdt0, hdt0, db, state, derivs);
  TIME_END("VerletPredict1");

  // Check if the timestep is still a good idea...
  if (dtcheck) {
    TIME_BEGIN("VerletDtCheck");
    const auto dtnew = this->selectDt(dtMin,
                                      dtMax,
                                      state,
                                      derivs);
    if (dtnew < dtcheckFrac*dt0) {
      this->currentTime(t);
      state.assign(*state0);
      return false;
      TIME_END("VerletDtCheck");
    }
    TIME_END("VerletDtCheck");
  }

  // Copy the mid-point state.
  TIME_BEGIN("VerletMidPointCopy");
  State<Dimension> state12(state);
  state12.copyState();
  TIME_END("VerletMidPointCopy");

  // Advance the position to the end of step using the half-step velocity.
  TIME_BEGIN("VerletPredict2");
  auto vel12 = state.fields(HydroFieldNames::velocity, Vector::zero());
  pos0 += dt0*vel12;

  // Predict state at the end-point, but override the positions with our time-centered prediction.
  state.update(derivs, hdt0, t, dt0);
  {
    auto pos = state.fields(HydroFieldNames::position, Vector::zero());
    pos.assignFields(pos0);
  }
  this->enforceBoundaries(state, derivs);
  this->applyGhostBoundaries(state, derivs);
  this->finalizeGhostBoundaries();
  this->postStateUpdate(t + dt0, dt0, db, state, derivs);
  TIME_END("VerletPredict2");

  // Evaluate the derivatives at the predicted end-point.
  TIME_BEGIN("VerletEvalDerivs2");
  this->currentTime(t + dt0);
  this->initializeDerivatives(t + dt0, dt0, state, derivs);
  derivs.Zero();
  this->evaluateDerivatives(t + dt0, dt0, db, state, derivs);
  this->finalizeDerivatives(t + dt0, dt0, db, state, derivs);
  TIME_END("VerletEvalDerivs2");

  // Check if the timestep is still a good idea...
  if (dtcheck) {
    TIME_BEGIN("VerletDtCheck");
    const auto dtnew = this->selectDt(dtMin,
                                      dtMax,
                                      state,
                                      derivs);
    if (dtnew < dtcheckFrac*dt0) {
      this->currentTime(t);
      state.assign(*state0);
      TIME_END("VerletDtCheck");
      return false;
    }
    TIME_END("VerletDtCheck");
  }

  // Correct the final state by the end-point derivatives.
  TIME_BEGIN("VerletUpdateState");
  state.assign(state12);
  // state.timeAdvanceOnly(false);
  state.update(derivs, hdt0, t + hdt0, dt0);
  {
    auto pos = state.fields(HydroFieldNames::position, Vector::zero());
    pos.assignFields(pos0);
  }
  this->enforceBoundaries(state, derivs);
  this->applyGhostBoundaries(state, derivs);
  this->finalizeGhostBoundaries();
  this->postStateUpdate(t + dt0, dt0, db, state, derivs);
  TIME_END("VerletUpdateState");

  // Apply any physics specific finalizations.
  TIME_BEGIN("VerletFinalize");
  this->postStepFinalize(t + dt0, dt0, state, derivs);

  // Set the new current time and last time step.
  this->currentCycle(this->currentCycle() + 1);
  this->lastDt(dt0);
  TIME_END("VerletFinalize");
  TIME_END("Verlet");

  return true;
}

}
