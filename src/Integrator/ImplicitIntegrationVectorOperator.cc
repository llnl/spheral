//---------------------------------Spheral++----------------------------------//
// ImplicitIntegrationVectorOperator
//
// Implements the VectorOperator interface for solving for implicit integration
//
// Created by JMO, Tue Feb 25 11:15:29 PST 2025
//----------------------------------------------------------------------------//
#include "Integrator/ImplicitIntegrationVectorOperator.hh"

using std::vector;

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename Dimension>
ImplicitIntegrationVectorOperator<Dimension>::
ImplicitIntegrationVectorOperator(const Scalar t,
                                  const Scalar dt,
                                  const Scalar beta,
                                  const State<Dimension>& state0,
                                  const StateDerivatives<Dimension>& derivs0,
                                  State<Dimension>& state,
                                  StateDerivatives<Dimension>& derivs,
                                  Integrator<Dimension>& integrator):
  SolverFunction(),
  mt(t),
  mdt(dt),
  mBeta(beta),
  mInputState0(),
  mDerivs0(),
  mStateRef(state),
  mDerivsRef(derivs),
  mIntegratorRef(integrator) {
  state0.serializeIndependentData(mInputState0);
  state0.serializeDerivatives(mDerivs0, derivs0);
  this->numUnknowns(mInputState0.size());
}

//------------------------------------------------------------------------------
// Apply the operator (compute the residual (output) given the state (input)
//------------------------------------------------------------------------------
template<typename Dimension>
void
ImplicitIntegrationVectorOperator<Dimension>::
operator()(std::vector<double>& outputResiduals,
           const std::vector<double>& inputState) const {
           
  const auto n = mInputState0.size();
  REQUIRE(inputState.size() == n);
  
  // Grab references to objects we need to coordinate
  auto& integrator = mIntegratorRef.get();
  auto& state = mStateRef.get();
  auto& derivs = mDerivsRef.get();
  auto& db = integrator.dataBase();

  // Unpack the inputState to our State independent variables, and update
  // all the following dependent state as well
  state.deserializeIndependentData(inputState);
  state.update(derivs, mdt, mt, mdt, true);
  integrator.applyGhostBoundaries(state, derivs);
  integrator.finalizeGhostBoundaries();
  integrator.postStateUpdate(mt + mdt, mdt, db, state, derivs);

  // Evaluate the new derivatives
  integrator.initializeDerivatives(mt, mdt, state, derivs);
  derivs.Zero();
  integrator.evaluateDerivatives(mt, mdt, db, state, derivs);
  integrator.finalizeDerivatives(mt, mdt, db, state, derivs);
  
  // Compute the current residual estimate, and we're done
  state.serializeDerivatives(outputResiduals, derivs);
  CHECK(outputResiduals.size() == n);
  for (auto i = 0u; i < n; ++i) {
    outputResiduals[i] = inputState[i] - mInputState0[i] - mdt*((1.0 - mBeta)*mDerivs0[i] + mBeta*outputResiduals[i]);
  }

  // cerr << "Input state:";
  // for (const auto x: inputState) cerr << " " << x;
  // cerr << endl
  //      << "Residuals:";
  // for (const auto x: outputResiduals) cerr << " " << x;
  // cerr << endl;
}

}
