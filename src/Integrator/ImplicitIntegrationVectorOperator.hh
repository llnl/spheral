//---------------------------------Spheral++----------------------------------//
// ImplicitIntegrationVectorOperator
//
// Implements the VectorOperator interface for solving for implicit integration
//
// Created by JMO, Tue Feb 25 11:15:29 PST 2025
//----------------------------------------------------------------------------//
#ifndef __Spheral_ImplicitIntegrationVectorOperator__
#define __Spheral_ImplicitIntegrationVectorOperator__

#include "Solvers/SolverFunction.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Integrator/Integrator.hh"

#include <vector>
#include <functional>

namespace Spheral {

template<typename Dimension>
class ImplicitIntegrationVectorOperator: public SolverFunction {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  // Constructor, destructor
  ImplicitIntegrationVectorOperator(const Scalar t,
                                    const Scalar dt,
                                    const Scalar beta,
                                    const State<Dimension>& state0,
                                    const StateDerivatives<Dimension>& derivs0,
                                    State<Dimension>& state,
                                    StateDerivatives<Dimension>& derivs,
                                    Integrator<Dimension>& integrator);
  ImplicitIntegrationVectorOperator(const ImplicitIntegrationVectorOperator& rhs) = default;
  ImplicitIntegrationVectorOperator& operator=(const ImplicitIntegrationVectorOperator& rhs) = default;
  virtual ~ImplicitIntegrationVectorOperator() = default;

  // Apply the operator (compute the residual (output) given the state (x)
  virtual void operator()(std::vector<double>& residuals,
                          const std::vector<double>& x) const override;

  // Handle internal data
  const std::vector<double>& inputState0()          const { return mInputState0; }
  const std::vector<double>& derivs0()              const { return mDerivs0; }
  State<Dimension>& state()                         const { return mStateRef.get(); }
  StateDerivatives<Dimension>& derivs()             const { return mDerivsRef.get(); }
  Integrator<Dimension>& integrator()               const { return mIntegratorRef.get(); }
  double t()                                        const { return mt; }
  double dt()                                       const { return mdt; }
  double beta()                                     const { return mBeta; }
  
  void state(State<Dimension>& x)                         { mStateRef = std::ref(x); }
  void derivs(StateDerivatives<Dimension>& x)             { mDerivsRef = std::ref(x); }
  void integrator(Integrator<Dimension>& x)               { mIntegratorRef = std::ref(x); }
  void t(const Scalar x)                                  { mt = x; }
  void dt(const Scalar x)                                 { mdt = x; }
  void beta(const Scalar x)                               { mBeta = x; }

  // Forbidden methods
  ImplicitIntegrationVectorOperator() = delete;

private:
  //--------------------------- Public Interface ---------------------------//
  double                                              mt, mdt, mBeta;
  std::vector<double>                                 mInputState0, mDerivs0;
  std::reference_wrapper<State<Dimension>>            mStateRef;
  std::reference_wrapper<StateDerivatives<Dimension>> mDerivsRef;
  std::reference_wrapper<Integrator<Dimension>>       mIntegratorRef;
};

}

#endif
