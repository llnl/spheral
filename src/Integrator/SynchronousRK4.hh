//---------------------------------Spheral++----------------------------------//
// SynchronousRK4 -- Advance the set of Physics packages in time using fourth
// order Runge-Kutta.  All packages are advanced at one timestep simultaneously
// each step, i.e., synchronously.
//
// Created by JMO, Thu Jun 14 10:22:47 PDT 2012
//----------------------------------------------------------------------------//
#ifndef SynchronousRK4_HH
#define SynchronousRK4_HH

#include "Integrator.hh"

#include <vector>

namespace Spheral {

template<typename Dimension>
class SynchronousRK4: public Integrator<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  // Constructors.
  SynchronousRK4(DataBase<Dimension>& dataBase,
                 const std::vector<Physics<Dimension>*>& physicsPackages);
  virtual ~SynchronousRK4() = default;
  SynchronousRK4& operator=(const SynchronousRK4& rhs) = default;

  // All Integrators are required to provide the single cycle method.
  virtual bool step(Scalar maxTime,
                    State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) override;

  // We need to make the simpler form of step visible!
  using Integrator<Dimension>::step;

  // Fobidden methods
  SynchronousRK4() = delete;
};

}

#endif
