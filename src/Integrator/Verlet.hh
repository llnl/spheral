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
#ifndef Verlet_HH
#define Verlet_HH

#include "Integrator.hh"

#include <vector>

namespace Spheral {

template<typename Dimension>
class Verlet: public Integrator<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  // Constructors.
  Verlet(DataBase<Dimension>& dataBase,
         const std::vector<Physics<Dimension>*>& physicsPackages);
  virtual ~Verlet() = default;
  Verlet& operator=(const Verlet& rhs) = default;

  // All Integrators are required to provide the single cycle method.
  virtual bool step(Scalar maxTime,
                    State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) override;

  // We need to make the simpler form of step visible!
  using Integrator<Dimension>::step;

  // Forbidden methods
  Verlet() = delete;

private:
  //--------------------------- Private Interface ---------------------------//
};

}

#endif
