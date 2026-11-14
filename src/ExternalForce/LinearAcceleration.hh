//---------------------------------Spheral++----------------------------------//
// LinearAcceleration -- Impose a linear acceleration as a function of position
// on a given set of nodes.
//
// Created by JMO, Wed Sep 29 15:45:23 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral_PhysicsSpace_LinearAcceleration__
#define __Spheral_PhysicsSpace_LinearAcceleration__

#include "Physics/GenericBodyForce.hh"

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class NodeList;
template<typename Dimension> class DataBase;

template<typename Dimension>
class LinearAcceleration: public GenericBodyForce<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  using TimeStepType = typename Physics<Dimension>::TimeStepType;

  // Constructors.
  LinearAcceleration(const Scalar a0,
                     const Scalar aslope);

  // Destructor.
  virtual ~LinearAcceleration() = default;

  // This is the derivative method that all BodyPotential classes must provide.
  virtual 
  void evaluateDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivs) const override;

  // Provide the timestep appropriate for this package.
  virtual TimeStepType dt(const DataBase<Dimension>& dataBase, 
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs,
                          const Scalar currentTime) const override;

  //! Required label for Physics interface.
  virtual std::string label() const override { return "LinearAcceleration"; }

  // Access the acceleration parameters.
  Scalar a0()     const { return ma0; }
  Scalar aslope() const { return maslope; }

  // No default constructor, copying, or assignment.
  LinearAcceleration() = delete;
  LinearAcceleration(const LinearAcceleration& rhs) = delete;
  LinearAcceleration& operator=(const LinearAcceleration& rhs) = delete;

private:
  //--------------------------- Public Interface ---------------------------//
  Scalar ma0;
  Scalar maslope;
};

}

#endif
