//---------------------------------Spheral++----------------------------------//
// IncrementASPHHtensor
//
// Specialized version of UpdatePolicy for time integrating the H tensor.
//
// Created by JMO, Mon Oct  7 13:31:02 PDT 2024
//----------------------------------------------------------------------------//
#ifndef __Spheral_IncrementASPHHtensor_hh__
#define __Spheral_IncrementASPHHtensor_hh__

#include "DataBase/UpdatePolicyBase.hh"
#include "SmoothingScale/ASPHRadialFunctor.hh"

#include <memory>  // std::shared_ptr

namespace Spheral {

// Forward declarations.
template<typename Dimension> class StateDerivatives;

template<typename Dimension>
class IncrementASPHHtensor: public UpdatePolicyBase<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using KeyType = typename UpdatePolicyBase<Dimension>::KeyType;
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using SymTensor = typename Dimension::SymTensor;
  using RadialFunctorType = ASPHRadialFunctor<Dimension>;

  // Constructors, destructor.
  IncrementASPHHtensor(const bool fixShape,
                       const bool radialOnly,
                       std::shared_ptr<RadialFunctorType> radialFunctorPtr);
  virtual ~IncrementASPHHtensor()   {}
  IncrementASPHHtensor(const IncrementASPHHtensor& rhs) = delete;
  IncrementASPHHtensor& operator=(const IncrementASPHHtensor& rhs) = delete;
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt) override;

  // Access the internal state
  bool fixShape() const             { return mFixShape; }
  bool radialOnly() const           { return mRadialOnly; }

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const override;

  static const std::string prefix() { return "delta "; }

private:
  //--------------------------- Private Interface ---------------------------//
  bool mFixShape, mRadialOnly;
  std::shared_ptr<RadialFunctorType> mRadialFunctorPtr;
};

}

#endif
