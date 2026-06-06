//---------------------------------Spheral++----------------------------------//
// ArtificialConductionPolicy -- Override the default energy policy in the
// presence of artificial conduction.
//
// Created by CDR, 9/30/2014
//----------------------------------------------------------------------------//

#ifndef __ArtificialConductionPolicy_hh__
#define __ArtificialConductionPolicy_hh__

#include <string>
#include "DataBase/FieldUpdatePolicy.hh"

namespace Spheral {
    
// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class FluidNodeList;
template<typename Dimension, typename DataType> class Field;

template<typename Dimension>
class ArtificialConductionPolicy: public FieldUpdatePolicy<Dimension, typename Dimension::Scalar> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using Scalar = typename Dimension::Scalar;
  using SymTensor = typename Dimension::SymTensor;
  using KeyType = typename UpdatePolicyBase<Dimension>::KeyType;
  using PolicyPointer = typename State<Dimension>::PolicyPointer;

  // Constructors, destructor.
  ArtificialConductionPolicy(PolicyPointer& energyPolicy);
  virtual ~ArtificialConductionPolicy() = default;

  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt) override;

  virtual void updateAsIncrement(const KeyType& key,
                                 State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivs,
                                 const double multiplier,
                                 const double t,
                                 const double dt) override;

  void conduct(const KeyType& key,
               State<Dimension>& state,
               StateDerivatives<Dimension>& derivs,
               const double multiplier,
               const double t,
               const double dt);

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const override;

  // Forbidden methods
  ArtificialConductionPolicy(const ArtificialConductionPolicy& rhs) = delete;
  ArtificialConductionPolicy& operator=(const ArtificialConductionPolicy& rhs) = delete;

private:
  //--------------------------- Private Interface ---------------------------//
  typename State<Dimension>::PolicyPointer mEnergyPolicy;
};

}

#endif
