//---------------------------------Spheral++----------------------------------//
// LongitudinalSoundSpeedPolicy -- An implementation of UpdatePolicyBase
// specialized for the updating the dependent longitudinal sound speed.
//
// Created by JMO, Sun Jan 23 14:14:24 PST 2005
//----------------------------------------------------------------------------//
#ifndef __Spheral_LongitudinalSoundSpeedPolicy_hh__
#define __Spheral_LongitudinalSoundSpeedPolicy_hh__

#include <string>

#include "DataBase/FieldUpdatePolicy.hh"
#include "NodeList/SolidNodeList.hh"

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class FluidNodeList;
template<typename Dimension, typename DataType> class Field;

template<typename Dimension>
class LongitudinalSoundSpeedPolicy: public FieldUpdatePolicy<Dimension, typename Dimension::Scalar> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using KeyType = typename UpdatePolicyBase<Dimension>::KeyType;

  // Constructors, destructor.
  LongitudinalSoundSpeedPolicy(const SolidNodeList<Dimension>& nodes);
  virtual ~LongitudinalSoundSpeedPolicy() = default;
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt);

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const;

  // Forbidden methods
  LongitudinalSoundSpeedPolicy(const LongitudinalSoundSpeedPolicy& rhs) = delete;
  LongitudinalSoundSpeedPolicy& operator=(const LongitudinalSoundSpeedPolicy& rhs) = delete;

private:
  //--------------------------- Private Interface ---------------------------//
  const SolidNodeList<Dimension>& mSolidNodeList;
};

}

#endif
