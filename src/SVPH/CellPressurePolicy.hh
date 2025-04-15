//---------------------------------Spheral++----------------------------------//
// CellPressurePolicy
//
// Created by JMO, Sun Aug 25 14:40:42 PDT 2013
//----------------------------------------------------------------------------//
#ifndef __Spheral_CellPressurePolicy_hh__
#define __Spheral_CellPressurePolicy_hh__

#include <string>

#include "DataBase/FieldUpdatePolicy.hh"

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class FluidNodeList;
template<typename Dimension, typename DataType> class Field;

template<typename Dimension>
class CellPressurePolicy: public FieldUpdatePolicy<Dimension, typename Dimension::Scalar> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using Scalar = typename Dimension::Scalar;
  using KeyType = typename FieldUpdatePolicy<Dimension, Scalar>::KeyType;

  // Constructors, destructor.
  CellPressurePolicy();
  virtual ~CellPressurePolicy() = default;
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt) override;

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const override;

  // Forbidden methods
  CellPressurePolicy(const CellPressurePolicy& rhs) = delete;
  CellPressurePolicy& operator=(const CellPressurePolicy& rhs) = delete;
};

}

#endif
