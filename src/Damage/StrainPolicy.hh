//---------------------------------Spheral++----------------------------------//
// StrainPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the strain.
//
// Created by JMO, Sun Sep 26 16:15:17 PDT 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral_StrainPolicy_hh__
#define __Spheral_StrainPolicy_hh__

#include "DataBase/FieldUpdatePolicy.hh"

#include <string>

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class FluidNodeList;
template<typename Dimension, typename DataType> class Field;
template<typename Dimension> class StrainModel;

template<typename Dimension>
class StrainPolicy: public FieldUpdatePolicy<Dimension, typename Dimension::SymTensor> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using FieldType = Field<Dimension, Scalar>;
  using KeyType = typename FieldUpdatePolicy<Dimension, SymTensor>::KeyType;

  // Constructors, destructor.
  StrainPolicy();
  virtual ~StrainPolicy() = default;
  
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
  StrainPolicy(const StrainPolicy& rhs) = delete;
  StrainPolicy& operator=(const StrainPolicy& rhs) = delete;
};

}

#endif
