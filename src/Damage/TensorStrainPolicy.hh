//---------------------------------Spheral++----------------------------------//
// TensorStrainPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the tensor strain.
//
// Created by JMO, Mon Oct 17 10:56:28 PDT 2005
//----------------------------------------------------------------------------//
#ifndef __Spheral_TensorStrainPolicy_hh__
#define __Spheral_TensorStrainPolicy_hh__

#include "DataBase/FieldUpdatePolicy.hh"
#include "TensorDamageModel.hh"

#include <string>

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class FluidNodeList;
template<typename Dimension, typename DataType> class Field;

template<typename Dimension>
class TensorStrainPolicy: public FieldUpdatePolicy<Dimension, typename Dimension::SymTensor> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using KeyType = typename FieldUpdatePolicy<Dimension, SymTensor>::KeyType;

  // Constructors, destructor.
  TensorStrainPolicy(const TensorStrainAlgorithm strainType);
  virtual ~TensorStrainPolicy() = default;
  
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
  TensorStrainPolicy() = delete;
  TensorStrainPolicy(const TensorStrainPolicy& rhs) = delete;
  TensorStrainPolicy& operator=(const TensorStrainPolicy& rhs) = delete;

private:
  //--------------------------- Private Interface ---------------------------//
  TensorStrainAlgorithm mStrainType;
};

}

#endif
