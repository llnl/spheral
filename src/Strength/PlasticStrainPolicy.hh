//---------------------------------Spheral++----------------------------------//
// PlasticStrainPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the dependent plastic strain state.
//
// Created by JMO, Wed Sep 15 23:07:21 PDT 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral_PlasticStrainPolicy_hh__
#define __Spheral_PlasticStrainPolicy_hh__

#include "DataBase/FieldUpdatePolicy.hh"

#include <string>

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class FluidNodeList;
template<typename Dimension, typename DataType> class Field;

template<typename Dimension>
class PlasticStrainPolicy: public FieldUpdatePolicy<Dimension, typename Dimension::Scalar> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using KeyType = typename FieldUpdatePolicy<Dimension, Scalar>::KeyType;

  // Constructors, destructor.
  PlasticStrainPolicy();
  virtual ~PlasticStrainPolicy() = default;
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt) override;

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const override;

  // Fobridden methjods
  PlasticStrainPolicy(const PlasticStrainPolicy& rhs) = delete;
  PlasticStrainPolicy& operator=(const PlasticStrainPolicy& rhs) = delete;
};

}

#endif
