//---------------------------------Spheral++----------------------------------//
// YoungsModulusPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the dependent Youngs modulus.
//
// Created by JMO, Sun Jan 23 14:14:24 PST 2005
//----------------------------------------------------------------------------//
#ifndef __Spheral_YoungsModulusPolicy_hh__
#define __Spheral_YoungsModulusPolicy_hh__

#include "DataBase/FieldUpdatePolicy.hh"
#include "NodeList/SolidNodeList.hh"

#include <string>

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class FluidNodeList;
template<typename Dimension, typename DataType> class Field;

template<typename Dimension>
class YoungsModulusPolicy: public FieldUpdatePolicy<Dimension, typename Dimension::Scalar> {

public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using KeyType = typename UpdatePolicyBase<Dimension>::KeyType;

  // Constructors, destructor.
  YoungsModulusPolicy(const SolidNodeList<Dimension>& nodes);
  virtual ~YoungsModulusPolicy() = default;
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt);

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const;

  // Provide a static method to compute Youngs modulus
  static double YoungsModulus(const double K,     // bulk modulus
                              const double mu);   // shear modulus

  // Forbidden methods
  YoungsModulusPolicy(const YoungsModulusPolicy& rhs) = delete;
  YoungsModulusPolicy& operator=(const YoungsModulusPolicy& rhs) = delete;

private:
  //--------------------------- Private Interface ---------------------------//
  const SolidNodeList<Dimension>& mSolidNodeList;
};

}

#endif
