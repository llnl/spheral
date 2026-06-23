//---------------------------------Spheral++----------------------------------//
// AxisymmetricMassDensityPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the 3D mass density in RZ calculations.
//
// Created by JMO, Thu May  7 15:20:58 PDT 2026
//----------------------------------------------------------------------------//
#ifndef __Spheral_AxisymmetricMassDensityPolicy_hh__
#define __Spheral_AxisymmetricMassDensityPolicy_hh__

#include "DataBase/FieldUpdatePolicy.hh"
#include "Geometry/Dimension.hh"

#include <string>

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class FluidNodeList;
template<typename Dimension, typename DataType> class FieldList;

class AxisymmetricMassDensityPolicy: public FieldUpdatePolicy<Dim<2>, Dim<2>::Scalar> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using Dimension = Dim<2>;
  using Scalar = Dimension::Scalar;
  using Vector = Dimension::Vector;
  using SymTensor = Dimension::SymTensor;
  using KeyType = UpdatePolicyBase<Dimension>::KeyType;

  // Constructors, destructor.
  AxisymmetricMassDensityPolicy(const Scalar rhoMin,
                                const Scalar rhoMax);
  virtual ~AxisymmetricMassDensityPolicy() = default;
  
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
  AxisymmetricMassDensityPolicy() = delete;
  AxisymmetricMassDensityPolicy(const AxisymmetricMassDensityPolicy& rhs) = delete;
  AxisymmetricMassDensityPolicy& operator=(const AxisymmetricMassDensityPolicy& rhs) = delete;

private:
  Scalar mRhoMin, mRhoMax;
};

}

#endif
