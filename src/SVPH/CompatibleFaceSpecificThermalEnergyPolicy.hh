//---------------------------------Spheral++----------------------------------//
// CompatibleFaceSpecificThermalEnergyPolicy
//
// Created by JMO, Fri Aug 16 14:47:57 PDT 2013
//----------------------------------------------------------------------------//
#ifndef __Spheral_CompatibleFaceSpecificThermalEnergyPolicy_hh__
#define __Spheral_CompatibleFaceSpecificThermalEnergyPolicy_hh__

#include <string>

#include "DataBase/IncrementState.hh"
#include "Physics/Physics.hh"
#include "Boundary/Boundary.hh"

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class FluidNodeList;
template<typename Dimension, typename DataType> class Field;
template<typename Dimension> class DataBase;
template<typename Dimension> class TableKernel;

template<typename Dimension>
class CompatibleFaceSpecificThermalEnergyPolicy: 
    public IncrementState<Dimension, typename Dimension::Scalar> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using KeyType = typename FieldUpdatePolicy<Dimension, Scalar>::KeyType;
  using ConstBoundaryIterator = typename Physics<Dimension>::ConstBoundaryIterator;

  // Constructors, destructor.
  CompatibleFaceSpecificThermalEnergyPolicy(const TableKernel<Dimension>& W,
                                            const DataBase<Dimension>& dataBase,
                                            ConstBoundaryIterator boundaryBegin,
                                            ConstBoundaryIterator boundaryEnd);
  virtual ~CompatibleFaceSpecificThermalEnergyPolicy() = default;
  
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
  CompatibleFaceSpecificThermalEnergyPolicy(const CompatibleFaceSpecificThermalEnergyPolicy& rhs) = delete;
  CompatibleFaceSpecificThermalEnergyPolicy& operator=(const CompatibleFaceSpecificThermalEnergyPolicy& rhs) = delete;

private:
  //--------------------------- Private Interface ---------------------------//
  const TableKernel<Dimension>& mW;
  const DataBase<Dimension>& mDataBase;
  ConstBoundaryIterator mBoundaryBegin, mBoundaryEnd;
  static bool mFired;
};

}

#endif
