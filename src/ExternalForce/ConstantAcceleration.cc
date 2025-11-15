//---------------------------------Spheral++----------------------------------//
// ConstantAcceleration -- Impose a constant acceleration on a given set of
// nodes.
//
// Created by JMO, Mon Sep 27 23:01:15 PDT 2004
//----------------------------------------------------------------------------//
#include "ConstantAcceleration.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/IncrementState.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "Utilities/DBC.hh"

using std::vector;
using std::string;
using std::pair;
using std::make_pair;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
ConstantAcceleration<Dimension>::
ConstantAcceleration(const Vector a0,
                     const NodeList<Dimension>& nodeList,
                     const vector<int>& indices):
  GenericBodyForce<Dimension>(),
  ma0(a0),
  mNodeListPtr(&nodeList),
  mFlagsPtr(std::make_shared<Field<Dimension, int>>("constant acceleration flags for " + nodeList.name(), nodeList, 0)) {
  for (const auto i: indices) {
    (*mFlagsPtr)(i) = 1;
  }
}

//------------------------------------------------------------------------------
// Constructor for all nodes in NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
ConstantAcceleration<Dimension>::
ConstantAcceleration(const Vector a0,
                     const NodeList<Dimension>& nodeList):
  GenericBodyForce<Dimension>(),
  ma0(a0),
  mNodeListPtr(&nodeList),
  mFlagsPtr() {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
ConstantAcceleration<Dimension>::
~ConstantAcceleration() {
}

//------------------------------------------------------------------------------
// Apply the acceleration to the derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConstantAcceleration<Dimension>::
evaluateDerivatives(const typename Dimension::Scalar /*time*/,
                    const typename Dimension::Scalar /*dt*/,
                    const DataBase<Dimension>& /*dataBase*/,
                    const State<Dimension>& /*state*/,
                    StateDerivatives<Dimension>& derivs) const {

  // Get the acceleration state we'll be modifying.
  const typename State<Dimension>::KeyType key = State<Dimension>::buildFieldKey(IncrementState<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::velocity,
                                                                                 mNodeListPtr->name());
  Field<Dimension, Vector>& DvDt = derivs.field(key, Vector::zero());

  // Increment the acceleration.
  const unsigned n = mNodeListPtr->numNodes();
  for (unsigned i = 0; i != n; ++i) {
    if (not mFlagsPtr or ((*mFlagsPtr)(i) == 1)) {
      DvDt(i) += ma0;
    }
  }
}

//------------------------------------------------------------------------------
// Calculate the timestep constraint.
//------------------------------------------------------------------------------
template<typename Dimension>
typename ConstantAcceleration<Dimension>::TimeStepType
ConstantAcceleration<Dimension>::
dt(const DataBase<Dimension>& /*dataBase*/, 
   const State<Dimension>& /*state*/,
   const StateDerivatives<Dimension>& /*derivs*/,
   const typename Dimension::Scalar /*currentTime*/) const {

  return TimeStepType(FLT_MAX, "No vote.");
}

}
