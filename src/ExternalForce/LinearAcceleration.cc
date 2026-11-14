//---------------------------------Spheral++----------------------------------//
// LinearAcceleration -- Impose a linear acceleration as a function of position
// on a given set of nodes.
//
// Created by JMO, Wed Sep 29 15:45:23 2004
//----------------------------------------------------------------------------//
#include "LinearAcceleration.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/IncrementState.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Field/FieldList.hh"
#include "Utilities/DBC.hh"

#include <limits>

namespace Spheral {


//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
LinearAcceleration<Dimension>::
LinearAcceleration(const Scalar a0,
                   const Scalar aslope):
  GenericBodyForce<Dimension>(),
  ma0(a0),
  maslope(aslope) {
}

//------------------------------------------------------------------------------
// Apply the acceleration to the derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
LinearAcceleration<Dimension>::
evaluateDerivatives(const typename Dimension::Scalar /*time*/,
                    const typename Dimension::Scalar /*dt*/,
                    const DataBase<Dimension>& /*dataBase*/,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) const {

  // Get the state.
  const auto r = state.fields(HydroFieldNames::position, Vector::zero());
  auto       DvDt = derivs.fields(IncrementState<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::velocity, Vector::zero());
  const auto numNodeLists = r.numFields();
  CHECK(DvDt.numFields() == numNodeLists);

  // Increment the acceleration.
  for (auto k = 0u; k < numNodeLists; ++k) {
    const auto n = r[k]->numInternalElements();
    for (auto i = 0u; i < n; ++i) {
      const auto x = r(k, i).x();
      DvDt(k, i)[0] += ma0 + maslope*x;
    }
  }
}

//------------------------------------------------------------------------------
// Calculate the timestep constraint.
//------------------------------------------------------------------------------
template<typename Dimension>
typename LinearAcceleration<Dimension>::TimeStepType
LinearAcceleration<Dimension>::
dt(const DataBase<Dimension>& /*dataBase*/, 
   const State<Dimension>& /*state*/,
   const StateDerivatives<Dimension>& /*derivs*/,
   const typename Dimension::Scalar /*currentTime*/) const {
  return TimeStepType(std::numeric_limits<Scalar>::max(), "No vote.");
}

}
