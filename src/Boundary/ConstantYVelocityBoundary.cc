//---------------------------------Spheral++----------------------------------//
// ConstantYVelocityBoundary -- A boundary condition to enforce a constant 
// velocity on a given set of nodes.
//----------------------------------------------------------------------------//
#include "ConstantYVelocityBoundary.hh"
#include "Field/FieldList.hh"
#include "Field/Field.hh"
#include "Hydro/HydroFieldNames.hh"

#include "Utilities/DBC.hh"

using std::vector;

namespace Spheral {

//------------------------------------------------------------------------------
// Construct with the given set of nodes.
//------------------------------------------------------------------------------
template<typename Dimension>
ConstantYVelocityBoundary<Dimension>::
ConstantYVelocityBoundary(const NodeList<Dimension>& nodeList,
                          const vector<size_t>& nodeIndices):
  ConstantVelocityBoundary<Dimension>(nodeList, nodeIndices) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
ConstantYVelocityBoundary<Dimension>::~ConstantYVelocityBoundary() {
}

//------------------------------------------------------------------------------
// Apply the boundary condition to fields of different DataTypes.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConstantYVelocityBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::Vector>& field) const {

  REQUIRE(this->valid());

  // Is this field the velocity on the NodeList we're watching?
  if (field.nodeListPtr() == &(this->nodeList()) &&
      field.name() == HydroFieldNames::velocity) {

    // This is the velocity field, so enforce the boundary.
    size_t k = 0u;
    const auto nodeIDs = this->nodeIndices();
    for (auto i: nodeIDs) {
      CHECK(i < field.numElements());
      CHECK(k < this->velocityCondition().size());
      field[i].y(this->velocityCondition()[k].y());
      ++k;
    }
    CHECK(k == nodeIDs.size());
  }
}

}
