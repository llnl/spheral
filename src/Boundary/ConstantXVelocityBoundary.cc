//---------------------------------Spheral++----------------------------------//
// ConstantXVelocityBoundary -- A boundary condition to enforce a constant 
// velocity on a given set of nodes.
//----------------------------------------------------------------------------//
#include "ConstantXVelocityBoundary.hh"
#include "Field/FieldList.hh"
#include "Field/Field.hh"
#include "Hydro/HydroFieldNames.hh"

#include "Utilities/DBC.hh"

using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {

//------------------------------------------------------------------------------
// Construct with the given set of nodes.
//------------------------------------------------------------------------------
template<typename Dimension>
ConstantXVelocityBoundary<Dimension>::
ConstantXVelocityBoundary(const NodeList<Dimension>& nodeList,
                          const vector<size_t>& nodeIndices):
  ConstantVelocityBoundary<Dimension>(nodeList, nodeIndices) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
ConstantXVelocityBoundary<Dimension>::~ConstantXVelocityBoundary() {
}

//------------------------------------------------------------------------------
// Apply the boundary condition to fields of different DataTypes.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConstantXVelocityBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::Vector>& field) const {

  REQUIRE(this->valid());

  // Is this field the velocity on the NodeList we're watching?
  if (field.nodeListPtr() == &(this->nodeList()) &&
      field.name() == HydroFieldNames::velocity) {

    // This is the velocity field, so enforce the boundary.
    size_t k = 0;
    const auto nodeIDs = this->nodeIndices();
    for (auto i: nodeIDs) {
      CHECK(i < field.numElements());
      CHECK(k < this->velocityCondition().size());
      field[i].x(this->velocityCondition()[k].x());
      ++k;
    }
    CHECK(k == nodeIDs.size());
  }
}

}
