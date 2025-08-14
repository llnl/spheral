//---------------------------------Spheral++----------------------------------//
// ConstantRVelocityBoundary -- A boundary condition to enforce a constant 
// velocity on a given set of nodes.
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
#include "Field/FieldList.hh"
#include "Field/Field.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Utilities/DBC.hh"

#include "ConstantRVelocityBoundary.hh"

using std::vector;
using std::string;
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
ConstantRVelocityBoundary<Dimension>::
ConstantRVelocityBoundary(const NodeList<Dimension>& nodeList,
                          const vector<size_t>& nodeIndices):
  ConstantVelocityBoundary<Dimension>(nodeList, nodeIndices),
  mRadialVelocity() {
  mRadialVelocity.reserve(nodeIndices.size());
  const Field<Dimension, Vector>& positions = nodeList.positions();
  const Field<Dimension, Vector>& velocities = nodeList.velocity();
  for (auto i: nodeIndices) {
    mRadialVelocity.push_back(velocities(i).dot(positions(i).unitVector()));
  }
  ENSURE(mRadialVelocity.size() == nodeIndices.size());
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
ConstantRVelocityBoundary<Dimension>::~ConstantRVelocityBoundary() {
}

//------------------------------------------------------------------------------
// Apply the boundary condition to fields of different DataTypes.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConstantRVelocityBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::Vector>& field) const {

  REQUIRE(this->valid());

  // Is this field the velocity on the NodeList we're watching?
  if (field.nodeListPtr() == &(this->nodeList()) &&
      field.name() == HydroFieldNames::velocity) {

    // This is the velocity field, so enforce the boundary.
    size_t k = 0u;
    const auto nodeIDs = this->nodeIndices();
    const auto& positions = this->nodeList().positions();
    CHECK(nodeIDs.size() == mRadialVelocity.size());
    for (auto i: nodeIDs) {
      CHECK(i < field.numElements());
      CHECK(k < mRadialVelocity.size());
      const auto runit = positions(i).unitVector();
      const auto vperp = field[i] - field[i].dot(runit)*runit;
      field[i] = mRadialVelocity[k]*runit + vperp;
      ++k;
    }
    CHECK(k == nodeIDs.size() and k == mRadialVelocity.size());
  }
}

//------------------------------------------------------------------------------
// Dump the state.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConstantRVelocityBoundary<Dimension>::
dumpState(FileIO& file, const string& pathName) const {
  ConstantVelocityBoundary<Dimension>::dumpState(file, pathName);
  file.write(mRadialVelocity, pathName + "/radialVelocities");
}

//------------------------------------------------------------------------------
// Read the state.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConstantRVelocityBoundary<Dimension>::
restoreState(const FileIO& file, const string& pathName) {
  ConstantVelocityBoundary<Dimension>::restoreState(file, pathName);
  file.read(mRadialVelocity, pathName + "/radialVelocities");
}

}
