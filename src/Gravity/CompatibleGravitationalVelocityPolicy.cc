//---------------------------------Spheral++----------------------------------//
// CompatibleGravitationalVelocityPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the gravitationally driven velocity in an exactly energy conservative manner.
// 
// Created by JMO, Mon Oct  2 13:44:59 PDT 2017
//----------------------------------------------------------------------------//
#include "CompatibleGravitationalVelocityPolicy.hh"
#include "Hydro/HydroFieldNames.hh"
#include "NodeList/NodeList.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/FieldUpdatePolicyBase.hh"
#include "DataBase/IncrementFieldList.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Utilities/DBC.hh"
#include "Utilities/SpheralFunctions.hh"

#include <vector>
#include <limits>
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
CompatibleGravitationalVelocityPolicy<Dimension>::
CompatibleGravitationalVelocityPolicy(const DataBase<Dimension>& dataBase,
                                      const Scalar G,
                                      const Scalar softeningLength):
  IncrementFieldList<Dimension, typename Dimension::Vector>(HydroFieldNames::position),
  mDataBasePtr(&dataBase),
  mG(G),
  mSofteningLength(softeningLength),
  mPositions0(dataBase.globalPosition()),
  mVelocity0(dataBase.globalVelocity()) {
  mPositions0.copyFields();
  mVelocity0.copyFields();
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
CompatibleGravitationalVelocityPolicy<Dimension>::
~CompatibleGravitationalVelocityPolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CompatibleGravitationalVelocityPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double /*t*/,
       const double /*dt*/) {

//   // HACK!
//   std::cerr.setf(std::ios::scientific, std::ios::floatfield);
//   std::cerr.precision(15);

  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == HydroFieldNames::velocity and 
          nodeListKey == UpdatePolicyBase<Dimension>::wildcard());
  auto vel = state.fields(fieldKey, Vector());
  const unsigned numFields = vel.numFields();

  // Get the state fields.
  const auto  mass = state.fields(HydroFieldNames::mass, Scalar());
  const auto  positions1 = state.fields(HydroFieldNames::position, Vector::zero());
  const auto  DvDt = derivs.fields(this->prefix() + HydroFieldNames::velocity, Vector::zero());

  // Prepare to accumulate velocity magnitude changes.
  auto dvel = mDataBasePtr->newGlobalFieldList(Scalar(0.0), "velocity mag change");

  // Walk all the NodeLists.
  for (size_t nodeListi = 0; nodeListi != numFields; ++nodeListi) {
    const auto ni = vel[nodeListi]->numInternalElements();

    // Iterate over the internal nodes of this NodeList.
    for (auto i = 0u; i < ni; ++i) {

      // State for node i.
      const auto  mi = mass(nodeListi, i);
      const auto& xi0 = mPositions0(nodeListi, i);
      const auto& xi1 = positions1(nodeListi, i);
      const auto& vi0 = mVelocity0(nodeListi, i);
      const auto  vi1 = vel(nodeListi, i) + multiplier*DvDt(nodeListi, i);
      auto&       dveli = dvel(nodeListi, i);

      // Iterate over the neighbor points.
      for (size_t nodeListj = 0; nodeListj != numFields; ++nodeListj) {
        const unsigned nj = vel[nodeListj]->numInternalElements();
        for (auto j = 0u; j < nj; ++j) {
          if (nodeListi != nodeListj or i != j) { // no self-interaction

            // State for node j.
            const auto  mj = mass(nodeListj, j);
            const auto& xj0 = mPositions0(nodeListj, j);
            const auto& xj1 = positions1(nodeListj, j);

            // Solve for the delta j->i that would give us the expected energy balance (gravitational->kinetic).
            const auto xji0 = xj0 - xi0;
            const auto xji1 = xj1 - xi1;
            const auto rij0 = sqrt(xji0.magnitude2() + mSofteningLength*mSofteningLength);
            const auto rij1 = sqrt(xji1.magnitude2() + mSofteningLength*mSofteningLength);
            const auto deltavij2 = 2.0*mG*mj/mi*(1.0/rij1 - 1.0/rij0);                       // Divide by dt^2 for aij2
            dveli += sqrt(std::abs(deltavij2))*sgn(deltavij2);
            cerr << "Diff : " << vi0.magnitude() << " " << dveli << " " << xji0 << " " << xji1 << " " << (xji1 - xji0) << endl;
          }
        }
      }

      // Now we can update the velocity.
      const auto vhati = vi1.unitVector();    // direction of time integrated velocity answer
      vel(nodeListi, i) = (vi0.magnitude() + dveli) * vhati;
    }
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
CompatibleGravitationalVelocityPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const CompatibleGravitationalVelocityPolicy<Dimension>* rhsPtr = dynamic_cast<const CompatibleGravitationalVelocityPolicy<Dimension>*>(&rhs);
  if (rhsPtr == 0) {
    return false;
  } else {
    return true;
  }
}

}

