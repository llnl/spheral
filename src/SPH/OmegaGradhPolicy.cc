//---------------------------------Spheral++----------------------------------//
// OmegaGradhPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the grad h correction terms.
//
// Created by JMO, Tue Oct 30 15:42:41 PDT 2007
//----------------------------------------------------------------------------//
#include "OmegaGradhPolicy.hh"
#include "Hydro/HydroFieldNames.hh"
#include "NodeList/NodeList.hh"
#include "NodeList/FluidNodeList.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Boundary/Boundary.hh"
#include "Geometry/GeometricUtilities.hh"
#include "Utilities/DBC.hh"
#include "Kernel/TableKernel.hh"

#include <algorithm>

namespace Spheral {


//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
OmegaGradhPolicy<Dimension>::
OmegaGradhPolicy(const DataBase<Dimension>& dataBase):
  UpdatePolicyBase<Dimension, Scalar>(HydroFieldNames::position,
                                      HydroFieldNames::H),
  mDataBase(dataBase) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
OmegaGradhPolicy<Dimension>::
~OmegaGradhPolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
OmegaGradhPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {
  REQUIRE(key.second == HydroFieldNames::omegaGradh);

  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Get the FluidNodeList.
  const FluidNodeList<Dimension>* nodeListPtr = dynamic_cast<const FluidNodeList<Dimension>*>(key.first);
  CHECK(nodeListPtr != 0);

  // Grab the state we need.
  const TableKernel<Dimension>& W = nodeListPtr->kernel();
  const ConnectivityMap<Dimension>& connectivityMap = mDataBase.connectivityMap();
  const vector<const NodeList<Dimension>*>& nodeLists = connectivityMap.nodeLists();
  const int numNodeLists = nodeLists.size();
  const int nodeListi = distance(nodeLists.begin(), find(nodeLists.begin(), nodeLists.end(), key.first));
  CHECK(nodeListi < numNodeLists);

  const FieldList<Dimension, Vector> position = state.vectorFields(HydroFieldNames::position);
  const FieldList<Dimension, SymTensor> H = state.symTensorFields(HydroFieldNames::H);
  const Field<Dimension, Vector>& positionThis = *position[nodeListi];
  const Field<Dimension, SymTensor>& HThis = *H[nodeListi];

  // Initialize the result field.
  Field<Dimension, Scalar>& omega = state.scalarField(key);

  // Iterate over the nodes in this NodeList.
  for (int i = 0; i != nodeListPtr->numInternalNodes(); ++i) {

    // State for node i.
    const vector< vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(nodeListPtr, i);
    const Vector& ri = positionThis(i);
    const SymTensor& Hi = HThis(i);
    const Scalar Hdeti = Hi.Determinant();
    Scalar& omegai = omega(i);
    Scalar ni = W(0.0, Hdeti);
    Scalar gradi = 0.0;

    // Iterate over all NodeLists.
    for (int nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {

      // Iterate over the neighbors in this NodeList.
      const vector<int>& connectivity = fullConnectivity[nodeListj];
      if (connectivity.size() > 0) {
        const Field<Dimension, Vector>& positionThem = *position[nodeListj];

        // Iterate over the neighbors.
        for (vector<int>::const_iterator jItr = connectivity.begin();
             jItr != connectivity.end();
             ++jItr) {

          // State of node j.
          const int j = *jItr;
          const Vector& rj = positionThem(j);

          // Add the contribution.
          const Vector rij = ri - rj;
          const Scalar etai = (Hi*rij).magnitude();
          ni += W(etai, Hdeti);
          gradi += etai * W.grad(etai, Hdeti);

        }
      }
    }

    // Finish off the correction for node i.
    CHECK(ni > 0.0);
    omegai = -gradi/(Dimension::nDim * ni);
  }

}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
OmegaGradhPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension, Scalar>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const OmegaGradhPolicy<Dimension>* rhsPtr = dynamic_cast<const OmegaGradhPolicy<Dimension>*>(&rhs);
  if (rhsPtr == 0) {
    return false;
  } else {
    return true;
  }
}

}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
namespace Spheral {
  template class OmegaGradhPolicy<Dim<1> >;
  template class OmegaGradhPolicy<Dim<2> >;
  template class OmegaGradhPolicy<Dim<3> >;
}
