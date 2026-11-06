//---------------------------------Spheral++----------------------------------//
// FieldListSecondDerivatives
// A set of experimental methods for evaluating the second derivative of a 
// FieldList.
//
// Created by JMO, Wed Dec 18 22:46:54 PST 2002
//----------------------------------------------------------------------------//

#include "FieldListSecondDerivatives.hh"
#include "FieldListFunctions.hh"
#include "Field/FieldList.hh"
#include "Field/Field.hh"
#include "NodeList/NodeList.hh"
#include "Kernel/TableKernel.hh"
#include "Boundary/Boundary.hh"

#include <vector>
using std::vector;

namespace Spheral {


//------------------------------------------------------------------------------
// Calculate the gradient of the divergence of a Vector FieldList.
// Explicit method that performs the div and grad operations as sequential first
// derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, typename Dimension::Vector>
gradDivVectorFieldList
(const FieldList<Dimension, typename Dimension::Vector>& fieldList,
 const FieldList<Dimension, typename Dimension::Vector>& position,
 const FieldList<Dimension, typename Dimension::Scalar>& weight,
 const FieldList<Dimension, typename Dimension::Scalar>& mass,
 const FieldList<Dimension, typename Dimension::Scalar>& density,
 const FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
 const TableKernel<Dimension>& kernel,
 const vector<Boundary<Dimension>*>&  boundaries) {

  // Some convenient typedefs.
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;

  // First evaluate the divergence of the input field list.
  FieldList<Dimension, Scalar> divField = divergence(fieldList,
                                                     position,
                                                     weight,
                                                     mass,
                                                     density,
                                                     Hfield,
                                                     kernel);

  // Apply boundary conditions to the divergence.
  for (typename vector<Boundary<Dimension>*>::const_iterator bcItr = boundaries.begin();
       bcItr < boundaries.end();
       ++bcItr) {
    (*bcItr)->applyFieldListGhostBoundary(divField);
  }
  for (typename vector<Boundary<Dimension>*>::const_iterator bcItr = boundaries.begin();
       bcItr < boundaries.end();
       ++bcItr) (*bcItr)->finalizeGhostBoundary();

  // Now take the gradient of this.
  FieldList<Dimension, Vector> result = gradient(divField,
                                                 position,
                                                 weight,
                                                 mass,
                                                 density,
                                                 Hfield,
                                                 kernel);

  // That's it.
  return result;
}

}
