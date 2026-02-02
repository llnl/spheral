//---------------------------------Spheral++----------------------------------//
// FieldListSecondDerivatives
// A set of experimental methods for evaluating the second derivative of a 
// FieldList.
//
// Created by JMO, Wed Dec 18 22:46:54 PST 2002
//----------------------------------------------------------------------------//

#include "FieldListSecondDerivatives.hh"
#include "Field/FieldList.hh"
#include "Field/Field.hh"
#include "Field/NodeIterators.hh"
#include "NodeList/NodeList.hh"
#include "Neighbor/Neighbor.hh"
#include "Kernel/TableKernel.hh"
#include "Boundary/Boundary.hh"
#include "Utilities/rotationMatrix.hh"

#include <vector>

namespace Spheral {

using std::vector;

//------------------------------------------------------------------------------
// Calculate the gradient of the divergence of a Vector FieldList.
// Simplest method, just use the second derivative of the kernel.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, typename Dimension::Vector>
gradDivVectorFieldListSimple
(const FieldList<Dimension, typename Dimension::Vector>& fieldList,
 const FieldList<Dimension, typename Dimension::Vector>& position,
 const FieldList<Dimension, typename Dimension::Scalar>& weight,
 const FieldList<Dimension, typename Dimension::Scalar>& /*mass*/,
 const FieldList<Dimension, typename Dimension::Scalar>& /*rho*/,
 const FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
 const TableKernel<Dimension>& kernel) {

  // Some convenient typedefs.
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Return FieldList.
  FieldList<Dimension, Vector> result;
  vector< vector<bool> > flagNodeDone(fieldList.numFields());
  result.copyFields();
  for (typename FieldList<Dimension, Vector>::const_iterator fieldItr = fieldList.begin();
       fieldItr < fieldList.end(); 
       ++fieldItr) {
    result.appendField(Field<Dimension, Vector>("grad div " + (*fieldItr)->name(), (*fieldItr)->nodeList()));
    flagNodeDone[fieldItr - fieldList.begin()].resize((*fieldItr)->nodeListPtr()->numInternalNodes(), false);
  }

  // Loop over all the elements in the input FieldList.
  for (InternalNodeIterator<Dimension> nodeItr = fieldList.internalNodeBegin();
       nodeItr < fieldList.internalNodeEnd();
       ++nodeItr) {

    // Check if this node has been done yet.
    if (!flagNodeDone[nodeItr.fieldID()][nodeItr.nodeID()]) {

      // We will do the batch of master nodes associated with this node together.
      // Set the neighbor information.
      vector<vector<int>> masterLists, coarseNeighbors, refineNeighbors;
      fieldList.setMasterNodeLists(position(nodeItr), Hfield(nodeItr), masterLists, coarseNeighbors);

      // Now loop over all the master nodes.
      for (MasterNodeIterator<Dimension> masterItr = fieldList.masterNodeBegin(masterLists);
           masterItr < fieldList.masterNodeEnd();
           ++masterItr) {
        CHECK(flagNodeDone[masterItr.fieldID()][masterItr.nodeID()] == false);

        // Set the refined neighbor information for this master node.
        fieldList.setRefineNodeLists(position(masterItr), Hfield(masterItr), coarseNeighbors, refineNeighbors);

        // State for node i.
        const Vector& ri = position(masterItr);
        const SymTensor& Hi = Hfield(masterItr);
        const Scalar Hdeti = Hi.Determinant();
        const Scalar& weighti = weight(masterItr);
        //const Scalar& mi = mass(masterItr);
        //const Scalar& rhoi = rho(masterItr);
        const Vector& fieldi = fieldList(masterItr);

        // Loop over the refined neighbors.
        Scalar Wnormalization = weighti*kernel(0.0, Hdeti);
        vector<Tensor> D2fDx2i(Dimension::nDim);
        for (RefineNodeIterator<Dimension> neighborItr = fieldList.refineNodeBegin(refineNeighbors);
             neighborItr < fieldList.refineNodeEnd();
             ++neighborItr) {
          if (neighborItr != masterItr) {

            const Vector& rj = position(neighborItr);
            const SymTensor& Hj = Hfield(neighborItr);
            const Scalar Hdetj = Hj.Determinant();
            const Scalar& weightj = weight(neighborItr);
            //const Scalar& mj = mass(neighborItr);
            //const Scalar& rhoj = rho(neighborItr);
            const Vector& fieldj = fieldList(neighborItr);

            const Vector rji = rj - ri;
            const Vector etai = Hi*rji;
            const Vector etaj = Hj*rji;
            const Scalar etaMagi = etai.magnitude();
            const Scalar etaMagj = etaj.magnitude();

            const Vector etaiNorm = etai.unitVector();
            const Vector etajNorm = etaj.unitVector();
            const Vector Hetai = Hi*etaiNorm;
            const Vector Hetaj = Hj*etajNorm;

            const Scalar Wi = kernel(etaMagi, Hdeti);
            const Scalar Wj = kernel(etaMagj, Hdetj);

            const Scalar getai = kernel.grad(etaMagi, Hdeti);
            const Scalar getaj = kernel.grad(etaMagj, Hdetj);

            const Vector gWi = Hetai*getai;
            const Vector gWj = Hetaj*getaj;

            // Get the symmetrized kernel weighting for this node pair.
            Scalar Wij;
            Vector gWij;
            switch((*fieldList.begin())->nodeListPtr()->neighbor().neighborSearchType()) {
            case NeighborSearchType::GatherScatter:
              Wij = 0.5*(Wi + Wj);
              gWij = 0.5*(gWi + gWj);
              break;

            case NeighborSearchType::Gather:
              Wij = Wi;
              gWij = gWi;
              break;

            case NeighborSearchType::Scatter:
              Wij = Wj;
              gWij = gWj;
              break;

            default:
              VERIFY2(false, "Unhandled neighbor search type.");
            }

            // Sum this pairs contribution to the elements.
            Wnormalization += weightj*Wij;
            const Vector fji = fieldj - fieldi;
            const Scalar frach = 0.01/Hi.Trace()/Dimension::nDim;
            const Vector rjiUnit = rji.unitVector();
            const Tensor R = rotationMatrix(rjiUnit);
            const Tensor Rinverse = R.Transpose();
            const Scalar dxp = rji.magnitude()/(rji.magnitude2() + frach*frach);

            // Find the local contribution to dfdx.
            const Vector dfdxcol = (R*fji)*dxp;

            // Take the second derivative in the prime frame.
            const Tensor d2fdx2 = dfdxcol*(R*gWij);

            // Now increment the total result in the prime frame, and then rotate
            // it back.
            result(masterItr) = R*result(masterItr);
            result(masterItr) -= weightj*(d2fdx2.diagonalElements());
            result(masterItr) = Rinverse*result(masterItr);
          }
        }

        // Put together the elements into the final grad div answer.
        CHECK(Wnormalization > 0.0);
        result(masterItr) *= 2.0/Wnormalization;

        // This master node is finished.
        flagNodeDone[masterItr.fieldID()][masterItr.nodeID()] = true;

      }
    }
  }

  // After we're done, all nodes in all NodeLists should be flagged as done.
  for (typename vector< vector<bool> >::const_iterator flagNodeItr = flagNodeDone.begin();
       flagNodeItr < flagNodeDone.end();
       ++flagNodeItr) {
    int checkcount = count(flagNodeItr->begin(), flagNodeItr->end(), false);
    if (checkcount > 0) {
      cerr << "Error in FieldList::smoothFields: Not all values determined on exit "
           << checkcount << endl;
    }
    CHECK(checkcount == 0);
  }

  return result;
}

}
