//---------------------------------Spheral++----------------------------------//
// divergencePairWise -- Calculate a low order estimate of the gradient based on
// the pairwise formulation.
//
// Created by JMO, Sat Apr 12 15:58:53 PDT 2003
//----------------------------------------------------------------------------//

#include "PairWiseFieldListFunctions.hh"
#include "Field/FieldList.hh"
#include "Field/Field.hh"
#include "Field/NodeIterators.hh"
#include "NodeList/NodeList.hh"
#include "Neighbor/Neighbor.hh"
#include "Kernel/TableKernel.hh"
#include "Geometry/MathTraits.hh"
#include "incrementPairWiseDelta.hh"

namespace Spheral {

using std::vector;

//------------------------------------------------------------------------------
// Calculate the divergence of a FieldList.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
FieldList<Dimension, typename MathTraits<Dimension, DataType>::DivergenceType>
divergencePairWise
(const FieldList<Dimension, DataType>& fieldList,
 const FieldList<Dimension, typename Dimension::Vector>& position,
 const FieldList<Dimension, typename Dimension::Scalar>& weight,
 const FieldList<Dimension, typename Dimension::Scalar>& /*mass*/,
 const FieldList<Dimension, typename Dimension::Scalar>& /*rho*/,
 const FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
 const TableKernel<Dimension>& kernel) {

  // Typedef's to ease typing/understandability.
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename MathTraits<Dimension, DataType>::DivergenceType DivergenceType;

  // Return FieldList.
  FieldList<Dimension, DivergenceType> result;
  vector< vector<bool> > flagNodeDone(fieldList.numFields());
  result.copyFields();
  for (typename FieldList<Dimension, DataType>::const_iterator
         fieldItr = fieldList.begin();
       fieldItr < fieldList.end(); 
       ++fieldItr) {
    result.appendField(Field<Dimension, DivergenceType>("div pair wise", (*fieldItr)->nodeList()));
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
        //const Scalar& rhoi = rho(masterItr);
        //const Scalar& weighti = weight(masterItr);
        const DataType& fieldi = fieldList(masterItr);

        // Loop over the refined neighbors.
        Tensor gradi(0.0);
        Tensor normalization(0.0);
        for (RefineNodeIterator<Dimension> neighborItr = fieldList.refineNodeBegin(refineNeighbors);
             neighborItr < fieldList.refineNodeEnd();
             ++neighborItr) {
          if (neighborItr != masterItr) {

            const Vector& rj = position(neighborItr);
            const SymTensor& Hj = Hfield(neighborItr);
            //const Scalar& mj = mass(neighborItr);
            const Scalar& weightj = weight(neighborItr);
            const DataType& fieldj = fieldList(neighborItr);

            const Vector rij = ri - rj;
            const Vector etai = Hi*rij;
            const Vector etaj = Hj*rij;

            // Get the symmetrized kernel weighting for this node pair.
            Scalar Wij;
            switch((*fieldList.begin())->nodeListPtr()->neighbor().neighborSearchType()) {
            case NeighborSearchType::GatherScatter:
              Wij = 0.5*(kernel(etai.magnitude(), Hi.Determinant()) + 
                         kernel(etaj.magnitude(), Hj.Determinant()));
              break;

            case NeighborSearchType::Gather:
              Wij = kernel(etai.magnitude(), Hi.Determinant());
              break;

            case NeighborSearchType::Scatter:
              Wij = kernel(etaj.magnitude(), Hj.Determinant());
              break;

            default:
              VERIFY2(false, "Unhandled neighbor search type.");
            }

            // Add this nodes contribution to the master value.
            const Scalar wij = weightj*Wij;
            const DataType fij = fieldi - fieldj;
            incrementPairWiseDelta<Dimension>(gradi,
                                              normalization,
                                              wij,
                                              fij,
                                              rij,
                                              Hi);

//             // Add this nodes contribution to the master value.
//             const Vector fji = fieldj - fieldi;
//             const Scalar wij = weightj*Wij;
//             result(masterItr) += wij*fji.dot(rji)/(rji.magnitude2() + 1.0e-10);
//             normalization += wij;
          }
        }
        renormalizeGradient<Dimension>(gradi,
                                       normalization);

        // Now contract to the divergence.
        result(masterItr) = gradi.Trace();

//         CHECK(normalization > 0.0);
//         result(masterItr) /= normalization;

        // This master node is finished.
        flagNodeDone[masterItr.fieldID()][masterItr.nodeID()] = true;
      }
    }
  }

  // After we're done, all nodes in all NodeLists should be flagged as done.
  for (typename vector< vector<bool> >::const_iterator
         flagNodeItr = flagNodeDone.begin();
       flagNodeItr < flagNodeDone.end();
       ++flagNodeItr) {
    int checkcount = count(flagNodeItr->begin(), flagNodeItr->end(), false);
    if (checkcount > 0) {
      cerr << "Error in FieldList::divergence: Not all values determined on exit "
           << checkcount << endl;
    }
    CHECK(checkcount == 0);
  }

  return result;
}

}
