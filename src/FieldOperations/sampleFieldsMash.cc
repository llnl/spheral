//---------------------------------Spheral++----------------------------------//
// FieldListFunctionsMash -- A set of global functions which can be applied to
// FieldLists using MASH prescriptions.
//
// Created by JMO, Wed Dec  6 21:09:29 PST 2000
//----------------------------------------------------------------------------//
#include "FieldListFunctionsMash.hh"
#include "Field/FieldList.hh"
#include "Field/Field.hh"
#include "Field/NodeIterators.hh"
#include "NodeList/NodeList.hh"
#include "Neighbor/Neighbor.hh"
#include "Kernel/TableKernel.hh"
#include "Geometry/MathTraits.hh"

namespace Spheral {

using std::vector;

//------------------------------------------------------------------------------
// Return a MASH sampled version of the given FieldList at the new positions.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
FieldList<Dimension, DataType>
sampleFieldsMash(const FieldList<Dimension, DataType>& fieldList,
                 const FieldList<Dimension, typename Dimension::Vector>& position,
                 const FieldList<Dimension, typename Dimension::Scalar>& weight,
                 const FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
                 const TableKernel<Dimension>& kernel,
                 const FieldList<Dimension, typename Dimension::Vector>& samplePositions,
                 const FieldList<Dimension, typename Dimension::Scalar>& sampleWeight,
                 const FieldList<Dimension, typename Dimension::SymTensor>& sampleHfield) {

  // Some convenient typedefs.
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;

  // Pre-conditions.
  CHECK(fieldList.numFields() == position.numFields());
  CHECK(fieldList.numFields() == weight.numFields());
  CHECK(fieldList.numFields() == Hfield.numFields());
  for (auto i = 0u; i < fieldList.numFields(); ++i) {
    CHECK(fieldList[i]->nodeListPtr() == position[i]->nodeListPtr());
    CHECK(fieldList[i]->nodeListPtr() == weight[i]->nodeListPtr());
    CHECK(fieldList[i]->nodeListPtr() == Hfield[i]->nodeListPtr());
  }
  CHECK(samplePositions.numFields() == sampleWeight.numFields());
  CHECK(samplePositions.numFields() == sampleHfield.numFields());
  for (auto i = 0u; i < samplePositions.numFields(); ++i) {
    CHECK(samplePositions[i]->nodeListPtr() == sampleWeight[i]->nodeListPtr());
    CHECK(samplePositions[i]->nodeListPtr() == sampleHfield[i]->nodeListPtr());
  }

  // Return FieldList.
  FieldList<Dimension, DataType> result;
  FieldList<Dimension, Scalar> normalization;
  result.copyFields();
  normalization.copyFields();
  for (typename FieldList<Dimension, Vector>::const_iterator
         fieldItr = samplePositions.begin();
       fieldItr < samplePositions.end(); 
       ++fieldItr) {
    result.appendField(Field<Dimension, DataType>("sample mash", (*fieldItr)->nodeList()));
    normalization.appendField(Field<Dimension, Scalar>("normalization", (*fieldItr)->nodeList(), 1.0e-30));
  }
  vector< vector<bool> > flagNodeDone(fieldList.numFields());
  for (typename FieldList<Dimension, DataType>::const_iterator
         fieldItr = fieldList.begin();
       fieldItr < fieldList.end(); 
       ++fieldItr) {
    flagNodeDone[fieldItr - fieldList.begin()].resize((*fieldItr)->nodeListPtr()->numInternalNodes(), false);
  }

  // Loop over all the positions in the fieldList we're sampling.
  for (InternalNodeIterator<Dimension> nodeItr = fieldList.internalNodeBegin();
       nodeItr < fieldList.internalNodeEnd();
       ++nodeItr) {

    // Check if this node has been done yet.
    if (!flagNodeDone[nodeItr.fieldID()][nodeItr.nodeID()]) {

      // Set the neighbor info over the positions we're sampling to.
      vector<vector<int>> masterLists, coarseNeighbors, refineNeighbors,
                          masterListsSample, coarseNeighborsSample, refineNeighborsSample;
      fieldList.setMasterNodeLists(position(nodeItr), Hfield(nodeItr), masterLists, coarseNeighbors);
      samplePositions.setMasterNodeLists(position(nodeItr), Hfield(nodeItr), masterListsSample, coarseNeighborsSample);

      // Loop over the set of master nodes in the FieldList we're sampling from.
      for (MasterNodeIterator<Dimension> masterItr = fieldList.masterNodeBegin(masterLists);
           masterItr < fieldList.masterNodeEnd();
           ++masterItr) {
        CHECK(flagNodeDone[masterItr.fieldID()][masterItr.nodeID()] == false);
      
        // Sample node (i) state.
        const Vector& ri = position(masterItr);
        const SymTensor& Hi = Hfield(masterItr);
        const Scalar& weighti = weight(masterItr);
        const DataType& fieldi = fieldList(masterItr);

        // Refine the set of nodes we're sampling to for this position.
        samplePositions.setRefineNodeLists(ri, Hi, masterListsSample, coarseNeighborsSample);

        // Loop over the refined neighbors.
        for (RefineNodeIterator<Dimension> neighborItr = samplePositions.refineNodeBegin(refineNeighborsSample);
             neighborItr < samplePositions.refineNodeEnd();
             ++neighborItr) {

          // Node j's state.
          const Vector& rj = samplePositions(neighborItr);
          const SymTensor& Hj = sampleHfield(neighborItr);
          const Scalar& weightj = sampleWeight(neighborItr);

          const Vector rij = ri - rj;
          const Vector etai = Hi*rij;
          const Vector etaj = Hj*rij;

          // Calculate the kernel estimates for each node.
          Scalar Wi = kernel(etai.magnitude(), 1.0);
          Scalar Wj = kernel(etaj.magnitude(), 1.0);

          // Get the symmetrized kernel weighting for this node pair.
          Scalar Wij, weightij;
          switch(neighborItr.nodeListPtr()->neighbor().neighborSearchType()) {
          case NeighborSearchType::GatherScatter:
            Wij = 0.5*(Wi + Wj);
            weightij = 0.5*(weighti + weightj);
            break;

          case NeighborSearchType::Gather:
            Wij = Wi;
            weightij = weighti;
            break;

          case NeighborSearchType::Scatter:
            Wij = Wj;
            weightij = weightj;
            break;

          default:
            VERIFY2(false, "Unhandled neighbor search type.");
          }

          // Add this nodes contribution to the master value.
          normalization(neighborItr) += weightij*Wij;
          result(neighborItr) += fieldi*weightij*Wij;
        }
      
        // Flag this master node as done.
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
      cerr << "Error in FieldList::sampleFieldsMash: Not all values determined on exit "
           << checkcount << endl;
    }
    CHECK(checkcount == 0);
  }

  // Normalize the estimate.
  result /= normalization;
  return result;
}

}
