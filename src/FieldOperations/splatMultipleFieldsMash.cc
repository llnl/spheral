//---------------------------------Spheral++----------------------------------//
// FieldListFunctionsMash -- A set of global functions which can be applied to
// FieldLists using MASH prescriptions.
//
// Created by JMO, Tue Jun  4 17:07:31 PDT 2002
//----------------------------------------------------------------------------//
#include "FieldListFunctionsMash.hh"
#include "Field/FieldList.hh"
#include "Field/FieldListSet.hh"
#include "Field/Field.hh"
#include "Field/NodeIterators.hh"
#include "NodeList/NodeList.hh"
#include "NodeList/FluidNodeList.hh"
#include "Neighbor/Neighbor.hh"
#include "Kernel/TableKernel.hh"
#include "Geometry/MathTraits.hh"
#include "Boundary/Boundary.hh"

#include "Utilities/DBC.hh"

namespace Spheral {

using std::vector;

//------------------------------------------------------------------------------
// Return a MASH donated version of the given FieldList at the new positions.
// This version uses the "splat" operation, guaranteeing that this operation is
// conservative in the sense that summing over the original Field values and
// the new Field values yeilds the same total.
// This version operates on a set of FieldLists at the same time, purely as an 
// optimzation over the single FieldList version, splatFieldsMash.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldListSet<Dimension>
splatMultipleFieldsMash(const FieldListSet<Dimension>& fieldListSet,
                        const FieldList<Dimension, typename Dimension::Vector>& position,
                        const FieldList<Dimension, typename Dimension::Scalar>& weight,
                        const FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
                        const TableKernel<Dimension>& kernel,
                        const FieldList<Dimension, typename Dimension::Vector>& samplePositions,
                        const FieldList<Dimension, typename Dimension::Scalar>& sampleWeight,
                        const FieldList<Dimension, typename Dimension::SymTensor>& sampleHfield,
                        const vector<Boundary<Dimension>*>& boundaries) {

  // Some convenient typedefs.
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Pre-conditions.
  for (typename vector< FieldList<Dimension, Scalar> >::const_iterator fieldListItr = 
         fieldListSet.ScalarFieldLists.begin();
       fieldListItr < fieldListSet.ScalarFieldLists.end();
       ++fieldListItr) {
    const FieldList<Dimension, Scalar>& fieldList = *fieldListItr;
    REQUIRE(fieldList.numFields() == position.numFields());
    REQUIRE(fieldList.numFields() == weight.numFields());
    REQUIRE(fieldList.numFields() == Hfield.numFields());
    for (auto i = 0u; i < fieldList.numFields(); ++i) {
      REQUIRE(fieldList[i]->nodeListPtr() == position[i]->nodeListPtr());
      REQUIRE(fieldList[i]->nodeListPtr() == weight[i]->nodeListPtr());
      REQUIRE(fieldList[i]->nodeListPtr() == Hfield[i]->nodeListPtr());
    }
  }
  for (typename vector< FieldList<Dimension, Vector> >::const_iterator fieldListItr = 
         fieldListSet.VectorFieldLists.begin();
       fieldListItr < fieldListSet.VectorFieldLists.end();
       ++fieldListItr) {
    const FieldList<Dimension, Vector>& fieldList = *fieldListItr;
    REQUIRE(fieldList.numFields() == position.numFields());
    REQUIRE(fieldList.numFields() == weight.numFields());
    REQUIRE(fieldList.numFields() == Hfield.numFields());
    for (auto i = 0u; i < fieldList.numFields(); ++i) {
      REQUIRE(fieldList[i]->nodeListPtr() == position[i]->nodeListPtr());
      REQUIRE(fieldList[i]->nodeListPtr() == weight[i]->nodeListPtr());
      REQUIRE(fieldList[i]->nodeListPtr() == Hfield[i]->nodeListPtr());
    }
  }
  for (typename vector< FieldList<Dimension, Tensor> >::const_iterator fieldListItr = 
         fieldListSet.TensorFieldLists.begin();
       fieldListItr < fieldListSet.TensorFieldLists.end();
       ++fieldListItr) {
    const FieldList<Dimension, Tensor>& fieldList = *fieldListItr;
    REQUIRE(fieldList.numFields() == position.numFields());
    REQUIRE(fieldList.numFields() == weight.numFields());
    REQUIRE(fieldList.numFields() == Hfield.numFields());
    for (auto i = 0u; i < fieldList.numFields(); ++i) {
      REQUIRE(fieldList[i]->nodeListPtr() == position[i]->nodeListPtr());
      REQUIRE(fieldList[i]->nodeListPtr() == weight[i]->nodeListPtr());
      REQUIRE(fieldList[i]->nodeListPtr() == Hfield[i]->nodeListPtr());
    }
  }
  for (typename vector< FieldList<Dimension, SymTensor> >::const_iterator fieldListItr = 
         fieldListSet.SymTensorFieldLists.begin();
       fieldListItr < fieldListSet.SymTensorFieldLists.end();
       ++fieldListItr) {
    const FieldList<Dimension, SymTensor>& fieldList = *fieldListItr;
    REQUIRE(fieldList.numFields() == position.numFields());
    REQUIRE(fieldList.numFields() == weight.numFields());
    REQUIRE(fieldList.numFields() == Hfield.numFields());
    for (auto i = 0u; i < fieldList.numFields(); ++i) {
      REQUIRE(fieldList[i]->nodeListPtr() == position[i]->nodeListPtr());
      REQUIRE(fieldList[i]->nodeListPtr() == weight[i]->nodeListPtr());
      REQUIRE(fieldList[i]->nodeListPtr() == Hfield[i]->nodeListPtr());
    }
  }
  REQUIRE(samplePositions.numFields() == sampleWeight.numFields());
  REQUIRE(samplePositions.numFields() == sampleHfield.numFields());
  for (auto i = 0u; i < samplePositions.numFields(); ++i) {
    REQUIRE(samplePositions[i]->nodeListPtr() == sampleWeight[i]->nodeListPtr());
    REQUIRE(samplePositions[i]->nodeListPtr() == sampleHfield[i]->nodeListPtr());
  }

  // Return FieldList.
  FieldListSet<Dimension> resultSet;
  resultSet.ScalarFieldLists.resize(fieldListSet.ScalarFieldLists.size());
  resultSet.VectorFieldLists.resize(fieldListSet.VectorFieldLists.size());
  resultSet.TensorFieldLists.resize(fieldListSet.TensorFieldLists.size());
  resultSet.SymTensorFieldLists.resize(fieldListSet.SymTensorFieldLists.size());
  for (typename vector<FieldList<Dimension, Scalar> >::iterator resultItr = 
         resultSet.ScalarFieldLists.begin();
       resultItr < resultSet.ScalarFieldLists.end();
       ++resultItr) resultItr->copyFields();
  for (typename vector<FieldList<Dimension, Vector> >::iterator resultItr = 
         resultSet.VectorFieldLists.begin();
       resultItr < resultSet.VectorFieldLists.end();
       ++resultItr) resultItr->copyFields();
  for (typename vector<FieldList<Dimension, Tensor> >::iterator resultItr = 
         resultSet.TensorFieldLists.begin();
       resultItr < resultSet.TensorFieldLists.end();
       ++resultItr) resultItr->copyFields();
  for (typename vector<FieldList<Dimension, SymTensor> >::iterator resultItr = 
         resultSet.SymTensorFieldLists.begin();
       resultItr < resultSet.SymTensorFieldLists.end();
       ++resultItr) resultItr->copyFields();
  for (typename FieldList<Dimension, Vector>::const_iterator fieldItr = samplePositions.begin();
       fieldItr < samplePositions.end(); 
       ++fieldItr) {
    for (typename vector<FieldList<Dimension, Scalar> >::iterator resultItr = 
           resultSet.ScalarFieldLists.begin();
         resultItr < resultSet.ScalarFieldLists.end();
         ++resultItr) 
      resultItr->appendField(Field<Dimension, Scalar>("splat" + (*fieldItr)->name(), (*fieldItr)->nodeList()));
    for (typename vector<FieldList<Dimension, Vector> >::iterator resultItr = 
           resultSet.VectorFieldLists.begin();
         resultItr < resultSet.VectorFieldLists.end();
         ++resultItr) 
      resultItr->appendField(Field<Dimension, Vector>("splat" + (*fieldItr)->name(), (*fieldItr)->nodeList()));
    for (typename vector<FieldList<Dimension, Tensor> >::iterator resultItr = 
           resultSet.TensorFieldLists.begin();
         resultItr < resultSet.TensorFieldLists.end();
         ++resultItr) 
      resultItr->appendField(Field<Dimension, Tensor>("splat" + (*fieldItr)->name(), (*fieldItr)->nodeList()));
    for (typename vector<FieldList<Dimension, SymTensor> >::iterator resultItr = 
           resultSet.SymTensorFieldLists.begin();
         resultItr < resultSet.SymTensorFieldLists.end();
         ++resultItr) 
      resultItr->appendField(Field<Dimension, SymTensor>("splat" + (*fieldItr)->name(), (*fieldItr)->nodeList()));
  }
  
  FieldList<Dimension, int> flagNodeDone(FieldStorageType::CopyFields);
  FieldList<Dimension, Scalar> normalization(FieldStorageType::CopyFields);
  for (typename FieldList<Dimension, Vector>::const_iterator fieldItr = position.begin();
       fieldItr < position.end(); 
       ++fieldItr) {
    flagNodeDone.appendNewField("flag nodes", (*fieldItr)->nodeList(), 0);
    normalization.appendNewField("normalization", (*fieldItr)->nodeList(), 1.0e-30);
  }

  // This first pass counts the total weight from each donor node, so we can compute the normalization correctly.
  // Loop over all the positions in the donor fieldList.
  for (InternalNodeIterator<Dimension> nodeItr = position.internalNodeBegin();
       nodeItr < position.internalNodeEnd();
       ++nodeItr) {

    // Check if this node has been done yet.
    if (flagNodeDone(nodeItr) == 0) {

      // Set the neighbor info over the positions we're sampling to.
      vector<vector<int>> masterLists, coarseNeighbors, refineNeighbors,
                          masterListsSample, coarseNeighborsSample, refineNeighborsSample;
      position.setMasterNodeLists(position(nodeItr), Hfield(nodeItr), masterLists, coarseNeighbors);
      samplePositions.setMasterNodeLists(position(nodeItr), Hfield(nodeItr), masterListsSample, coarseNeighborsSample);

      // Loop over the set of master nodes in the FieldList we're sampling from.
      for (MasterNodeIterator<Dimension> masterItr = position.masterNodeBegin(masterLists);
           masterItr < position.masterNodeEnd();
           ++masterItr) {
        CHECK(flagNodeDone(masterItr) == 0);
   
        // Sample node (i) state.
        const Vector& ri = position(masterItr);
        const SymTensor& Hi = Hfield(masterItr);
        const Scalar weighti = weight(masterItr);

        // Refine the set of nodes we're sampling to for this position.
        samplePositions.setRefineNodeLists(ri, Hi, coarseNeighborsSample, refineNeighborsSample);

        // Loop over the refined neighbors, and determine the normalization
        // constant.
        for (RefineNodeIterator<Dimension> neighborItr = samplePositions.refineNodeBegin(refineNeighborsSample);
             neighborItr < samplePositions.refineNodeEnd();
             ++neighborItr) {

          // Node j's state.
          const Vector& rj = samplePositions(neighborItr);
          const SymTensor& Hj = sampleHfield(neighborItr);
          const Scalar weightj = sampleWeight(neighborItr);

          const Vector rij = ri - rj;
          const Scalar etai = (Hi*rij).magnitude();
          const Scalar etaj = (Hj*rij).magnitude();

          // Calculate the kernel estimates for each node.
          Scalar Wi = kernel(etai, 1.0);
          Scalar Wj = kernel(etaj, 1.0);

          const Scalar Wij = max(Wi, Wj);
          const Scalar weightij = 0.5*(weighti + weightj);

          // // Get the symmetrized kernel weighting for this node pair.
          // Scalar Wij, weightij;
          // switch(neighborItr.nodeListPtr()->neighbor().neighborSearchType()) {
          // case GatherScatter:
          //   Wij = 0.5*(Wi + Wj);
          //   weightij = 0.5*(weighti + weightj);
          //   break;

          // case Gather:
          //   Wij = Wi;
          //   weightij = weighti;
          //   break;

          // case Scatter:
          //   Wij = Wj;
          //   weightij = weightj;
          //   break;

          // default:
          //   VERIFY2(false, "Unhandled neighbor search type.");
          // }

          // Add this nodes contribution to the master value.
          normalization(masterItr) += weightij*Wij;
        }
        CHECK(normalization(masterItr) > 0.0);
        normalization(masterItr) = 1.0/normalization(masterItr);

        // Flag this master node as done.
        flagNodeDone(masterItr) = 1;
      }
    }
  }

  // After we're done, all nodes in the sampling from NodeLists should be flagged as done.
  BEGIN_CONTRACT_SCOPE
  {
    for (InternalNodeIterator<Dimension> nodeItr = position.internalNodeBegin();
         nodeItr != position.internalNodeEnd();
         ++nodeItr) CHECK2(flagNodeDone(nodeItr) == 1, nodeItr.fieldID() << " " << nodeItr.nodeID());
  }
  END_CONTRACT_SCOPE

  // Apply boundaries to the donor information.
  for (typename vector<Boundary<Dimension>*>::const_iterator bcItr = boundaries.begin();
       bcItr != boundaries.end();
       ++bcItr) {
    (*bcItr)->applyFieldListGhostBoundary(normalization);
    (*bcItr)->finalizeGhostBoundary();
  }

  // Iterate over the masters again and do the splatting.
  flagNodeDone = 0;
  for (AllNodeIterator<Dimension> nodeItr = position.nodeBegin();
       nodeItr < position.nodeEnd();
       ++nodeItr) {
    CHECK(flagNodeDone(nodeItr) == 0);

    // Set the neighbor info over the positions we're sampling to.
    vector<vector<int>> masterLists, coarseNeighbors, refineNeighbors,
                        masterListsSample, coarseNeighborsSample, refineNeighborsSample;
    position.setMasterNodeLists(position(nodeItr), Hfield(nodeItr), masterLists, coarseNeighbors);
    samplePositions.setMasterNodeLists(position(nodeItr), Hfield(nodeItr), masterListsSample, coarseNeighborsSample);

    // Sample node (i) state.
    const Vector& ri = position(nodeItr);
    const SymTensor& Hi = Hfield(nodeItr);
    const Scalar weighti = weight(nodeItr);

    // Refine the set of nodes we're donating from to for this position.
    samplePositions.setRefineNodeLists(ri, Hi, coarseNeighborsSample, refineNeighborsSample);

    // Loop over the refined neighbors again, and do the splat of the donor node
    // values to each of the sample nodes.
    for (RefineNodeIterator<Dimension> neighborItr = samplePositions.refineNodeBegin(refineNeighborsSample);
         neighborItr < samplePositions.refineNodeEnd();
         ++neighborItr) {

      // Node j's state.
      const Vector& rj = samplePositions(neighborItr);
      const SymTensor& Hj = sampleHfield(neighborItr);
      const Scalar weightj = sampleWeight(neighborItr);

      const Vector rij = ri - rj;
      const Scalar etai = (Hi*rij).magnitude();
      const Scalar etaj = (Hj*rij).magnitude();
       
      // Calculate the kernel estimates for each node.
      Scalar Wi = kernel(etai, 1.0);
      Scalar Wj = kernel(etaj, 1.0);

      const Scalar Wij = max(Wi, Wj);
      const Scalar weightij = 0.5*(weighti + weightj);

      // // Get the symmetrized kernel weighting for this node pair.
      // Scalar Wij, weightij;
      // switch(neighborItr.nodeListPtr()->neighbor().neighborSearchType()) {
      // case GatherScatter:
      //   Wij = 0.5*(Wi + Wj);
      //   weightij = 0.5*(weighti + weightj);
      //   break;

      // case Gather:
      //   Wij = Wj;
      //   weightij = weightj;
      //   break;
         
      // case Scatter:
      //   Wij = Wi;
      //   weightij = weighti;
      //   break;

      // default:
      //   VERIFY2(false, "Unhandled neighbor search type.");
      // }

      // Loop over all the FieldLists we're sampling from, and add their contributions
      // to their correspoding result FieldList.
      const Scalar localWeight = weightij*Wij*normalization(nodeItr);
      for (auto i = 0u; i < fieldListSet.ScalarFieldLists.size(); ++i) {
        const FieldList<Dimension, Scalar>& fieldList = fieldListSet.ScalarFieldLists[i];
        FieldList<Dimension, Scalar>& result = resultSet.ScalarFieldLists[i];
        result(neighborItr) += fieldList(nodeItr)*localWeight;
      }
      for (auto i = 0u; i < fieldListSet.VectorFieldLists.size(); ++i) {
        const FieldList<Dimension, Vector>& fieldList = fieldListSet.VectorFieldLists[i];
        FieldList<Dimension, Vector>& result = resultSet.VectorFieldLists[i];
        result(neighborItr) += fieldList(nodeItr)*localWeight;
      }
      for (auto i = 0u; i < fieldListSet.TensorFieldLists.size(); ++i) {
        const FieldList<Dimension, Tensor>& fieldList = fieldListSet.TensorFieldLists[i];
        FieldList<Dimension, Tensor>& result = resultSet.TensorFieldLists[i];
        result(neighborItr) += fieldList(nodeItr)*localWeight;
      }
      for (auto i = 0u; i < fieldListSet.SymTensorFieldLists.size(); ++i) {
        const FieldList<Dimension, SymTensor>& fieldList = fieldListSet.SymTensorFieldLists[i];
        FieldList<Dimension, SymTensor>& result = resultSet.SymTensorFieldLists[i];
        result(neighborItr) += fieldList(nodeItr)*localWeight;
      }
    }
   
    // Flag this master node as done.
    flagNodeDone(nodeItr) = 1;
  }

  // After we're done, all nodes in the sampling from NodeLists should be flagged as done.
  BEGIN_CONTRACT_SCOPE
  {
    for (InternalNodeIterator<Dimension> nodeItr = position.internalNodeBegin();
         nodeItr != position.internalNodeEnd();
         ++nodeItr) CHECK2(flagNodeDone(nodeItr) == 1, nodeItr.fieldID() << " " << nodeItr.nodeID());
  }
  END_CONTRACT_SCOPE

  return resultSet;
}

}
