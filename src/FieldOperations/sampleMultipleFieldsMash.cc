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
#include "Neighbor/Neighbor.hh"
#include "Kernel/TableKernel.hh"
#include "Geometry/MathTraits.hh"

#include "Utilities/DBC.hh"

namespace Spheral {

using std::vector;

//------------------------------------------------------------------------------
// Return a MASH sampled version of the given FieldList at the new positions.
// This version operates on a set of FieldLists at the same time, purely as an 
// optimzation over the single FieldList version, sampleFieldsMash.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldListSet<Dimension>
sampleMultipleFieldsMash(const FieldListSet<Dimension>& fieldListSet,
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

  // Return FieldListSet.
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
      resultItr->appendField(Field<Dimension, Scalar>("sample" + (*fieldItr)->name(), (*fieldItr)->nodeList()));
    for (typename vector<FieldList<Dimension, Vector> >::iterator resultItr = 
           resultSet.VectorFieldLists.begin();
         resultItr < resultSet.VectorFieldLists.end();
         ++resultItr) 
      resultItr->appendField(Field<Dimension, Vector>("sample" + (*fieldItr)->name(), (*fieldItr)->nodeList()));
    for (typename vector<FieldList<Dimension, Tensor> >::iterator resultItr = 
           resultSet.TensorFieldLists.begin();
         resultItr < resultSet.TensorFieldLists.end();
         ++resultItr) 
      resultItr->appendField(Field<Dimension, Tensor>("sample" + (*fieldItr)->name(), (*fieldItr)->nodeList()));
    for (typename vector<FieldList<Dimension, SymTensor> >::iterator resultItr = 
           resultSet.SymTensorFieldLists.begin();
         resultItr < resultSet.SymTensorFieldLists.end();
         ++resultItr) 
      resultItr->appendField(Field<Dimension, SymTensor>("sample" + (*fieldItr)->name(), (*fieldItr)->nodeList()));
  }
  vector< vector<bool> > flagNodeDone(samplePositions.numFields());
  for (typename FieldList<Dimension, Vector>::const_iterator fieldItr = samplePositions.begin();
       fieldItr < samplePositions.end(); 
       ++fieldItr) {
    flagNodeDone[fieldItr - samplePositions.begin()].resize((*fieldItr)->nodeListPtr()->numInternalNodes(), false);
  }

  // Loop over all the positions in the fieldList we're sampling to.
  for (InternalNodeIterator<Dimension> nodeItr = samplePositions.internalNodeBegin();
       nodeItr < samplePositions.internalNodeEnd();
       ++nodeItr) {

    // Check if this node has been done yet.
    if (!flagNodeDone[nodeItr.fieldID()][nodeItr.nodeID()]) {

      // Set the neighbor info over the positions we're sampling to.
      vector<vector<int>> masterLists, coarseNeighbors, refineNeighbors,
                          masterListsSample, coarseNeighborsSample, refineNeighborsSample;
      position.setMasterNodeLists(samplePositions(nodeItr), sampleHfield(nodeItr), masterLists, coarseNeighbors);
      samplePositions.setMasterNodeLists(samplePositions(nodeItr), sampleHfield(nodeItr), masterListsSample, coarseNeighborsSample);

      // Loop over the set of master nodes in the FieldList we're sampling to.
      for (MasterNodeIterator<Dimension> masterItr = samplePositions.masterNodeBegin(masterListsSample);
           masterItr < samplePositions.masterNodeEnd();
           ++masterItr) {
        CHECK(flagNodeDone[masterItr.fieldID()][masterItr.nodeID()] == false);
      
        // Sample node (i) state.
        const Vector& ri = samplePositions(masterItr);
        const SymTensor& Hi = sampleHfield(masterItr);
        const Scalar& weighti = sampleWeight(masterItr);

        // Refine the set of nodes we're sampling from for this position.
        position.setRefineNodeLists(ri, Hi, coarseNeighbors, refineNeighbors);

        // Loop over the refined neighbors, and calculate the local values.
        Scalar normalization = 1.0e-30;
        for (RefineNodeIterator<Dimension> neighborItr = position.refineNodeBegin(refineNeighbors);
             neighborItr < position.refineNodeEnd();
             ++neighborItr) {

          // Node j's state.
          const Vector& rj = position(neighborItr);
          const SymTensor& Hj = Hfield(neighborItr);
          const Scalar& weightj = weight(neighborItr);

          const Vector rij = ri - rj;
          const Scalar etai = (Hi*rij).magnitude();
          const Scalar etaj = (Hj*rij).magnitude();

          // Calculate the kernel estimates for each node.
          Scalar Wi = kernel(etai, 1.0);
          Scalar Wj = kernel(etaj, 1.0);

          // Get the symmetrized kernel weighting for this node pair.
          Scalar Wij, weightij;
          switch(masterItr.nodeListPtr()->neighbor().neighborSearchType()) {
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

          // Increment the MASH normalization.
          const Scalar localWeight = weightij*Wij;
          normalization += localWeight;

          // Loop over all the result FieldLists, and increment their values 
          // for this pair interaction.
          for (auto i = 0u; i < resultSet.ScalarFieldLists.size(); ++i) {
            FieldList<Dimension, Scalar>& result = resultSet.ScalarFieldLists[i];
            const FieldList<Dimension, Scalar>& fieldList = fieldListSet.ScalarFieldLists[i];
            result(masterItr) += fieldList(neighborItr)*localWeight;
          }
          for (auto i = 0u; i < resultSet.VectorFieldLists.size(); ++i) {
            FieldList<Dimension, Vector>& result = resultSet.VectorFieldLists[i];
            const FieldList<Dimension, Vector>& fieldList = fieldListSet.VectorFieldLists[i];
            result(masterItr) += fieldList(neighborItr)*localWeight;
          }
          for (auto i = 0u; i < resultSet.TensorFieldLists.size(); ++i) {
            FieldList<Dimension, Tensor>& result = resultSet.TensorFieldLists[i];
            const FieldList<Dimension, Tensor>& fieldList = fieldListSet.TensorFieldLists[i];
            result(masterItr) += fieldList(neighborItr)*localWeight;
          }
          for (auto i = 0u; i < resultSet.SymTensorFieldLists.size(); ++i) {
            FieldList<Dimension, SymTensor>& result = resultSet.SymTensorFieldLists[i];
            const FieldList<Dimension, SymTensor>& fieldList = fieldListSet.SymTensorFieldLists[i];
            result(masterItr) += fieldList(neighborItr)*localWeight;
          }
        }

        // Normalize the local estimates.
        CHECK(normalization > 0.0);
        for (auto i = 0u; i < resultSet.ScalarFieldLists.size(); ++i) {
          FieldList<Dimension, Scalar>& result = resultSet.ScalarFieldLists[i];
          result(masterItr) /= normalization;
        }

        // Flag this master node as done.
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
      cerr << "Error in FieldList::sampleFieldsMash: Not all values determined on exit "
           << checkcount << endl;
    }
    CHECK(checkcount == 0);
  }

  return resultSet;
}

}
