//---------------------------------Spheral++----------------------------------//
// FieldListFunctionsMash -- A set of global functions which can be applied to
// FieldLists using MASH prescriptions.
//
// Created by JMO, Tue Jun  4 17:07:31 PDT 2002
//----------------------------------------------------------------------------//
#include "FieldListFunctionsMash.hh"
#include "FieldOperations/FieldListConnectivity.hh"
#include "Field/FieldList.hh"
#include "Field/FieldListSet.hh"
#include "Field/Field.hh"
#include "NodeList/NodeList.hh"
#include "NodeList/FluidNodeList.hh"
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
  
  FieldList<Dimension, Scalar> normalization(FieldStorageType::CopyFields);
  for (typename FieldList<Dimension, Vector>::const_iterator fieldItr = position.begin();
       fieldItr < position.end(); 
       ++fieldItr) {
    normalization.appendNewField("normalization", (*fieldItr)->nodeList(), 1.0e-30);
  }

  FieldListPairConnectivity<Dimension, Vector, Vector>
    connectivity(position, samplePositions, kernel.kernelExtent(), true);

  // First normalize each internal donor.
  for (auto fieldi = 0u; fieldi < position.numFields(); ++fieldi) {
    const auto numInternal = position[fieldi]->nodeListPtr()->numInternalNodes();
    for (auto i = 0u; i < numInternal; ++i) {
        const Vector& ri = position(fieldi, i);
        const SymTensor& Hi = Hfield(fieldi, i);
        const Scalar weighti = weight(fieldi, i);
        connectivity.forEachSecondNeighbor(fieldi, i, [&](const size_t samplej, const size_t j) {
          const Vector& rj = samplePositions(samplej, j);
          const SymTensor& Hj = sampleHfield(samplej, j);
          const Scalar weightj = sampleWeight(samplej, j);

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
          // Historical Gather/Scatter alternatives intentionally replaced by max(Wi, Wj).
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
          normalization(fieldi, i) += weightij*Wij;
        });
        CHECK(normalization(fieldi, i) > 0.0);
        normalization(fieldi, i) = 1.0/normalization(fieldi, i);
    }
  }

  // Apply boundaries to the donor information.
  for (typename vector<Boundary<Dimension>*>::const_iterator bcItr = boundaries.begin();
       bcItr != boundaries.end();
       ++bcItr) {
    (*bcItr)->applyFieldListGhostBoundary(normalization);
    (*bcItr)->finalizeGhostBoundary();
  }

  // Splat internal and ghost donors after the boundary normalization update.
  for (auto fieldi = 0u; fieldi < position.numFields(); ++fieldi) {
    const auto numNodes = position[fieldi]->nodeListPtr()->numNodes();
    for (auto i = 0u; i < numNodes; ++i) {
      const Vector& ri = position(fieldi, i);
      const SymTensor& Hi = Hfield(fieldi, i);
      const Scalar weighti = weight(fieldi, i);
      connectivity.forEachSecondNeighbor(fieldi, i, [&](const size_t samplej, const size_t j) {
      const Vector& rj = samplePositions(samplej, j);
      const SymTensor& Hj = sampleHfield(samplej, j);
      const Scalar weightj = sampleWeight(samplej, j);

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
      // Historical Gather/Scatter alternatives intentionally replaced by max(Wi, Wj).
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
      const Scalar localWeight = weightij*Wij*normalization(fieldi, i);
      for (auto k = 0u; k < fieldListSet.ScalarFieldLists.size(); ++k) {
        const FieldList<Dimension, Scalar>& fieldList = fieldListSet.ScalarFieldLists[k];
        FieldList<Dimension, Scalar>& result = resultSet.ScalarFieldLists[k];
        result(samplej, j) += fieldList(fieldi, i)*localWeight;
      }
      for (auto k = 0u; k < fieldListSet.VectorFieldLists.size(); ++k) {
        const FieldList<Dimension, Vector>& fieldList = fieldListSet.VectorFieldLists[k];
        FieldList<Dimension, Vector>& result = resultSet.VectorFieldLists[k];
        result(samplej, j) += fieldList(fieldi, i)*localWeight;
      }
      for (auto k = 0u; k < fieldListSet.TensorFieldLists.size(); ++k) {
        const FieldList<Dimension, Tensor>& fieldList = fieldListSet.TensorFieldLists[k];
        FieldList<Dimension, Tensor>& result = resultSet.TensorFieldLists[k];
        result(samplej, j) += fieldList(fieldi, i)*localWeight;
      }
      for (auto k = 0u; k < fieldListSet.SymTensorFieldLists.size(); ++k) {
        const FieldList<Dimension, SymTensor>& fieldList = fieldListSet.SymTensorFieldLists[k];
        FieldList<Dimension, SymTensor>& result = resultSet.SymTensorFieldLists[k];
        result(samplej, j) += fieldList(fieldi, i)*localWeight;
      }
      });
    }
  }

  return resultSet;
}

}
