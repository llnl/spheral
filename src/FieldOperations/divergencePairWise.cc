//---------------------------------Spheral++----------------------------------//
// divergencePairWise -- Calculate a low order estimate of the gradient based on
// the pairwise formulation.
//
// Created by JMO, Sat Apr 12 15:58:53 PDT 2003
//----------------------------------------------------------------------------//

#include "PairWiseFieldListFunctions.hh"
#include "FieldOperations/FieldListConnectivity.hh"
#include "Field/FieldList.hh"
#include "Field/Field.hh"
#include "NodeList/NodeList.hh"
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
  result.copyFields();
  for (typename FieldList<Dimension, DataType>::const_iterator
         fieldItr = fieldList.begin();
       fieldItr < fieldList.end(); 
       ++fieldItr) {
    result.appendField(Field<Dimension, DivergenceType>("div pair wise", (*fieldItr)->nodeList()));
  }
  FieldListConnectivity<Dimension, DataType> connectivity(fieldList, kernel.kernelExtent());
  const auto searchType = connectivity.searchType();
  for (auto fieldi = 0u; fieldi < fieldList.numFields(); ++fieldi) {
    const auto numInternal = fieldList[fieldi]->nodeListPtr()->numInternalNodes();
    for (auto i = 0u; i < numInternal; ++i) {
      const Vector& ri = position(fieldi, i);
      const SymTensor& Hi = Hfield(fieldi, i);
      const DataType& fi = fieldList(fieldi, i);
      Tensor gradi(0.0);
      Tensor normalization(0.0);
      connectivity.forEachNeighbor(fieldi, i, [&](const size_t fieldj, const size_t j) {
            const Vector& rj = position(fieldj, j);
            const SymTensor& Hj = Hfield(fieldj, j);
            const Scalar& weightj = weight(fieldj, j);
            const DataType& fj = fieldList(fieldj, j);

            const Vector rij = ri - rj;
            const Vector etai = Hi*rij;
            const Vector etaj = Hj*rij;

            // Get the symmetrized kernel weighting for this node pair.
            Scalar Wij;
            switch(searchType) {
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
            const DataType fij = fi - fj;
            incrementPairWiseDelta<Dimension>(gradi,
                                              normalization,
                                              wij,
                                              fij,
                                              rij,
                                              Hi);

      });
        renormalizeGradient<Dimension>(gradi,
                                       normalization);

        // Now contract to the divergence.
        result(fieldi, i) = gradi.Trace();

    }
  }

  return result;
}

}
