//---------------------------------Spheral++----------------------------------//
// FieldListFunctionsMash -- A set of global functions which can be applied to
// FieldLists using MASH prescriptions.
//
// Created by JMO, Wed Dec  6 21:09:29 PST 2000
//----------------------------------------------------------------------------//

#include "FieldListFunctionsMash.hh"
#include "FieldOperations/FieldListConnectivity.hh"
#include "Field/FieldList.hh"
#include "Field/Field.hh"
#include "NodeList/NodeList.hh"
#include "Kernel/TableKernel.hh"
#include "Geometry/MathTraits.hh"

namespace Spheral {

using std::vector;

//------------------------------------------------------------------------------
// Calculate the divergence of a FieldList.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
FieldList<Dimension, typename MathTraits<Dimension, DataType>::DivergenceType>
divergenceMash(const FieldList<Dimension, DataType>& fieldList,
               const FieldList<Dimension, typename Dimension::Vector>& position,
               const FieldList<Dimension, typename Dimension::Scalar>& weight,
               const FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
               const TableKernel<Dimension>& kernel) {

  // Typedef's to ease typing/understandability.
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename MathTraits<Dimension, DataType>::DivergenceType DivergenceType;

  // Return FieldList.
  FieldList<Dimension, DivergenceType> result;
  result.copyFields();
  for (typename FieldList<Dimension, DataType>::const_iterator fieldItr = fieldList.begin();
       fieldItr < fieldList.end(); 
       ++fieldItr) {
    result.appendField(Field<Dimension, DivergenceType>("div mash", (*fieldItr)->nodeList()));
  }
  FieldListConnectivity<Dimension, DataType> connectivity(fieldList, kernel.kernelExtent());
  const auto searchType = connectivity.searchType();
  for (auto fieldi = 0u; fieldi < fieldList.numFields(); ++fieldi) {
    const auto numInternal = fieldList[fieldi]->nodeListPtr()->numInternalNodes();
    for (auto i = 0u; i < numInternal; ++i) {
      const Vector& ri = position(fieldi, i);
      const SymTensor& Hi = Hfield(fieldi, i);
      const DataType& fi = fieldList(fieldi, i);
      Scalar normalization = weight(fieldi, i)*kernel(0.0, 1.0);
      connectivity.forEachNeighbor(fieldi, i, [&](const size_t fieldj, const size_t j) {
          const Vector& rj = position(fieldj, j);
          const SymTensor& Hj = Hfield(fieldj, j);
          const Scalar& weightj = weight(fieldj, j);
          const DataType& fj = fieldList(fieldj, j);

          const Vector rij = ri - rj;
          const Vector etai = Hi*rij;
          const Vector etaj = Hj*rij;
          const Vector etaiNorm = etai.unitVector();
          const Vector etajNorm = etaj.unitVector();
          const Scalar etaiMag = etai.magnitude();
          const Scalar etajMag = etaj.magnitude();

          // Get the symmetrized kernel gradient for this node pair.
          Scalar Wij;
          Vector gradWij;
          switch(searchType) {
          case NeighborSearchType::GatherScatter:
            Wij = 0.5*(kernel(etaiMag, 1.0) + kernel(etajMag, 1.0));
            gradWij = 0.5*(Hi*etaiNorm*kernel.grad(etaiMag, 1.0) + 
                           Hj*etajNorm*kernel.grad(etajMag, 1.0));
            break;

          case NeighborSearchType::Gather:
            Wij = kernel(etaiMag, 1.0);
            gradWij = Hi*etaiNorm*kernel.grad(etaiMag, 1.0);
            break;

          case NeighborSearchType::Scatter:
            Wij = kernel(etajMag, 1.0);
            gradWij = Hj*etajNorm*kernel.grad(etajMag, 1.0);
            break;

          default:
            VERIFY2(false, "Unhandled neighbor search type.");
          }

          // Add this nodes contribution to the master value.
          normalization += weightj*Wij;
          result(fieldi, i) += weightj*(fj - fi).dot(gradWij);

      });

        // Normalize the divergence value.
        CHECK(normalization > 0.0);
      result(fieldi, i) /= normalization;
    }
  }

  return result;
}

}
