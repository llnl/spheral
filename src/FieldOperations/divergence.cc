//---------------------------------Spheral++----------------------------------//
// FieldListFunctions -- A set of global functions which can be applied to
// FieldLists.
//
// Created by JMO, Wed Dec  6 21:09:29 PST 2000
//----------------------------------------------------------------------------//

#include "FieldListFunctions.hh"
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
divergence(const FieldList<Dimension, DataType>& fieldList,
           const FieldList<Dimension, typename Dimension::Vector>& position,
           const FieldList<Dimension, typename Dimension::Scalar>& /*weight*/,
           const FieldList<Dimension, typename Dimension::Scalar>& mass,
           const FieldList<Dimension, typename Dimension::Scalar>& rho,
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
  for (typename FieldList<Dimension, DataType>::const_iterator
         fieldItr = fieldList.begin();
       fieldItr < fieldList.end(); 
       ++fieldItr) {
    result.appendField(Field<Dimension, DivergenceType>("div", (*fieldItr)->nodeList()));
  }
  FieldListConnectivity<Dimension, DataType> connectivity(fieldList, kernel.kernelExtent());
  const auto searchType = connectivity.searchType();
  for (auto fieldi = 0u; fieldi < fieldList.numFields(); ++fieldi) {
    const auto numInternal = fieldList[fieldi]->nodeListPtr()->numInternalNodes();
    for (auto i = 0u; i < numInternal; ++i) {
      const Vector& ri = position(fieldi, i);
      const SymTensor& Hi = Hfield(fieldi, i);
      const Scalar& rhoi = rho(fieldi, i);
      const DataType& fi = fieldList(fieldi, i);
      connectivity.forEachNeighbor(fieldi, i, [&](const size_t fieldj, const size_t j) {
          const Vector& rj = position(fieldj, j);
          const SymTensor& Hj = Hfield(fieldj, j);
          const Scalar& mj = mass(fieldj, j);
          const DataType& fj = fieldList(fieldj, j);

          const Vector rij = ri - rj;
          const Vector etai = Hi*rij;
          const Vector etaj = Hj*rij;
          const Vector etaiNorm = etai.unitVector();
          const Vector etajNorm = etaj.unitVector();

          // Get the symmetrized kernel gradient for this node pair.
          Vector gradWij;
          switch(searchType) {
          case NeighborSearchType::GatherScatter:
            gradWij = 0.5*(Hi*etaiNorm*kernel.grad(etai.magnitude(), Hi.Determinant()) + 
                           Hj*etajNorm*kernel.grad(etaj.magnitude(), Hj.Determinant()));
            break;

          case NeighborSearchType::Gather:
            gradWij = Hi*etaiNorm*kernel.grad(etai.magnitude(), Hi.Determinant());
            break;

          case NeighborSearchType::Scatter:
            gradWij = Hj*etajNorm*kernel.grad(etaj.magnitude(), Hj.Determinant());
            break;

          default:
            VERIFY2(false, "Unhandled neighbor search type.");
          }

          // Add this nodes contribution to the master value.
          result(fieldi, i) += mj*(fj - fi).dot(gradWij);
      });

        // Normalize by density.
        CHECK(rhoi > 0.0);
      result(fieldi, i) /= rhoi;
    }
  }

  return result;
}

}
