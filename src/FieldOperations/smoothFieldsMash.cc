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
// Return a monotonic smoothed estimate of the given FieldList by using a
// MASH like prescription.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
FieldList<Dimension, DataType>
smoothFieldsMash(const FieldList<Dimension, DataType>& fieldList,
                 const FieldList<Dimension, typename Dimension::Vector>& position,
                 const FieldList<Dimension, typename Dimension::Scalar>& weight,
                 const FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
                 const TableKernel<Dimension>& kernel) {

  // Some convenient typedefs.
  typedef typename Dimension::Scalar Scalar;

  // Return FieldList.
  FieldList<Dimension, DataType> result;
  result.copyFields();
  for (auto fieldi = 0u; fieldi < fieldList.numFields(); ++fieldi) {
    result.appendField(Field<Dimension, DataType>("smooth", fieldList[fieldi]->nodeList()));
  }

  FieldListConnectivity<Dimension, DataType> connectivity(fieldList, kernel.kernelExtent());
  const auto searchType = connectivity.searchType();
  for (auto fieldi = 0u; fieldi < fieldList.numFields(); ++fieldi) {
    const auto numInternal = fieldList[fieldi]->nodeListPtr()->numInternalNodes();
    for (auto i = 0u; i < numInternal; ++i) {
      const auto& ri = position(fieldi, i);
      const auto& Hi = Hfield(fieldi, i);
      Scalar totalWeight = 0.0;
      const auto addContribution = [&](const size_t fieldj, const size_t j) {
          const auto rij = ri - position(fieldj, j);
          const auto Wi = kernel((Hi*rij).magnitude(), 1.0);
          const auto Wj = kernel((Hfield(fieldj, j)*rij).magnitude(), 1.0);
          Scalar Wij = 0.0;
          switch(searchType) {
          case NeighborSearchType::GatherScatter:
            Wij = 0.5*(Wi + Wj);
            break;
          case NeighborSearchType::Gather:
            Wij = Wi;
            break;
          case NeighborSearchType::Scatter:
            Wij = Wj;
            break;
          default:
            VERIFY2(false, "Unhandled neighbor search type.");
          }
          totalWeight += weight(fieldj, j)*Wij;
          result(fieldi, i) += fieldList(fieldj, j)*weight(fieldj, j)*Wij;
      };
      addContribution(fieldi, i);
      connectivity.forEachNeighbor(fieldi, i, addContribution);
      CHECK(totalWeight > 0.0);
      result(fieldi, i) /= totalWeight;
    }
  }

  return result;
}

}
