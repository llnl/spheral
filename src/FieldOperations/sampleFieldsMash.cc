#include "FieldListFunctionsMash.hh"
#include "FieldOperations/FieldListConnectivity.hh"
#include "Field/FieldList.hh"
#include "Field/Field.hh"
#include "NodeList/NodeList.hh"
#include "Kernel/TableKernel.hh"

namespace Spheral {

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
  using Scalar = typename Dimension::Scalar;
  CHECK(fieldList.numFields() == position.numFields());
  CHECK(fieldList.numFields() == weight.numFields());
  CHECK(fieldList.numFields() == Hfield.numFields());
  for (auto fieldi = 0u; fieldi < fieldList.numFields(); ++fieldi) {
    CHECK(fieldList[fieldi]->nodeListPtr() == position[fieldi]->nodeListPtr());
    CHECK(fieldList[fieldi]->nodeListPtr() == weight[fieldi]->nodeListPtr());
    CHECK(fieldList[fieldi]->nodeListPtr() == Hfield[fieldi]->nodeListPtr());
  }
  CHECK(samplePositions.numFields() == sampleWeight.numFields());
  CHECK(samplePositions.numFields() == sampleHfield.numFields());
  for (auto samplei = 0u; samplei < samplePositions.numFields(); ++samplei) {
    CHECK(samplePositions[samplei]->nodeListPtr() == sampleWeight[samplei]->nodeListPtr());
    CHECK(samplePositions[samplei]->nodeListPtr() == sampleHfield[samplei]->nodeListPtr());
  }

  FieldList<Dimension, DataType> result;
  FieldList<Dimension, Scalar> normalization;
  result.copyFields();
  normalization.copyFields();
  for (auto samplei = 0u; samplei < samplePositions.numFields(); ++samplei) {
    auto& nodes = samplePositions[samplei]->nodeList();
    result.appendField(Field<Dimension, DataType>("sample mash", nodes));
    normalization.appendField(Field<Dimension, Scalar>("normalization", nodes, 1.0e-30));
  }

  FieldListPairConnectivity<Dimension, DataType, typename Dimension::Vector>
    connectivity(fieldList, samplePositions, kernel.kernelExtent());
  const auto searchType = connectivity.searchType();
  for (auto fieldi = 0u; fieldi < fieldList.numFields(); ++fieldi) {
    const auto numInternal = fieldList[fieldi]->nodeListPtr()->numInternalNodes();
    for (auto i = 0u; i < numInternal; ++i) {
      const auto& ri = position(fieldi, i);
      const auto& Hi = Hfield(fieldi, i);
      const auto weighti = weight(fieldi, i);
      connectivity.forEachSecondNeighbor(fieldi, i, [&](const size_t samplej, const size_t j) {
        const auto rij = ri - samplePositions(samplej, j);
        const auto Wi = kernel((Hi*rij).magnitude(), 1.0);
        const auto Wj = kernel((sampleHfield(samplej, j)*rij).magnitude(), 1.0);
        const auto weightj = sampleWeight(samplej, j);
        Scalar Wij, weightij;
        switch (searchType) {
        case NeighborSearchType::GatherScatter: Wij = 0.5*(Wi + Wj); weightij = 0.5*(weighti + weightj); break;
        case NeighborSearchType::Gather: Wij = Wi; weightij = weighti; break;
        case NeighborSearchType::Scatter: Wij = Wj; weightij = weightj; break;
        default: VERIFY2(false, "Unhandled neighbor search type.");
        }
        normalization(samplej, j) += weightij*Wij;
        result(samplej, j) += fieldList(fieldi, i)*weightij*Wij;
      });
    }
  }
  result /= normalization;
  return result;
}

}
