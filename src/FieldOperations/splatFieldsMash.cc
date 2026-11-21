#include "FieldListFunctionsMash.hh"
#include "FieldOperations/FieldListConnectivity.hh"
#include "Field/FieldList.hh"
#include "Field/Field.hh"
#include "NodeList/NodeList.hh"
#include "Kernel/TableKernel.hh"

namespace Spheral {

template<typename Dimension, typename DataType>
FieldList<Dimension, DataType>
splatFieldsMash(const FieldList<Dimension, DataType>& fieldList,
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
  result.copyFields();
  for (auto samplei = 0u; samplei < samplePositions.numFields(); ++samplei) {
    result.appendField(Field<Dimension, DataType>("splat mash", samplePositions[samplei]->nodeList()));
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
      Scalar totalWeight = 1.0e-30;
      connectivity.forEachSecondNeighbor(fieldi, i, [&](const size_t samplej, const size_t j) {
        const auto rij = ri - samplePositions(samplej, j);
        const auto Wi = kernel((Hi*rij).magnitude(), 1.0);
        const auto Wj = kernel((sampleHfield(samplej, j)*rij).magnitude(), 1.0);
        const auto weightj = sampleWeight(samplej, j);
        switch (searchType) {
        case NeighborSearchType::GatherScatter: totalWeight += 0.25*(weighti + weightj)*(Wi + Wj); break;
        case NeighborSearchType::Gather: totalWeight += weighti*Wi; break;
        case NeighborSearchType::Scatter: totalWeight += weightj*Wj; break;
        default: VERIFY2(false, "Unhandled neighbor search type.");
        }
      });
      CHECK(totalWeight > 0.0);
      connectivity.forEachSecondNeighbor(fieldi, i, [&](const size_t samplej, const size_t j) {
        const auto rij = ri - samplePositions(samplej, j);
        const auto Wi = kernel((Hi*rij).magnitude(), 1.0);
        const auto Wj = kernel((sampleHfield(samplej, j)*rij).magnitude(), 1.0);
        const auto weightj = sampleWeight(samplej, j);
        Scalar pairWeight;
        switch (searchType) {
        case NeighborSearchType::GatherScatter: pairWeight = 0.25*(weighti + weightj)*(Wi + Wj); break;
        case NeighborSearchType::Gather: pairWeight = weighti*Wi; break;
        case NeighborSearchType::Scatter: pairWeight = weightj*Wj; break;
        default: VERIFY2(false, "Unhandled neighbor search type.");
        }
        result(samplej, j) += fieldList(fieldi, i)*pairWeight/totalWeight;
      });
    }
  }
  return result;
}

}
