#include "FieldListFunctionsMash.hh"
#include "FieldOperations/FieldListConnectivity.hh"
#include "Field/FieldList.hh"
#include "Field/Field.hh"
#include "NodeList/NodeList.hh"
#include "Kernel/TableKernel.hh"

namespace Spheral {

template<typename Dimension, typename DataType>
FieldList<Dimension, DataType>
smoothFieldsMash2(const FieldList<Dimension, DataType>& fieldList,
                  const FieldList<Dimension, typename Dimension::Vector>& position,
                  const FieldList<Dimension, typename Dimension::Scalar>& weight,
                  const FieldList<Dimension, typename Dimension::Scalar>& weightDensity,
                  const FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
                  const TableKernel<Dimension>& kernel) {
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;

  FieldList<Dimension, DataType> result;
  FieldList<Dimension, Scalar> a, b;
  result.copyFields();
  a.copyFields();
  b.copyFields();
  for (auto fieldi = 0u; fieldi < fieldList.numFields(); ++fieldi) {
    auto& nodes = fieldList[fieldi]->nodeList();
    result.appendField(Field<Dimension, DataType>("smooth", nodes));
    a.appendField(Field<Dimension, Scalar>("a", nodes));
    b.appendField(Field<Dimension, Scalar>("b", nodes));
  }

  FieldListConnectivity<Dimension, DataType> connectivity(fieldList, kernel.kernelExtent());
  const auto searchType = connectivity.searchType();
  const auto kernelWeight = [&](const auto& Hi, const auto& Hj, const auto& rij) {
    const auto Wi = kernel((Hi*rij).magnitude(), 1.0);
    const auto Wj = kernel((Hj*rij).magnitude(), 1.0);
    switch (searchType) {
    case NeighborSearchType::GatherScatter: return 0.5*(Wi + Wj);
    case NeighborSearchType::Gather: return Wi;
    case NeighborSearchType::Scatter: return Wj;
    default: VERIFY2(false, "Unhandled neighbor search type."); return 0.0;
    }
  };

  for (auto fieldi = 0u; fieldi < fieldList.numFields(); ++fieldi) {
    const auto numInternal = fieldList[fieldi]->nodeListPtr()->numInternalNodes();
    for (auto i = 0u; i < numInternal; ++i) {
      const auto& ri = position(fieldi, i);
      const auto& Hi = Hfield(fieldi, i);
      Scalar totalWeight = weight(fieldi, i)*kernel(0.0, 1.0);
      Vector rCenter;
      connectivity.forEachNeighbor(fieldi, i, [&](const size_t fieldj, const size_t j) {
        const auto rij = ri - position(fieldj, j);
        const auto w = weight(fieldj, j)*kernelWeight(Hi, Hfield(fieldj, j), rij);
        totalWeight += w;
        rCenter -= w*rij;
      });
      CHECK(totalWeight > 0.0);
      rCenter /= totalWeight;
      a(fieldi, i) = weightDensity(fieldi, i);
      b(fieldi, i) = 2.0*weightDensity(fieldi, i)*rCenter.x()*Hi.xx()*Hi.xx();
    }
  }

  for (auto fieldi = 0u; fieldi < fieldList.numFields(); ++fieldi) {
    const auto numInternal = fieldList[fieldi]->nodeListPtr()->numInternalNodes();
    for (auto i = 0u; i < numInternal; ++i) {
      const auto& ri = position(fieldi, i);
      const auto& Hi = Hfield(fieldi, i);
      Scalar totalWeight = weight(fieldi, i)*kernel(0.0, 1.0)*a(fieldi, i);
      result(fieldi, i) = fieldList(fieldi, i)*totalWeight;
      connectivity.forEachNeighbor(fieldi, i, [&](const size_t fieldj, const size_t j) {
        const auto rij = ri - position(fieldj, j);
        const auto Wij = kernelWeight(Hi, Hfield(fieldj, j), rij)*
                         (a(fieldi, i) + b(fieldi, i)*rij.x());
        const auto w = weight(fieldj, j)*Wij;
        totalWeight += w;
        result(fieldi, i) += fieldList(fieldj, j)*w;
      });
      CHECK(totalWeight > 0.0);
      result(fieldi, i) /= totalWeight;
    }
  }
  return result;
}

}
