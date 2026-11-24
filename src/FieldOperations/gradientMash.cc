#include "FieldListFunctionsMash.hh"
#include "FieldOperations/FieldListConnectivity.hh"
#include "Field/FieldList.hh"
#include "Field/Field.hh"
#include "NodeList/NodeList.hh"
#include "Kernel/TableKernel.hh"
#include "Geometry/MathTraits.hh"

namespace Spheral {

template<typename Dimension, typename DataType>
FieldList<Dimension, typename MathTraits<Dimension, DataType>::GradientType>
gradientMash(const FieldList<Dimension, DataType>& fieldList,
             const FieldList<Dimension, typename Dimension::Vector>& position,
             const FieldList<Dimension, typename Dimension::Scalar>& weight,
             const FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
             const TableKernel<Dimension>& kernel) {
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using GradientType = typename MathTraits<Dimension, DataType>::GradientType;
  FieldList<Dimension, GradientType> result;
  result.copyFields();
  for (auto fieldi = 0u; fieldi < fieldList.numFields(); ++fieldi) {
    result.appendField(Field<Dimension, GradientType>("grad mash", fieldList[fieldi]->nodeList()));
  }
  FieldListConnectivity<Dimension, DataType> connectivity(fieldList, kernel.kernelExtent());
  const auto searchType = connectivity.searchType();
  for (auto fieldi = 0u; fieldi < fieldList.numFields(); ++fieldi) {
    const auto numInternal = fieldList[fieldi]->nodeListPtr()->numInternalNodes();
    for (auto i = 0u; i < numInternal; ++i) {
      const auto& ri = position(fieldi, i);
      const auto& Hi = Hfield(fieldi, i);
      const auto& Fi = fieldList(fieldi, i);
      Scalar normalization = weight(fieldi, i)*kernel(0.0, 1.0);
      connectivity.forEachNeighbor(fieldi, i, [&](const size_t fieldj, const size_t j) {
        const auto rij = ri - position(fieldj, j);
        const auto& Hj = Hfield(fieldj, j);
        const auto etai = Hi*rij;
        const auto etaj = Hj*rij;
        const auto etaiMag = etai.magnitude();
        const auto etajMag = etaj.magnitude();
        Scalar Wij;
        Vector gradWij;
        switch (searchType) {
        case NeighborSearchType::GatherScatter:
          Wij = 0.5*(kernel(etaiMag, 1.0) + kernel(etajMag, 1.0));
          gradWij = 0.5*(Hi*etai.unitVector()*kernel.grad(etaiMag, 1.0) +
                         Hj*etaj.unitVector()*kernel.grad(etajMag, 1.0));
          break;
        case NeighborSearchType::Gather:
          Wij = kernel(etaiMag, 1.0);
          gradWij = Hi*etai.unitVector()*kernel.grad(etaiMag, 1.0);
          break;
        case NeighborSearchType::Scatter:
          Wij = kernel(etajMag, 1.0);
          gradWij = Hj*etaj.unitVector()*kernel.grad(etajMag, 1.0);
          break;
        default:
          VERIFY2(false, "Unhandled neighbor search type.");
        }
        const auto weightj = weight(fieldj, j);
        normalization += weightj*Wij;
        result(fieldi, i) += weightj*(fieldList(fieldj, j) - Fi)*gradWij;
      });
      CHECK(normalization > 0.0);
      result(fieldi, i) /= normalization;
    }
  }
  return result;
}

template<typename Dimension, typename DataType>
FieldList<Dimension, typename MathTraits<Dimension, DataType>::GradientType>
gradientMash2(const FieldList<Dimension, DataType>& fieldList,
              const FieldList<Dimension, typename Dimension::Vector>& position,
              const FieldList<Dimension, typename Dimension::Scalar>& weight,
              const FieldList<Dimension, typename Dimension::Scalar>& weightDensity,
              const FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
              const TableKernel<Dimension>& kernel) {
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using GradientType = typename MathTraits<Dimension, DataType>::GradientType;
  FieldList<Dimension, GradientType> result;
  FieldList<Dimension, DataType> smoothFieldList;
  FieldList<Dimension, Scalar> a, b;
  result.copyFields();
  smoothFieldList.copyFields();
  a.copyFields();
  b.copyFields();
  for (auto fieldi = 0u; fieldi < fieldList.numFields(); ++fieldi) {
    auto& nodes = fieldList[fieldi]->nodeList();
    result.appendField(Field<Dimension, GradientType>("gradient", nodes));
    smoothFieldList.appendField(Field<Dimension, DataType>("smooth", nodes));
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
      const auto Wself = kernel(0.0, 1.0);
      const auto weighti = weight(fieldi, i);
      Scalar normalization = weighti*Wself*a(fieldi, i);
      smoothFieldList(fieldi, i) = fieldList(fieldi, i)*normalization;
      connectivity.forEachNeighbor(fieldi, i, [&](const size_t fieldj, const size_t j) {
        const auto rij = ri - position(fieldj, j);
        const auto correction = a(fieldi, i) + b(fieldi, i)*rij.x();
        const auto w = weight(fieldj, j)*kernelWeight(Hi, Hfield(fieldj, j), rij)*correction;
        normalization += w;
        smoothFieldList(fieldi, i) += fieldList(fieldj, j)*w;
      });
      CHECK(normalization > 0.0);
      smoothFieldList(fieldi, i) /= normalization;
    }
  }

  for (auto fieldi = 0u; fieldi < fieldList.numFields(); ++fieldi) {
    const auto numInternal = fieldList[fieldi]->nodeListPtr()->numInternalNodes();
    for (auto i = 0u; i < numInternal; ++i) {
      const auto& ri = position(fieldi, i);
      const auto& Hi = Hfield(fieldi, i);
      const auto Wself = kernel(0.0, 1.0);
      const auto weighti = weight(fieldi, i);
      Scalar normalization = weighti*Wself*a(fieldi, i);
      result(fieldi, i) += weighti*(fieldList(fieldi, i) - smoothFieldList(fieldi, i))*
                           Vector(-b(fieldi, i)*Wself*a(fieldi, i));
      connectivity.forEachNeighbor(fieldi, i, [&](const size_t fieldj, const size_t j) {
        const auto rij = ri - position(fieldj, j);
        const auto& Hj = Hfield(fieldj, j);
        const auto etai = Hi*rij;
        const auto etaj = Hj*rij;
        const auto etaiMag = etai.magnitude();
        const auto etajMag = etaj.magnitude();
        Scalar Wij;
        Vector gradWij;
        switch (searchType) {
        case NeighborSearchType::GatherScatter:
          Wij = 0.5*(kernel(etaiMag, 1.0) + kernel(etajMag, 1.0));
          gradWij = 0.5*(Hi*etai.unitVector()*kernel.grad(etaiMag, 1.0) +
                         Hj*etaj.unitVector()*kernel.grad(etajMag, 1.0));
          break;
        case NeighborSearchType::Gather:
          Wij = kernel(etaiMag, 1.0);
          gradWij = Hi*etai.unitVector()*kernel.grad(etaiMag, 1.0);
          break;
        case NeighborSearchType::Scatter:
          Wij = kernel(etajMag, 1.0);
          gradWij = Hj*etaj.unitVector()*kernel.grad(etajMag, 1.0);
          break;
        default:
          VERIFY2(false, "Unhandled neighbor search type.");
        }
        const auto correction = a(fieldi, i) + b(fieldi, i)*rij.x();
        Wij *= correction;
        gradWij = correction*gradWij - Vector(b(fieldi, i)*Wij);
        const auto weightj = weight(fieldj, j);
        normalization += weightj*Wij;
        result(fieldi, i) += weightj*(fieldList(fieldj, j) - smoothFieldList(fieldi, i))*gradWij;
      });
      CHECK(normalization > 0.0);
      result(fieldi, i) /= normalization;
    }
  }
  return result;
}

}
