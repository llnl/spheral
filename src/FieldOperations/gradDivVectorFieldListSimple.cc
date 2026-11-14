#include "FieldListSecondDerivatives.hh"
#include "FieldOperations/FieldListConnectivity.hh"
#include "Field/FieldList.hh"
#include "Field/Field.hh"
#include "NodeList/NodeList.hh"
#include "Kernel/TableKernel.hh"
#include "Utilities/rotationMatrix.hh"

namespace Spheral {

template<typename Dimension>
FieldList<Dimension, typename Dimension::Vector>
gradDivVectorFieldListSimple
(const FieldList<Dimension, typename Dimension::Vector>& fieldList,
 const FieldList<Dimension, typename Dimension::Vector>& position,
 const FieldList<Dimension, typename Dimension::Scalar>& weight,
 const FieldList<Dimension, typename Dimension::Scalar>&,
 const FieldList<Dimension, typename Dimension::Scalar>&,
 const FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
 const TableKernel<Dimension>& kernel) {
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;

  FieldList<Dimension, Vector> result;
  result.copyFields();
  for (auto fieldi = 0u; fieldi < fieldList.numFields(); ++fieldi) {
    result.appendField(Field<Dimension, Vector>("grad div " + fieldList[fieldi]->name(), fieldList[fieldi]->nodeList()));
  }
  FieldListConnectivity<Dimension, Vector> connectivity(fieldList, kernel.kernelExtent());
  const auto searchType = connectivity.searchType();
  for (auto fieldi = 0u; fieldi < fieldList.numFields(); ++fieldi) {
    const auto numInternal = fieldList[fieldi]->nodeListPtr()->numInternalNodes();
    for (auto i = 0u; i < numInternal; ++i) {
      const auto& ri = position(fieldi, i);
      const auto& Hi = Hfield(fieldi, i);
      const auto Hdeti = Hi.Determinant();
      const auto& fi = fieldList(fieldi, i);
      Scalar Wnormalization = weight(fieldi, i)*kernel(0.0, Hdeti);
      connectivity.forEachNeighbor(fieldi, i, [&](const size_t fieldj, const size_t j) {
        const auto& Hj = Hfield(fieldj, j);
        const auto Hdetj = Hj.Determinant();
        const auto rji = position(fieldj, j) - ri;
        const auto etai = Hi*rji;
        const auto etaj = Hj*rji;
        const auto Wi = kernel(etai.magnitude(), Hdeti);
        const auto Wj = kernel(etaj.magnitude(), Hdetj);
        const auto gWi = Hi*etai.unitVector()*kernel.grad(etai.magnitude(), Hdeti);
        const auto gWj = Hj*etaj.unitVector()*kernel.grad(etaj.magnitude(), Hdetj);
        Scalar Wij;
        Vector gWij;
        switch (searchType) {
        case NeighborSearchType::GatherScatter: Wij = 0.5*(Wi + Wj); gWij = 0.5*(gWi + gWj); break;
        case NeighborSearchType::Gather: Wij = Wi; gWij = gWi; break;
        case NeighborSearchType::Scatter: Wij = Wj; gWij = gWj; break;
        default: VERIFY2(false, "Unhandled neighbor search type.");
        }
        const auto weightj = weight(fieldj, j);
        Wnormalization += weightj*Wij;
        const auto frach = 0.01/Hi.Trace()/Dimension::nDim;
        const auto R = rotationMatrix(rji.unitVector());
        const auto dxp = rji.magnitude()/(rji.magnitude2() + frach*frach);
        const auto dfdxcol = (R*(fieldList(fieldj, j) - fi))*dxp;
        const Tensor d2fdx2 = dfdxcol*(R*gWij);
        result(fieldi, i) = R*result(fieldi, i);
        result(fieldi, i) -= weightj*d2fdx2.diagonalElements();
        result(fieldi, i) = R.Transpose()*result(fieldi, i);
      });
      CHECK(Wnormalization > 0.0);
      result(fieldi, i) *= 2.0/Wnormalization;
    }
  }
  return result;
}

}
