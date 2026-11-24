#include "FieldListSecondDerivatives.hh"
#include "FieldOperations/FieldListConnectivity.hh"
#include "Field/FieldList.hh"
#include "Field/Field.hh"
#include "NodeList/NodeList.hh"
#include "Kernel/TableKernel.hh"
#include "Utilities/rotationMatrix.hh"
#include <cmath>
#include <vector>

namespace Spheral {

template<typename Dimension>
FieldList<Dimension, typename Dimension::Vector>
gradDivVectorFieldListPairWise
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
      std::vector<Tensor> grad2i(Dimension::nDim);
      Tensor normalization(0.0);
      connectivity.forEachNeighbor(fieldi, i, [&](const size_t fieldj, const size_t j) {
        const auto& Hj = Hfield(fieldj, j);
        const auto rij = ri - position(fieldj, j);
        const auto etai = Hi*rij;
        const auto etaj = Hj*rij;
        const auto Wi = kernel(etai.magnitude(), Hdeti);
        const auto Wj = kernel(etaj.magnitude(), Hj.Determinant());
        const auto gWi = Hi*etai.unitVector()*kernel.grad(etai.magnitude(), Hdeti);
        const auto gWj = Hj*etaj.unitVector()*kernel.grad(etaj.magnitude(), Hj.Determinant());
        Scalar Wij;
        Vector gWij;
        switch (searchType) {
        case NeighborSearchType::GatherScatter: Wij = 0.5*(Wi + Wj); gWij = 0.5*(gWi + gWj); break;
        case NeighborSearchType::Gather: Wij = Wi; gWij = gWi; break;
        case NeighborSearchType::Scatter: Wij = Wj; gWij = gWj; break;
        default: VERIFY2(false, "Unhandled neighbor search type.");
        }
        const auto weightj = weight(fieldj, j);
        const auto wij = weightj*Wij;
        const auto gwij = weightj*gWij;
        const auto frach = 0.01/Hi.Trace()/Dimension::nDim;
        const auto R = rotationMatrix(-rij.unitVector());
        const auto dxp = rij.magnitude()/(rij.magnitude2() + frach*frach);
        const auto dfdxpcol = -dxp*R*(fi - fieldList(fieldj, j));
        Tensor dfdxp;
        dfdxp.setColumn(0, dfdxpcol);
        dfdxp.rotationalTransform(R.Transpose());
        for (auto k = 0; k < Dimension::nDim; ++k) grad2i[k] += dfdxp*gwij(k);
        Tensor tweight(0.0);
        for (auto k = 0; k < Dimension::nDim; ++k) tweight(k, 0) = wij;
        tweight.rotationalTransform(R.Transpose());
        for (auto k = 0; k < Dimension::nDim; ++k) {
          for (auto m = 0; m < Dimension::nDim; ++m) normalization(k, m) += std::abs(tweight(k, m));
        }
      });
      for (auto k = 0; k < Dimension::nDim; ++k) {
        for (auto m = 0; m < Dimension::nDim; ++m) {
          for (auto n = 0; n < Dimension::nDim; ++n) grad2i[k](m, n) /= normalization(m, n) + 1.0e-10;
        }
        result(fieldi, i)(k) = grad2i[k].Trace();
      }
    }
  }
  return result;
}

}
