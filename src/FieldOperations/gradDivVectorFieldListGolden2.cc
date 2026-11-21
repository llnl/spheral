#include "FieldListSecondDerivatives.hh"
#include "FieldOperations/FieldListConnectivity.hh"
#include "Field/FieldList.hh"
#include "Field/Field.hh"
#include "NodeList/NodeList.hh"
#include "Kernel/TableKernel.hh"
#include "Utilities/rotationMatrix.hh"
#include <vector>

namespace Spheral {

template<typename Dimension>
FieldList<Dimension, typename Dimension::Vector>
gradDivVectorFieldListGolden2
(const FieldList<Dimension, typename Dimension::Vector>& fieldList,
 const FieldList<Dimension, typename Dimension::Vector>& position,
 const FieldList<Dimension, typename Dimension::Scalar>&,
 const FieldList<Dimension, typename Dimension::Scalar>& mass,
 const FieldList<Dimension, typename Dimension::Scalar>& rho,
 const FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
 const TableKernel<Dimension>& kernel) {
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
      const auto rhoi = rho(fieldi, i);
      const auto& fi = fieldList(fieldi, i);
      std::vector<Tensor> grad2elements(Dimension::nDim);
      Tensor gradF;
      Vector gradRho;
      connectivity.forEachNeighbor(fieldi, i, [&](const size_t fieldj, const size_t j) {
        const auto& Hj = Hfield(fieldj, j);
        const auto rji = position(fieldj, j) - ri;
        const auto etai = Hi*rji;
        const auto etaj = Hj*rji;
        const auto gWi = Hi*etai.unitVector()*kernel.grad(etai.magnitude(), Hi.Determinant());
        const auto gWj = Hj*etaj.unitVector()*kernel.grad(etaj.magnitude(), Hj.Determinant());
        Vector gWij;
        switch (searchType) {
        case NeighborSearchType::GatherScatter: gWij = 0.5*(gWi + gWj); break;
        case NeighborSearchType::Gather: gWij = gWi; break;
        case NeighborSearchType::Scatter: gWij = gWj; break;
        default: VERIFY2(false, "Unhandled neighbor search type.");
        }
        const auto mj = mass(fieldj, j);
        const auto rhoj = rho(fieldj, j);
        const auto fji = fieldList(fieldj, j) - fi;
        const auto frach = 0.01/Hi.Trace()/Dimension::nDim;
        const auto R = rotationMatrix(rji.unitVector());
        const auto dxp = rji.magnitude()/(rji.magnitude2() + frach*frach);
        Tensor dfdxp;
        dfdxp.setColumn(0, R*fji*dxp);
        const auto dfdx = R.Transpose()*dfdxp*R;
        for (auto gamma = 0; gamma < Dimension::nDim; ++gamma) {
          for (auto beta = 0; beta < Dimension::nDim; ++beta) {
            for (auto alpha = 0; alpha < Dimension::nDim; ++alpha) {
              grad2elements[gamma](beta, alpha) -= 2.0*mj*dfdx(gamma, alpha)*gWij(beta);
            }
            gradF(gamma, beta) += mj*fji(gamma)*gWij(beta);
          }
          gradRho(gamma) += mj*(rhoj - rhoi)*gWij(gamma);
        }
      });
      CHECK(rhoi > 0.0);
      gradF /= rhoi;
      gradRho /= rhoi;
      for (auto beta = 0; beta < Dimension::nDim; ++beta) {
        for (auto alpha = 0; alpha < Dimension::nDim; ++alpha) {
          result(fieldi, i)(beta) += grad2elements[alpha](beta, alpha) -
                                     gradRho(beta)*gradF(alpha, alpha) -
                                     gradRho(alpha)*gradF(alpha, beta);
        }
      }
      result(fieldi, i) /= rhoi;
    }
  }
  return result;
}

}
