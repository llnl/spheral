#include "FieldListSecondDerivatives.hh"
#include "FieldOperations/FieldListConnectivity.hh"
#include "Field/FieldList.hh"
#include "Field/Field.hh"
#include "NodeList/NodeList.hh"
#include "Kernel/TableKernel.hh"
#include <vector>

namespace Spheral {

template<typename Dimension>
FieldList<Dimension, typename Dimension::Vector>
gradDivVectorFieldListGolden
(const FieldList<Dimension, typename Dimension::Vector>& fieldList,
 const FieldList<Dimension, typename Dimension::Vector>& position,
 const FieldList<Dimension, typename Dimension::Scalar>&,
 const FieldList<Dimension, typename Dimension::Scalar>& mass,
 const FieldList<Dimension, typename Dimension::Scalar>& rho,
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
  const Scalar cutoff2 = kernel.kernelExtent()*kernel.kernelExtent();
  FieldListConnectivity<Dimension, Vector> connectivity(fieldList, kernel.kernelExtent());
  const auto searchType = connectivity.searchType();
  for (auto fieldi = 0u; fieldi < fieldList.numFields(); ++fieldi) {
    const auto numInternal = fieldList[fieldi]->nodeListPtr()->numInternalNodes();
    for (auto i = 0u; i < numInternal; ++i) {
      const auto& ri = position(fieldi, i);
      const auto& Hi = Hfield(fieldi, i);
      const auto rhoi = rho(fieldi, i);
      const auto& fi = fieldList(fieldi, i);
      CHECK(Hi.Determinant() > 0.0 and mass(fieldi, i) > 0.0 and rhoi > 0.0);
      std::vector<Tensor> grad2elements(Dimension::nDim);
      Vector gradRho;
      Tensor gradF;
      connectivity.forEachNeighbor(fieldi, i, [&](const size_t fieldj, const size_t j) {
        const auto& Hj = Hfield(fieldj, j);
        const auto mj = mass(fieldj, j);
        const auto rhoj = rho(fieldj, j);
        const auto rij = ri - position(fieldj, j);
        const auto etai = Hi*rij;
        const auto etaj = Hj*rij;
        if (etai.magnitude2() < cutoff2 or etaj.magnitude2() < cutoff2) {
          const auto Hetai = Hi*etai.unitVector();
          const auto Hetaj = Hj*etaj.unitVector();
          const auto getai = kernel.grad(etai.magnitude(), Hi.Determinant());
          const auto getaj = kernel.grad(etaj.magnitude(), Hj.Determinant());
          const auto gWi = Hetai*getai;
          const auto gWj = Hetaj*getaj;
          const auto acki = Hetai.dyad(Hetai);
          const auto ackj = Hetaj.dyad(Hetaj);
          const auto g2Wi = acki*kernel.grad2(etai.magnitude(), Hi.Determinant()) +
                            ((Hi*Hi).Transpose() - acki)/(etai.magnitude() + 1.0e-30)*getai;
          const auto g2Wj = ackj*kernel.grad2(etaj.magnitude(), Hj.Determinant()) +
                            ((Hj*Hj).Transpose() - ackj)/(etaj.magnitude() + 1.0e-30)*getaj;
          Vector gWij;
          Tensor g2Wij;
          switch (searchType) {
          case NeighborSearchType::GatherScatter: gWij = 0.5*(gWi + gWj); g2Wij = 0.5*(g2Wi + g2Wj); break;
          case NeighborSearchType::Gather: gWij = gWi; g2Wij = g2Wi; break;
          case NeighborSearchType::Scatter: gWij = gWj; g2Wij = g2Wj; break;
          default: VERIFY2(false, "Unhandled neighbor search type.");
          }
          const auto fij = fi - fieldList(fieldj, j);
          for (auto gamma = 0; gamma < Dimension::nDim; ++gamma) {
            for (auto beta = 0; beta < Dimension::nDim; ++beta) {
              for (auto alpha = 0; alpha < Dimension::nDim; ++alpha) {
                grad2elements[gamma](beta, alpha) -= mj*fij(gamma)*g2Wij(beta, alpha);
              }
              gradF(gamma, beta) -= mj*fij(gamma)*gWij(beta);
            }
            gradRho(gamma) -= mj*(rhoi - rhoj)*gWij(gamma);
          }
        }
      });
      gradRho /= rhoi;
      gradF /= rhoi;
      for (auto beta = 0; beta < Dimension::nDim; ++beta) {
        for (auto alpha = 0; alpha < Dimension::nDim; ++alpha) {
          result(fieldi, i)(beta) += grad2elements[alpha](beta, alpha) -
                                     gradRho(alpha)*gradF(alpha, beta) -
                                     gradRho(beta)*gradF(alpha, alpha);
        }
      }
      result(fieldi, i) /= rhoi;
    }
  }
  return result;
}

}
