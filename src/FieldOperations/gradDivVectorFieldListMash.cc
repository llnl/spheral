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
gradDivVectorFieldListMash
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
  constexpr Scalar epsilon2 = 1.0e-4;
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
      const auto& fi = fieldList(fieldi, i);
      CHECK(Hi.Determinant() > 0.0 and mass(fieldi, i) > 0.0 and rho(fieldi, i) > 0.0);
      Scalar normalization = 1.0e-20;
      std::vector<Tensor> grad2elements(Dimension::nDim);
      Tensor gradF;
      Vector gradW;
      connectivity.forEachNeighbor(fieldi, i, [&](const size_t fieldj, const size_t j) {
        const auto& Hj = Hfield(fieldj, j);
        const auto mj = mass(fieldj, j);
        CHECK(Hj.Determinant() > 0.0 and mj > 0.0 and rho(fieldj, j) > 0.0);
        const auto rji = position(fieldj, j) - ri;
        const auto etai = Hi*rji;
        const auto etaj = Hj*rji;
        if (etai.magnitude2() < cutoff2 or etaj.magnitude2() < cutoff2) {
          const auto Hetai = Hi*etai.unitVector();
          const auto Hetaj = Hj*etaj.unitVector();
          auto getai = kernel.grad(etai.magnitude(), Hi.Determinant());
          auto getaj = kernel.grad(etaj.magnitude(), Hj.Determinant());
          auto g2etai = kernel.grad2(etai.magnitude(), Hi.Determinant());
          auto g2etaj = kernel.grad2(etaj.magnitude(), Hj.Determinant());
          auto Wi = kernel(etai.magnitude(), Hi.Determinant());
          auto Wj = kernel(etaj.magnitude(), Hj.Determinant());
          const auto thpti = 1.0/(etai.magnitude2() + epsilon2);
          g2etai = g2etai*thpti - 4.0*etai.magnitude()*thpti*thpti*getai +
                    (6.0*etai.magnitude2() - epsilon2)*thpti*thpti*thpti*Wi;
          getai = getai*thpti - 2.0*etai.magnitude()*thpti*thpti*Wi;
          Wi *= thpti;
          const auto thptj = 1.0/(etaj.magnitude2() + epsilon2);
          g2etaj = g2etaj*thptj - 4.0*etaj.magnitude()*thptj*thptj*getaj +
                    (6.0*etaj.magnitude2() - epsilon2)*thptj*thptj*thptj*Wj;
          getaj = getaj*thptj - 2.0*etaj.magnitude()*thptj*thptj*Wj;
          Wj *= thptj;
          const auto gWi = Hetai*getai;
          const auto gWj = Hetaj*getaj;
          const auto acki = Hetai.dyad(Hetai);
          const auto ackj = Hetaj.dyad(Hetaj);
          const auto g2Wi = acki*g2etai + ((Hi*Hi).Transpose() - acki)/(etai.magnitude() + 1.0e-30)*getai;
          const auto g2Wj = ackj*g2etaj + ((Hj*Hj).Transpose() - ackj)/(etaj.magnitude() + 1.0e-30)*getaj;
          Scalar Wij;
          Vector gWij;
          Tensor g2Wij;
          switch (searchType) {
          case NeighborSearchType::GatherScatter: Wij = 0.5*(Wi + Wj); gWij = 0.5*(gWi + gWj); g2Wij = 0.5*(g2Wi + g2Wj); break;
          case NeighborSearchType::Gather: Wij = Wi; gWij = gWi; g2Wij = g2Wi; break;
          case NeighborSearchType::Scatter: Wij = Wj; gWij = gWj; g2Wij = g2Wj; break;
          default: VERIFY2(false, "Unhandled neighbor search type.");
          }
          const auto fji = fieldList(fieldj, j) - fi;
          normalization += mj*Wij;
          for (auto gamma = 0; gamma < Dimension::nDim; ++gamma) {
            for (auto beta = 0; beta < Dimension::nDim; ++beta) {
              for (auto alpha = 0; alpha < Dimension::nDim; ++alpha) {
                grad2elements[gamma](beta, alpha) += mj*fji(gamma)*g2Wij(beta, alpha);
              }
              gradF(gamma, beta) += mj*fji(gamma)*gWij(beta);
            }
            gradW(gamma) += mj*gWij(gamma);
          }
        }
      });
      CHECK(normalization > 0.0);
      gradF /= normalization;
      gradW /= normalization;
      for (auto beta = 0; beta < Dimension::nDim; ++beta) {
        for (auto alpha = 0; alpha < Dimension::nDim; ++alpha) {
          result(fieldi, i)(beta) += grad2elements[alpha](beta, alpha)/normalization -
                                     gradF(alpha, alpha)*gradW(beta);
        }
      }
    }
  }
  return result;
}

}
