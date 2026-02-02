//------------------------------------------------------------------------------
// Compute the SPH mass density summation.
//------------------------------------------------------------------------------

#include "computeSPHSumMassDensity.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "NodeList/NodeList.hh"
#include "Hydro/HydroFieldNames.hh"

namespace Spheral {

template<typename Dimension, typename KernelType>
void
computeSPHSumMassDensity(const ConnectivityMap<Dimension>& connectivityMap,
                         const KernelType& W,
                         const bool sumOverAllNodeLists,
                         const FieldList<Dimension, typename Dimension::Vector>& position,
                         const FieldList<Dimension, typename Dimension::Scalar>& mass,
                         const FieldList<Dimension, typename Dimension::SymTensor>& H,
                         FieldList<Dimension, typename Dimension::Scalar>& massDensity) {

  // Pre-conditions.
  const auto numNodeLists = massDensity.size();
  CONTRACT_VAR(numNodeLists);
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(mass.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);

  // The set of interacting node pairs.
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();

  // First the self contribution.
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto n = massDensity[nodeListi]->numInternalElements();
#pragma omp parallel for
    for (auto i = 0u; i < n; ++i) {
      const auto& posi = position(nodeListi, i);
      const auto  mi = mass(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();
      const auto  etaii = Hi*posi;
      massDensity(nodeListi, i) = mi*W(etaii, etaii, Hdeti);
    }
  }

  // Walk the pairs and sum their contributions.
#pragma omp parallel
  {
    size_t i, j, nodeListi, nodeListj;
    auto massDensity_thread = massDensity.threadCopy();

#pragma omp for
    for (auto k = 0u; k < npairs; ++k) {
      i = pairs[k].i_node;
      j = pairs[k].j_node;
      nodeListi = pairs[k].i_list;
      nodeListj = pairs[k].j_list;

      // State for node i
      const auto& ri = position(nodeListi, i);
      const auto  mi = mass(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();

      // State for node j
      const auto& rj = position(nodeListj, j);
      const auto  mj = mass(nodeListj, j);
      const auto& Hj = H(nodeListj, j);
      const auto  Hdetj = Hj.Determinant();

      // Kernel weighting and gradient.
      const auto Wi = W(Hi*rj, Hi*ri, Hdeti);
      const auto Wj = W(Hj*rj, Hj*ri, Hdetj);

      // Sum the pair-wise contributions.
      massDensity_thread(nodeListi, i) += (nodeListi == nodeListj ? mj : mi)*Wj;
      massDensity_thread(nodeListj, j) += (nodeListi == nodeListj ? mi : mj)*Wi;
    }

#pragma omp critical
    {
      massDensity_thread.threadReduce();
    }
  }
}

}
