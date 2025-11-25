//------------------------------------------------------------------------------
// Compute the SPH mass density summation.
//------------------------------------------------------------------------------
#include "computeSPHSumMassDensity.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"
#include "NodeList/NodeList.hh"
#include "Hydro/HydroFieldNames.hh"

using std::vector;
using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {

template<typename Dimension, typename KernelType>
void
correctSPHSumMassDensity(const ConnectivityMap<Dimension>& connectivityMap,
                         const KernelType& W,
                         const bool sumOverAllNodeLists,
                         const FieldList<Dimension, typename Dimension::Vector>& position,
                         const FieldList<Dimension, typename Dimension::Scalar>& mass,
                         const FieldList<Dimension, typename Dimension::SymTensor>& H,
                         FieldList<Dimension, typename Dimension::Scalar>& massDensity) {

  // Pre-conditions.
  const size_t numNodeLists = massDensity.size();
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(mass.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);

  typedef typename Dimension::Scalar Scalar;

  // The set of interacting node pairs (includes self interaction)
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();

  // Prepare the kernel sum correction field.
  FieldList<Dimension, Scalar> sumUnity(FieldStorageType::CopyFields);
  for (auto nodeListi = 0u; nodeListi != numNodeLists; ++nodeListi) {
    sumUnity.appendNewField("SPH sum unity check", massDensity[nodeListi]->nodeList(), 0.0);
  }

  // Now the pair contributions.
#pragma omp parallel
  {
    int i, j, nodeListi, nodeListj;
    auto sumUnity_thread = sumUnity.threadCopy();

#pragma omp for
    for (auto k = 0u; k < npairs; ++k) {
      i = pairs[k].i_node;
      j = pairs[k].j_node;
      nodeListi = pairs[k].i_list;
      nodeListj = pairs[k].j_list;

      // State for node i
      const auto& ri = position(nodeListi, i);
      const auto  mi = mass(nodeListi, i);
      const auto  rhoi = massDensity(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();
      CHECK(rhoi > 0.0);

      // State for node j
      const auto& rj = position(nodeListj, j);
      const auto  mj = mass(nodeListj, j);
      const auto  rhoj = massDensity(nodeListj, j);
      const auto& Hj = H(nodeListj, j);
      const auto  Hdetj = Hj.Determinant();
      CHECK(rhoj > 0.0);

      // Kernel weighting and gradient.
      const auto Wi = W(Hi*rj, Hi*ri, Hdeti);
      const auto Wj = W(Hj*rj, Hj*ri, Hdetj);

      // Sum the pair-wise contributions.
      sumUnity_thread(nodeListi, i) += mj/rhoj*Wj;
      sumUnity_thread(nodeListj, j) += mi/rhoi*Wi;
    }

#pragma omp critical
    {
      sumUnity_thread.threadReduce();
    } // OMP critical
  }   // OMP parallel

  // Apply the correction.
  massDensity /= sumUnity;
}

}
