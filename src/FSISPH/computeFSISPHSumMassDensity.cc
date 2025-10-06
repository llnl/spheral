#include "FSISPH/computeFSISPHSumMassDensity.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"
#include "NodeList/NodeList.hh"

namespace Spheral{

template<typename Dimension>
void
computeFSISPHSumMassDensity(const ConnectivityMap<Dimension>& connectivityMap,
                            const TableKernel<Dimension>& W,
                            const std::vector<int>& sumDensityNodeLists,
                            const FieldList<Dimension, typename Dimension::Vector>& position,
                            const FieldList<Dimension, typename Dimension::Scalar>& mass,
                            const FieldList<Dimension, typename Dimension::SymTensor>& H,
                            const bool consistentSum,
                            FieldList<Dimension, typename Dimension::Scalar>& massDensity) {

  // Pre-conditions.
  const auto numNodeLists = massDensity.size();
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(mass.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);

  // Some useful variables.
  const auto W0 = W.kernelValue(0.0, 1.0);

  // The set of interacting node pairs.
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();

  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto n = massDensity[nodeListi]->numInternalElements();
    if (sumDensityNodeLists[nodeListi]==1){
#pragma omp parallel for
      for (auto i = 0u; i < n; ++i) {
        const auto  mi = mass(nodeListi, i);
        const auto& Hi = H(nodeListi, i);
        const auto  Hdeti = Hi.Determinant();
        massDensity(nodeListi,i) =  mi*Hdeti*W0;
      }
    }
  }


  // Now the pair contributions.
#pragma omp parallel
  {
    int i, j, nodeListi, nodeListj;
    auto massDensity_thread = massDensity.threadCopy();

#pragma omp for
    for (auto k = 0u; k < npairs; ++k) {
      i = pairs[k].i_node;
      j = pairs[k].j_node;
      nodeListi = pairs[k].i_list;
      nodeListj = pairs[k].j_list;

      const auto  mi = mass(nodeListi, i);
      const auto  mj = mass(nodeListj, j);
      const auto& ri = position(nodeListi, i);
      const auto& rj = position(nodeListj, j);
      const auto  rij = ri - rj;
      const auto normalSum = (nodeListi == nodeListj) || consistentSum;

      if(sumDensityNodeLists[nodeListi]==1){
        const auto& Hi = H(nodeListi, i);
        const auto  Hdeti = Hi.Determinant();
        const auto  etai = (Hi*rij).magnitude();
        const auto  Wi = W.kernelValue(etai, Hdeti);

        massDensity_thread(nodeListi, i) += (normalSum ? mj : mi)*Wi; 
      } 
      
      if(sumDensityNodeLists[nodeListj]==1){
        const auto& Hj = H(nodeListj, j);
        const auto  Hdetj = Hj.Determinant();
        const auto  etaj = (Hj*rij).magnitude();
        const auto  Wj = W.kernelValue(etaj, Hdetj);
    
        massDensity_thread(nodeListj, j) += (normalSum ? mi : mj)*Wj;
      }
    }

#pragma omp critical
    {
      massDensity_thread.threadReduce();
    }
  }


} // function


} //spheral namespace
