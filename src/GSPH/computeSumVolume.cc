//---------------------------------Spheral++----------------------------------//
// Compute volume from inverse of the kernel summation.
//
//   Hopkins P.F. (2015) "A New Class of Accurate, Mesh-Free Hydrodynamic 
//   Simulation Methods," MNRAS, 450(1):53-110
//
// J.M. Pearl 2022
//----------------------------------------------------------------------------//

#include "GSPH/computeSumVolume.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"
#include "NodeList/NodeList.hh"

namespace Spheral{

template<typename Dimension>
void
computeSumVolume(const ConnectivityMap<Dimension>& connectivityMap,
                 const TableKernel<Dimension>& W,
                 const FieldList<Dimension, typename Dimension::Vector>& position,
                 const FieldList<Dimension, typename Dimension::SymTensor>& H,
                       FieldList<Dimension, typename Dimension::Scalar>& volume) {

  // Pre-conditions.
  const auto numNodeLists = volume.size();
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(volume.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);

  // start fresh
  volume.Zero();

  // Some useful variables.
  const auto W0 = W.kernelValue(0.0, 1.0);

  // The set of interacting node pairs.
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();

  // Now the pair contributions.
#pragma omp parallel
  {
    int i, j, nodeListi, nodeListj;
    auto volume_thread = volume.threadCopy();

#pragma omp for
    for (auto k = 0u; k < npairs; ++k) {
      i = pairs[k].i_node;
      j = pairs[k].j_node;
      nodeListi = pairs[k].i_list;
      nodeListj = pairs[k].j_list;

      const auto& ri = position(nodeListi, i);
      const auto& Hi = H(nodeListi, i);  

      const auto& rj = position(nodeListj, j);
      const auto& Hj = H(nodeListj, j);

      const auto  rij = ri - rj;

      const auto  Hdeti = Hi.Determinant();
      const auto  etai = (Hi*rij).magnitude();
      const auto  Wi = W.kernelValue(etai, Hdeti);

      const auto  Hdetj = Hj.Determinant();
      const auto  etaj = (Hj*rij).magnitude();
      const auto  Wj = W.kernelValue(etaj, Hdetj);

      volume_thread(nodeListi, i) += Wi; 
      volume_thread(nodeListj, j) += Wj;

    }

#pragma omp critical
    {
      volume_thread.threadReduce();
    }
  }    // omp critical region

  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto n = volume[nodeListi]->numInternalElements();
#pragma omp parallel for
    for (auto i = 0u; i < n; ++i) {
      const auto& Hi = H(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();
      volume(nodeListi,i) =  1.0/(volume(nodeListi,i)+Hdeti*W0);
    }
  }   
}     // function
}     // spheral namespace
