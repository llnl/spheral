//------------------------------------------------------------------------------
// Compute the RK volume summation.
//------------------------------------------------------------------------------
#include "computeOccupancyVolume.hh"
#include "NodeList/NodeList.hh"
#include "Hydro/HydroFieldNames.hh"

namespace Spheral {

template<typename Dimension>
void
computeOccupancyVolume(const ConnectivityMap<Dimension>& connectivityMap,
                       const TableKernel<Dimension>& W,
                       const FieldList<Dimension, typename Dimension::Vector>& position,
                       const FieldList<Dimension, typename Dimension::SymTensor>& H,
                       FieldList<Dimension, typename Dimension::Scalar>& vol) {

  // Pre-conditions.
  const auto numNodeLists = vol.size();
  CONTRACT_VAR(position);
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);

  // Zero it out.
  vol = 1.0;

  // Extent of the kernel.
  const auto kernelExtent = W.kernelExtent();
  const auto volFactor = 2.0*kernelExtent;

  // The set of interacting node pairs.
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();

#pragma omp parallel
  {
    // Some scratch variables.
    int i, j, nodeListi, nodeListj;
    auto vol_thread = vol.threadCopy();

#pragma omp for
    for (auto k = 0u; k < npairs; ++k) {
      i = pairs[k].i_node;
      j = pairs[k].j_node;
      nodeListi = pairs[k].i_list;
      nodeListj = pairs[k].j_list;

      vol_thread(nodeListi, i) += 1;
      vol_thread(nodeListj, j) += 1;
    }

#pragma omp critical
    {
      vol_thread.threadReduce();
    } // OMP critical
  }   // OMP parallel

  // Finish the volume
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto ni = vol[nodeListi]->numInternalElements();
#pragma omp parallel for
    for (auto i = 0u; i < ni; ++i) {
      const auto& Hi = H(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();
      vol(nodeListi, i) = volFactor/(Hdeti*vol(nodeListi, i));
    }
  }
}

}

