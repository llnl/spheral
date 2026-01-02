#include "computeRKSumVolume.hh"
#include "NodeList/NodeList.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Geometry/Dimension.hh"

using std::vector;

namespace Spheral {

//------------------------------------------------------------------------------
// Function to compute the per dimension volume multiplier.
//------------------------------------------------------------------------------
namespace {

constexpr double volumeElement(const Dim<1>&) {
  return 2.0;
}
  
constexpr double volumeElement(const Dim<2>&) {
  return M_PI;
}
  
constexpr double volumeElement(const Dim<3>&) {
  return 4.0/3.0*M_PI;
}

}

//------------------------------------------------------------------------------
// Compute the RK volume summation.
//------------------------------------------------------------------------------
template<typename Dimension>
void
computeRKSumVolume(const ConnectivityMap<Dimension>& connectivityMap,
                   const TableKernel<Dimension>& W,
                   const FieldList<Dimension, typename Dimension::Vector>& position,
                   const FieldList<Dimension, typename Dimension::Scalar>& /*mass*/,
                   const FieldList<Dimension, typename Dimension::SymTensor>& H,
                   FieldList<Dimension, typename Dimension::Scalar>& vol) {

  // Pre-conditions.
  const auto numNodeLists = vol.size();
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);

  // Get the maximum allowed volume in eta space.
  const auto etaVolMax = Dimension::pownu(0.5*W.kernelExtent()) * volumeElement(Dimension());

  // Zero it out.
  vol = 0.0;

  // Some useful variables.
  const auto W0 = W.kernelValue(0.0, 1.0);

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

      // Get the state for node i.
      const auto& ri = position(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();

      // And node j
      const auto& rj = position(nodeListj, j);
      const auto& Hj = H(nodeListj, j);
      const auto  Hdetj = Hj.Determinant();

      // Kernel weighting and gradient.
      const auto rij = ri - rj;
      const auto etai = (Hi*rij).magnitude();
      const auto etaj = (Hj*rij).magnitude();
      const auto Wi = W.kernelValue(etai, Hdeti);
      const auto Wj = W.kernelValue(etaj, Hdetj);

      // Sum the pair-wise contributions.
      vol_thread(nodeListi, i) += Wi;
      vol_thread(nodeListj, j) += Wj;
    }

#pragma omp critical
    {
      vol_thread.threadReduce();
    }
  }

  // The self contribution.
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto ni = vol[nodeListi]->numInternalElements();

#pragma omp parallel for
    for (auto i = 0u; i < ni; ++i) {
      const auto& Hi = H(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();
  
      // Add the self-contribution.
      vol(nodeListi, i) += Hdeti*W0;
      CHECK(vol(nodeListi, i) > 0.0);
      vol(nodeListi, i) = min(etaVolMax/Hdeti, 1.0/vol(nodeListi, i));
    }
  }
}

}
