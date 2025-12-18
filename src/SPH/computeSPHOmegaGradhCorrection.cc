//------------------------------------------------------------------------------
// Compute the SPH grad h correction due to Springel et al.
//------------------------------------------------------------------------------
#include "computeSPHOmegaGradhCorrection.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"
#include "NodeList/NodeList.hh"
#include "Hydro/HydroFieldNames.hh"

namespace Spheral {

template<typename Dimension>
void
computeSPHOmegaGradhCorrection(const ConnectivityMap<Dimension>& connectivityMap,
                               const TableKernel<Dimension>& W,
                               const FieldList<Dimension, typename Dimension::Vector>& position,
                               const FieldList<Dimension, typename Dimension::SymTensor>& H,
                               FieldList<Dimension, typename Dimension::Scalar>& omegaGradh) {

  // Pre-conditions.
  const size_t numNodeLists = omegaGradh.size();
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);

  typedef typename Dimension::Scalar Scalar;

  // Some useful variables.
  const auto W0 = W.kernelValue(0.0, 1.0);

  // Zero out the result.
  omegaGradh = 0.0;

  // Prepare a FieldList to hold the sum gradient.
  FieldList<Dimension, Scalar> gradsum(FieldStorageType::CopyFields);
  for (auto nodeListi = 0u; nodeListi != numNodeLists; ++nodeListi) {
    const auto& nodeList = omegaGradh[nodeListi]->nodeList();
    gradsum.appendNewField("sum of the gradient", nodeList, 0.0);
  }

  // The set of interacting node pairs.
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();

  // Walk all the interacting pairs.
#pragma omp parallel
  {
    // Thread private scratch variables
    int i, j, nodeListi, nodeListj;
    Scalar Wi, gWi, Wj, gWj;

    typename SpheralThreads<Dimension>::FieldListStack threadStack;
    auto omegaGradh_thread = omegaGradh.threadCopy(threadStack);
    auto gradsum_thread = gradsum.threadCopy(threadStack);

#pragma omp for
    for (auto kk = 0u; kk < npairs; ++kk) {
      i = pairs[kk].i_node;
      j = pairs[kk].j_node;
      nodeListi = pairs[kk].i_list;
      nodeListj = pairs[kk].j_list;

      const auto& ri = position(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();

      const auto& rj = position(nodeListj, j);
      const auto& Hj = H(nodeListj, j);
      const auto  Hdetj = Hj.Determinant();

      const auto rij = ri - rj;
      const auto etai = (Hi*rij).magnitude();
      const auto etaj = (Hj*rij).magnitude();

      // Kernel weighting and gradient.
      W.kernelAndGradValue(etai, Hdeti, Wi, gWi);
      W.kernelAndGradValue(etaj, Hdetj, Wj, gWj);

      // Sum the pair-wise contributions.
      omegaGradh_thread(nodeListi, i) += Wi;
      omegaGradh_thread(nodeListj, j) += Wj;

      gradsum_thread(nodeListi, i) += etai*gWi;
      gradsum_thread(nodeListj, j) += etaj*gWj;
    }

    // Reduce the thread values to the master.
    threadReduceFieldLists<Dimension>(threadStack);

  }    // OMP parallel
  
  // Finish up for each point.
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto ni = omegaGradh[nodeListi]->numInternalElements();
#pragma omp parallel for
    for (auto i = 0u; i < ni; ++i) {

      // If this point is isolated, we punt to unity.
      if (connectivityMap.numNeighborsForNode(nodeListi, i) == 0) {
        omegaGradh(nodeListi, i) = 1.0;

      } else {
        const Scalar Hdeti = H(nodeListi, i).Determinant();
        omegaGradh(nodeListi, i) += Hdeti*W0;
        CHECK(omegaGradh(nodeListi, i) > 0.0);
        omegaGradh(nodeListi, i) = std::max(1.0e-30, -gradsum(nodeListi, i)/(Dimension::nDim * omegaGradh(nodeListi, i)));

      }

      // Post-conditions.
      ENSURE2(omegaGradh(nodeListi, i) >= 0.0, 
              nodeListi << " " << i << " " << omegaGradh(nodeListi, i));
    }
  }
}

}
