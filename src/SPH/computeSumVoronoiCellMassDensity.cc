//------------------------------------------------------------------------------
// Compute the Voronoi cell mass density summation.
//------------------------------------------------------------------------------
#include "computeSumVoronoiCellMassDensity.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"
#include "NodeList/NodeList.hh"
#include "NodeList/FluidNodeList.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Utilities/safeInv.hh"

using std::vector;
using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {

template<typename Dimension>
void
computeSumVoronoiCellMassDensity(const ConnectivityMap<Dimension>& connectivityMap,
                                 const TableKernel<Dimension>& W,
                                 const FieldList<Dimension, typename Dimension::Vector>& position,
                                 const FieldList<Dimension, typename Dimension::Scalar>& mass,
                                 const FieldList<Dimension, typename Dimension::Scalar>& volume,
                                 const FieldList<Dimension, typename Dimension::SymTensor>& H,
                                 FieldList<Dimension, typename Dimension::Scalar>& massDensity) {

  // Pre-conditions.
  const size_t numNodeLists = massDensity.size();
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(mass.size() == numNodeLists);
  REQUIRE(volume.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);

  typedef typename Dimension::Scalar Scalar;

  // Some useful variables.
  const auto W0 = W.kernelValue(0.0, 1.0);
  // Zero out the result, and prepare a FieldList to hold the effective volume.
  massDensity = 0.0;
  FieldList<Dimension, Scalar> Veff(FieldStorageType::CopyFields);
  for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const NodeList<Dimension>& nodeList = massDensity[nodeListi]->nodeList();
    Veff.appendNewField("effective volume", nodeList, 0.0);
  }

  // The set of interacting node pairs.
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();

  // Walk all the interacting pairs.
#pragma omp parallel
  {
    // Thread private scratch variables
    int i, j, nodeListi;

    typename SpheralThreads<Dimension>::FieldListStack threadStack;
    auto Veff_thread = Veff.threadCopy(threadStack);
    auto massDensity_thread = massDensity.threadCopy(threadStack);

#pragma omp for
    for (auto kk = 0u; kk < npairs; ++kk) {
      i = pairs[kk].i_node;
      j = pairs[kk].j_node;
      nodeListi = pairs[kk].i_list;
      //nodeListj = pairs[kk].j_list;

      // State for node i
      const auto& ri = position(nodeListi, i);
      const auto  mi = mass(nodeListi, i);
      const auto  Vi = volume(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();

      // State for node j
      const auto& rj = position(nodeListi, j);
      const auto  mj = mass(nodeListi, j);
      const auto  Vj = volume(nodeListi, j);
      const auto& Hj = H(nodeListi, j);
      const auto  Hdetj = Hj.Determinant();

      // Kernel weighting and gradient.
      const auto rij = ri - rj;
      const auto etai = (Hi*rij).magnitude();
      const auto etaj = (Hj*rij).magnitude();
      const auto Wi = W.kernelValue(etai, Hdeti);
      const auto Wj = W.kernelValue(etaj, Hdetj);

      // Sum the pair-wise contributions.
      Veff_thread(nodeListi, i) += Vj*Wi;
      massDensity_thread(nodeListi, i) += mj*Wi;

      Veff_thread(nodeListi, j) += Vi*Wj;
      massDensity_thread(nodeListi, j) += mi*Wj;
    }

    // Reduce the thread values to the master.
    threadReduceFieldLists<Dimension>(threadStack);

  }   // OMP parallel


  // Finalize the density for each point.
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto& nodeList = dynamic_cast<const FluidNodeList<Dimension>&>(massDensity[0]->nodeList());
    const auto ni = nodeList.numInternalNodes();
    const auto rhoMin = nodeList.rhoMin();
    const auto rhoMax = nodeList.rhoMax();

#pragma omp parallel for
    for (auto i = 0u; i < ni; ++i) {
      const auto  mi = mass(nodeListi, i);
      const auto  Vi = volume(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();
      massDensity(nodeListi, i) = max(rhoMin, 
                                      min(rhoMax,
                                          (massDensity(nodeListi, i) + mi*Hdeti*W0) * 
                                          safeInv(Veff(nodeListi, i) + Vi*Hdeti*W0)));
      CHECK(massDensity(nodeListi, i) > 0.0);
    }
  }
}

}
