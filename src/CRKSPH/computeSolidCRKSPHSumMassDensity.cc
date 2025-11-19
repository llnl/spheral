//------------------------------------------------------------------------------
// Compute the CRKSPH mass density summation.
//------------------------------------------------------------------------------
#include "computeSolidCRKSPHSumMassDensity.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"
#include "NodeList/NodeList.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Utilities/NodeCoupling.hh"

using std::vector;
using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {

template<typename Dimension>
void
computeSolidCRKSPHSumMassDensity(const ConnectivityMap<Dimension>& connectivityMap,
                                 const TableKernel<Dimension>& W,
                                 const FieldList<Dimension, typename Dimension::Vector>& position,
                                 const FieldList<Dimension, typename Dimension::Scalar>& mass,
                                 const FieldList<Dimension, typename Dimension::SymTensor>& H,
                                 const FieldList<Dimension, typename Dimension::Scalar>& massDensity0,
                                 const NodeCoupling& nodeCoupling,
                                 FieldList<Dimension, typename Dimension::Scalar>& massDensity) {

  // Pre-conditions.
  const size_t numNodeLists = massDensity.size();
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(mass.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;

  const auto W0 = W.kernelValue(0.0, 1.0);

  // Prepare to sum the correction.
  FieldList<Dimension, Scalar> m0(FieldStorageType::CopyFields);
  for (auto nodeListi = 0u; nodeListi != numNodeLists; ++nodeListi) {
    m0.appendNewField("zeroth correction", position[nodeListi]->nodeList(), 0.0);
  }
  massDensity = 0.0;

  // The set of interacting node pairs.
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();

  // Walk the FluidNodeLists and sum the new mass density.
#pragma omp parallel
  {
    // Some scratch variables.
    int i, j, nodeListi, nodeListj;
    Vector rij, etai, etaj;

    typename SpheralThreads<Dimension>::FieldListStack threadStack;
    auto massDensity_thread = massDensity.threadCopy(threadStack);
    auto m0_thread = m0.threadCopy(threadStack);

#pragma omp for
    for (auto k = 0u; k < npairs; ++k) {
      i = pairs[k].i_node;
      j = pairs[k].j_node;
      nodeListi = pairs[k].i_list;
      nodeListj = pairs[k].j_list;

      // Check the coupling of these points.
      const auto fij = nodeCoupling(pairs[k]);
      if (fij > 0.0) {

        // Get the state for node i.
        const auto& ri = position(nodeListi, i);
        const auto  mi = mass(nodeListi, i);
        const auto& Hi = H(nodeListi, i);
        const auto  Hdeti = Hi.Determinant();
        const auto  rho0i = massDensity0(nodeListi, i);
        const auto  wi = mi/rho0i;

        // State for node j.
        const auto& rj = position(nodeListj, j);
        const auto  mj = mass(nodeListj, j);
        const auto& Hj = H(nodeListj, j);
        const auto  Hdetj = Hj.Determinant();
        const auto  rho0j = massDensity0(nodeListj, j);
        const auto  wj = mj/rho0j;

        // Kernel weighting and gradient.
        const auto rij = ri - rj;
        const auto etai = (Hi*rij).magnitude();
        const auto etaj = (Hj*rij).magnitude();
        const auto Wi = W.kernelValue(etai, Hdeti);
        const auto Wj = W.kernelValue(etaj, Hdetj);

        // Sum the pair-wise contributions.
        massDensity_thread(nodeListi, i) += fij*wj*Wj*mi;
        massDensity_thread(nodeListj, j) += fij*wi*Wi*mj;
        m0_thread(nodeListi, i) += fij*wj*Wj*wi;
        m0_thread(nodeListj, j) += fij*wi*Wi*wj;
      }
    }
      
    // Reduce the thread values to the master.
    threadReduceFieldLists<Dimension>(threadStack);

  } // OMP parallel
  
  // Finish for each point.
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto& nodeList = dynamic_cast<const FluidNodeList<Dimension>&>(massDensity[nodeListi]->nodeList());
    const auto ni = nodeList.numInternalNodes();
    const auto rhoMin = nodeList.rhoMin();
    const auto rhoMax = nodeList.rhoMax();

#pragma omp parallel for
    for (auto i = 0u; i < ni; ++i) {

      // Get the state for node i.
      const auto  mi = mass(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();
      const auto  rho0i = massDensity0(nodeListi, i);
      const auto  wi = mi/rho0i;
      massDensity(nodeListi, i) = max(rhoMin, 
                                      min(rhoMax,
                                          (massDensity(nodeListi, i) + wi*Hdeti*W0*mi)/
                                          (m0(nodeListi, i) + wi*Hdeti*W0*wi)));
      CHECK(massDensity(nodeListi, i) > 0.0);
    }
  }
}

}

