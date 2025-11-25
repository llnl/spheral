//---------------------------------Spheral++----------------------------------//
// initializes the pressure and velocity gradients for Riemann solver - based
// SPH varients
//
// J.M. Pearl 2023
//----------------------------------------------------------------------------//

#include "GSPH/initializeGradients.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"
#include "NodeList/NodeList.hh"

namespace Spheral{

template<typename Dimension>
void
initializeGradients(const ConnectivityMap<Dimension>& connectivityMap,
                    const TableKernel<Dimension>& W,
                    const FieldList<Dimension, typename Dimension::Vector>& position,
                    const FieldList<Dimension, typename Dimension::SymTensor>& H,
                    const FieldList<Dimension, typename Dimension::Scalar>& volume,
                    const FieldList<Dimension, typename Dimension::Scalar>& pressure,
                    const FieldList<Dimension, typename Dimension::Vector>& velocity,
                          FieldList<Dimension, typename Dimension::Tensor>& M,
                          FieldList<Dimension, typename Dimension::Vector>& DpDx,
                          FieldList<Dimension, typename Dimension::Tensor>& DvDx) {

  typedef typename Dimension::Tensor Tensor;

  const auto& nodeLists = connectivityMap.nodeLists();
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();
  const auto  numNodeLists = nodeLists.size();

  REQUIRE(volume.size() == numNodeLists);
  REQUIRE(velocity.size() == numNodeLists);
  REQUIRE(pressure.size() == numNodeLists);
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);

  REQUIRE(DpDx.size() == numNodeLists);
  REQUIRE(DvDx.size() == numNodeLists);
  REQUIRE(M.size() == numNodeLists);

#pragma omp parallel
  {
    // Thread private scratch variables
    int i, j, nodeListi, nodeListj;

    typename SpheralThreads<Dimension>::FieldListStack threadStack;
    auto M_thread = M.threadCopy(threadStack);
    auto DpDx_thread = DpDx.threadCopy(threadStack);
    auto DvDx_thread = DvDx.threadCopy(threadStack);

#pragma omp for
    for (auto kk = 0u; kk < npairs; ++kk) {
      i = pairs[kk].i_node;
      j = pairs[kk].j_node;
      nodeListi = pairs[kk].i_list;
      nodeListj = pairs[kk].j_list;
      
      // Get the state for node i.
      const auto& vi = velocity(nodeListi, i);
      const auto& Pi = pressure(nodeListi, i);
      const auto& ri = position(nodeListi, i);
      const auto& voli = volume(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();

      CHECK(voli > 0.0);
      CHECK(Hdeti > 0.0);

      auto& DpDxi = DpDx(nodeListi, i);
      auto& DvDxi = DvDx(nodeListi, i);
      auto& Mi = M(nodeListi, i);

      // Get the state for node j
      const auto& vj = velocity(nodeListj, j);
      const auto& Pj = pressure(nodeListj, j);
      const auto& rj = position(nodeListj, j);
      const auto& volj = volume(nodeListj, j);
      const auto& Hj = H(nodeListj, j);
      const auto  Hdetj = Hj.Determinant();

      CHECK(volj > 0.0);
      CHECK(Hdetj > 0.0);

      auto& DpDxj = DpDx(nodeListj, j);
      auto& DvDxj = DvDx(nodeListj, j);
      auto& Mj = M(nodeListj, j);

      const auto rij = ri - rj;
      const auto vij = vi - vj;
      const auto Pij = Pi - Pj;

      const auto etai = Hi*rij;
      const auto etaj = Hj*rij;
      const auto etaMagi = etai.magnitude();
      const auto etaMagj = etaj.magnitude();

      CHECK(etaMagi >= 0.0);
      CHECK(etaMagj >= 0.0);

      const auto gWi = W.gradValue(etaMagi, Hdeti);
      const auto Hetai = Hi*etai.unitVector();
      const auto gradWi = gWi*Hetai;

      const auto gWj = W.gradValue(etaMagj, Hdetj);
      const auto Hetaj = Hj*etaj.unitVector();
      const auto gradWj = gWj*Hetaj;

      const auto gradPsii = volj*gradWi;
      const auto gradPsij = voli*gradWj;

      Mi -= rij.dyad(gradPsii);
      Mj -= rij.dyad(gradPsij);
      
      DpDxi -= Pij*gradPsii;
      DpDxj -= Pij*gradPsij;

      DvDxi -= vij.dyad(gradPsii);
      DvDxj -= vij.dyad(gradPsij);

    } // loop over pairs
    // Reduce the thread values to the master.
    threadReduceFieldLists<Dimension>(threadStack);
  }   // OpenMP parallel region
  
  // Finish up the spatial gradient calculation
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto& nodeList = volume[nodeListi]->nodeList();
    const auto ni = nodeList.numInternalNodes();
#pragma omp parallel for
    for (auto i = 0u; i < ni; ++i) {
      const auto  numNeighborsi = connectivityMap.numNeighborsForNode(nodeListi, i);
      
      auto& Mi = M(nodeListi, i);
      auto& DpDxi = DpDx(nodeListi, i);
      auto& DvDxi = DvDx(nodeListi, i);

      const auto Mdeti = std::abs(Mi.Determinant());

      const auto enoughNeighbors =  numNeighborsi > Dimension::pownu(2);
      const auto goodM =  (Mdeti > 1e-2 and enoughNeighbors);                   

      Mi = ( goodM ? Mi.Inverse() : Tensor::one());

      DpDxi = Mi.Transpose()*DpDxi;
      DvDxi = DvDxi*Mi;

    } // loop nodes
  }   // loop nodelists
}     // function
}     // spheral namespace
