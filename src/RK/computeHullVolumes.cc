//---------------------------------Spheral++------------------------------------
// Compute the volume per point based on convex hulls.
//------------------------------------------------------------------------------
#include "computeHullVolumes.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "NodeList/NodeList.hh"
#include "Hydro/HydroFieldNames.hh"

using std::vector;
using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {

template<typename Dimension>
void
computeHullVolumes(const ConnectivityMap<Dimension>& connectivityMap,
                   const typename Dimension::Scalar kernelExtent,
                   const FieldList<Dimension, typename Dimension::Vector>& position,
                   const FieldList<Dimension, typename Dimension::SymTensor>& H,
                   FieldList<Dimension, typename Dimension::Scalar>& volume) {

  // Pre-conditions.
  const auto numNodeLists = volume.size();
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);
  REQUIRE(kernelExtent > 0.0);

  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::FacetedVolume FacetedVolume;

  FieldList<Dimension, vector<Vector>> etaInv(FieldStorageType::CopyFields);
  for (const auto& fieldPtr: position) {
    etaInv.appendNewField("eta inv", fieldPtr->nodeList(), vector<Vector>());
  }

  // The set of interacting node pairs.
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();

#pragma omp parallel
  {
    // Some scratch variables.
    int i, j, nodeListi, nodeListj;

    auto etaInv_thread = etaInv.threadCopy();

    // Collect the half-way positions of all neighbors
#pragma omp for
    for (auto k = 0u; k < npairs; ++k) {
      i = pairs[k].i_node;
      j = pairs[k].j_node;
      nodeListi = pairs[k].i_list;
      nodeListj = pairs[k].j_list;

      // Get the state for node i.
      const auto& ri = position(nodeListi, i);
      const auto& Hi = H(nodeListi, i);

      // State for j
      const auto& rj = position(nodeListj, j);
      const auto& Hj = H(nodeListj, j);

      const auto  rji = 0.5*(rj - ri);
      const auto  etai =  Hi*rji;
      const auto  etaj = -Hj*rji;
      const auto  etaiMag = etai.magnitude();
      const auto  etajMag = etaj.magnitude();

      if (etaiMag < kernelExtent) {
        const auto etaiHat = etai.unitVector();
        etaInv_thread(nodeListi, i).push_back(1.0/max(etaiMag, 1.0e-30) * etaiHat);
      }

      if (etajMag < kernelExtent) {
        const auto etajHat = etaj.unitVector();
        etaInv_thread(nodeListj, j).push_back(1.0/max(etajMag, 1.0e-30) * etajHat);
      }
    }

#pragma omp critical
    {
      etaInv_thread.threadReduce();
    } // OMP critical
  }   // OMP parallel

    // Now we can do each node independently.
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto ni = position[nodeListi]->numInternalElements();
#pragma omp parallel for
    for (auto i = 0u; i < ni; ++i) {

      // Get the state for node i.
      const auto& Hi = H(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();

      // Build the hull of the inverse.
      const FacetedVolume hullInv(etaInv(nodeListi, i));

      // Use the vertices selected by the inverse hull to construct the
      // volume of the node.
      vector<Vector> eta;
      const vector<Vector>& vertsInv = hullInv.vertices();
      CHECK((Dimension::nDim == 1 and vertsInv.size() == 2) or
            (Dimension::nDim == 2 and vertsInv.size() >= 3) or
            (Dimension::nDim == 3 and vertsInv.size() >= 4));
      for (typename std::vector<Vector>::const_iterator itr = vertsInv.begin();
           itr != vertsInv.end();
           ++itr) {
        if (itr->magnitude2() < 1.0e-30) {
          eta.push_back(Vector::zero());
        } else {
          eta.push_back(1.0/sqrt(itr->magnitude2()) * itr->unitVector());
        }
      }

      // And we have it.
      const FacetedVolume polyeta = FacetedVolume(eta);
      volume(nodeListi, i) = polyeta.volume()/Hdeti;
    }
  }
}

}

