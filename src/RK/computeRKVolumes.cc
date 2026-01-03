//------------------------------------------------------------------------------
// Compute the RK volumes
//------------------------------------------------------------------------------
#include "computeRKVolumes.hh"

#include "VoronoiCells/computeVoronoiVolume.hh"
#include "computeHullVolumes.hh"
#include "computeRKSumVolume.hh"
#include "computeHVolumes.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"
#include "NodeList/NodeList.hh"

#include <vector>
#include <iostream>

namespace Spheral {

//------------------------------------------------------------------------------
// Compute the volumes
//------------------------------------------------------------------------------
template<typename Dimension>
void
computeRKVolumes(const ConnectivityMap<Dimension>& connectivityMap,
                 const TableKernel<Dimension>& W,
                 const FieldList<Dimension, typename Dimension::Vector>& position,
                 const FieldList<Dimension, typename Dimension::Scalar>& mass,
                 const FieldList<Dimension, typename Dimension::Scalar>& massDensity,
                 const FieldList<Dimension, typename Dimension::SymTensor>& H,
                 const FieldList<Dimension, typename Dimension::SymTensor>& damage,
                 const std::vector<typename Dimension::FacetedVolume>& facetedBoundaries,
                 const std::vector<std::vector<typename Dimension::FacetedVolume>>& facetedHoles,
                 const std::vector<Boundary<Dimension>*>& boundaryConditions,
                 const RKVolumeType volumeType,
                 FieldList<Dimension, int>& surfacePoint,
                 FieldList<Dimension, typename Dimension::Vector>& deltaCentroid,
                 FieldList<Dimension, std::vector<typename Dimension::Vector>>& etaVoidPoints,
                 FieldList<Dimension, typename Dimension::FacetedVolume>& cells,
                 FieldList<Dimension, std::vector<CellFaceFlag>>& cellFaceFlags,
                 FieldList<Dimension, typename Dimension::Scalar>& volume) {
  typedef typename Dimension::Scalar Scalar;
  
  switch(volumeType) {
  case RKVolumeType::RKMassOverDensity:
    CHECK(mass.size() == massDensity.size());
    volume.assignFields(mass/massDensity);
    break;
    
  case RKVolumeType::RKSumVolume:
    computeRKSumVolume(connectivityMap, W, position, mass, H, volume);
    break;

  case RKVolumeType::RKVoronoiVolume:
    {
      FieldList<Dimension, typename Dimension::Scalar> weights;
      volume.assignFields(mass/massDensity);
      computeVoronoiVolume(position, H, connectivityMap, damage,
                           facetedBoundaries, facetedHoles, boundaryConditions, weights,
                           surfacePoint, volume, deltaCentroid, etaVoidPoints,
                           cells, cellFaceFlags); 
    }
    break;

  case RKVolumeType::RKHullVolume:
    computeHullVolumes(connectivityMap, W.kernelExtent(), position, H, volume);
    break;

  case RKVolumeType::HVolume:
    {
      const Scalar nPerh = volume.nodeListPtrs()[0]->nodesPerSmoothingScale();
      computeHVolumes(nPerh, H, volume);
    }
    break;

  default:
    VERIFY2(false, "Unknown RK volume weighting.");
  }
}

} // end namespace Spheral
