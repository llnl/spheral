//---------------------------------Spheral++------------------------------------
// Compute the volume per point based on the Voronoi tessellation.
//------------------------------------------------------------------------------
#ifndef __Spheral__computeVoronoiVolume__
#define __Spheral__computeVoronoiVolume__

#include "Geometry/CellFaceFlag.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Boundary/Boundary.hh"

#include <vector>

namespace Spheral {

template<typename Dimension>
void
computeVoronoiVolume(const FieldList<Dimension, typename Dimension::Vector>& position,
                     const FieldList<Dimension, typename Dimension::SymTensor>& H,
                     const ConnectivityMap<Dimension >& connectivityMap,
                     const FieldList<Dimension, typename Dimension::SymTensor>& damage,
                     const std::vector<typename Dimension::FacetedVolume>& facetedBoundaries,
                     const std::vector<std::vector<typename Dimension::FacetedVolume> >& holes,
                     const std::vector<Boundary<Dimension>*>& boundaries,
                     const FieldList<Dimension, typename Dimension::Scalar>& weight,
                     FieldList<Dimension, int>& surfacePoint,
                     FieldList<Dimension, typename Dimension::Scalar>& vol,
                     FieldList<Dimension, typename Dimension::Vector>& deltaMedian,
                     FieldList<Dimension, std::vector<typename Dimension::Vector>>& etaVoidPoints,
                     FieldList<Dimension, typename Dimension::FacetedVolume>& cells,
                     FieldList<Dimension, std::vector<CellFaceFlag>>& cellFaceFlags);

#ifdef SPHERAL1D
// 1D is specialized
template<>
void
computeVoronoiVolume(const FieldList<Dim<1>, Dim<1>::Vector>& position,
                     const FieldList<Dim<1>, Dim<1>::SymTensor>& H,
                     const ConnectivityMap<Dim<1> >& connectivityMap,
                     const FieldList<Dim<1>, Dim<1>::SymTensor>& damage,
                     const std::vector<Dim<1>::FacetedVolume>& facetedBoundaries,
                     const std::vector<std::vector<Dim<1>::FacetedVolume> >& holes,
                     const std::vector<Boundary<Dim<1>>*>& boundaries,
                     const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                     FieldList<Dim<1>, int>& surfacePoint,
                     FieldList<Dim<1>, Dim<1>::Scalar>& vol,
                     FieldList<Dim<1>, Dim<1>::Vector>& deltaMedian,
                     FieldList<Dim<1>, std::vector<Dim<1>::Vector>>& etaVoidPoints,
                     FieldList<Dim<1>, Dim<1>::FacetedVolume>& cells,
                     FieldList<Dim<1>, std::vector<CellFaceFlag>>& cellFaceFlags);
#endif

}

#endif
