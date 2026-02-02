//---------------------------------Spheral++------------------------------------
// Compute the volume per point using an inverse convex hull.
// Optionally this volume can then be clipped to the Voronoi.
//------------------------------------------------------------------------------
#ifndef __Spheral__computeHullVolume_v2__
#define __Spheral__computeHullVolume_v2__

#include "Geometry/CellFaceFlag.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"

namespace Spheral {

template<typename Dimension>
void
computeHullVolume(const FieldList<Dimension, typename Dimension::Vector>& position,
                  const FieldList<Dimension, typename Dimension::SymTensor>& H,
                  const ConnectivityMap<Dimension >& connectivityMap,
                  const bool clipToVoronoi,
                  FieldList<Dimension, int>& surfacePoint,
                  FieldList<Dimension, typename Dimension::Scalar>& vol,
                  FieldList<Dimension, typename Dimension::FacetedVolume>& cells);

}

#endif
