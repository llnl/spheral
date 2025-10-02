//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "CRKSPH/editMultimaterialSurfaceTopology.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template void editMultimaterialSurfaceTopology(FieldList<Dim<1>, int>& surfacePoint,
                                               ConnectivityMap<Dim<1>>& connectivityMap);
#endif

#if defined(SPHERAL_ENABLE_2D)
template void editMultimaterialSurfaceTopology(FieldList<Dim<2>, int>& surfacePoint,
                                               ConnectivityMap<Dim<2>>& connectivityMap);
#endif

#if defined(SPHERAL_ENABLE_3D)
template void editMultimaterialSurfaceTopology(FieldList<Dim<3>, int>& surfacePoint,
                                               ConnectivityMap<Dim<3>>& connectivityMap);
#endif

}
