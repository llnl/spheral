//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "GSPH/Policies/CompatibleMFVSpecificThermalEnergyPolicy.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class CompatibleMFVSpecificThermalEnergyPolicy<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class CompatibleMFVSpecificThermalEnergyPolicy<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class CompatibleMFVSpecificThermalEnergyPolicy<Dim<3>>;
#endif

}
