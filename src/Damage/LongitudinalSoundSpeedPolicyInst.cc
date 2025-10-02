//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "Damage/LongitudinalSoundSpeedPolicy.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class LongitudinalSoundSpeedPolicy<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class LongitudinalSoundSpeedPolicy<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class LongitudinalSoundSpeedPolicy<Dim<3>>;
#endif

}
