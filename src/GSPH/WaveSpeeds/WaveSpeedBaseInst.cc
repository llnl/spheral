//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "GSPH/WaveSpeeds/WaveSpeedBase.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class WaveSpeedBase<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class WaveSpeedBase<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class WaveSpeedBase<Dim<3>>;
#endif

}
