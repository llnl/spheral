//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "GSPH/WaveSpeeds/EinfeldtWaveSpeed.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class EinfeldtWaveSpeed<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class EinfeldtWaveSpeed<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class EinfeldtWaveSpeed<Dim<3>>;
#endif

}
