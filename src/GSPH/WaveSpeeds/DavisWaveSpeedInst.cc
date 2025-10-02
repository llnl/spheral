//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "GSPH/WaveSpeeds/DavisWaveSpeed.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class DavisWaveSpeed<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class DavisWaveSpeed<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class DavisWaveSpeed<Dim<3>>;
#endif

}
