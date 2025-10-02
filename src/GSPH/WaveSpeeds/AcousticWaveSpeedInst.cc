//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "GSPH/WaveSpeeds/AcousticWaveSpeed.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class AcousticWaveSpeed<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class AcousticWaveSpeed<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class AcousticWaveSpeed<Dim<3>>;
#endif

}
