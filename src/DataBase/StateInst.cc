//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "DataBase/State.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class State<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class State<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class State<Dim<3>>;
#endif

}
