//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "Field/InternalNodeIterator.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class InternalNodeIterator<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class InternalNodeIterator<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class InternalNodeIterator<Dim<3>>;
#endif

}
