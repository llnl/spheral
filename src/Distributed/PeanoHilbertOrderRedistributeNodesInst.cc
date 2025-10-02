//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "Distributed/PeanoHilbertOrderRedistributeNodes.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class PeanoHilbertOrderRedistributeNodes<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class PeanoHilbertOrderRedistributeNodes<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class PeanoHilbertOrderRedistributeNodes<Dim<3>>;
#endif

}
