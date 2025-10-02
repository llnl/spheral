//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "GSPH/RiemannSolvers/RiemannSolverBase.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template class RiemannSolverBase<Dim<1>>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class RiemannSolverBase<Dim<2>>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class RiemannSolverBase<Dim<3>>;
#endif

}
