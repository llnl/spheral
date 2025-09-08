text = """
//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------
#include "Gravity/NBodyGravity.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {

  template class NBodyGravity<Dim< %(ndim)s > >;

} // end namespace Spheral

"""
