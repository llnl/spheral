text = """
//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------
#include "SPHGravity/SPHGravity.cc"

namespace Spheral {
  template class SPHGravity<Dim< %(ndim)s > >;
} // end namespace Spheral

"""
