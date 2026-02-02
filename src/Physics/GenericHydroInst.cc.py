text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Physics/GenericHydro.cc"

namespace Spheral {
  template class GenericHydro< Dim< %(ndim)s > >;
}
"""
