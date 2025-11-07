text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "ArtificialViscosity/TensorMonaghanGingoldViscosityView.cc"

namespace Spheral {
  template class TensorMonaghanGingoldViscosityView< Dim< %(ndim)s > >;
}
"""
