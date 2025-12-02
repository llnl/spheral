text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "ArtificialViscosity/TensorMonaghanGingoldViscosityView.hh"

namespace Spheral {
  template class TensorMonaghanGingoldViscosityView< Dim< %(ndim)s > >;
}
"""
