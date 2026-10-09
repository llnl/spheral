text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "ArtificialViscosity/TensorMonaghanGingoldViscosity.hh"

namespace Spheral {
  template class TensorMonaghanGingoldViscosity< Dim< %(ndim)s > >;
  template class TensorMonaghanGingoldViscosityView< Dim< %(ndim)s > >;
}
"""
