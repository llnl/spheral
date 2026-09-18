text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "ArtificialViscosity/PythonScalarArtificialViscosity.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class PythonScalarArtificialViscosity<Dim<%(ndim)s>>;
}
"""
