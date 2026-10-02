text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "ArtificialViscosity/PythonArtificialViscosity.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class PythonArtificialViscosity<Dim<%(ndim)s>, Dim<%(ndim)s>::Scalar>;
  template class PythonArtificialViscosity<Dim<%(ndim)s>, Dim<%(ndim)s>::Tensor>;
}
"""
