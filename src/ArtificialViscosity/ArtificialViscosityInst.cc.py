text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "ArtificialViscosity/ArtificialViscosity.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class ArtificialViscosity<Dim<%(ndim)s>>;
  //template class ArtificialViscosityView<Dim<%(ndim)s>, Dim<%(ndim)s>::Scalar>;
  //template class ArtificialViscosityView<Dim<%(ndim)s>, Dim<%(ndim)s>::Tensor>;
}
"""
