text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "Damage/TensorStrainPolicy.cc"

namespace Spheral {
  template class TensorStrainPolicy<Dim<%(ndim)s>>;
}
"""
