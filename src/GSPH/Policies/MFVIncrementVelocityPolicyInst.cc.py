text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "GSPH/Policies/MFVIncrementVelocityPolicy.cc"

namespace Spheral {
  template class MFVIncrementVelocityPolicy<Dim< %(ndim)s >>;
}
"""
