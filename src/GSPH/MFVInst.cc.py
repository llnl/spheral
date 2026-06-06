text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "GSPH/MFV.cc"
#include "GSPH/MFVEvaluateDerivatives.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class MFV<Dim<%(ndim)s>>;
}
"""
