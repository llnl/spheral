text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "GSPH/MFM.cc"
#include "GSPH/MFMEvaluateDerivatives.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class MFM<Dim<%(ndim)s>>;
}
"""
