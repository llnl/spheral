text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "GSPH/GSPH.cc"
#include "GSPH/GSPHEvaluateDerivatives.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class GSPH<Dim<%(ndim)s>>;
}
"""
