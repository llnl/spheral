text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "FSISPH/SolidFSISPH.cc"
#include "FSISPH/SolidFSISPHEvaluateDerivatives.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class SolidFSISPH<Dim<%(ndim)s>>;
}
"""
