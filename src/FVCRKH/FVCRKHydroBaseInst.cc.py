text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "FVCRKHydroBase.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
template class FVCRKHydroBase< Dim< %(ndim)s > >;
}
"""
