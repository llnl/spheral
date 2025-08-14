text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "CRKSPH/SolidCRKSPH.cc"

namespace Spheral {
template class SolidCRKSPH< Dim< %(ndim)s > >;
}
"""
