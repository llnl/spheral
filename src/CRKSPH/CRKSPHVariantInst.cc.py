text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "CRKSPH/CRKSPHVariant.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
template class CRKSPHVariant< Dim< %(ndim)s > >;
}
"""
