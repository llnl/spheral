text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "CRKSPH/CRKSPHBase.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
template class CRKSPHBase<Dim<%(ndim)s>>;
}
"""
