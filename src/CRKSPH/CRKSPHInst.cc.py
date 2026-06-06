text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "CRKSPH/CRKSPH.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
template class CRKSPH<Dim<%(ndim)s>>;
}
"""
