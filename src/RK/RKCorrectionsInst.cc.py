text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "RK/RKCorrections.cc"
#include "Geometry/Dimension.hh"
namespace Spheral {
template class RKCorrections<Dim<%(ndim)s>>;
}
"""
