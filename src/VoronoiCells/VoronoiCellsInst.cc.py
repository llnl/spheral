text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "VoronoiCells/VoronoiCells.cc"
#include "Geometry/Dimension.hh"
namespace Spheral {
template class VoronoiCells<Dim<%(ndim)s>>;
}
"""
