text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "VoronoiCells/IncrementVoronoiCells.cc"
#include "Geometry/Dimension.hh"
namespace Spheral {
template class IncrementVoronoiCells<Dim<%(ndim)s>>;
}
"""
