text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "VoronoiCells/UpdateVoronoiCells.cc"
#include "Geometry/Dimension.hh"
namespace Spheral {
template class UpdateVoronoiCells<Dim<%(ndim)s>>;
}
"""
