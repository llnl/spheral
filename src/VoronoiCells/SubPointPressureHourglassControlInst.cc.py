text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "VoronoiCells/SubPointPressureHourglassControl.cc"
#include "Geometry/Dimension.hh"
namespace Spheral {
template class SubPointPressureHourglassControl<Dim<%(ndim)s>>;
}
"""
