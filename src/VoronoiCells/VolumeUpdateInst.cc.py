text = """
#include "VoronoiCells/VolumeUpdate.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class VolumeUpdate<Dim<%(ndim)s>>;
}
"""
