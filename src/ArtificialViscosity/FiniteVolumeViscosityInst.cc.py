text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "ArtificialViscosity/FiniteVolumeViscosity.cc"
#include "ArtificialViscosity/FiniteVolumeViscosityView.cc"

namespace Spheral {
  template class FiniteVolumeViscosity< Dim< %(ndim)s > >;
  template class FiniteVolumeViscosityView< Dim< %(ndim)s > >;
}
"""
