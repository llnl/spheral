text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "GSPH/Policies/CompatibleMFVSpecificThermalEnergyPolicy.cc"

namespace Spheral {
  template class CompatibleMFVSpecificThermalEnergyPolicy<Dim< %(ndim)s > >;
}
"""
