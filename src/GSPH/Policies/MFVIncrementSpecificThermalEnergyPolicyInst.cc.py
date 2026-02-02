text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "GSPH/Policies/MFVIncrementSpecificThermalEnergyPolicy.cc"

namespace Spheral {
  template class MFVIncrementSpecificThermalEnergyPolicy<Dim< %(ndim)s >>;
}
"""
