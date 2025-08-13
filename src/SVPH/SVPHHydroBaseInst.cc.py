//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SVPH/SVPHHydroBase.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class SVPHHydroBase< Dim<1> >;
  template class SVPHHydroBase< Dim<2> >;
  template class SVPHHydroBase< Dim<3> >;
}
