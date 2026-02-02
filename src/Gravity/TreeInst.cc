//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------
#include "Tree.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {

#ifdef SPHERAL2D
  template class Tree<Dim<2>>;
#endif
#ifdef SPHERAL3D
  template class Tree<Dim<3>>;
#endif

}
