text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "copy2polytope.cc"

namespace Spheral {
  template void copy2polytope(const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::FacetedVolume>& cells,
                              polytope::Tessellation<%(ndim)s, double>& mesh);
}
"""

# Only supports 2D and 3D
dimensions = [2, 3]
