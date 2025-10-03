dimensions = [2, 3]

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
