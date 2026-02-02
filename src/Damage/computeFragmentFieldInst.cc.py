text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Damage/computeFragmentField.cc"

namespace Spheral {
  template Field<Dim< %(ndim)s >, int> computeFragmentField(const NodeList<Dim< %(ndim)s > >&, const double, const Field<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&, const Field<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>&, const Field<Dim< %(ndim)s >, int>&, const double, const double, const bool);
}
"""
