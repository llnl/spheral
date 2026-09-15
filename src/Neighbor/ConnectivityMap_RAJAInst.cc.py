text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Neighbor/ConnectivityMap_RAJA.cc"
#include "Geometry/Dimension.hh"

template class Spheral::ConnectivityMap_RAJA<Spheral::Dim< %(ndim)s > >;
"""
