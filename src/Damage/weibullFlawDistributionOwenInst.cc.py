text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "Damage/weibullFlawDistributionOwen.cc"

template Spheral::Field<Spheral::Dim< %(ndim)s >, std::vector<double> > 
Spheral::weibullFlawDistributionOwen<Spheral::Dim< %(ndim)s > >(const unsigned,
                                                                const double,
                                                                const double,
                                                                const Spheral::FluidNodeList<Spheral::Dim< %(ndim)s > >&,
                                                                const Spheral::State<Spheral::Dim< %(ndim)s > >&,
                                                                const size_t,
                                                                const double,
                                                                const Spheral::Field<Spheral::Dim<%(ndim)s>, int>&);
"""
