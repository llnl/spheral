text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "NodeGenerators/relaxNodeDistribution.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template void relaxNodeDistribution(DataBase<Dim< %(ndim)s > >& dataBase,
                                      const Dim< %(ndim)s >::FacetedVolume& boundary,
                                      const std::vector<Boundary<Dim< %(ndim)s > >*>& boundaries,
                                      const TableKernel<Dim< %(ndim)s > >& W,
                                      const WeightingFunctor<Dim< %(ndim)s > >& weightingFunctor,
                                      const WeightingFunctor<Dim< %(ndim)s > >& massDensityFunctor,
                                      const double targetMass,
                                      const int maxIterations,
                                      const double tolerance);
}
"""
