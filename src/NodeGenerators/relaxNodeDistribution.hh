//------------------------------------------------------------------------------
// Centroidally relax a node distribution in a boundary.
// Optionally the user can specify a weighting function for the nodes.
//------------------------------------------------------------------------------
#ifndef __Spheral_relaxNodeDistribution__
#define __Spheral_relaxNodeDistribution__

#include "DataBase/DataBase.hh"
#include "Boundary/Boundary.hh"
#include "Kernel/TableKernel.hh"
#include "Geometry/Dimension.hh"

#include <vector>

namespace Spheral {

//------------------------------------------------------------------------------
// Base class for weighting the nodes.
//------------------------------------------------------------------------------
template<typename Dimension>
struct WeightingFunctor {
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using FacetedVolume = typename Dimension::FacetedVolume;

  WeightingFunctor() {}
  virtual ~WeightingFunctor() {}
  virtual double operator()(const Vector& pos, const FacetedVolume& boundary) const {
    return this->__call__(pos, boundary);
  }
  virtual double __call__(const Vector& /*pos*/, const FacetedVolume& /*boundary*/) const { return 1.0; }
};

//------------------------------------------------------------------------------
// The main method.
//------------------------------------------------------------------------------
template<typename Dimension>
void
relaxNodeDistribution(DataBase<Dimension>& dataBase,
                      const typename Dimension::FacetedVolume& boundary,
                      const std::vector<Boundary<Dimension>*>& boundaries,
                      const TableKernel<Dimension>& W,
                      const WeightingFunctor<Dimension>& weightingFunctor,
                      const WeightingFunctor<Dimension>& massDensityFunctor,
                      const double targetMass,
                      const int maxIterations,
                      const double tolerance);
}

#endif
