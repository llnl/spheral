//---------------------------------Spheral++----------------------------------//
// IntegrationKernel
//----------------------------------------------------------------------------//
#ifndef __Spheral_IntegrationKernel__
#define __Spheral_IntegrationKernel__

#include "Field/FieldList.hh"
#include "Geometry/Dimension.hh"

#include <utility>
#include <vector>

namespace Spheral {

template<typename Dimension>
class IntegrationKernel {
public:
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;
  
  IntegrationKernel() { }

  virtual ~IntegrationKernel() { };

  virtual double extent(const Scalar Hmult) const = 0;
  
  virtual void evaluate(const Vector& xp,
                        const std::vector<std::pair<int, int>>& indices,
                        const FieldList<Dimension, Vector>& position,
                        const FieldList<Dimension, SymTensor>& H,
                        const FieldList<Dimension, Scalar>& volume,
                        const Scalar Hmult,
                        std::vector<Scalar>& values,
                        std::vector<Vector>& dvalues) const = 0;
};

} // end namespace Spheral

#endif
