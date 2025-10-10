//---------------------------------Spheral++------------------------------------
// Compute the RK interpolation at each point.
//------------------------------------------------------------------------------
#ifndef __Spheral__interpolateRK__
#define __Spheral__interpolateRK__

#include "Utilities/NodeCoupling.hh"
#include "RK/RKCorrectionParams.hh"
#include "RK/ReproducingKernel.hh"

#include <vector>
#include <variant>

namespace Spheral {

// Forward declarations.
template<typename Dimension> class ConnectivityMap;
template<typename Dimension, typename DataType> class FieldList;

template<typename Dimension>
std::vector<std::variant<FieldList<Dimension, typename Dimension::Scalar>,
                         FieldList<Dimension, typename Dimension::Vector>,
                         FieldList<Dimension, typename Dimension::Tensor>,
                         FieldList<Dimension, typename Dimension::SymTensor>,
                         FieldList<Dimension, typename Dimension::ThirdRankTensor>>>
interpolateRK(const std::vector<std::variant<FieldList<Dimension, typename Dimension::Scalar>,
                                             FieldList<Dimension, typename Dimension::Vector>,
                                             FieldList<Dimension, typename Dimension::Tensor>,
                                             FieldList<Dimension, typename Dimension::SymTensor>,
                                             FieldList<Dimension, typename Dimension::ThirdRankTensor>>>& fieldLists,
                  const FieldList<Dimension, typename Dimension::Vector>& position,
                  const FieldList<Dimension, typename Dimension::Scalar>& weight,
                  const FieldList<Dimension, typename Dimension::SymTensor>& H,
                  const ConnectivityMap<Dimension>& connectivityMap,
                  const ReproducingKernel<Dimension>& WR,
                  const FieldList<Dimension, RKCoefficients<Dimension>>& corrections,
                  const NodeCoupling& nodeCoupling = NodeCoupling());

}

#endif
