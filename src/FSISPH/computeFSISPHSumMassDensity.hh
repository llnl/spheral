//---------------------------------Spheral++------------------------------------
// Compute normals at an interface
//------------------------------------------------------------------------------
#ifndef __Spheral__computeFSISPHSumMassDensity__
#define __Spheral__computeFSISPHSumMassDensity__

#include <vector>

namespace Spheral {

  // Forward declarations.
  template<typename Dimension> class ConnectivityMap;
  template<typename Dimension> class TableKernel;
  template<typename Dimension, typename DataType> class FieldList;


template<typename Dimension>
void
computeFSISPHSumMassDensity(const ConnectivityMap<Dimension>& connectivityMap,
                            const TableKernel<Dimension>& W,
                            const std::vector<int>& sumDensityNodeLists,
                            const FieldList<Dimension, typename Dimension::Vector>& position,
                            const FieldList<Dimension, typename Dimension::Scalar>& mass,
                            const FieldList<Dimension, typename Dimension::SymTensor>& H,
                            const bool consistentSum,
                            FieldList<Dimension, typename Dimension::Scalar>& massDensity);

}

 
 #endif
