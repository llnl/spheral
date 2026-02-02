//---------------------------------Spheral++----------------------------------//
// initializes the pressure and velocity gradients for Riemann solver - based
// SPH varients
//
// J.M. Pearl 2023
//----------------------------------------------------------------------------//

#ifndef __Spheral__initializeGradients__
#define __Spheral__initializeGradients__

#include <vector>

namespace Spheral {

  // Forward declarations.
  template<typename Dimension> class ConnectivityMap;
  template<typename Dimension> class TableKernel;
  template<typename Dimension, typename DataType> class FieldList;


template<typename Dimension>
void
initializeGradients(const ConnectivityMap<Dimension>& connectivityMap,
                    const TableKernel<Dimension>& W,
                    const FieldList<Dimension, typename Dimension::Vector>& position,
                    const FieldList<Dimension, typename Dimension::SymTensor>& H,
                    const FieldList<Dimension, typename Dimension::Scalar>& volume,
                    const FieldList<Dimension, typename Dimension::Scalar>& pressure,
                    const FieldList<Dimension, typename Dimension::Vector>& velocity,
                          FieldList<Dimension, typename Dimension::Tensor>& M,
                          FieldList<Dimension, typename Dimension::Vector>& DpDx,
                          FieldList<Dimension, typename Dimension::Tensor>& DvDx);

}

 
 #endif