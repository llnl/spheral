//---------------------------------Spheral++----------------------------------//
// Compute volume from inverse of the kernel summation.
//
//   Hopkins P.F. (2015) "A New Class of Accurate, Mesh-Free Hydrodynamic 
//   Simulation Methods," MNRAS, 450(1):53-110
//
// J.M. Pearl 2022
//----------------------------------------------------------------------------//

#ifndef __Spheral__computeSumVolume__
#define __Spheral__computeSumVolume__

#include <vector>

namespace Spheral {

  // Forward declarations.
  template<typename Dimension> class ConnectivityMap;
  template<typename Dimension> class TableKernel;
  template<typename Dimension, typename DataType> class FieldList;


template<typename Dimension>
void
computeSumVolume(const ConnectivityMap<Dimension>& connectivityMap,
                 const TableKernel<Dimension>& W,
                 const FieldList<Dimension, typename Dimension::Vector>& position,
                 const FieldList<Dimension, typename Dimension::SymTensor>& H,
                       FieldList<Dimension, typename Dimension::Scalar>& volume);

}

 
 #endif