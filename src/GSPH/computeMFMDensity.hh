//---------------------------------Spheral++----------------------------------//
// Compute the density from m/V
//
// J.M. Pearl 2022
//----------------------------------------------------------------------------//

#ifndef __Spheral__computeMFMDensity__
#define __Spheral__computeMFMDensity__

#include <vector>

namespace Spheral {

  // Forward declarations.
  template<typename Dimension> class ConnectivityMap;
  template<typename Dimension> class TableKernel;
  template<typename Dimension, typename DataType> class FieldList;


template<typename Dimension>
void
computeMFMDensity(const FieldList<Dimension, typename Dimension::Scalar>& mass,
                  const FieldList<Dimension, typename Dimension::Scalar>& massDensity,
                        FieldList<Dimension, typename Dimension::Scalar>& volume);

}

 
 #endif