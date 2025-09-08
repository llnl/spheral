//---------------------------------Spheral++----------------------------------//
//
// J.M. Pearl 2023
//----------------------------------------------------------------------------//
#ifndef __Spheral_setUniqueNodeIDs_hh__
#define __Spheral_setUniqueNodeIDs_hh__

namespace Spheral {

template<typename Dimension, typename DataType> class FieldList;

template<typename Dimension>
void
setUniqueNodeIDs(FieldList<Dimension,size_t>& uniqueIndex);

}

#include "setUniqueNodeIDsInline.hh"

#endif
