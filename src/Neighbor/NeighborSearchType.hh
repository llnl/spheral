//---------------------------------Spheral++----------------------------------//
// NeighborSearchType
//
// The enum for the different neighbor search choices
//----------------------------------------------------------------------------//
#ifndef __Spheral_NeighborSearchType__
#define __Spheral_NeighborSearchType__

namespace Spheral {

enum class NeighborSearchType {
  Gather = 1,
  Scatter = 2,
  GatherScatter = 3
};

}

#endif
