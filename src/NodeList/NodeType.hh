//---------------------------------Spheral++----------------------------------//
// NodeType -- tags used to distinguish internal and ghost nodes.
//----------------------------------------------------------------------------//
#ifndef __Spheral_NodeType_hh__
#define __Spheral_NodeType_hh__

namespace Spheral {

enum class NodeType {
  InternalNode = 0,
  GhostNode = 1
};

}

#endif
