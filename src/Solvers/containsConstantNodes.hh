//---------------------------------Spheral++----------------------------------//
// containsConstantNodes
//
// Checks whether boundaries contain constant nodes
//----------------------------------------------------------------------------//
#ifndef __Spheral_containsConstantNodes_hh__
#define __Spheral_containsConstantNodes_hh__

#include <vector>

namespace Spheral {

template<typename Dimension> class Boundary;

template<typename Dimension>
inline bool containsConstantNodes(const Boundary<Dimension>* boundary);

template<typename Dimension>
inline bool containsConstantNodes(const std::vector<Boundary<Dimension>*>& boundaries);

} // end namespace Spheral

#include "containsConstantNodesInline.hh"

#endif
