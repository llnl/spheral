#include <typeinfo>

namespace Spheral {

//------------------------------------------------------------------------------
// Return the number of nodes this Boundary is going to create.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
size_t
ConstantBoundary<Dimension>::
numConstantNodes() const {
  return mNumConstantNodes;
}

//------------------------------------------------------------------------------
// Return the NodeList this Boundary is defined on.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const NodeList<Dimension>&
ConstantBoundary<Dimension>::
nodeList() const {
  CHECK(mNodeListPtr != 0);
  return *mNodeListPtr;
}

//------------------------------------------------------------------------------
// The reflection operator.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Tensor
ConstantBoundary<Dimension>::
reflectOperator() const {
  return mReflectOperator;
}

}
