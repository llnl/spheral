#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Access the Hextent variable.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
NestedGridRedistributeNodes<Dimension>::Hextent() const {
  return mHextent;
}

template<typename Dimension>
inline
void
NestedGridRedistributeNodes<Dimension>::Hextent(const double val) {
  REQUIRE(distinctlyGreaterThan(val, 0.0));
  mHextent = val;
}

}
