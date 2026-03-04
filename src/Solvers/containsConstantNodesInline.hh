#include "Boundary/ConstantBoundary.hh"
#include "Boundary/InflowOutflowBoundary.hh"

namespace Spheral {

template<typename Dimension>
inline
bool
containsConstantNodes(const Boundary<Dimension>* boundary) {
  // This is bad. It would probably be better to add a bool to the boundaries
  // that is true for boundaries that contain constant nodes.
  return (dynamic_cast<const ConstantBoundary<Dimension>*>(boundary) != nullptr
          || dynamic_cast<const InflowOutflowBoundary<Dimension>*>(boundary) != nullptr);
}

template<typename Dimension>
inline
bool
containsConstantNodes(const std::vector<Boundary<Dimension>*>& boundaries) {
  for (auto& boundary : boundaries) {
    if (containsConstantNodes(boundary)) {
      return true;
    }
  }
  return false;
}

} // end namespace Spheral
