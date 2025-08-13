#ifndef _Spheral_NodePairIdxType_hh_
#define _Spheral_NodePairIdxType_hh_

#include "Utilities/size_t_bits.hh"
#include "Utilities/DBC.hh"

#include <iostream>
#include <iostream>
#include <functional>   // hash

// These are based on what we get from size_t_bits
#define MAX_NODE_INDEX (size_t(1u) << ((SIZE_T_BITS - 10)/2))
#define MAX_NODELIST_INDEX (size_t(1u) << 5)

namespace Spheral {

//------------------------------------------------------------------------------
struct NodePairIdxType {
  using HashType = size_t;
  size_t i_node, i_list, j_node, j_list;
  double f_couple;                       // An arbitrary fraction in [0,1] to hold the effective coupling of the pair

  NodePairIdxType(size_t i_n,
                  size_t i_l,
                  size_t j_n,
                  size_t j_l,
                  double f = 1.0): i_node(i_n), i_list(i_l), j_node(j_n), j_list(j_l), f_couple(f) {}

  HashType hash() const {
    // We do this with simple bit shifting, requiring max values for the integer
    // components.  We assume the
    //    i_list, j_list < 32 (2^5)
    //    i_node, j_node < 134217728 (2^27) (on 64 bit machines)
    REQUIRE(i_node < MAX_NODE_INDEX);
    REQUIRE(j_node < MAX_NODE_INDEX);
    REQUIRE(i_list < MAX_NODELIST_INDEX);
    REQUIRE(j_list < MAX_NODELIST_INDEX);
    const auto flip = (i_list > j_list or
                       (i_list == j_list and i_node > j_node));
    const size_t i_l = flip ? j_list : i_list;
    const size_t i_n = flip ? j_node : i_node;
    const size_t j_l = flip ? i_list : j_list;
    const size_t j_n = flip ? i_node : j_node;
    return ((i_l << (SIZE_T_BITS - 5)) +
            (i_n << (SIZE_T_BITS/2)) +
            (j_l << (SIZE_T_BITS/2 - 5)) +
            j_n);
  }

  // Comparisons
  bool operator==(const NodePairIdxType& val) const { return (this->hash() == val.hash()); }
  bool operator< (const NodePairIdxType& val) const { return (this->hash() <  val.hash()); }
};

//------------------------------------------------------------------------------
// Output for NodePairIdxType
inline
std::ostream& operator<<(std::ostream& os, const NodePairIdxType& x) {
  os << "[(" << x.i_list << " " << x.i_node << ") (" << x.j_list << " " << x.j_node << ")]";
  return os;
}

} //namespace Spheral

//------------------------------------------------------------------------------
// Provide a method of hashing NodePairIdxType
namespace std {
  template<>
  struct hash<Spheral::NodePairIdxType> {
    size_t operator()(const Spheral::NodePairIdxType& x) const {
      return x.hash();
    }
  };
} // namespace std

#endif // _Spheral_NodePairIdxType_hh_
