#ifndef _Spheral_NodePairList_hh_
#define _Spheral_NodePairList_hh_

#include "Neighbor/NodePairIdxType.hh"

#include <vector>
#include <unordered_map>

namespace Spheral {

//------------------------------------------------------------------------------
class NodePairList {
public:
  using ContainerType = std::vector<NodePairIdxType>;
  using value_type = typename ContainerType::value_type;
  using reference = typename ContainerType::reference;
  using const_reference = typename ContainerType::const_reference;
  using iterator = typename ContainerType::iterator;
  using const_iterator = typename ContainerType::const_iterator;
  using reverse_iterator = typename ContainerType::reverse_iterator;
  using const_reverse_iterator = typename ContainerType::const_reverse_iterator;

  NodePairList()                                                               = default;
  NodePairList(const NodePairList& rhs)                                        = default;
  ~NodePairList()                                                              = default;
  NodePairList& operator=(const NodePairList& rhs)                             = default;
  void push_back(NodePairIdxType nodePair)                                     { mNodePairList.push_back(nodePair); }
  void clear()                                                                 { mNodePairList.clear(); mPair2Index.clear(); }
  void reserve(const size_t n)                                                 { mNodePairList.reserve(n); }
  void resize(const size_t n) { mNodePairList.resize(n, NodePairIdxType(-1,-1,-1,-1)); }
  size_t size() const                                                          { return mNodePairList.size(); }
  void makeUnique();

  // Iterators
  iterator begin()                                                             { return mNodePairList.begin(); }
  iterator end()                                                               { return mNodePairList.end(); }
  const_iterator begin() const                                                 { return mNodePairList.begin(); }
  const_iterator end() const                                                   { return mNodePairList.end(); }

  // Reverse iterators
  reverse_iterator rbegin()                                                    { return mNodePairList.rbegin(); }
  reverse_iterator rend()                                                      { return mNodePairList.rend(); }
  const_reverse_iterator rbegin() const                                        { return mNodePairList.rbegin(); }
  const_reverse_iterator rend() const                                          { return mNodePairList.rend(); }

  // Indexing
  reference operator[](const size_t i)                                         { return mNodePairList[i]; }
  reference operator()(const NodePairIdxType& x)                               { return mNodePairList[index(x)]; }
  reference operator()(const size_t i_node,
                       const size_t i_list,
                       const size_t j_node,
                       const size_t j_list)                                    { return mNodePairList[index(NodePairIdxType(i_node, i_list, j_node, j_list))]; }

  const_reference operator[](const size_t i) const                             { return mNodePairList[i]; }
  const_reference operator()(const NodePairIdxType& x) const                   { return mNodePairList[index(x)]; }
  const_reference operator()(const size_t i_node,
                             const size_t i_list,
                             const size_t j_node,
                             const size_t j_list) const                        { return mNodePairList[index(NodePairIdxType(i_node, i_list, j_node, j_list))]; }

  // Inserting
  template<typename InputIterator>
  iterator insert(const_iterator pos, InputIterator first, InputIterator last) { return mNodePairList.insert(pos, first, last); }

  // Find the index corresponding to the given pair
  size_t index(const NodePairIdxType& x) const;

  // Compute the lookup table for Pair->index
  void computeLookup() const;

private:
  ContainerType mNodePairList;
  mutable std::unordered_map<NodePairIdxType, size_t> mPair2Index;  // mutable for lazy evaluation in index
};

}

#endif // _Spheral_NeighbourSpace_NodePairList_hh_
