#ifndef Spheral_NodePairList_hh
#define Spheral_NodePairList_hh

#include "Neighbor/NodePairIdxType.hh"
#include "NodePairListView.hh"
#include "config.hh"
#include "chai/ManagedArray.hpp"
#include "chai/ExecutionSpaces.hpp"
#include "Utilities/CHAI_MA_wrapper.hh"
#include "chai/config.hpp"

#include <vector>
#include <unordered_map>

namespace Spheral {

//------------------------------------------------------------------------------
class NodePairList : public NodePairListView {
public:
  using ContainerType = std::vector<NodePairIdxType>;
  using value_type = typename ContainerType::value_type;
  using reference = typename ContainerType::reference;
  using const_reference = typename ContainerType::const_reference;
  using iterator = typename ContainerType::iterator;
  using const_iterator = typename ContainerType::const_iterator;
  using reverse_iterator = typename ContainerType::reverse_iterator;
  using const_reverse_iterator = typename ContainerType::const_reverse_iterator;

  NodePairList()                                             = default;

  // Constructor: copies underlying data
  NodePairList(const ContainerType& vals);

  // Constructor: moves underlying data
  NodePairList(ContainerType&& vals) noexcept;

  NodePairList(const NodePairList& rhs);
  NodePairList& operator=(const NodePairList& rhs);

  ~NodePairList()                                            { mData.free(); }

  void fill(const ContainerType& vals);
  void clear();

  // Iterators
  iterator begin()                                           { return mNodePairList.begin(); }
  iterator end()                                             { return mNodePairList.end(); }
  const_iterator begin() const                               { return mNodePairList.begin(); }
  const_iterator end() const                                 { return mNodePairList.end(); }

  // Reverse iterators
  reverse_iterator rbegin()                                  { return mNodePairList.rbegin(); }
  reverse_iterator rend()                                    { return mNodePairList.rend(); }
  const_reverse_iterator rbegin() const                      { return mNodePairList.rbegin(); }
  const_reverse_iterator rend() const                        { return mNodePairList.rend(); }

  // Indexing
  reference operator()(const NodePairIdxType& x)             { return mNodePairList[index(x)]; }
  reference operator()(const size_t i_node,
                       const size_t i_list,
                       const size_t j_node,
                       const size_t j_list)                  { return mNodePairList[index(NodePairIdxType(i_node, i_list, j_node, j_list))]; }

  const_reference operator()(const NodePairIdxType& x) const { return mNodePairList[index(x)]; }
  const_reference operator()(const size_t i_node,
                             const size_t i_list,
                             const size_t j_node,
                             const size_t j_list) const      { return mNodePairList[index(NodePairIdxType(i_node, i_list, j_node, j_list))]; }

  // Inserting (not performant, avoid if possible)
  template<typename InputIterator>
  iterator insert(const_iterator pos, InputIterator first, InputIterator last) {
    iterator n = mNodePairList.insert(pos, first, last);
    initView();
    return n;
  }

  // Find the index corresponding to the given pair
  size_t index(const NodePairIdxType& x) const;

  // Compute the lookup table for Pair->index
  void computeLookup() const;

  inline NodePairListView view() {
    return static_cast<NodePairListView>(*this);
  }

  void initView() {
    initMAView(mData, mNodePairList);
  }

#ifndef CHAI_DISABLE_RM
  template<typename F> inline
  void setUserCallback(F&& extension) {
#ifndef CHAI_DISABLE_RM
    mData.setUserCallback(getNPLCallback(std::forward<F>(extension)));
#endif
  }
#endif

protected:
  template<typename F>
  auto getNPLCallback(F callback) {
    return [callback](
      const chai::PointerRecord * record,
      chai::Action action,
      chai::ExecutionSpace space) {
             callback(record, action, space);
           };
  }
private:
  ContainerType mNodePairList;
  mutable std::unordered_map<NodePairIdxType, size_t> mPair2Index;  // mutable for lazy evaluation in index
};

}

#endif // Spheral_NodePairList_hh
