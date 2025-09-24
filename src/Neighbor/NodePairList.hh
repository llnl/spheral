#ifndef Spheral_NodePairList_hh
#define Spheral_NodePairList_hh

#include "Neighbor/NodePairIdxType.hh"

#include <vector>
#include <unordered_map>
#include "config.hh"
#include "chai/ManagedArray.hpp"
#include "chai/ExecutionSpaces.hpp"
#include "Utilities/CHAI_MA_wrapper.hh"

namespace Spheral {

class NodePairListView : public chai::CHAICopyable {
  using MAContainer = typename chai::ManagedArray<NodePairIdxType>;

public:
  SPHERAL_HOST_DEVICE NodePairListView() = default;
  SPHERAL_HOST_DEVICE ~NodePairListView() = default;
  SPHERAL_HOST NodePairListView(MAContainer const &d) : mData(d) {}

  SPHERAL_HOST_DEVICE
  NodePairIdxType& operator[](const size_t i) { return mData[i]; }

  SPHERAL_HOST_DEVICE
  NodePairIdxType& operator[](const size_t i) const { return mData[i]; }

  SPHERAL_HOST_DEVICE
  size_t size() const { return mData.size(); }
  SPHERAL_HOST_DEVICE
  const NodePairIdxType* data() const { return mData.data(); }

  void move(chai::ExecutionSpace space) { mData.move(space); }

  SPHERAL_HOST
  void touch(chai::ExecutionSpace space) { mData.registerTouch(space); }

protected:
  MAContainer mData;
};

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
    initMA();
    return n;
  }

  // Find the index corresponding to the given pair
  size_t index(const NodePairIdxType& x) const;

  // Compute the lookup table for Pair->index
  void computeLookup() const;

  inline NodePairListView view() {
    return static_cast<NodePairListView>(*this);
  }

  void initMA() {
    initializeManagedArray(mData, mNodePairList);
  }

  template<typename F> inline
  void setUserCallback(F&& extension) {
    mData.setUserCallback(getNPLCallback(std::forward<F>(extension)));
  }

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
