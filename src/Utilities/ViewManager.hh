//------------------------------------------------------------------------------
// ViewManager
//
// Manage the movement for a collection of View objects between different
// memory spaces (GPU->CPU, CPU->GPU, etc)
//------------------------------------------------------------------------------
#ifndef __Spheral_ViewManager__
#define __Spheral_ViewManager__

#include "Field/FieldViewBase.hh"
#include "Field/FieldListViewBase.hh"
#include "Neighbor/PairwiseFieldViewBase.hh"

#include "chai/ManagedArray.hpp"
#include "chai/ExecutionSpaces.hpp"

#include <functional>
#include <variant>
#include <vector>
#include <any>

namespace Spheral {
  
template<typename Dimension>
class ViewManager {
public:

  // A variant type of all the View types we can manage
  using ViewType = std::variant<std::reference_wrapper<FieldViewBase<Dimension>>,
                                std::reference_wrapper<FieldListViewBase<Dimension>>,
                                std::reference_wrapper<PairwiseFieldViewBase>>;
  using KeyType = typename StateBase<Dimension>::KeyType;

  // Constructors
  ViewManager() = default;
  ViewManager(const StateBase<Dimension>& state);
  ~ViewManager() = default;

  // Add View types for us to manage
  template<typename T> void enroll(T& view);

  // Get FieldViews from the StateBase
  template<typename Value> FieldView<Dimension, Value> field(const KeyType& key);
  template<typename Value> FieldView<Dimension, Value> field(const KeyType& key,
                                                             const Value& dummy);

  // Get FieldListViews from the StateBase
  template<typename Value> FieldListView<Dimension, Value> fields(const KeyType& key,
                                                                  bool allowNone = false);
  template<typename Value> FieldListView<Dimension, Value> fields(const KeyType& key,
                                                                  const Value& dummy,
                                                                  bool allowNone = false);

  // Get an arbitrary type from the StateBase
  template<typename Value, typename View = typename Value::ViewType> View get(const KeyType& key);

  // Generic methods we can call on all stored View types
  void move(chai::ExecutionSpace space);
  void touch(chai::ExecutionSpace space);

  // No copying or assignment
  ViewManager(const ViewManager&) = delete;
  ViewManager& operator=(const ViewManager&) = delete;

private:
  std::vector<ViewType> mViews;
  std::vector<std::any> mValueCache;
  const StateBase<Dimension>* mStateBasePtr = nullptr;
};

}

#include "ViewManagerInline.hh"

#endif
