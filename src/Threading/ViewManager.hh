//------------------------------------------------------------------------------
// ViewManager
//
// Manage the movement for a collection of View objects between different
// memory spaces (GPU->CPU, CPU->GPU, etc)
//------------------------------------------------------------------------------
#ifndef __Spheral_ViewManager__
#define __Spheral_ViewManager__

#include "Threading/ManagedView.hh"

#include <memory>
#include <vector>
#include <any>

namespace Spheral {
  
template<typename Dimension>
class ViewManager {
public:

  using KeyType = typename StateBase<Dimension>::KeyType;

  // Constructors
  ViewManager() = default;
  ViewManager(const StateBase<Dimension>& state);
  ~ViewManager() = default;

  // Add View types for us to manage
  template<typename ViewType> void enroll(ViewType& view);

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
  std::vector<std::unique_ptr<ManagedViewBase>> mViewPtrs;
  std::vector<std::any> mValueCache;
  const StateBase<Dimension>* mStateBasePtr = nullptr;
};

}

#include "ViewManagerInline.hh"

#endif
