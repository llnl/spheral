//------------------------------------------------------------------------------
// ViewManager
//
// Manage the movement for a collection of View objects between different
// memory spaces (GPU->CPU, CPU->GPU, etc)
//------------------------------------------------------------------------------

namespace Spheral {

//------------------------------------------------------------------------------
// Construct with a StateBase
//------------------------------------------------------------------------------
template<typename Dimension>
inline
ViewManager<Dimension>::ViewManager(const StateBase<Dimension>& state):
  mViews(),
  mValueCache(),
  mStateBasePtr(&state) {
}

//------------------------------------------------------------------------------
// Enroll a View
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename T>
inline
void
ViewManager<Dimension>::enroll(T& view) {
  mViews.emplace_back(std::ref(view));
}

//------------------------------------------------------------------------------
// Get a FieldView from the StateBase
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename Value>
inline
FieldView<Dimension, Value> 
ViewManager<Dimension>::field(const KeyType& key) {
  REQUIRE2(mStateBasePtr != nullptr,
           "ViewManager ERROR: attempt to extract field " << key << " from null State pointer");
  auto& f = mStateBasePtr->template field<Value>(key);
  auto fv = f.view();
  this->enroll(fv);
  return fv;
}

template<typename Dimension>
template<typename Value>
inline
FieldView<Dimension, Value> 
ViewManager<Dimension>::field(const KeyType& key,
                              const Value& dummy) {
  return this->template field<Value>(key);
}

//------------------------------------------------------------------------------
// Get a FieldListView from the StateBase
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename Value>
inline
FieldListView<Dimension, Value> 
ViewManager<Dimension>::fields(const KeyType& key,
                               bool allowNone) {
  REQUIRE2(mStateBasePtr != nullptr,
           "ViewManager ERROR: attempt to extract fields " << key << " from null State pointer");
  auto fl0 = mStateBasePtr->template fields<Value>(key, allowNone);
  mValueCache.push_back(fl0);
  auto& fl = std::any_cast<FieldList<Dimension, Value>&>(mValueCache.back());
  auto flv = fl.view();
  this->enroll(flv);
  return flv;
}

template<typename Dimension>
template<typename Value>
inline
FieldListView<Dimension, Value> 
ViewManager<Dimension>::fields(const KeyType& key,
                               const Value& dummy,
                               bool allowNone) {
  return this->template fields<Value>(key, allowNone);
}

//------------------------------------------------------------------------------
// Get an arbitrary type from the StateBase
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename Value, typename View>
inline
View
ViewManager<Dimension>::get(const KeyType& key) {
  REQUIRE2(mStateBasePtr != nullptr,
           "ViewManager ERROR: attempt to extract type for key " << key << " from null State pointer");
  auto& val = mStateBasePtr->template get<Value>(key);
  auto view = val.view();
  this->enroll(view);
  return view;
}

//------------------------------------------------------------------------------
// move all views
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
ViewManager<Dimension>::move(chai::ExecutionSpace space) {
  for (auto& v: mViews) {
    std::visit([&](auto&& v_ref_wrapper) { v_ref_wrapper.get().move(space); }, v);
  }
}

//------------------------------------------------------------------------------
// touch all views
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
ViewManager<Dimension>::touch(chai::ExecutionSpace space) {
  for (auto& v: mViews) {
    std::visit([&](auto&& v_ref_wrapper) { v_ref_wrapper.get().touch(space); }, v);
  }
}

}
