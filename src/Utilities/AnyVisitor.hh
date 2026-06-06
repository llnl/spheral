//---------------------------------Spheral++----------------------------------//
// Collect visitor methods to apply to std::any object holders
//
// This allows us to use the visitor pattern with containers of std::any
// obfuscated objects similarly to the std::variant pattern.
//----------------------------------------------------------------------------//
#ifndef __Spheral_AnyVisitor__
#define __Spheral_AnyVisitor__

#include <typeindex>
#include <unordered_map>

namespace Spheral {

template<typename RETURNT, typename... ARGS>
class AnyVisitor {
public:
  using VisitorFunc = std::function<RETURNT (ARGS...)>;

  template<typename T, typename... EXTRAARGS>
  RETURNT visit(T value, EXTRAARGS&&... extraargs) const  {
    auto it = mVisitors.find(std::type_index(value.type()));
    if (it == mVisitors.end()) {
      VERIFY2(false, "AnyVisitor ERROR: unable to process unknown data type " << std::quoted(value.type().name()));
    }
    return it->second(value, extraargs...);
  }

  template<typename T>
  void addVisitor(VisitorFunc visitor) {
    mVisitors[std::type_index(typeid(T))] = visitor;
  }

private:
  std::unordered_map<std::type_index, VisitorFunc> mVisitors;
};

}

#endif
