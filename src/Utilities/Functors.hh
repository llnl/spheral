//------------------------------------------------------------------------------
// This is kinda silly, but we provide an overridable interface to some functors
// here to assist passing functions from Python to C++.
//
// We can probably do something clever with variadic arguments to generalize
// this...
//------------------------------------------------------------------------------
#ifndef __Spheral_Overridable_Functors__
#define __Spheral_Overridable_Functors__

#include "Geometry/Dimension.hh"

namespace Spheral {
namespace PythonBoundFunctors {

// retT F(argT) 
template<typename argT, typename retT>
class SpheralFunctor {
public:
  SpheralFunctor() {};
  virtual ~SpheralFunctor() {};
  virtual retT operator()(const argT& x1) const = 0;
};

// retT F(argT1, argT2) 
template<typename argT1, typename argT2, typename retT>
class Spheral2ArgFunctor {
public:
  Spheral2ArgFunctor() {};
  virtual ~Spheral2ArgFunctor() {};
  virtual retT operator()(const argT1& x1, const argT2& x2) const = 0;
};

// retT F(argT1, argT2, argT3) 
template<typename argT1, typename argT2, typename argT3, typename retT>
class Spheral3ArgFunctor {
public:
  Spheral3ArgFunctor() {};
  virtual ~Spheral3ArgFunctor() {};
  virtual retT operator()(const argT1& x1, const argT2& x2, const argT3& x3) const = 0;
};

// retT F(argT1, argT2, argT3, argT4) 
template<typename argT1, typename argT2, typename argT3, typename argT4, typename retT>
class Spheral4ArgFunctor {
public:
  Spheral4ArgFunctor() {};
  virtual ~Spheral4ArgFunctor() {};
  virtual retT operator()(const argT1& x1, const argT2& x2, const argT3& x3, const argT3& x4) const = 0;
};

}
}

#endif
