#ifndef __Spheral_Silo_vectorstringutils__
#define __Spheral_Silo_vectorstringutils__

#include <string>
#include <sstream>
#include <vector>
#include "Geometry/Dimension.hh"

#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Generic vector -> string conversion.
//------------------------------------------------------------------------------
template<typename T>
inline
std::string
vector2string(const std::vector<T>& val,
              const int precision) {
  std::ostringstream ss;
  ss.precision(precision);
  ss << val.size() << " ";
  for (typename std::vector<T>::const_iterator itr = val.begin();
       itr != val.end();
       ++itr) ss << *itr << " ";
  ss << std::ends;
  return ss.str();
}

//------------------------------------------------------------------------------
// Generic string -> vector conversion.
//------------------------------------------------------------------------------
template<typename T>
inline
std::vector<T>
string2vector(const std::string& val) {
  std::istringstream ss(val);
  int size;
  ss >> size;
  std::vector<T> result;
  result.reserve(size);
  T x;
  while (ss >> x) result.push_back(x);
  VERIFY((int)result.size() == size);
  return result;
}

//------------------------------------------------------------------------------
// Declare specializations for strings.
//------------------------------------------------------------------------------
template<> std::string vector2string(const std::vector<std::string>& val, const int precision);
template<> std::vector<std::string> string2vector<std::string>(const std::string& val);

}

#endif
