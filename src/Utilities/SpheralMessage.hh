//---------------------------------Spheral++----------------------------------//
// SpheralMessage
//
// Utilities to facilitate messages and warnings.
//
// Created by JMO, Mon Dec 16 15:20:42 PST 2024
//----------------------------------------------------------------------------//
#ifndef __Spheral_SpheralMessage__
#define __Spheral_SpheralMessage__

#include "Distributed/Process.hh"

#include <iostream>
#include <sstream>

namespace Spheral{
namespace Detail {

inline void emit_message(std::ostream& ss,
                         const char* level,
                         const char* file,
                         const int line,
                         const std::ostringstream& os) {
  ss << level << " [" << file << ":" << line << "]: " << os.str() << '\n';
}

inline void emit_message(std::ostream& ss,
                         const char* level,
                         const std::ostringstream& os) {
  ss << level << ": " << os.str() << '\n';
}

}
}

#define BUILD_SPHERAL_MSG_STREAM(expr)                               \
    ([&]() {                                                         \
       std::ostringstream os;                                        \
       os << expr;                                                   \
       return os;                                                    \
    }())

#define SpheralMessage(expr) do {                                       \
  if (Spheral::Process::getRank() == 0)  {                              \
    Spheral::Detail::emit_message(std::cout, "INFO",                    \
                                  BUILD_SPHERAL_MSG_STREAM(expr));      \
  }                                                                     \
} while(0)
  
#define SpheralWarning(expr) do {                                       \
  if (Spheral::Process::getRank() == 0)  {                              \
    Spheral::Detail::emit_message(std::cerr, "WARNING",                 \
                                  __FILE__, __LINE__,                   \
                                  BUILD_SPHERAL_MSG_STREAM(expr));      \
  }                                                                     \
} while(0)

#define SpheralError(expr) do {                                         \
  if (Spheral::Process::getRank() == 0)  {                              \
    Spheral::Detail::emit_message(std::cout, "ERROR",                   \
                                  __FILE__, __LINE__,                   \
                                  BUILD_SPHERAL_MSG_STREAM(expr));      \
  }                                                                     \
} while(0)
  
#define SpheralDeprecationWarning(expr) do {                            \
  if (Spheral::Process::getRank() == 0)  {                              \
    Spheral::Detail::emit_message(std::cout, "DEPRECATION WARNING",     \
                                  BUILD_SPHERAL_MSG_STREAM(expr));      \
  }                                                                     \
} while(0)

#endif

