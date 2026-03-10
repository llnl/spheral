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

// namespace Detail {

// class Message : public std::ostream {
// public:
//   std::reference_wrapper<std::ostream> os;
//   enum LogLevel { INFO, WARNING, ERROR };
//   Message(LogLevel level = INFO): os(level > INFO ? std::cerr : std::cout) {}
//   template<typename T> std::ostream& operator<< (const T& t) {
//     return (Process::getRank() == 0 ? os.get() << t : os.get());
//   }
// };
// } // Detail

// template<typename T>
// std::ostream& SpheralMessage(T& msg) {
//   Detail::Message result(Detail::Message::INFO);
//   result << msg;
//   return result.os.get();
// }

// template<typename T>
// std::ostream& DeprecationWarning(T& msg) {
//   Detail::Message result(Detail::Message::WARNING);
//   result << "DEPRECATION Warning: ";
//   result << msg;
//   return result.os.get();
// }

#define SpheralMessage(msg)                                          \
  if (Spheral::Process::getRank() == 0)  {                           \
    std::cout << "INFO: " << msg << std::endl;                       \
  }
  
#define SpheralError(msg)                                            \
  if (Spheral::Process::getRank() == 0)  {                           \
    std::cerr << "ERROR: " << msg << std::endl;                      \
  }
  
#define DeprecationWarning(msg)                                      \
  if (Spheral::Process::getRank() == 0)  {                           \
    std::cerr << "DEPRECATION Warning: " << msg << std::endl;        \
  }

#endif

