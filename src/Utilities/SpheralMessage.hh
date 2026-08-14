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

class InfoStream {
public:
  explicit InfoStream(std::ostream& ss):
    mss(ss) {
  }

  template<typename MsgType>
  InfoStream(std::ostream& ss, MsgType msg):
    InfoStream(ss) {
    if (Process::getRank() == 0) mss << msg;
  }

  template<typename T>
  InfoStream& operator<<(const T& stuff) {
    if (Process::getRank() == 0) mss << stuff;
    return *this;
  }

  InfoStream& operator<<(std::ostream& (*manipulator)(std::ostream&)) {
    if (Process::getRank() == 0) manipulator(mss);
    return *this;
  }

  InfoStream& operator<<(std::ios_base& (*manipulator)(std::ios_base&)) {
    if (Process::getRank() == 0) manipulator(mss);
    return *this;
  }

private:
  std::ostream& mss;
};

}
}

#define SpheralMessage            Spheral::Detail::InfoStream(std::cout) << "INFO: "

#define SpheralWarning            Spheral::Detail::InfoStream(std::cerr) << "WARNING: " << __FILE__ << " " << __LINE__ << " : "

#define SpheralError              Spheral::Detail::InfoStream(std::cerr) << "ERROR: "   << __FILE__ << " " << __LINE__ << " : "

#define SpheralDeprecationWarning Spheral::Detail::InfoStream(std::cout) << "DEPRECATION WARNING: "

#endif

