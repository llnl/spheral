//------------------------------------------------------------------------------
// timingUtilities
//
// A set of inline helper methods to encapsulate how we do timing.
//
// JMO:  Tue Dec  9 10:31:14 PST 2008
//------------------------------------------------------------------------------
#ifndef __Spheral__timingUtilities__
#define __Spheral__timingUtilities__

#include <chrono>   // C++11

namespace Spheral {
//------------------------------------------------------------------------------
// Get the current clock time.
//------------------------------------------------------------------------------
struct Timing {
  typedef std::chrono::time_point<std::chrono::high_resolution_clock> Time;
  typedef std::chrono::microseconds duration;
  // typedef std::chrono::duration<double, std::chrono::high_resolution_clock> duration;

  static Time currentTime() { return std::chrono::high_resolution_clock::now(); }
  static double difference(const Time& t1, const Time& t2) { return (std::chrono::duration_cast<duration>(t2 - t1)).count(); }

  static duration zero() { return duration(0); }
};

} // namespace Spheral

#endif
