// This is a simple little method we use to initialize the Tau profiling package.
#include "TAU.h"

#include "Distributed/Process.hh"

namespace Spheral {
inline
void
initializeTau() {
#ifdef PROFILING_ON
  TAU_PROFILE("initializeTau", "", TAU_USER);
  int myid = Process::getRank();
  TAU_PROFILE_SET_NODE(myid);
#endif
}

}
