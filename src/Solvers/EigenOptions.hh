//---------------------------------Spheral++----------------------------------//
// EigenOptions
//
// Holds the options for Eigen solvers and preconditioners
//----------------------------------------------------------------------------//
#ifndef __Spheral_EigenOptions_hh__
#define __Spheral_EigenOptions_hh__

namespace Spheral {

struct EigenOptions {
  EigenOptions() { }

  bool qr = false;
  
}; // end struct EigenOptions

} // end namespace Spheral

#endif
