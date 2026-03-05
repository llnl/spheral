//---------------------------------Spheral++----------------------------------//
// Initialize Spheral's use of HYPRE constructs
//----------------------------------------------------------------------------//
#ifndef Spheral_initializeHypre
#define Spheral_initializeHypre

namespace Spheral {

void initializeHypre();   // Set up HYPRE at the start of a Spheral run
void finalizeHypre();     // Finalize HYPRE at the end of a Spheral run

}

#endif
