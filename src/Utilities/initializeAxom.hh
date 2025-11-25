//---------------------------------Spheral++----------------------------------//
// Initialize Spheral's use of Axom constructs
//----------------------------------------------------------------------------//
#ifndef Spheral_initializeAxom
#define Spheral_initializeAxom

namespace Spheral {

void initializeAxom();   // Set up Axom at the start of a Spheral run
void finalizeAxom();     // Finalize Axom at the end of a Spheral run

}

#endif
