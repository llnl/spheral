//---------------------------------Spheral++----------------------------------//
// Initialize Spheral's use of Axom constructs
//----------------------------------------------------------------------------//
#include "Utilities/initializeAxom.hh"
#include "axom/slic.hpp"

namespace Spheral {

//------------------------------------------------------------------------------
// Set up Axom at the start of Spheral
//------------------------------------------------------------------------------
void initializeAxom() {
  axom::slic::initialize();
}


//------------------------------------------------------------------------------
// Finalize Axom at the termination of a Spheral session
//------------------------------------------------------------------------------
void finalizeAxom() {
  axom::slic::finalize();
}

}
