//---------------------------------Spheral++----------------------------------//
// Initialize Spheral's use of HYPRE constructs.
// HYPRE 3.0.0 requires HYPRE_Initialize() before any API call.
//----------------------------------------------------------------------------//
#include "Utilities/initializeHypre.hh"

#ifdef SPHERAL_ENABLE_SOLVERS
#include "HYPRE_utilities.h"
#endif

namespace Spheral {

//------------------------------------------------------------------------------
// Set up HYPRE at the start of Spheral
//------------------------------------------------------------------------------
void initializeHypre() {
#ifdef SPHERAL_ENABLE_SOLVERS
  HYPRE_Initialize();
#endif
}

//------------------------------------------------------------------------------
// Finalize HYPRE at the termination of a Spheral session
//------------------------------------------------------------------------------
void finalizeHypre() {
#ifdef SPHERAL_ENABLE_SOLVERS
  if (HYPRE_Initialized() && !HYPRE_Finalized()) {
    HYPRE_Finalize();
  }
#endif
}

}
