//---------------------------------Spheral++----------------------------------//
// Communicator
// A singleton object which holds the MPI communicator Spheral uses.
//
// Created by JMO, Wed Jan 30 14:52:05 PST 2013
//----------------------------------------------------------------------------//
#include "Communicator.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Default constructor (private).
//------------------------------------------------------------------------------
Communicator::Communicator() {
#ifdef SPHERAL_ENABLE_MPI
  mCommunicator = MPI_COMM_WORLD;
#else
  mCommunicator = 0;
#endif
}

//------------------------------------------------------------------------------
// Destructor (private).
//------------------------------------------------------------------------------
Communicator::~Communicator() {
}

//------------------------------------------------------------------------------
// Public routines
//------------------------------------------------------------------------------

MPI_Comm* Communicator::comm_ptr() {
#ifdef SPHERAL_ENABLE_MPI
  return &(instance().mCommunicator);
#else
  return nullptr;
#endif
}

void Communicator::finalize() {
#ifdef SPHERAL_ENABLE_MPI
  int finalized = 0;
  MPI_Finalized(&finalized);
  if (finalized != 0) {
    int finalize = MPI_Finalize();
    if (finalize != 0) {
      char* string = nullptr;
      int resultlen = 0;
      MPI_Error_string(finalize, string, &resultlen);
      VERIFY2(finalize, string);
    }
  }
#endif
}
}
