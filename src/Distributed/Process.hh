//---------------------------------------------------------------------------
// \author $Author: mikeowen $
// \version $Revision: 1588 $
// \date $Date: 2005-04-26 17:52:41 -0700 (Tue, 26 Apr 2005) $
//---------------------------------------------------------------------------

// The Process "class" is a set of static methods that give you information 
// about the running process.  So far, most of these methods have to do with 
// the number and rank of running processes in parallel runs.

#ifndef PROCESS_HH
#define PROCESS_HH

#include "config.hh"

namespace Spheral {
class Process
{
   public:

   //! Returns the MPI rank of this process.
   static int getRank();

   //! Returns the total number of MPI processes in this run.
   static int getTotalNumberOfProcesses();

   //! Terminate all processes with prejudice
   static void haltAll (const char* msg);

   private:

#ifdef ENABLE_MPI
   // MPI rank of this process.
   static int sRank;

   // Total number of MPI processes in this run.
   static int sTotalProcs;
#endif

}; // end class Process
}

#endif
