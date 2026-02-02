//---------------------------------Spheral++----------------------------------//
// waitAllWithDeadlockDetection
// 
// Replacement for MPI_Waitall, which automatically spits out after a 
// predetermined time interval any processes that are stuck at the wait.
// Useful for detecting, you know, deadlocks.
//
// Directly cribbed from Jeff Johnson's original code, and cast here as a 
// standalone dohickey.
//----------------------------------------------------------------------------//
#include "Utilities/DBC.hh"

#include <mpi.h>
#ifndef WIN32
#include <unistd.h>
#endif

#include <time.h>
#include <vector>
#include <sstream>
#include <algorithm>

using std::vector;
using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {

void 
waitallWithDeadlockDetection(const string label,
                             const vector<int>& sendProcIDs,
                             const vector<int>& recvProcIDs, 
                             vector<MPI_Request>& sendRequests, 
                             vector<MPI_Request>& recvRequests, 
                             vector<MPI_Status>& sendStatuses,
                             vector<MPI_Status>& recvStatuses,
                             MPI_Comm comm) {
  REQUIRE(sendProcIDs.size() == sendRequests.size());
  REQUIRE(sendProcIDs.size() == sendStatuses.size());
  REQUIRE(recvProcIDs.size() == recvRequests.size());
  REQUIRE(recvProcIDs.size() == recvStatuses.size());

  // Here we can use our own hair-brained deadlock detection if we 
  // want!  This will report deadlocks and may come in handy.
  // If we are caught here for a number of seconds, there's a good chance
  // that we're deadlocked.  So we start a timer here and count up to 
  // N seconds, reporting the state of the communication if we reach
  // that threshhold.  We leave it to the user/developer to decide 
  // what to do from there.
  time_t deadlockInterval = 5; // 5 seconds == possible deadlock.
  time_t startTime = time(0);
  vector<int> sendDone(sendRequests.size(), 0);
  vector<int> recvDone(recvRequests.size(), 0);
  while (count(sendDone.begin(), sendDone.end(), 0) + 
         count(recvDone.begin(), recvDone.end(), 0) != 0) {
    for (size_t i = 0; i != sendRequests.size(); ++i) {
      if (!sendDone[i]) MPI_Test(&sendRequests[i], &sendDone[i], &sendStatuses[i]);
    }
    for (size_t i = 0; i != recvRequests.size(); ++i) {
      if (!recvDone[i]) MPI_Test(&recvRequests[i], &recvDone[i], &recvStatuses[i]);
    }

    // Measure the time that has passed so far.
    time_t timePassed = time(0) - startTime;

    // If it's greater than our deadlock interval, print out a 
    // deadlock warning.
    if (timePassed > deadlockInterval) {
      int rank;
      MPI_Comm_rank(comm, &rank);
      vector<size_t> outstandingRecvs, outstandingSends, completedRecvs, completedSends;
      for (size_t i = 0; i != sendDone.size(); ++i) {
        if (sendDone[i] == 0) {
          outstandingSends.push_back(sendProcIDs[i]);
        } else {
          completedSends.push_back(sendProcIDs[i]);
        }
      }
      for (size_t i = 0; i != recvDone.size(); ++i) {
        if (recvDone[i] == 0) {
          outstandingRecvs.push_back(recvProcIDs[i]);
        } else {
          completedRecvs.push_back(recvProcIDs[i]);
        }
      }

      // Generate the report.  We generate a huge string representing the 
      // report and then print everything out in process order.
      std::stringstream report;

      // Report header: only printed by the 0th process.
      MPI_Barrier(comm);
      if (rank == 0) {
        report << label << "\n";
        report << "\n\n---------------------\n"
               << "MPI DEADLOCK WARNING:\n"
               << "---------------------\n"
               << "Possible deadlock due to the following:\n\n";
      }

      if (!completedSends.empty()) {
        report << rank << ": Successfully sent data to: ";
        for (size_t i = 0; i < completedSends.size(); ++i)
          report << completedSends[i] << ' ';
        report << "\n";
      }

      if (!completedRecvs.empty()) {
        report << rank << ": Successfully received data from: ";
        for (size_t i = 0; i < completedRecvs.size(); ++i)
          report << completedRecvs[i] << ' ';
        report << "\n";
      } // end if

      if (!outstandingSends.empty()) {
        report << rank << ": Still sending data to: ";
        for (size_t i = 0; i < outstandingSends.size(); ++i)
          report << outstandingSends[i] << ' ';
        report << "\n";
      }

      if (!outstandingRecvs.empty()) {
        report << rank << ": Still expecting data from: ";
        for (size_t i = 0; i < outstandingRecvs.size(); ++i)
          report << outstandingRecvs[i] << ' ';
        report << "\n\0";
      }

      // Print stuff out in order.
      int numProcs;
      MPI_Comm_size(comm, &numProcs);
      for (int i = 0; i != numProcs; ++i) {
        if (rank == i) {
          std::cout << report.str();
        } else {
          // Sleep for 1 ms (1000 us).
#if WIN32
           _sleep(1000);
#else
          usleep(1000);
#endif
        }
      }

      // Reset the start time.
      startTime = time(0);
    }
  }
}

}
