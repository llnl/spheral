//---------------------------------Spheral++----------------------------------//
// Initialize Spheral's use of Adiak
//----------------------------------------------------------------------------//
#ifndef Spheral_initializeAdiak
#define Spheral_initializeAdiak

#include "config.hh"
#include "Distributed/Communicator.hh"
#include "GPUUtils.hh"
#include "adiak.hpp"

namespace Spheral {

void initializeAdiak() {
  adiak::init((void*) Communicator::comm_ptr());
  // Always collect some curated default adiak information
  adiak::adiakversion();
  adiak::user();
  adiak::uid();
  adiak::launchdate();
  adiak::workdir();
  adiak::hostname();
  adiak::clustername();
  adiak::walltime();
  adiak::cputime();
  adiak::jobsize();
  adiak::numhosts();
  adiak::hostlist();
  adiak::mpi_library_version();
  adiak::value("gpus_per_rank", GPUUtils::deviceCount());
}

enum adiak_categories
  {
   unset = 0,
   all,
   general,
   performance,
   control
  };

}

#endif
