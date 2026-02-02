#include "Utilities/DomainNode.hh"
#include "Communicator.hh"
#include "Utilities/DBC.hh"
#include "Process.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Get the domain ID.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
RedistributeNodes<Dimension>::domainID() const {
  int domainID = Process::getRank();
  ENSURE(domainID >= 0 && domainID < numDomains());
  return domainID;
}

//------------------------------------------------------------------------------
// Get the number of domains.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
RedistributeNodes<Dimension>::numDomains() const {
  int nProcs = Process::getTotalNumberOfProcesses();
  return nProcs;
}

//------------------------------------------------------------------------------
// Flag controlling how we compute the work.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
RedistributeNodes<Dimension>::computeWork() const {
  return mComputeWork;
}

template<typename Dimension>
inline
void
RedistributeNodes<Dimension>::computeWork(bool x) {
  mComputeWork = x;
}

}
