//---------------------------------Spheral++----------------------------------//
//
// J.M. Pearl 2023
//----------------------------------------------------------------------------//
#include "NodeList/NodeList.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"

#include <vector>
#include <tuple>

#ifdef USE_MPI
#include <mpi.h>
#include "Distributed/Communicator.hh"
#endif

namespace Spheral {


//------------------------------------------------------------------------------
// Compute a unique set of global node IDs for all nodes in the given set of
// NodeLists, returning the result as a FieldList<int>.
//------------------------------------------------------------------------------
template<typename Dimension>
// inline
void
setUniqueNodeIDs(FieldList<Dimension,size_t>& uniqueIndex) {

  // Prepare the result.
  const auto numFields = uniqueIndex.numFields();
  const auto maxUnique = uniqueIndex.max();

#ifdef USE_MPI

  // This processors domain id.
  const int procID = Process::getRank();
  const int numProcs = Process::getTotalNumberOfProcesses();

  // Count up how many nodes have uninitialized unique indices
  size_t numDomainNodes = 0u;
  for (auto fieldi = 0u; fieldi < numFields; ++fieldi) {
    const auto& nodeList = uniqueIndex[fieldi]->nodeList();
    const auto ni = nodeList.numInternalNodes();
    for (auto i = 0u; i < ni; ++i) {
      if(uniqueIndex(fieldi,i)==0u) ++numDomainNodes;
    }
  }
  
  // loop processors and send cumulative number to the next
  auto beginID = maxUnique;
  for (auto sendProc = 0; sendProc < numProcs - 1; ++sendProc) {
    size_t sendProcDomainNodes = numDomainNodes;
    MPI_Bcast(&sendProcDomainNodes, 1, DataTypeTraits<size_t>::MpiDataType(), sendProc, Communicator::communicator());
    if (procID > sendProc) beginID += sendProcDomainNodes;
  }

  // initialize the previously uninitialized
  size_t k = 1u;
  for (auto fieldi = 0u; fieldi < numFields; ++fieldi) {
    const auto& nodeList = uniqueIndex[fieldi]->nodeList();
    const auto ni = nodeList.numInternalNodes();

    for (auto i = 0u; i < ni; ++i){
      if(uniqueIndex(fieldi,i)==0u){
        uniqueIndex(fieldi,i) = beginID + k;
        k++;
      }
    }
  }

#else

  size_t k = 1u;
  for (auto fieldi = 0u; fieldi < numFields; ++fieldi) {
    const auto& nodeList = uniqueIndex[fieldi]->nodeList();
    const auto ni = nodeList.numInternalNodes();

    for (auto i = 0u; i < ni; ++i){
      if(uniqueIndex(fieldi,i)==0u){
        uniqueIndex(fieldi,i) = maxUnique + k;
        k++;
      }
    }
  }


#endif

}

}
