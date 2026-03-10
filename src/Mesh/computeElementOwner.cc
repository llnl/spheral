//---------------------------------Spheral++----------------------------------//
// computeElementOwner
//
// Utility method to reduce the shared elements to the lowest common processor.
//
// Created by JMO, Mon Oct 29 22:01:30 PDT 2012
//----------------------------------------------------------------------------//
#include "Utilities/DBC.hh"
#include "Utilities/SpheralFunctions.hh"

#include "computeElementOwner.hh"
#include <algorithm>

using std::vector;
using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {

vector<unsigned>
computeElementOwner(const unsigned numElements,
                    const vector<unsigned>& neighborDomains,
                    const vector<vector<unsigned> >& sharedElements) {
  REQUIRE(neighborDomains.size() == sharedElements.size());

  const unsigned rank = Process::getRank();
  const unsigned numProcs = Process::getTotalNumberOfProcesses();

  // Prepare the result.
  vector<unsigned> result(numElements, rank);

  // There's really only work if we have more than one domain.
  if (numProcs > 1) {
    const unsigned numDomains = neighborDomains.size();
    for (unsigned idomain = 0; idomain != numDomains; ++idomain) {
      const unsigned neighbor = neighborDomains[idomain];
      for (const unsigned element: sharedElements[idomain]) {
        CHECK(element < numElements);
        result[element] = std::min(result[element], neighbor);
      }
    }
  }

  return result;
}

}
