//------------------------------------------------------------------------------
// Construct a flaw distribution for a set of nodes according to a Weibull 
// (power-law) distribution.
//------------------------------------------------------------------------------
#include <set>
#include <algorithm>
#include <limits>
#include <unordered_map>
#include <random>

#include "weibullFlawDistributionBenzAsphaug.hh"
#include "Utilities/globalNodeIDs.hh"
#include "Utilities/mortonOrderIndices.hh"
#include "Utilities/nodeOrdering.hh"
#include "NodeList/FluidNodeList.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Strength/SolidFieldNames.hh"
#include "DataBase/State.hh"
#include "Distributed/Communicator.hh"
#include "Distributed/allReduce.hh"

using std::unordered_map;
using std::vector;
using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {

//------------------------------------------------------------------------------
// This version uses the Benz-Asphaug algorithm, stepping up deterministically
// in flaw energy based on a minimum chosen from the simulation volume.
//------------------------------------------------------------------------------
template<typename Dimension>
Field<Dimension, vector<double> >
weibullFlawDistributionBenzAsphaug(double volume,
                                   const double volumeStretchFactor,
                                   const unsigned seed,
                                   const double kWeibull,
                                   const double mWeibull,
                                   const FluidNodeList<Dimension>& nodeList,
                                   const State<Dimension>& state,
                                   const size_t minFlawsPerNode,
                                   const size_t minTotalFlaws,
                                   const Field<Dimension, int>& mask) {

  // Pre-conditions.
  REQUIRE(volume >= 0.0);
  REQUIRE(volumeStretchFactor >= 1.0);
  REQUIRE(kWeibull >= 0.0);
  REQUIRE(mWeibull > 0.0);
  REQUIRE(minFlawsPerNode > 0);
  REQUIRE(minTotalFlaws > 0);
  REQUIRE(mask.nodeListPtr() == &nodeList);

  // Prepare the result.
  Field<Dimension, vector<double> > flaws("Weibull flaw distribution",
                                          nodeList);

  // Construct unique global IDs for all nodes in the NodeList.
  const auto n = max(size_t(1u), numGlobalNodes(nodeList));
  const Field<Dimension, size_t> globalIDs = globalNodeIDs(nodeList);

  // Prepare a table to faciliate looking local IDs from global.
  unordered_map<size_t, size_t> global2local;
  const auto nlocal = nodeList.numInternalNodes();
#pragma omp parallel for
  for (auto i = 0u; i < nlocal; i++) global2local[globalIDs(i)] = i;
  CHECK(global2local.size() == nodeList.numInternalNodes());

  // Prepare an int per *each* node, so that each process can keep track of how many
  // flaws are seeded globally (avoiding communication at the expense of memory).
  vector<size_t> numFlawsPerNode((size_t) n, 0);

  // Identify the rank and number of domains.
  const auto procID = Process::getRank();

  // Only proceed if there are nodes to initialize!
  if (n > 0) {

    // If the user did not speicify a volume, we compute it from the information
    // in the State.
    if (volume == 0.0) {
      auto buildKey = [&](const std::string& fkey) -> std::string { return StateBase<Dimension>::buildFieldKey(fkey, nodeList.name()); };
      const auto& mass = state.field(buildKey(HydroFieldNames::mass), 0.0);
      const auto& rho = (state.registered(buildKey(SolidFieldNames::porositySolidDensity)) ?
                         state.field(buildKey(SolidFieldNames::porositySolidDensity), 0.0) :
                         state.field(buildKey(HydroFieldNames::massDensity), 0.0));
#pragma omp parallel for
      for (auto i = 0u; i < nlocal; i++) {
        CHECK(rho(i) > 0.0);
        volume += mass(i)/rho(i);
      }
      volume = allReduce(volume, SPHERAL_OP_SUM);
    }
    volume = std::max(volume, 1e-100);
    CHECK(volume > 0.0);

    // Compute the minimum (starting) failure strain.
    const double mInv = 1.0/(mWeibull + 1.0e-50);
    const double epsMin = pow(kWeibull*volume*volumeStretchFactor, -mInv);
    CHECK(epsMin > 0.0);

    // Construct a random number generator.
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> uniform01(0.0, 1.0);

    // Loop and initialize flaws until:
    // a) every node has the minimum number of flaws per node, and
    // b) we meet the minimum number of total flaws.
    size_t numCompletedNodes = 0;
    size_t ienergy = 1;
    while ((numCompletedNodes < n) || (ienergy <= minTotalFlaws)) {

      // Randomly select a global node.
      const auto iglobal = size_t(uniform01(gen) * n);
      CHECK(iglobal < n);

      // Increment the number of flaws for this node, and check if this
      // completes this node.
      ++numFlawsPerNode[iglobal];
      if (numFlawsPerNode[iglobal] == minFlawsPerNode) ++numCompletedNodes;

      // Is this node one of ours?
      const auto itr = global2local.find(iglobal);
      if (itr != global2local.end()) {

        const unsigned i = itr->second;
        CHECK(i < nodeList.numInternalNodes());
        if (mask(i) == 1) {

          // The activation energy.
          const double epsij = epsMin * pow(ienergy*volumeStretchFactor, mInv);

          // Add a flaw with this activation energy to this node.
          flaws(i).push_back(epsij);
        }
      }

      // Increment the energy multiplier.
      ++ienergy;
    }

    // Sort the flaws on each node by energy.
    unsigned minNumFlaws = INT_MAX;
    unsigned maxNumFlaws = 0;
    unsigned totalNumFlaws = 0;
    double epsMax = 0.0;
    double sumFlaws = 0.0;
    for (auto i = 0u; i < nodeList.numInternalNodes(); ++i) {
      minNumFlaws = min(minNumFlaws, unsigned(flaws(i).size()));
      maxNumFlaws = max(maxNumFlaws, unsigned(flaws(i).size()));
      totalNumFlaws += flaws(i).size();
      if (mask(i) == 1) {
        sort(flaws(i).begin(), flaws(i).end());
        epsMax = max(epsMax, flaws(i).back());
        for (auto j = 0u; j != flaws(i).size(); ++j) sumFlaws += flaws(i)[j];
      }
    }

    // Prepare some diagnostic output.
    const auto nused = std::max(1, mask.sumElements());
    minNumFlaws = allReduce(minNumFlaws, SPHERAL_OP_MIN);
    maxNumFlaws = allReduce(maxNumFlaws, SPHERAL_OP_MAX);
    totalNumFlaws = allReduce(totalNumFlaws, SPHERAL_OP_SUM);
    epsMax = allReduce(epsMax, SPHERAL_OP_MAX);
    sumFlaws = allReduce(sumFlaws, SPHERAL_OP_SUM);
    if (procID == 0) {
      cerr << "weibullFlawDistributionBenzAsphaug: Min num flaws per node: " << minNumFlaws << endl
           << "                                    Max num flaws per node: " << maxNumFlaws << endl
           << "                                    Total num flaws       : " << totalNumFlaws << endl
           << "                                    Avg flaws per node    : " << totalNumFlaws / nused << endl
           << "                                    Min flaw strain       : " << epsMin << endl
           << "                                    Max flaw strain       : " << epsMax << endl
           << "                                    Avg node failure      : " << sumFlaws / nused << endl;
    }
  }

  // That's it.
  BEGIN_CONTRACT_SCOPE
  {
    for (int i = 0; i != (int)nodeList.numInternalNodes(); ++i) {
      if (mask(i) == 1) {
        ENSURE(flaws(i).size() >= minFlawsPerNode);
        for (auto itr = flaws(i).begin() + 1;
             itr != flaws(i).end();
             ++itr) ENSURE(*itr >= *(itr - 1));
      }
    }
  }
  END_CONTRACT_SCOPE

  return flaws;
}

}

