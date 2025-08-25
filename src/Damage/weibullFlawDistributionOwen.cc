//------------------------------------------------------------------------------
// Construct a flaw distribution for a set of nodes according to a Weibull 
// (power-law) distribution.
//------------------------------------------------------------------------------
#include "weibullFlawDistributionOwen.hh"
#include "Utilities/globalNodeIDs.hh"
#include "Utilities/mortonOrderIndices.hh"
#include "NodeList/FluidNodeList.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Strength/SolidFieldNames.hh"
#include "DataBase/State.hh"
#include "Distributed/Communicator.hh"
#include "Distributed/allReduce.hh"

#include <boost/functional/hash.hpp>  // hash_combine

#include <set>
#include <algorithm>
#include <limits>
#include <unordered_map>
#include <random>

using std::unordered_map;
using std::vector;
using std::string;
using std::pair;
using std::make_pair;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {

//------------------------------------------------------------------------------
// This version uses my own algorithm, stochastically seeding flaws in the range
// [0, epsmax] where epsmax is chosen per node based on the nodal volume.
//------------------------------------------------------------------------------
template<typename Dimension>
Field<Dimension, vector<double> >
weibullFlawDistributionOwen(const unsigned seed,
                            const double kWeibull,
                            const double mWeibull,
                            const FluidNodeList<Dimension>& nodeList,
                            const State<Dimension>& state,
                            const size_t minFlawsPerNode,
                            const double volumeMultiplier,
                            const Field<Dimension, int>& mask) {

  // Pre-conditions.
  REQUIRE(kWeibull >= 0.0);
  REQUIRE(mWeibull > 0.0);
  REQUIRE(minFlawsPerNode > 0);
  REQUIRE(mask.nodeListPtr() == &nodeList);

  typedef KeyTraits::Key Key;

  // Prepare the result.
  Field<Dimension, vector<double> > flaws("Weibull flaw distribution",
                                          nodeList);

  // Find the unique ordering to the nodes.
  DataBase<Dimension> db;
  db.appendNodeList(const_cast<FluidNodeList<Dimension>&>(nodeList));
  FieldList<Dimension, Key> keyList = mortonOrderIndices(db);
  const auto nglobal = db.globalNumInternalNodes();
  const auto nlocal = nodeList.numInternalNodes();

  // Is there anything to do?
  if (nglobal > 0) {

    // Identify the rank and number of domains.
    const auto procID = Process::getRank();

    // State for this NodeList.
    auto buildKey = [&](const std::string& fkey) -> std::string { return StateBase<Dimension>::buildFieldKey(fkey, nodeList.name()); };
    const auto& mass = state.field(buildKey(HydroFieldNames::mass), 0.0);
    const auto& rho = (state.registered(buildKey(SolidFieldNames::porositySolidDensity)) ?
                       state.field(buildKey(SolidFieldNames::porositySolidDensity), 0.0) :
                       state.field(buildKey(HydroFieldNames::massDensity), 0.0));

    // Construct a random number generator for each point.
    // Note we hash with the ordering key to generate a unique but reproducible sequence for each point.
    vector<std::mt19937> gens(nlocal);
#pragma omp parallel for
    for (auto i = 0u; i < nlocal; ++i) {
#ifdef __APPLE__
      std::size_t seedi = seed;
#else
      Key seedi = seed;
#endif
      boost::hash_combine(seedi, keyList(0,i));
      gens[i].seed(seedi);
    }
    vector<std::uniform_real_distribution<double>> uniform01(nlocal, std::uniform_real_distribution<double>(0.0, 1.0));

    // Find the minimum and maximum node volumes.
    double Vmin = std::numeric_limits<double>::max(), 
           Vmax = std::numeric_limits<double>::min();
    for (auto i = 0u; i < nlocal; i++) {
      if (mask(i) == 1) {
        const double Vi = mass(i)/rho(i);
        Vmin = min(Vmin, Vi);
        Vmax = max(Vmax, Vi);
      }
    }
    Vmin = allReduce(Vmin*volumeMultiplier, SPHERAL_OP_MIN);
    Vmax = allReduce(Vmax*volumeMultiplier, SPHERAL_OP_MAX);
    CHECK(Vmin > 0.0);
    CHECK(Vmax >= Vmin);

    // Compute the maximum strain we expect for the minimum volume.
    const double epsMax2m = minFlawsPerNode/(kWeibull*Vmin);  // epsmax ** m

    // Based on this compute the maximum number of flaws any node will have.  We'll use this to
    // spin the random number generator without extra communiction.
    const auto maxFlawsPerNode = std::max(1u, unsigned(kWeibull*Vmax*epsMax2m + 0.5));

    // Generate the flaws on each node indepedently.
    const double mInv = 1.0/mWeibull;
    unsigned minNumFlaws = std::numeric_limits<int>::max();
    unsigned maxNumFlaws = 0u;
    unsigned totalNumFlaws = 0u;
    double epsMin = std::numeric_limits<double>::max();
    double epsMax = std::numeric_limits<double>::min();
    double sumFlaws = 0.0;
#pragma omp parallel for
    for (auto i = 0u; i < nlocal; ++i) {
      if (mask(i) == 1) {
        CHECK(rho(i) > 0.0);
        const auto Vi = mass(i)/rho(i) * volumeMultiplier;
        CHECK(Vi > 0.0);
        const auto numFlawsi = std::max(1u, std::min(maxFlawsPerNode, unsigned(kWeibull*Vi*epsMax2m + 0.5)));
        const auto Ai = numFlawsi/(kWeibull*Vi);
        CHECK(Ai > 0.0);
        for (auto j = 0u; j < numFlawsi; ++j) {
          flaws(i).push_back(pow(Ai * uniform01[i](gens[i]), mInv));
        }
        // Sort the flaws on each node by energy.
        sort(flaws(i).begin(), flaws(i).end());
#pragma omp critical
        {
          minNumFlaws = min(minNumFlaws, unsigned(flaws(i).size()));
          maxNumFlaws = max(maxNumFlaws, unsigned(flaws(i).size()));
          totalNumFlaws += flaws(i).size();
          epsMin = min(epsMin, flaws(i).front());
          epsMax = max(epsMax, flaws(i).back());
          for (auto j = 0u; j < flaws(i).size(); ++j) sumFlaws += flaws(i)[j];
        }
      }
    }

    // Some diagnostic output.
    const auto nused = std::max(1, mask.sumElements());
    if (nglobal > 0) {
      minNumFlaws = allReduce(minNumFlaws, SPHERAL_OP_MIN);
      maxNumFlaws = allReduce(maxNumFlaws, SPHERAL_OP_MAX);
      totalNumFlaws = allReduce(totalNumFlaws, SPHERAL_OP_SUM);
      epsMin = allReduce(epsMin, SPHERAL_OP_MIN);
      epsMax = allReduce(epsMax, SPHERAL_OP_MAX);
      sumFlaws = allReduce(sumFlaws, SPHERAL_OP_SUM);
    }
    if (procID == 0) {
      cerr << "weibullFlawDistributionOwen: Min num flaws per node: " << minNumFlaws << endl
           << "                             Max num flaws per node: " << maxNumFlaws << endl
           << "                             Total num flaws       : " << totalNumFlaws << endl
           << "                             Avg flaws per node    : " << totalNumFlaws / nused << endl
           << "                             Min flaw strain       : " << epsMin << endl
           << "                             Max flaw strain       : " << epsMax << endl
           << "                             Avg node failure      : " << sumFlaws / nused << endl;
    }

    // That's it.
    BEGIN_CONTRACT_SCOPE
    {
      for (auto i = 0u; i < nlocal; ++i) {
        if (mask(i) == 1) {
          for (vector<double>::const_iterator itr = flaws(i).begin() + 1;
               itr != flaws(i).end();
               ++itr) ENSURE(*itr >= *(itr - 1));
        }
      }
    }
    END_CONTRACT_SCOPE
  }

  return flaws;
}

}

