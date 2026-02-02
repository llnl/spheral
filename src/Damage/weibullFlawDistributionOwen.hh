//------------------------------------------------------------------------------
// Construct a flaw distribution for a set of nodes according to a Weibull 
// (power-law) distribution.
//------------------------------------------------------------------------------
#ifndef __Spheral_weibullFlawDistribution__
#define __Spheral_weibullFlawDistribution__

#include <vector>

// Foward declarations.
namespace Spheral {
  template<typename Dimension> class FluidNodeList;
  template<typename Dimension, typename DataType> class Field;
  template<typename Dimension> class State;
}

namespace Spheral {

//------------------------------------------------------------------------------
// Implements the Owen algorithm, stochastically seeding flaws with a maximum
// value per node chosen based on the volume of the node.
//------------------------------------------------------------------------------
template<typename Dimension>
Field<Dimension, std::vector<double> >
weibullFlawDistributionOwen(const unsigned seed,
                            const double kWeibull,
                            const double mWeibull,
                            const FluidNodeList<Dimension>& nodeList,
                            const State<Dimension>& state,
                            const size_t minFlawsPerNode,
                            const double volumeMultiplier,
                            const Field<Dimension, int>& mask);

}

#endif
