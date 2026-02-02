//------------------------------------------------------------------------------
// Construct a flaw distribution for a set of nodes according to a Weibull 
// (power-law) distribution.
//------------------------------------------------------------------------------
#ifndef __Spheral_weibullFlawDistributionBenzAsphaug__
#define __Spheral_weibullFlawDistributionBenzAsphaug__

#include <vector>

// Foward declarations.
namespace Spheral {
  template<typename Dimension> class FluidNodeList;
  template<typename Dimension, typename DataType> class Field;
  template<typename Dimension> class State;
}

namespace Spheral {

//------------------------------------------------------------------------------
// Implements the Benz-Asphaug algorithm, starting with a minimum based
// on the volume of the simulation.
//------------------------------------------------------------------------------
template<typename Dimension>
Field<Dimension, std::vector<double> >
weibullFlawDistributionBenzAsphaug(double volume,
                                   const double volumeStretchFactor,
                                   const unsigned seed,
                                   const double kWeibull,
                                   const double mWeibull,
                                   const FluidNodeList<Dimension>& nodeList,
                                   const State<Dimension>& state,
                                   const size_t minFlawsPerNode,
                                   const size_t minTotalFlaws,
                                   const Field<Dimension, int>& mask);

}

#endif
