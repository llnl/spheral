//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "Damage/weibullFlawDistributionBenzAsphaug.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template Field<Dim<1>, std::vector<double>> 
weibullFlawDistributionBenzAsphaug<Dim<1>>(double,
                                           const double,
                                           const unsigned,
                                           const double,
                                           const double,
                                           const FluidNodeList<Dim<1>>&,
                                           const State<Dim<1>>&,
                                           const size_t,
                                           const size_t,
                                           const Field<Dim<1>, int>&);
#endif

#if defined(SPHERAL_ENABLE_2D)
template Field<Dim<2>, std::vector<double>> 
weibullFlawDistributionBenzAsphaug<Dim<2>>(double,
                                           const double,
                                           const unsigned,
                                           const double,
                                           const double,
                                           const FluidNodeList<Dim<2>>&,
                                           const State<Dim<2>>&,
                                           const size_t,
                                           const size_t,
                                           const Field<Dim<2>, int>&);
#endif

#if defined(SPHERAL_ENABLE_3D)
template Field<Dim<3>, std::vector<double>> 
weibullFlawDistributionBenzAsphaug<Dim<3>>(double,
                                           const double,
                                           const unsigned,
                                           const double,
                                           const double,
                                           const FluidNodeList<Dim<3>>&,
                                           const State<Dim<3>>&,
                                           const size_t,
                                           const size_t,
                                           const Field<Dim<3>, int>&);
#endif

}
