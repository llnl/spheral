//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "Damage/computeFragmentField.cc"

namespace Spheral {

#if defined(SPHERAL_ENABLE_1D)
template Field<Dim<1>, int> computeFragmentField(const NodeList<Dim<1>>&,
                                                 const double,
                                                 const Field<Dim<1>, Dim<1>::Scalar>&,
                                                 const Field<Dim<1>, Dim<1>::SymTensor>&,
                                                 const Field<Dim<1>, int>&,
                                                 const double,
                                                 const double,
                                                 const bool);
#endif

#if defined(SPHERAL_ENABLE_2D)
template Field<Dim<2>, int> computeFragmentField(const NodeList<Dim<2>>&,
                                                 const double,
                                                 const Field<Dim<2>, Dim<2>::Scalar>&,
                                                 const Field<Dim<2>, Dim<2>::SymTensor>&,
                                                 const Field<Dim<2>, int>&,
                                                 const double,
                                                 const double,
                                                 const bool);
#endif

#if defined(SPHERAL_ENABLE_3D)
template Field<Dim<3>, int> computeFragmentField(const NodeList<Dim<3>>&,
                                                 const double,
                                                 const Field<Dim<3>, Dim<3>::Scalar>&,
                                                 const Field<Dim<3>, Dim<3>::SymTensor>&,
                                                 const Field<Dim<3>, int>&,
                                                 const double,
                                                 const double,
                                                 const bool);
#endif

}
