specializations = {
    1:
"""
    template<> const double SuperGaussianKernel<Dim<1> >::mKW = 0.5*3.0;
    template<> const double SuperGaussianKernel<Dim<1> >::mKGW = 0.5*1.0;
""",
    2:
"""
    template<> const double SuperGaussianKernel<Dim<2> >::mKW = 0.5*4.0;
    template<> const double SuperGaussianKernel<Dim<2> >::mKGW = 0.5*2.0;
""",
    3:
"""
    template<> const double SuperGaussianKernel<Dim<3> >::mKW = 0.5*5.0;
    template<> const double SuperGaussianKernel<Dim<3> >::mKGW = 0.5*3.0;
"""
}

text = """
//---------------------------------Spheral++----------------------------------//
// SuperGaussianKernel -- The super gaussian interpolation kernel.
//
// Created by JMO, Wed Dec  1 22:38:21 PST 1999
//----------------------------------------------------------------------------//
#include "Kernel/SuperGaussianKernel.hh"

#include <math.h>

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
namespace Spheral {
    template class SuperGaussianKernel<Dim< %(ndim)s > >;
}
"""
