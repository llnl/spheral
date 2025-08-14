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
"""

text += """
    template class SuperGaussianKernel<Dim< %(ndim)s > >;
}
"""
