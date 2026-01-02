text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Kernel/TableKernel.cc"
#include "Geometry/Dimension.hh"

#include "Kernel/BSplineKernel.hh"
#include "Kernel/W4SplineKernel.hh"
#include "Kernel/GaussianKernel.hh"
#include "Kernel/SuperGaussianKernel.hh"
#include "Kernel/PiGaussianKernel.hh"
#include "Kernel/SincKernel.hh"
#include "Kernel/NSincPolynomialKernel.hh"
#include "Kernel/NBSplineKernel.hh"
#include "Kernel/HatKernel.hh"
#include "Kernel/QuarticSplineKernel.hh"
#include "Kernel/QuinticSplineKernel.hh"
#include "Kernel/WendlandC2Kernel.hh"
#include "Kernel/WendlandC4Kernel.hh"
#include "Kernel/WendlandC6Kernel.hh"
#include "Kernel/ExpInvKernel.hh"

namespace Spheral {
  template class TableKernel< Dim< %(ndim)s > >;
  template class TableKernelView< Dim< %(ndim)s > >;
"""

for Wname in ("BSplineKernel",
              "W4SplineKernel",
              "GaussianKernel",
              "SuperGaussianKernel",
              "PiGaussianKernel",
              "SincKernel",
              "NSincPolynomialKernel",
              "NBSplineKernel",
              "HatKernel",
              "QuarticSplineKernel",
              "QuinticSplineKernel",
              "WendlandC2Kernel",
              "WendlandC4Kernel",
              "WendlandC6Kernel",
              "ExpInvKernel"):
    text += """
  template TableKernel<Dim<%%(ndim)s>>::TableKernel(const %(Wname)s<Dim<%%(ndim)s>>&, const unsigned, const double, const double);
""" % {"Wname" : Wname}

text += """
}
"""
