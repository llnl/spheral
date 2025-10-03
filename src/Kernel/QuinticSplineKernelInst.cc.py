specializations = {
    1:
"""
template<>
QuinticSplineKernel< Dim<1> >::QuinticSplineKernel():
  Kernel<Dim<1>, QuinticSplineKernel< Dim<1> > >() {
  setVolumeNormalization(FastMath::pow5(3.0)/40.0);
  setKernelExtent(1.0);
  setInflectionPoint(0.342037); // (2.0/15.0*(7.0 - pow(2.0, 1.0/3.0) - pow(22.0, 2.0/3.0)));
}
""",
    2:
"""
template<>
QuinticSplineKernel< Dim<2> >::QuinticSplineKernel():
  Kernel<Dim<2>, QuinticSplineKernel< Dim<2> > >() {
  setVolumeNormalization(FastMath::pow7(3.0)*7.0/(478.0*M_PI));
  setKernelExtent(1.0);
  setInflectionPoint(0.342037); // (2.0/15.0*(7.0 - pow(2.0, 1.0/3.0) - pow(22.0, 2.0/3.0)));
}
""",
    3:
"""
template<>
QuinticSplineKernel< Dim<3> >::QuinticSplineKernel():
  Kernel<Dim<3>, QuinticSplineKernel< Dim<3> > >() {
  setVolumeNormalization(FastMath::pow7(3.0)/(40.0*M_PI));
  setKernelExtent(1.0);
  setInflectionPoint(0.342037); // (2.0/15.0*(7.0 - pow(2.0, 1.0/3.0) - pow(22.0, 2.0/3.0)));
}
"""
}

text = """
//---------------------------------Spheral++----------------------------------//
// QuinticSplineKernel -- A quintic spline, as described in
// Bonet & Kulasegaruam 2002, Appl. Math. Comput., 126, 133-155.
//
// Kernel extent: 2.0
//
// Created by JMO, Wed Jul  9 16:24:25 PDT 2008
//----------------------------------------------------------------------------//
#include "Kernel/QuinticSplineKernel.cc"

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

namespace Spheral {
template class QuinticSplineKernel< Dim< %(ndim)s > >;
}
"""
