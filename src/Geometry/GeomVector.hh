//---------------------------------Spheral++----------------------------------//
// GeomVector -- Geometric Vector Class.
//----------------------------------------------------------------------------//
#ifndef __Spheral_GeomVector_hh__
#define __Spheral_GeomVector_hh__

#if defined GEOMMEM_EIGEN
#include "Geometry/GeomVector_eigen.hh"
#elif defined GEOMMEM_ARRAY
#include "Geometry/GeomVector_array.hh"
#else
#include "Geometry/GeomVector_default.hh"
#endif
#endif
