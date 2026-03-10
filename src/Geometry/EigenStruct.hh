//---------------------------------Spheral++----------------------------------//
// EigenStruct, a data structure to hold eigen values and their associated
// eigen vectors.
//
// Created by J. Michael Owen, Thu Dec 28 16:10:05 PST 2000
//----------------------------------------------------------------------------//
#ifndef __Spheral_EigenStruct_hh__
#define __Spheral_EigenStruct_hh__

#include <iostream>

#include "Geometry/GeomVector_fwd.hh"
#include "Geometry/GeomTensor_fwd.hh"
#include "Geometry/GeomSymmetricTensor_fwd.hh"

namespace Spheral {

// EigenStruct
template<int nDim>
struct EigenStruct {
  GeomVector<nDim> eigenValues;
  GeomTensor<nDim> eigenVectors;
};

// Forward declarations.
template<int nDim> std::istream& operator>>(::std::istream is, Spheral::EigenStruct<nDim>& eigen);
template<int nDim> std::ostream& operator<<(::std::ostream os, Spheral::EigenStruct<nDim>& eigen);

}

// Include these here to avoid compiler complaints incomplete types.
#include "GeomVector.hh"
#include "GeomTensor.hh"
#include "EigenStructInline.hh"

#endif
