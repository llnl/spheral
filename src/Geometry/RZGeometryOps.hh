//------------------------------------------------------------------------------
// Provide operations with Geometry primitives useful for the RZ approximation
//
// Created by JMO, Tue May 19 11:28:13 PDT 2026
//------------------------------------------------------------------------------
#ifndef __Spheral_RZGeometryOps__
#define __Spheral_RZGeometryOps__

#include "Geometry/Dimension.hh"

namespace Spheral {

inline
Dim<3>::Tensor operator*(const Dim<3>::SymTensor& A,
                         const Dim<2>::Tensor& B) {
  REQUIRE(A[2] == 0.0 and A[4] == 0.0);  // Check we're consistent with RZ
  return Dim<3>::Tensor(A[0]*B[0] + A[1]*B[2], A[0]*B[1] + A[1]*B[3], 0.0,
                        A[1]*B[0] + A[3]*B[2], A[1]*B[1] + A[3]*B[3], 0.0,
                        0.0,                   0.0,                   0.0);

}

inline
Dim<3>::Tensor operator*(const Dim<2>::Tensor& B,
                         const Dim<3>::SymTensor& A) {
  REQUIRE(A[2] == 0.0 and A[4] == 0.0);  // Check we're consistent with RZ
  return Dim<3>::Tensor(A[0]*B[0] + A[1]*B[1], A[1]*B[0] + A[3]*B[1], 0.0,
                        A[0]*B[2] + A[1]*B[3], A[1]*B[2] + A[3]*B[3], 0.0,
                        0.0,                   0.0,                   0.0);

}

}

#endif
