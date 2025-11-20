//---------------------------------Spheral++----------------------------------//
// GeomTensor -- the full tensor class.
//----------------------------------------------------------------------------//
#include <cmath>
#include <limits.h>
#include <float.h>

#include "GeomTensor.hh"
#include "findEigenValues3.hh"
#include "Utilities/SpheralFunctions.hh"
#include "Utilities/rotationMatrix.hh"


namespace Spheral {

//------------------------------------------------------------------------------
// Find the eigenvalues of a tensor.
//------------------------------------------------------------------------------
template<>
GeomVector<1>
GeomTensor<1>::eigenValues() const {
  return GeomVector<1>(this->mxx);
}

//----------------------------------------------------------------------
template<>
GeomVector<2>
GeomTensor<2>::eigenValues() const {
  const double b = Trace();
  const double c = Determinant();
  const double q = 0.5*(b + sgn(b)*sqrt(std::max(0.0, b*b - 4.0*c))) + 1.0e-10*sgn(b);
  CHECK(q != 0.0);
  return GeomVector<2>(q, c/q);
}

//----------------------------------------------------------------------
template<>
GeomVector<3>
GeomTensor<3>::eigenValues() const {
  return findEigenValues3<GeomTensor<3> >(*this);
}

}
