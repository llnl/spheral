// Helper method to compute the eigen vector corresponding to a given 
// eigen value for a 3x3 symmetric matrix.
#include "buildEigenVector.hh"
#include "Utilities/SpheralFunctions.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

using std::min;
using std::max;
using std::abs;

//------------------------------------------------------------------------------
// This version is the code cooked up by Doug Miller and JMO.
//------------------------------------------------------------------------------
Dim<3>::Vector
buildEigenVector(const Dim<3>::SymTensor& A,
                 const double lambda) {
  typedef Dim<3>::Vector Vector;
  typedef Dim<3>::SymTensor Tensor;

  // First build the matrix we're going to solve.
  const Tensor B(A.xx() - lambda, A.xy(), A.xz(),
                 A.yx(), A.yy() - lambda, A.yz(),
                 A.zx(), A.zy(), A.zz() - lambda);

  BEGIN_CONTRACT_SCOPE
  // Compute an appropriate tolerance for "zero" on the matrix.
  const double tol = 1.0e-3*std::max(1.0, A.diagonalElements().magnitude());
  CONTRACT_VAR(tol);
  REQUIRE(fuzzyEqual(B.Determinant(), 0.0, tol));
  END_CONTRACT_SCOPE

  const double B0 = B.xx();
  const double B1 = B.xy();
  const double B2 = B.xz();
  const double B3 = B.yx();
  const double B4 = B.yy();
  const double B5 = B.yz();
  const double B6 = B.zx();
  const double B7 = B.zy();
  const double B8 = B.zz();

  // Get the three cofactors for the diagonal elements.
  const Vector cofVector(B4*B8 - B5*B7,
                         B0*B8 - B2*B6,
                         B0*B4 - B1*B3);

  // We choose which element to normalize to 1 by the maximum cofactor.
  int i = -1;
  double cof = 0.0;
  for (int j = 0; j < 3; ++j) {
    if (abs(cofVector(j)) > abs(cof)) {
      cof = cofVector(j);
      i = j;
    }
  }
  CHECK(i >= 0 && i < 3);
  CHECK(abs(cof) > 0.0);

  Vector result;
  if (i == 2) {
    const double alpha = (B1*B5 - B2*B4)/cof;
    const double beta = (B2*B3 - B0*B5)/cof;
    const double x3 = 1.0/sqrt(alpha*alpha + beta*beta + 1.0);
    result.x(alpha*x3);
    result.y(beta*x3);
    result.z(x3);
  } else if (i == 1) {
    const double alpha = (B2*B7 - B1*B8)/cof;
    const double beta = (B1*B6 - B0*B7)/cof;
    const double x2 = 1.0/sqrt(alpha*alpha + beta*beta + 1.0);
    result.x(alpha*x2);
    result.y(x2);
    result.z(beta*x2);
  } else {
    const double alpha = (B5*B6 - B3*B8)/cof;
    const double beta = (B3*B7 - B4*B6)/cof;
    const double x1 = 1.0/sqrt(alpha*alpha + beta*beta + 1.0);
    result.x(x1);
    result.y(alpha*x1);
    result.z(beta*x1);
  }
  ENSURE(fuzzyEqual(result.magnitude2(), 1.0, 1.0e-5));
  ENSURE(fuzzyEqual((B*result).magnitude2(), 0.0, 1.0e-5));
  return result;
}

//------------------------------------------------------------------------------
// This version is due to David Eberly: www.geometrictools.com
//------------------------------------------------------------------------------
Dim<3>::Vector
buildUniqueEigenVector(const Dim<3>::SymTensor& A,
                       const double lambda,
                       const Dim<3>::Vector& U0,
                       const Dim<3>::Vector& U1) {

  const double degenerate = 1.0e-20;
  const double tolerance = 1.0e-5;

  typedef Dim<3>::Vector Vector;

  const double p00 = (A*U0).dot(U0) - lambda;
  const double p01 = (A*U0).dot(U1);
  const double p11 = (A*U1).dot(U1) - lambda;
    
  const double x1 = p00*p00 + p01*p01;
  const double x2 = p01*p01 + p11*p11;

  double c0, c1;
  if (x1 >= x2 && x1 > degenerate) {
    const double thpt = 1.0/std::sqrt(x1);
    c0 = thpt*p01;
    c1 = -thpt*p00;
  } else if (x2 > x1 && x2 > degenerate) {
    const double thpt = 1.0/std::sqrt(x2);
    c0 = -thpt*p11;
    c1 = thpt*p01;
  } else {
    // This darn thing is degenerate!
    CHECK(x1 < degenerate && x2 < degenerate);
    if (x1 > x2) {
      c0 = 0.0;
      c1 = 1.0;
    } else {
      c0 = 1.0;
      c1 = 0.0;
    }
  }
  CHECK(fuzzyEqual(c0*c0 + c1*c1, 1.0));
  
  const Vector V0 = c0*U0 + c1*U1;
  CONTRACT_VAR(tolerance);
  ENSURE(fuzzyEqual(V0.magnitude2(), 1.0, tolerance));
  ENSURE(fuzzyEqual(((A - lambda*Dim<3>::SymTensor::one())*V0).maxAbsElement(), 0.0, tolerance));
  return V0;
}

}
