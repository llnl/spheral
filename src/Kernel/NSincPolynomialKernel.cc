//---------------------------------Spheral++----------------------------------//
// NSincPolynomialKernel -- The sinc interpolation kernel: W = sin(pi*eta)/(pi*eta).
//
// Created by JMO, Tue Jan  7 15:01:13 PST 2003
//----------------------------------------------------------------------------//

#include "Kernel.hh"
#include "NSincPolynomialKernel.hh"

#include <math.h>
#include <vector>
using std::vector;

namespace Spheral {

//------------------------------------------------------------------------------
// Initialize the given vector<vector> with the the appropriate coefficients.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NSincPolynomialKernel<Dimension>::
setPolynomialCoefficients(const int order,
                          vector< vector<double> >& Aij) const {

  // Check the order, and size the coefficient matrix approppriately.
  REQUIRE(order == 1 || order == 3 || order == 5 || order == 7 || order == 9);
  const int numPolynomials = (order + 1)/2;
  Aij.resize(numPolynomials);
  for (int i = 0; i < numPolynomials; ++i) Aij[i].resize(order + 1);

  // Now, we explicitly set the coefficients for each allowed order.  Kind of 
  // ugly, but since this is seldom called I just don't care.
  double alpha;
  switch(order) {

  case 1:
    // Linear kernel.
    Aij[0][0] = 1.0;
    Aij[0][1] = -1.0;
    break;

  case 3:
    // Cubic kernel.
    alpha = -0.5;
    Aij[0][0] = 1.0;
    Aij[0][1] = 0.0;
    Aij[0][2] = -(alpha + 3.0);
    Aij[0][3] = alpha + 2.0;

    Aij[1][0] = -4.0*alpha;
    Aij[1][1] = 8.0*alpha;
    Aij[1][2] = -5.0*alpha;
    Aij[1][3] = alpha;
    break;

  case 5:
    // Quintic kernel.
    alpha = 3.0/64.0;
    Aij[0][0] = 1.0;
    Aij[0][1] = 0.0;
    Aij[0][2] = 8.0*alpha - 5.0/2.0;
    Aij[0][3] = 0.0;
    Aij[0][4] = -18.0*alpha + 45.0/16.0;
    Aij[0][5] = 10.*alpha - 21.0/16.0;

    Aij[1][0] = -66.0*alpha + 5.0;
    Aij[1][1] = 265.0*alpha - 15.0;
    Aij[1][2] = -392.0*alpha + 35.0/2.0;
    Aij[1][3] = 270.0*alpha - 10.0;
    Aij[1][4] = -88.0*alpha + 45.0/16.0;
    Aij[1][5] = 11.0*alpha - 5.0/16.0;

    Aij[2][0] = -162.0*alpha;
    Aij[2][1] = 297.0*alpha;
    Aij[2][2] = -216.0*alpha;
    Aij[2][3] = 78.0*alpha;
    Aij[2][4] = -14.0*alpha;
    Aij[2][5] = alpha;
    break;

  case 7:
    // Septic kernel.
    alpha = -71.0/83232.0;
    Aij[0][0] = 1.0;
    Aij[0][1] = 0.0;
    Aij[0][2] = -384.0*alpha - 1393.0/578.0;
    Aij[0][3] = 0.0;
    Aij[0][4] = 760.0*alpha + 1960/867.0;
    Aij[0][5] = 0.0;
    Aij[0][6] = -621.0*alpha - 1148.0/867.0;
    Aij[0][7] = 245.0*alpha + 821.0/1734.0;

    Aij[1][0] = -2352.0*alpha - 2233.0/1156.0;
    Aij[1][1] = 14168.0*alpha + 120407.0/6936.0;
    Aij[1][2] = -36000.0*alpha - 13006.0/289.0;
    Aij[1][3] = 47880.0*alpha + 127575.0/2312.0;
    Aij[1][4] = -35640.0*alpha - 128695.0/3468.0;
    Aij[1][5] = 14952.0*alpha + 32683.0/2312.0;
    Aij[1][6] = -3309.0*alpha - 2492.0/867.0;
    Aij[1][7] = 301.0*alpha + 1687.0/6936.0;

    Aij[2][0] = -47280.0*alpha - 8505.0/1156.0;
    Aij[2][1] = 133336.0*alpha + 42525.0/2312.0;
    Aij[2][2] = -157632.0*alpha - 5670.0/289.0;
    Aij[2][3] = 101640.0*alpha + 1575.0/136.0;
    Aij[2][4] = -38720.0*alpha - 4725.0/1156.0;
    Aij[2][5] = 8736.0*alpha + 1995.0/2312.0;
    Aij[2][6] = -1083.0*alpha - 175.0/1734.0;
    Aij[2][7] = 57.0*alpha + 35.0/6936.0;

    Aij[3][0] = -12288.0*alpha;
    Aij[3][1] = 22528.0*alpha;
    Aij[3][2] = -17664.0*alpha;
    Aij[3][3] = 7680.0*alpha;
    Aij[3][4] = -2000.0*alpha;
    Aij[3][5] = 312.0*alpha;
    Aij[3][6] = -27.0*alpha;
    Aij[3][7] = alpha;
    break;
  }
}

}

