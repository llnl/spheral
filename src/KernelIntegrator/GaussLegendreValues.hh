//---------------------------------Spheral++----------------------------------//
// GaussLegendreValues
//
// Positive ordinates for the Gauss-Legendre quadrature, accurate to 20 digits
//----------------------------------------------------------------------------//
#ifndef __Spheral_GaussLegendreValues_hh__
#define __Spheral_GaussLegendreValues_hh__

#include <vector>

#include "Geometry/Dimension.hh"

namespace Spheral {

struct GaussLegendreValues {
  static const std::vector<double> values2;
  static const std::vector<double> values4;
  static const std::vector<double> values8;
  static const std::vector<double> values16;
  static const std::vector<double> values32;
  static const std::vector<double> values64;
  static const std::vector<double> values128;
  static const std::vector<double> values256;
  static const std::vector<double> values512;
  static const std::vector<double> values1024;
  // static const std::vector<double> values2048;
  // static const std::vector<double> values4096;

  // Get order and number of ordinates in terms of one another
  static int numOrdinatesForOrder(const int order);
  static int orderForNumOrdinates(const int numOrdinates);
  
  // Get values from above list
  static const std::vector<double>& getValues(const int numOrdinates);
  
  // Get quadrature from -1 to 1 in the first dimension of the vector
  static void getQuadrature(const int numOrdinates,
                            std::vector<double>& weights,
                            std::vector<Dim<1>::Vector>& ordinates);
  static void getQuadrature(const int numOrdinates,
                            std::vector<double>& weights,
                            std::vector<Dim<2>::Vector>& ordinates);
  static void getQuadrature(const int numOrdinates,
                            std::vector<double>& weights,
                            std::vector<double>& ordinates);
};

} // end namespace Spheral

#endif
