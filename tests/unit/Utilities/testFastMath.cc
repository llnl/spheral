#include <cmath>
#include <limits>
#include <iostream>
#include <random>
#include <time.h>
#include <vector>

#include "Utilities/FastMath.hh"
#include "Utilities/DBC.hh"

using namespace Spheral;

int main(int argc, char* argv[]) {

  // Make some random number generators.
  // Doubles.
  std::random_device rd;
  std::mt19937 generator(rd());
  std::uniform_real_distribution<double> unitDoubleDistribution(0.0, 1.0);
  std::uniform_real_distribution<double> tinyDoubleDistribution(0.0, 1.0e-20);
  std::uniform_real_distribution<double> hugeDoubleDistribution(0.0, sqrt(numeric_limits<double>::max()));

  // Make a big vector of doubles.
  const int n = 10000;
  cerr << "Generating random input..." << endl;
  vector<double> xinput(3*n);
  for (int i = 0; i != n; ++i) xinput[i] = unitDoubleDistribution(generator);
  for (int i = n; i != 2*n; ++i) xinput[i] = tinyDoubleDistribution(generator);
  for (int i = 2*n; i != 3*n; ++i) xinput[i] = hugeDoubleDistribution(generator);

  //===========================================================================
  // Check the correctness of sqrt.
  //===========================================================================
  cerr << "Testing sqrt correctness." << endl;
  const double tol = 1.0e-7;
  for (int i = 0; i != xinput.size(); ++i) {
    if (not (fuzzyEqual(FastMath::SqrtHalley2(xinput[i]),
                        sqrt(xinput[i]),
                        tol))) {
      const double a = FastMath::nth_root<2>(xinput[i]);
      const double a2 = a*a;
      const double R = xinput[i];
      cerr << i << " " << xinput[i] << " : "
           << FastMath::nth_root<2>(xinput[i]) << " "
           << FastMath::sqrta_halley(FastMath::nth_root<2>(xinput[i]), xinput[i]) << " : "
           << a2 << " "
           << a2 + R + R + R << " " 
           << a*(a2 + R + R + R) << " "
           << a2 + a2 + a2 + R + 1.0e-50 << " "
           << a*(a2 + R + R + R)/(a2 + a2 + a2 + R + 1.0e-50) << " : "
           << FastMath::SqrtHalley2(xinput[i]) << " != " << sqrt(xinput[i]) << endl;
    }
    VERIFY(fuzzyEqual(FastMath::SqrtHalley2(xinput[i]),
                      sqrt(xinput[i]),
                      tol));
  }

  //===========================================================================
  // Check the correctness of Sqrt
  //===========================================================================
//   cerr << "Testing Sqrt correctness." << endl;
//   for (int i = 0; i != xinput.size(); ++i) {
//     if (not fuzzyEqual(double(FastMath::Sqrt(xinput[i])),
//                        sqrt(xinput[i]),
//                        tol)) {
//       cerr << FastMath::Sqrt(xinput[i]) << " "
//            << sqrt(xinput[i]) << " "
//            << endl;
//     }
//     VERIFY(fuzzyEqual(double(FastMath::Sqrt(xinput[i])),
//                       sqrt(xinput[i]),
//                       tol));
//   }

  //===========================================================================
  // Check the correctness of cube root.
  //===========================================================================
  cerr << "Testing cube root correctness." << endl;
  const double onethird = 1.0/3.0;
  for (int i = 0; i != xinput.size(); ++i) {
    if (not (fuzzyEqual(FastMath::CubeRootHalley2(xinput[i]),
                        pow(xinput[i], onethird),
                        tol)))
      cerr << i << " " << xinput[i] << " : "
           << FastMath::CubeRootHalley2(xinput[i]) << " != " << pow(xinput[i], onethird) << endl;
    VERIFY(fuzzyEqual(FastMath::CubeRootHalley2(xinput[i]),
                      pow(xinput[i], onethird),
                      tol));
  }

  //===========================================================================
  // Time the sqrt.
  //===========================================================================
  cerr << "Timing standard sqrt." << endl;
  const unsigned ntime = 1000000000;
  clock_t tsqrt;
  {
    const clock_t t0 = clock();
    double sum = 0.0;
    for (unsigned i = 0; i != ntime; ++i) sum += sqrt(double(i));
    const clock_t t1 = clock();
    tsqrt = t1 - t0;
    cerr << t0 << " " << t1 << " " << sum << endl;
  }

  // Time FastMath::SqrtHalley2
  clock_t tFMSqrtHalley2;
  {
    const clock_t t0 = clock();
    double sum = 0.0;
    for (unsigned i = 0; i != ntime; ++i) sum += FastMath::SqrtHalley2(double(i));
    const clock_t t1 = clock();
    tFMSqrtHalley2 = t1 - t0;
    cerr << t0 << " " << t1 << " " << sum << endl;
  }

  // Time FastMath::Sqrt
  clock_t tFMSqrt;
  {
    const clock_t t0 = clock();
    double sum = 0.0;
    for (unsigned i = 0; i != ntime; ++i) sum += FastMath::Sqrt(i);
    const clock_t t1 = clock();
    tFMSqrt = t1 - t0;
    cerr << t0 << " " << t1 << " " << sum << endl;
  }

  // Report sqrt timing.
  cerr << "Sqrt time (sqrt, FastMath::SqrtHalley2, ratio) : "
       << double(tsqrt)/CLOCKS_PER_SEC << " sec, "
       << double(tFMSqrtHalley2)/CLOCKS_PER_SEC << " sec, "
       << double(tsqrt)/double(tFMSqrtHalley2)
       << endl;
  cerr << "Sqrt time (sqrt, FastMath::Sqrt, ratio) : "
       << double(tsqrt)/CLOCKS_PER_SEC << " sec, "
       << double(tFMSqrt)/CLOCKS_PER_SEC << " sec, "
       << double(tsqrt)/double(tFMSqrt)
       << endl;

  // We are assuming that the sqrt function with all modern compilers is already optimized
  // better than our simple method here, so check that it beat our local version.
//   VERIFY(tsqrtratio < 1.0);

  //===========================================================================
  // Time the cube root.
  //===========================================================================
  cerr << "Timing standard pow(1/3)." << endl;
  clock_t tpow13;
  {
    const clock_t t0 = clock();
    double sum = 0.0;
    for (unsigned i = 0; i != ntime; ++i) sum += pow(double(i), onethird);
    const clock_t t1 = clock();
    tpow13 = t1 - t0;
    cerr << t0 << " " << t1 << " " << sum << endl;
  }

  // Time FastMath::CubeRootHalley2
  clock_t tFMcuberoot;
  {
    const clock_t t0 = clock();
    double sum = 0.0;
    for (unsigned i = 0; i != ntime; ++i) sum += FastMath::CubeRootHalley2(double(i));
    const clock_t t1 = clock();
    tFMcuberoot = t1 - t0;
    cerr << t0 << " " << t1 << " " << sum << endl;
  }

  // Report cube root timing.
  const double tcuberootratio = double(tpow13)/double(tFMcuberoot);
  cerr << "Cube root time (pow(1/3), FastMath::CubeRootHalley2, ratio) : "
       << double(tpow13)/CLOCKS_PER_SEC << " sec, "
       << double(tFMcuberoot)/CLOCKS_PER_SEC << " sec, "
       << tcuberootratio 
       << endl;

  VERIFY(tcuberootratio > 1.0);

  return 0;
}
