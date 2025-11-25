#include "findMatchingVertex.hh"

using std::vector;
using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {

//------------------------------------------------------------------------------
// The the matching vertex between two lists (assumed to be in opposite order.)
//------------------------------------------------------------------------------
unsigned
findMatchingVertex(const unsigned i, 
                   const vector<Dim<2>::Vector>& verticesi,
                   const vector<Dim<2>::Vector>& verticesj) {
  REQUIRE(i < verticesi.size());
  const unsigned ni = verticesi.size(), nj = verticesj.size();
  const unsigned ii = (i + 1) % ni;
  vector<double> chi2;
  chi2.reserve(nj);
  unsigned j, jj;
  for (j = 0; j != nj; ++j) {
    jj = (j == 0 ? (nj - 1) : (j - 1));
    chi2.push_back((verticesi[i] - verticesj[j]).magnitude2() +
                   (verticesi[ii] - verticesj[jj]).magnitude2());
  }
  j = distance(chi2.begin(), min_element(chi2.begin(), chi2.end()));
//   if (!(j < nj)) {
//     cerr << i << endl;
//     for (jj = 0; jj != ni; ++jj) cerr << " " << verticesi[jj];
//     cerr << endl;
//     for (jj = 0; jj != nj; ++jj) cerr << " " << verticesj[jj];
//     cerr << endl;
//   }
  ENSURE2(j < nj, "Unable to find matching vertex");
  return j;
}

//------------------------------------------------------------------------------
// Find the vertex in the second list that best lines up with the first in 
// the first list.  Assumes that the two lists are the same vertices in opposite
// order with arbitrary starting points.
//------------------------------------------------------------------------------
unsigned
findMatchingVertex(const std::vector<Dim<3>::Vector>& verticesi,
                   const std::vector<Dim<3>::Vector>& verticesj,
                   const std::vector<unsigned>& indicesi,
                   const std::vector<unsigned>& indicesj) {
  REQUIRE(indicesi.size() == indicesj.size());
  const unsigned n = indicesi.size();
  unsigned i, j, ji;
  vector<double> chi2;
  chi2.reserve(n);
  for (i = 0; i != n; ++i) {
    chi2.push_back(0.0);
    for (j = 0; j != n; ++j) {
      ji = (j <= i ? (i - j) : (n - (j - i)));
      CHECK(ji < n);
      CHECK(indicesi[j] < verticesi.size());
      CHECK(indicesj[ji] < verticesj.size());
      chi2.back() += (verticesi[indicesi[j]] -
                      verticesj[indicesj[ji]]).magnitude2();
    }
  }
  j = distance(chi2.begin(), min_element(chi2.begin(), chi2.end()));
  ENSURE2(j < n, "Unable to find matching vertex: " << j << " " << n);
  ENSURE2(chi2[j] < 1.0e-5, "Suspiciously large chi2 in findMatchingVertex!  " << chi2[j] << " : " << j);
  return j;
}

}
