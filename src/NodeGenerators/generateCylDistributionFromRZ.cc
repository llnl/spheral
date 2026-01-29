//------------------------------------------------------------------------------
// Helper method for the GenerateCylindricalNodeDistribution3d node generator
// to generate the spun node distribution.
//------------------------------------------------------------------------------

#include "Boundary/CylindricalBoundary.hh"
#include "Utilities/DBC.hh"
#include "Geometry/Dimension.hh"
#include "Distributed/allReduce.hh"

#include <vector>
#include <algorithm>

using std::vector;

namespace Spheral {

void
generateCylDistributionFromRZ(vector<double>& x,
                              vector<double>& y,
                              vector<double>& z,
                              vector<double>& m,
                              vector<Dim<3>::SymTensor>& H,
                              vector<int>& globalIDs,
                              vector<vector<double> >& extraFields,
                              const double nNodePerh,
                              const double kernelExtent,
                              const double phi,
                              const int procID,
                              const int nProcs) {

  using Vector = Dim<3>::Vector;
  using SymTensor = Dim<3>::SymTensor;

  // Pre-conditions.
  const auto n = x.size();
  const auto nextra = extraFields.size();
  VERIFY(y.size() == n and
         z.size() == n and
         m.size() == n and
         H.size() == n);
  for (auto i = 0u; i < nextra; ++i) VERIFY(extraFields[i].size() == n);
  VERIFY(size_t(std::count(z.begin(), z.end(), 0.0)) == n);
  VERIFY(nNodePerh > 0.0);
  VERIFY(kernelExtent > 0.0);
  VERIFY(phi > 0.0);
  VERIFY(nProcs >= 1);
  VERIFY(procID >= 0 and procID < nProcs);

  // Make an initial pass to determine how many nodes we're going
  // to generate.
  size_t ntot = 0u;
  for (auto i = 0u; i < n; ++i) {
    const auto& Hi = H[i];
    // const double hzi = 1.0/(Hi*Vector(0.0, 0.0, 1.0)).magnitude();
    // const double hzi = Hi.Inverse().Trace()/3.0;
    const auto hzi = Hi.Inverse().eigenValues().maxElement();
    const auto yi = y[i];
    // const double dphi = CylindricalBoundary::angularSpacing(yi, hzi, nNodePerh, kernelExtent);
    const auto nhoopsegment = max(1u, unsigned(phi*yi/(hzi/nNodePerh) + 0.5));
    const auto dphi = phi/nhoopsegment;
    CHECK(distinctlyGreaterThan(dphi, 0.0));
    const auto nsegment = max(1u, unsigned(phi/dphi + 0.5));
    ntot += nsegment;
  }

  // Determine how the global IDs should be partitioned between processors.
  // This could actually fail if we have more processors than points, so if that comes up we need to generalize...
  const size_t ndomain0 = ntot/nProcs;
  const size_t remainder = ntot % nProcs;
  VERIFY(remainder < size_t(nProcs));
  const size_t ndomain = ndomain0 + (size_t(procID) < remainder ? 1u : 0u);
  const size_t minGlobalID = procID*ndomain0 + min(size_t(procID), remainder);
  const size_t maxGlobalID = minGlobalID + ndomain - 1u;
  VERIFY(unsigned(procID) < (nProcs - 1u) || maxGlobalID == (ntot - 1u));
  VERIFY2(ntot < std::numeric_limits<size_t>::max(), "generateCylDistributionFromRZ ERROR: requested configuation requires " << ntot << " points be generated, which exceeds the maximum possible value of " << std::numeric_limits<size_t>::max());
  
  // Copy the input.
  vector<double> xrz(x), yrz(y), zrz(z), mrz(m);
  vector<SymTensor> Hrz(H);
  vector<vector<double> > extrasrz(extraFields);

  // Prepare the lists we're going to rebuild.
  x.clear();
  y.clear();
  z.clear();
  m.clear();
  H.clear();
  globalIDs.clear();
  extraFields.clear();
  extraFields.resize(nextra);

  // Iterate over the plane of input nodes, and rotate it out for the full 3-D 
  // distribution.
  size_t globalID = 0u;
  for (auto i = 0u; i < n; ++i) {
    const auto& Hi = Hrz[i];
    // const double hzi = 1.0/(Hi*Vector(0.0, 0.0, 1.0)).magnitude();
    // const double hzi = Hi.Inverse().Trace()/3.0;
    const auto hzi = Hi.Inverse().eigenValues().maxElement();
    const auto xi = xrz[i];
    const auto yi = yrz[i];
    const auto mi = mrz[i];
    // const int nhoopsegment = max(1, int(phi/CylindricalBoundary::angularSpacing(yi, hzi, nNodePerh, kernelExtent) + 0.5));
    const size_t nhoopsegment = max(1u, unsigned(phi*yi/(hzi/nNodePerh) + 0.5));
    const auto dphi = phi/nhoopsegment;
    const auto posi = Vector(xi, yi, 0.0);
    for (auto ihoop = 0u; ihoop < nhoopsegment; ++ihoop) {
      const auto phii = (double(ihoop) + 0.5)*dphi;
      if (size_t(globalID) >= minGlobalID and size_t(globalID) <= maxGlobalID) {
        const auto xj = xi;
        const auto yj = yi*cos(phii);
        const auto zj = yi*sin(phii);
        x.push_back(xj);
        y.push_back(yj);
        z.push_back(zj);
        m.push_back(mi/nhoopsegment * phi/(2.0*M_PI));
        globalIDs.push_back(globalID);
        const auto posj = Vector(xj, yj, zj);
        const auto R = CylindricalBoundary::reflectOperator(posi, posj);
        H.push_back((R*Hi*R).Symmetric());
        for (auto ikey = 0u; ikey < nextra; ++ikey) extraFields[ikey].push_back(extrasrz[ikey][i]);
      }
      ++globalID;
    }
  }

  // Post-conditions.
  VERIFY2(x.size() == ndomain and
          y.size() == ndomain and
          z.size() == ndomain and
          m.size() == ndomain and
          globalIDs.size() == ndomain and
          H.size() == ndomain,
          "Something wrong with the final array sizes: " << ndomain << " != " << x.size() << " " << y.size() << " " << z.size() << " " << globalIDs.size() << " " << H.size());
  for (auto ikey = 0u; ikey < nextra; ++ikey) VERIFY(extraFields[ikey].size() == ndomain);
  auto nglobal = x.size();
  if (nProcs > 1) {
    nglobal = allReduce(x.size(), SPHERAL_OP_SUM);
  }
  VERIFY(nglobal == ntot);

}

}
