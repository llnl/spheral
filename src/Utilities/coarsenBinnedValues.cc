//---------------------------------Spheral++----------------------------------//
// coarsenBinnedValues
//
// Given a lattice of Values, return the multi-level result of progressively
// coarsening the distribution a requeted number of levels.
//
// Created by JMO, Fri Feb 19 09:38:29 PST 2010
//----------------------------------------------------------------------------//
#include "coarsenBinnedValues.hh"

using std::vector;
using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {

template<typename Dimension>
void
coarsenWorkBins(vector<vector<double> >& workBins,
                const vector<unsigned>& ncells);

//------------------------------------------------------------------------------
// 1-D
//------------------------------------------------------------------------------
template<typename Value>
void
coarsenBinnedValues(vector<vector<Value> >& values,
                    const unsigned nxFine) {
  const int numLevels = values.size();
  REQUIRE(values.back().size() == nxFine);
  REQUIRE(nxFine / (1 << (numLevels - 1)) > 0);
  if (numLevels > 1) {
    for (int level = numLevels - 2; level >= 0; --level) {
      const unsigned nx = nxFine/(1 << (numLevels - level - 1));
      const size_t numLevelBins = nx;
      values[level] = vector<Value>(numLevelBins, 0.0);
      for (unsigned ix = 0; ix != nx; ++ix) {
        const unsigned ix0 = 2*ix;
        values[level][ix] = (values[level + 1][ix0] +
                             values[level + 1][ix0 + 1]);
      }
    }
  }
}

//------------------------------------------------------------------------------
// 2-D
//------------------------------------------------------------------------------
template<typename Value>
void
coarsenBinnedValues(vector<vector<Value> >& values,
                    const unsigned nxFine,
                    const unsigned nyFine) {
  const int numLevels = values.size();
  REQUIRE(values.back().size() == nxFine*nyFine);
  REQUIRE(nxFine / (1 << (numLevels - 1)) > 0);
  REQUIRE(nyFine / (1 << (numLevels - 1)) > 0);
  if (numLevels > 1) {
    for (int level = numLevels - 2; level >= 0; --level) {
      const unsigned nx = nxFine/(1 << (numLevels - level - 1));
      const unsigned ny = nyFine/(1 << (numLevels - level - 1));
      const unsigned nx0 = 2*nx;
      const size_t numLevelBins = nx*ny;
      values[level] = vector<Value>(numLevelBins, 0.0);
      for (unsigned iy = 0; iy != ny; ++iy) {
        const unsigned iy0 = 2*iy;
        for (unsigned ix = 0; ix != nx; ++ix) {
          const unsigned ix0 = 2*ix;
          values[level][ix + iy*nx] = (values[level + 1][ix0     + iy0*      nx0] +
                                       values[level + 1][ix0 + 1 + iy0*      nx0] +
                                       values[level + 1][ix0     + (iy0 + 1)*nx0] +
                                       values[level + 1][ix0 + 1 + (iy0 + 1)*nx0]);
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
// 3-D
//------------------------------------------------------------------------------
template<typename Value>
void
coarsenBinnedValues(vector<vector<Value> >& values,
                    const unsigned nxFine,
                    const unsigned nyFine,
                    const unsigned nzFine) {
  const int numLevels = values.size();
  REQUIRE(values.back().size() == nxFine*nyFine*nzFine);
  REQUIRE(nxFine / (1 << (numLevels - 1)) > 0);
  REQUIRE(nyFine / (1 << (numLevels - 1)) > 0);
  REQUIRE(nzFine / (1 << (numLevels - 1)) > 0);
  if (numLevels > 1) {
    for (int level = numLevels - 2; level >= 0; --level) {
      const unsigned nx = nxFine/(1 << (numLevels - level - 1));
      const unsigned ny = nyFine/(1 << (numLevels - level - 1));
      const unsigned nz = nzFine/(1 << (numLevels - level - 1));
      const unsigned nxy = nx*ny;
      const unsigned nx0 = 2*nx;
      const unsigned ny0 = 2*ny;
      const unsigned nxy0 = nx0*ny0;
      const size_t numLevelBins = nx*ny*nz;
      values[level] = vector<Value>(numLevelBins, 0.0);
      for (unsigned iz = 0; iz != nz; ++iz) {
        const unsigned iz0 = 2*iz;
        for (unsigned iy = 0; iy != ny; ++iy) {
          const unsigned iy0 = 2*iy;
          for (unsigned ix = 0; ix != nx; ++ix) {
            const unsigned ix0 = 2*ix;
            values[level][ix + iy*nx + iz*nxy] = (values[level + 1][ix0     + iy0*      nx0 + iz0*      nxy0] +
                                                  values[level + 1][ix0 + 1 + iy0*      nx0 + iz0*      nxy0] +
                                                  values[level + 1][ix0     + (iy0 + 1)*nx0 + iz0*      nxy0] +
                                                  values[level + 1][ix0 + 1 + (iy0 + 1)*nx0 + iz0*      nxy0] +
                                                  values[level + 1][ix0     + iy0*      nx0 + (iz0 + 1)*nxy0] +
                                                  values[level + 1][ix0 + 1 + iy0*      nx0 + (iz0 + 1)*nxy0] +
                                                  values[level + 1][ix0     + (iy0 + 1)*nx0 + (iz0 + 1)*nxy0] +
                                                  values[level + 1][ix0 + 1 + (iy0 + 1)*nx0 + (iz0 + 1)*nxy0]);
          }
        }
      }
    }
  }

  // Post-conditions.
  BEGIN_CONTRACT_SCOPE
  {
    for (size_t level = 1; (int)level < numLevels; ++level) {
      const unsigned nx = nxFine/(1U << (numLevels - level - 1));
      const unsigned ny = nyFine/(1U << (numLevels - level - 1));
      const unsigned nz = nzFine/(1U << (numLevels - level - 1));
      const unsigned nxparent = nx/2;
      const unsigned nyparent = ny/2;
      const unsigned nzparent = nz/2;
      for (unsigned izparent = 0; izparent != nzparent; ++izparent) {
        const unsigned iz0 = 2*izparent;
        for (unsigned iyparent = 0; iyparent != nyparent; ++iyparent) {
          const unsigned iy0 = 2*iyparent;
          for (unsigned ixparent = 0; ixparent != nxparent; ++ixparent) {
            const unsigned ix0 = 2*ixparent;
            const Value sum = (values[level][(ix0    ) + (iy0    )*nx + (iz0    )*nx*ny] +
                               values[level][(ix0 + 1) + (iy0    )*nx + (iz0    )*nx*ny] +
                               values[level][(ix0    ) + (iy0 + 1)*nx + (iz0    )*nx*ny] +
                               values[level][(ix0 + 1) + (iy0 + 1)*nx + (iz0    )*nx*ny] +
                               values[level][(ix0    ) + (iy0    )*nx + (iz0 + 1)*nx*ny] +
                               values[level][(ix0 + 1) + (iy0    )*nx + (iz0 + 1)*nx*ny] +
                               values[level][(ix0    ) + (iy0 + 1)*nx + (iz0 + 1)*nx*ny] +
                               values[level][(ix0 + 1) + (iy0 + 1)*nx + (iz0 + 1)*nx*ny]);
            CONTRACT_VAR(sum);
            ENSURE(fuzzyEqual(sum, values[level - 1][ixparent + iyparent*nxparent + izparent*nxparent*nyparent], 1.0e-10));
          }
        }
      }
    }
  }
  END_CONTRACT_SCOPE
}

}
