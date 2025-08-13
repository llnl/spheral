//---------------------------------Spheral++----------------------------------//
// peanoHilbertOrderIndices
//
// Compute the PeanoHilbert ordered hashed indices for the given set of NodeLists.
// 
// Algorithm described in
// Warren & Salmon (1995), Computer Physics Communications, 87, 266-290.
//
// Created by JMO, Sat Dec 20 22:36:58 PST 2008
//----------------------------------------------------------------------------//
#include "peanoHilbertOrderIndices.hh"
#include "PeanoHilbertTransform.hh"
#include "globalBoundingVolumes.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"
#include "Utilities/DBC.hh"

#include <limits>

using std::vector;
using std::string;
using std::pair;
using std::make_pair;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {

//------------------------------------------------------------------------------
// The specialized function (per dimension) to dive down the recursive levels
// of the Peano-Hilbert path and generate the encoded key for a given position.
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// 1-D
//------------------------------------------------------------------------------
inline
KeyTraits::Key
hashPosition(const Dim<1>::Vector& position,
             const Dim<1>::Vector& boxmin,
             const Dim<1>::Vector& boxmax) {
  const double boxsize = (1.0 + 1.0e-10)*(boxmax.x() - boxmin.x());
  CHECK(boxsize > 0.0);
  return KeyTraits::Key(int((position.x() - boxmin.x())/(boxsize/KeyTraits::numbits1d)));
}

//------------------------------------------------------------------------------
// 2-D
//------------------------------------------------------------------------------
inline
KeyTraits::Key
hashPosition(const Dim<2>::Vector& position,
             const Dim<2>::Vector& boxmin,
             const Dim<2>::Vector& boxmax) {

  typedef KeyTraits::Key Key;
  typedef Dim<2>::Vector Vector;
  typedef Dim<2>::Scalar Scalar;

  const Scalar tiny = std::numeric_limits<Scalar>::epsilon();

  // Pre-conditions.
  REQUIRE(boxmin.x() <= boxmax.x() and boxmin.y() <= boxmax.y());
  REQUIRE(position.x() >= boxmin.x() and position.x() <= boxmax.x());
  REQUIRE(position.y() >= boxmin.y() and position.y() <= boxmax.y());

  // The set of transforms that define the Peano-Hilbert fractal pattern.
  vector<PeanoHilbertTransform2d> transforms = {
    PeanoHilbertTransform2d(0, 1,
                            1, 0),
    PeanoHilbertTransform2d(1, 0,
                            0, 1),
    PeanoHilbertTransform2d(1, 0,
                            0, 1),
    PeanoHilbertTransform2d(0, -1,
                            -1, 0)
  };
  CHECK(transforms.size() == 4);

  // Get the box dimensions.
  const Vector boxsize = (1.0 + 1.0e-10)*(boxmax - boxmin);
  CHECK(boxsize > 0.0);

  // Prepare our result.
  Key result = KeyTraits::zero;

  // Initialize the starting transform.
  PeanoHilbertTransform2d T = transforms[1];

  // Determine the integer position on the finest level.
  const Vector deltar = position - boxmin;
  CHECK(deltar >= 0.0);
  const Key ncells = KeyTraits::one << KeyTraits::numbits1d;
  const Vector cellsize = boxsize/ncells;
  const int xfine = int(deltar.x()/std::max(cellsize.x(),tiny));
  const int yfine = int(deltar.y()/std::max(cellsize.y(),tiny));
  CHECK(xfine >= 0 and xfine <= (int)ncells);
  CHECK(yfine >= 0 and yfine <= (int)ncells);

  // Recursively quadrant the position, until we get to the desired level.
  for (auto level = 0u; level != KeyTraits::numbits1d; ++level) {

    // Compute the integer coordinates on this level.
    const int ncellDelta = KeyTraits::one << (KeyTraits::numbits1d - level);
    const int x = xfine / ncellDelta;
    const int y = yfine / ncellDelta;

    const Key ncells = KeyTraits::two << level;
    CONTRACT_VAR(ncells);
    CHECK(x >= 0 and x <= (int)ncells);
    CHECK(y >= 0 and y <= (int)ncells);

    // Decide which (lab frame) quadrant this position represents.
    const int xx = x % 2;
    const int yy = y % 2;
    CHECK(xx == 0 or xx == 1);
    CHECK(yy == 0 or yy == 1);

    // Extract which quad along the local Peano-Hilbert path the
    // lab quad translates to.
    const unsigned phquad = T(xx, yy);
    CHECK(phquad < 4U);

    // Increment the result.
    Key delta = phquad;
    delta <<= (2*(KeyTraits::numbits1d - level));
    result += delta;

    // Transform down to the next level.
    T = T(transforms[phquad]);
  }

  // That's it.
  ENSURE(result <= KeyTraits::maxKey);
  return result;
}

//------------------------------------------------------------------------------
// 3-D
//------------------------------------------------------------------------------
inline
KeyTraits::Key
hashPosition(const Dim<3>::Vector& position,
             const Dim<3>::Vector& boxmin,
             const Dim<3>::Vector& boxmax) {

  typedef KeyTraits::Key Key;
  typedef Dim<3>::Vector Vector;
  typedef Dim<3>::Scalar Scalar;
  
  const Scalar tiny = std::numeric_limits<Scalar>::epsilon();
  
  // Pre-conditions.
  REQUIRE(boxmin.x() <= boxmax.x() and 
          boxmin.y() <= boxmax.y() and
          boxmin.z() <= boxmax.z());
  REQUIRE(position.x() >= boxmin.x() and position.x() <= boxmax.x());
  REQUIRE(position.y() >= boxmin.y() and position.y() <= boxmax.y());
  REQUIRE(position.z() >= boxmin.z() and position.z() <= boxmax.z());

  // The set of transforms that define the Peano-Hilbert fractal pattern.
  vector<PeanoHilbertTransform3d> transforms = {
    PeanoHilbertTransform3d(1, 0, 0,
                            0, 0, 1,
                            0, 1, 0),
    PeanoHilbertTransform3d(0, 0, 1,
                            0, 1, 0,
                            1, 0, 0),
    PeanoHilbertTransform3d(0, 0, 1,    // Same as T1
                            0, 1, 0,
                            1, 0, 0),
    PeanoHilbertTransform3d(-1, 0, 0,
                            0, -1, 0,
                            0, 0, 1),
    PeanoHilbertTransform3d(-1, 0, 0,   // Same as T3
                            0, -1, 0,
                            0, 0, 1),
    PeanoHilbertTransform3d(0, 0, -1,
                            0, 1, 0,
                            -1, 0, 0),
    PeanoHilbertTransform3d(0, 0, -1,   // Same as T5
                            0, 1, 0,
                            -1, 0, 0),
    PeanoHilbertTransform3d(1, 0, 0,
                            0, 0, -1,
                            0, -1, 0)
  };
  CHECK(transforms.size() == 8);

  // Get the box dimensions.
  const Vector boxsize = (1.0 + 1.0e-10)*(boxmax - boxmin);
  CHECK(boxsize > 0.0);

  // Prepare our result.
  Key result = KeyTraits::zero;

  // Initialize the starting transform.
  PeanoHilbertTransform3d T = transforms[1];

  // Determine the integer position on the finest level.
  const Vector deltar = position - boxmin;
  CHECK(deltar >= 0.0);
  const Key ncells = KeyTraits::one << KeyTraits::numbits1d;
  const Vector cellsize = boxsize/ncells;
  const int xfine = int(deltar.x()/std::max(cellsize.x(),tiny));
  const int yfine = int(deltar.y()/std::max(cellsize.y(),tiny));
  const int zfine = int(deltar.z()/std::max(cellsize.z(),tiny));
  CHECK(xfine >= 0 and xfine <= (int)ncells);
  CHECK(yfine >= 0 and yfine <= (int)ncells);
  CHECK(zfine >= 0 and zfine <= (int)ncells);

  // Recursively quadrant the position, until we get to the desired level.
  for (auto level = 0u; level != KeyTraits::numbits1d; ++level) {

    // Compute the integer coordinates on this level.
    const int ncellDelta = KeyTraits::one << (KeyTraits::numbits1d - level);
    const int x = xfine / ncellDelta;
    const int y = yfine / ncellDelta;
    const int z = zfine / ncellDelta;

    const Key ncells = KeyTraits::two << level;
    CONTRACT_VAR(ncells);
    CHECK(x >= 0 and x <= (int)ncells);
    CHECK(y >= 0 and y <= (int)ncells);
    CHECK(z >= 0 and z <= (int)ncells);

    // Decide which (lab frame) quadrant this position represents.
    const int xx = x % 2;
    const int yy = y % 2;
    const int zz = z % 2;
    CHECK(xx == 0 or xx == 1);
    CHECK(yy == 0 or yy == 1);
    CHECK(zz == 0 or zz == 1);

    // Extract which quad along the local Peano-Hilbert path the
    // lab quad translates to.
    const unsigned phquad = T(xx, yy, zz);
    CHECK(phquad < 8U);

    // Increment the result.
    Key delta = phquad;
    delta <<= (3*(KeyTraits::numbits1d - level));
    result += delta;

    // Transform down to the next level.
    T = T(transforms[phquad]);
  }

  // That's it.
  ENSURE(result <= KeyTraits::maxKey);
  return result;
}

//------------------------------------------------------------------------------
// The Lawder-Butz algorithm for generating the encoded keys.
// If this all looks like C, it is!  I've lifted it from Lawders paper
// nearly intact, just changing enough to use our ideas of what the dimension
// and Key size is.
//------------------------------------------------------------------------------
// template<unsigned nDim>
// struct Hcode {
//   unsigned hcode[nDim];
// };

// //----------------------------------------
// // g_mask
// //----------------------------------------
// template<unsigned nDim>
// unsigned g_mask(unsigned i) {
//   REQUIRE(i < nDim);
//   return 1 << (nDim - i);
// }

// //----------------------------------------
// // calc_P2
// //----------------------------------------
// template<unsigned nDim>
// unsigned calc_P2(unsigned S) {
//   int i;
//   unsigned P;

//   P = S & g_mask<nDim>(0);
//   for (i = 1; i < nDim; i++) {
//     if (S & g_mask<nDim>(i) ^ (P >> 1) & g_mask<nDim>(i)) {
//       P |= g_mask<nDim>(i);
//     }
//   }
// }

// //----------------------------------------
// // calc_J
// //----------------------------------------
// template<unsigned nDim>
// unsigned calc_J(unsigned P) {
//   int i;
//   unsigned J;

//   J = nDim;
//   for (i = 1; i < nDim; i++) {
//     if ((P >> i & 1) == (P & 1)) {
//       continue;
//     } else {
//       break;
//     }
//   }
//   if (i != nDim) J -= i;
//   return J;
// }

// //----------------------------------------
// // calc_T
// //----------------------------------------
// template<unsigned nDim>
// unsigned calc_T(unsigned P) {
//   if (P < 3) return 0;
//   if (P % 2) return (P - 1)^(P - 1)/2;
//   return (P - 2)^(P - 2)/2;
// }

// //----------------------------------------
// // calc_tS_tT
// //----------------------------------------
// template<unsigned nDim>
// unsigned calc_tS_tT(unsigned xJ, unsigned val) {
//   unsigned retval, temp1, temp2;

//   retval = val;
//   if (xJ % nDim != 0) {
//     temp1 = val >> xJ % nDim;
//     temp2 = val << nDim - xJ % nDim;
//     retval = temp1 | temp2;
//     retval &= (1U << nDim) - 1;
//   }
//   return retval;
// }

// //----------------------------------------
// // H_encode
// // Map from an nDim dimensional point to
// // the Hilbert 1-D index.
// //----------------------------------------
// template<unsigned nDim>
// KeyTraits::Key H_encode(Hcode<nDim> pt) {
//   const int ORDER = 32;   // Assuming 32 bit unsigned ints.
//   unsigned 
//     mask = 1U << ORDER - 1, 
//     element, 
//     A,
//     W = 0,
//     S,
//     tS,
//     T,
//     tT,
//     J,
//     P = 0,
//     xJ;
//   Hcode h = {0};
//   int i = ORDER * nDim - nDim, j;

//   for (j = A = 0; j < nDim; j++) {
//     if (pt.hcode[j] & mask) A |= g_mask<nDim>(j);
//   }

//   S = tS = A;
//   P = calc_P2(S);
  
//   // Add in nDim bits to hcode.
//   element = i/ORDER;
//   if (i % ORDER > ORDER - nDim) {
//     h.hcode[element]     |= P << i % ORDER;
//     h.hcode[element + 1] |= P >> ORDER - i % ORDER;
//   } else {
//     h.hcode[element] |= P << i - element * ORDER;
//   }

//   J = calc_J(P);
//   xJ = J - 1;
//   T = calc_T(P);
//   tT = T;

//   for (i -= nDim, mask >>= 1; i >= 0; i -= nDim, mask >>= 1) {
//     for (j = A = 0; j < nDim; j++) {
//       if (pt.hcode[j] & mask) A |= g_mask<nDim>(j);
//     }
//     W ^= tT;
//     tS = A ^ W;
//     S = calc_tS_tT(xJ, tS);
//     P = calc_P2(S);
      
//     // Add in nDim bits to hcode.
//     element = i / ORDER;
//     if (i % ORDER > ORDER - nDim) {
//       h.hcode[element]     |= P << i % ORDER;
//       h.hcode[element + 1] |= P >> ORDER - i % ORDER;
//     } else {
//       h.hcode[element] |= P << i - element * ORDER;
//     }

//     if (i > 0) {
//       T = calc_T(P);
//       tT = calc_tS_tT(xJ, T);
//       J = calc_J(P);
//       xJ += J - 1;
//     }
//   }

//   // Convert the 3 unsigned int result to our single key.
//   typedef KeyTraits::Key Key;
//   Key result = KeyTraits::zero;
//   for (int i = 0; i != nDim; ++i) {
//     Key thpt = h[i] << (i * ORDER);
//     result &= thpt;
//   }
//   return result;

// }

//------------------------------------------------------------------------------
// Hash the node positions into their tree ordered indices.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, KeyTraits::Key>
peanoHilbertOrderIndices(const FieldList<Dimension, typename Dimension::Vector>& positions) {

  typedef typename KeyTraits::Key Key;
  typedef typename Dimension::Vector Vector;

  // Prepare the result.
  FieldList<Dimension, Key> result(FieldStorageType::CopyFields);
  const vector<NodeList<Dimension>*>& nodeListPtrs = positions.nodeListPtrs();
  for (const NodeList<Dimension>* nodeListPtr: nodeListPtrs) {
    result.appendNewField("hashed indices", *nodeListPtr, KeyTraits::zero);
  }

  // Get the bounding box and step sizes.
  Vector xmin, xmax;
  globalBoundingBox(positions, xmin, xmax, true);

  // Go over all nodes and hash each position.
  for (AllNodeIterator<Dimension> nodeItr = positions.nodeBegin();
       nodeItr != positions.nodeEnd();
       ++nodeItr) {
    result(nodeItr) = hashPosition(positions(nodeItr), xmin, xmax);
  }

  return result;
}

//------------------------------------------------------------------------------
// Hash the node positions into their tree ordered indices.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, KeyTraits::Key>
peanoHilbertOrderIndices(const DataBase<Dimension>& dataBase) {
  return peanoHilbertOrderIndices(dataBase.globalPosition());
}

}

