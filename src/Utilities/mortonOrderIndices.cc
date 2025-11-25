//---------------------------------Spheral++----------------------------------//
// mortonOrderIndices
//
// Compute the Morton ordered hashed indices for the given set of NodeLists.
// 
// Algorithm described in
// Warren & Salmon (1995), Computer Physics Communications, 87, 266-290.
//
// Created by JMO, Fri Dec 19 14:58:23 PST 2008
//----------------------------------------------------------------------------//
#include "mortonOrderIndices.hh"
#include "globalBoundingVolumes.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"
#include "Utilities/DBC.hh"

using std::vector;
using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {

//------------------------------------------------------------------------------
// Stand alone functions to do the interleaved bit hashing of positions.
//------------------------------------------------------------------------------
inline
KeyTraits::Key
hashx(const double x, const double dx) {
  REQUIRE(x >= 0.0);
  REQUIRE(dx > 0.0);
  const double sresult = x/dx;
  KeyTraits::Key iresult = static_cast<KeyTraits::Key>(sresult);
  const double rem = sresult - iresult*dx;
  if (fuzzyEqual(std::abs(rem), 0.0, 1.0e-10) and rem > 0.0) iresult -= KeyTraits::one;
  // const KeyTraits::Key result = static_cast<KeyTraits::Key>(x/dx);
  ENSURE(iresult <= KeyTraits::maxKey1d);
  return iresult;
}

// 1-D
inline
KeyTraits::Key
hashPosition(const Dim<1>::Vector& xoff,
             const Dim<1>::Vector& stepSize) {
  const KeyTraits::Key result = hashx(xoff.x(), stepSize.x());
  ENSURE(result <= KeyTraits::maxKey);
  return result;
}

// 2-D
inline
KeyTraits::Key
hashPosition(const Dim<2>::Vector& xoff,
             const Dim<2>::Vector& stepSize) {

  typedef KeyTraits::Key Key;

  // Get the bits for each dimension.
  const Key xx = hashx(xoff.x(), stepSize.x());
  const Key yy = hashx(xoff.y(), stepSize.y());

  // Interleave the bits.
  Key result = KeyTraits::zero;
  Key mask = KeyTraits::one;
  for (size_t i = 0; i != KeyTraits::numbits1d; ++i, mask <<= 1) {
    result += (((xx & mask) << i) +
               ((yy & mask) << (i + 1)));
  }

  // That's it.
  ENSURE(result <= KeyTraits::maxKey);
  return result;
}

// 3-D
inline
KeyTraits::Key
hashPosition(const Dim<3>::Vector& xoff,
             const Dim<3>::Vector& stepSize) {

  typedef KeyTraits::Key Key;

  // Get the bits for each dimension.
  const Key xx = hashx(xoff.x(), stepSize.x());
  const Key yy = hashx(xoff.y(), stepSize.y());
  const Key zz = hashx(xoff.z(), stepSize.z());

  // Interleave the bits.
  Key result = KeyTraits::zero;
  Key mask = KeyTraits::one;
  for (size_t i = 0; i != KeyTraits::numbits1d; ++i, mask <<= 1) {
    result += (((xx & mask) << i) +
               ((yy & mask) << (i + 1)) +
               ((zz & mask) << (i + 2)));
  }

  // That's it.
  ENSURE(result <= KeyTraits::maxKey);
  return result;
}

//------------------------------------------------------------------------------
// Hash the node positions into their tree ordered indices.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, KeyTraits::Key>
mortonOrderIndices(const FieldList<Dimension, typename Dimension::Vector>& positions) {

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
  const Vector stepSize = (xmax - xmin)/KeyTraits::maxKey1d;

  // Go over all nodes.
  for (AllNodeIterator<Dimension> nodeItr = positions.nodeBegin();
       nodeItr != positions.nodeEnd();
       ++nodeItr) {

    // Find the offset from the minimum coordinates.
    const Vector xoff = positions(nodeItr) - xmin;

    // Hash that sucker.
    result(nodeItr) = hashPosition(xoff, stepSize);
  }

  return result;
}

//------------------------------------------------------------------------------
// Hash the node positions into their tree ordered indices.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, KeyTraits::Key>
mortonOrderIndices(const DataBase<Dimension>& dataBase) {
  return mortonOrderIndices(dataBase.globalPosition());
}

//------------------------------------------------------------------------------
// Same as above except allowing a mask to be applied (ignoring nodes with 
// mask = 0).
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, KeyTraits::Key>
mortonOrderIndices(const DataBase<Dimension>& dataBase,
                    const FieldList<Dimension, int>& mask) {

  typedef typename KeyTraits::Key Key;
  typedef typename Dimension::Vector Vector;

  // Prepare the result.
  FieldList<Dimension, Key> result = dataBase.newGlobalFieldList(KeyTraits::zero, "hashed indices");

  // Get the bounding box and step sizes.
  Vector xmin, xmax;
  dataBase.boundingBox(xmin, xmax, mask, true);
  const Vector stepSize = (xmax - xmin)/KeyTraits::maxKey1d;

  // Go over all nodes.
  const FieldList<Dimension, Vector> positions = dataBase.globalPosition();
  for (AllNodeIterator<Dimension> nodeItr = dataBase.nodeBegin();
       nodeItr != dataBase.nodeEnd();
       ++nodeItr) {
    if (mask(nodeItr) != 0) {

      // Find the offset from the minimum coordinates.
      const Vector xoff = positions(nodeItr) - xmin;

      // Hash that sucker.
      result(nodeItr) = hashPosition(xoff, stepSize);
    }
  }

  return result;
}

}

