//---------------------------------Spheral++----------------------------------//
// NeighborBase -- Base class for Neighbor and NeighborView.  This class
// maintains the common data and interface between the two.
//
// Created by J. Michael Owen, Mon Sep 21 14:46:15 PDT 2026
//----------------------------------------------------------------------------//
#ifndef __Spheral_NeighborBase__
#define __Spheral_NeighborBase__

#include "config.hh"
#include "Geometry/Dimension.hh"
#include "Utilities/DBC.hh"

#include <span>

namespace Spheral {

template<typename Dimension> class NodeList;

enum class NeighborSearchType {
  Gather = 1,
  Scatter = 2,
  GatherScatter = 3
};

template<typename Dimension>
class NeighborBase {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  // Constructors and destructors
  SPHERAL_HOST_DEVICE NeighborBase(const NeighborSearchType searchType,
                                   const Scalar kernelExtent):
    mSearchType(searchType),
    mKernelExtent(kernelExtent) {}
  SPHERAL_HOST_DEVICE virtual ~NeighborBase() = default;

  // The type of neighbor search we wish to use.
  SPHERAL_HOST_DEVICE NeighborSearchType neighborSearchType() const                { return mSearchType; }
  // SPEHRAL_HOST_DEVICE void neighborSearchType(const NeighborSearchType searchType) { mSearchType = searchType; }

  // All neighboring classes need to now how far to sample.
  SPHERAL_HOST_DEVICE Scalar kernelExtent() const                                  { return mKernelExtent; }
  // SPHERAL_HOST_DEVICE void kernelExtent(const kernelExtent)                        { REQUIRE(kernelExtent > 0.0; mKernelExtent = kernelExtent; }

  // Determine the maximum extent of a given H smoothing scale along the
  // Cartesian axes.
  SPHERAL_HOST_DEVICE static Vector Hextent(const Scalar& H, const Scalar& kernelExtent);
  SPHERAL_HOST_DEVICE static Vector Hextent(const SymTensor& H, const Scalar& kernelExtent);

protected:
  //-------------------------- Protected Interface --------------------------//
  NeighborSearchType mSearchType;
  Scalar mKernelExtent;
};

// We explicitly specialize the Hextent method for 1, 2, & 3 dimensions.
template<> Dim<1>::Vector NeighborBase<Dim<1>>::Hextent(const Dim<1>::SymTensor&, const Scalar& kernelExtent);
template<> Dim<2>::Vector NeighborBase<Dim<2>>::Hextent(const Dim<2>::SymTensor&, const Scalar& kernelExtent);
template<> Dim<3>::Vector NeighborBase<Dim<3>>::Hextent(const Dim<3>::SymTensor&, const Scalar& kernelExtent);

}

#include "NeighborBaseInline.hh"

#endif
