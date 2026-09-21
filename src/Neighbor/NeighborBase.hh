//---------------------------------Spheral++----------------------------------//
// NeighborBase -- Base class for Neighbor and NeighborView.  This class
// maintains the common data and interface between the two.
//
// Created by J. Michael Owen, Mon Sep 21 14:46:15 PDT 2026
//----------------------------------------------------------------------------//
#ifndef __Spheral_NeighborBase__
#define __Spheral_NeighborBase__

#include "Geometry/Dimension.hh"

#include <vector>

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
                                   const Scalar kernelExtent);
  SPHERAL_HOST_DEVICE virtual ~NeighborBase() = default;

  // The type of neighbor search we wish to use.
  SPHERAL_HOST_DEVICE NeighborSearchType neighborSearchType() const { return mSearchType; }

  // All neighboring classes need to now how far to sample.
  SPHERAL_HOST_DEVICE Scalar kernelExtent() const                   { return mKernelExtent; }


  // Determine the maximum extent of a given H smoothing scale along the
  // Cartesian axes.
  SPHERAL_HOST_DEVICE static Vector HExtent(const Scalar& H, const Scalar& kernelExtent);
  SPHERAL_HOST_DEVICE static Vector HExtent(const SymTensor& H, const Scalar& kernelExtent);

private:
  //--------------------------- Private Interface ---------------------------//
  NeighborSearchType mSearchType;
  Scalar mKernelExtent;
};

// We explicitly specialize the HExtent method for 1, 2, & 3 dimensions.
template<> Dim<1>::Vector Neighbor< Dim<1> >::HExtent(const Dim<1>::SymTensor&, const Scalar kernelExtent);
template<> Dim<2>::Vector Neighbor< Dim<2> >::HExtent(const Dim<2>::SymTensor&, const Scalar kernelExtent);
template<> Dim<3>::Vector Neighbor< Dim<3> >::HExtent(const Dim<3>::SymTensor&, const Scalar kernelExtent);

}

#include "NeighborBaseInline.hh"

#endif
