//---------------------------------Spheral++----------------------------------//
// GeomThirdRankTensor -- Geometric Tensor (rank 3) Class.
//
// This is a very simple, limited functionality rank 3 tensor.  Assumes
// the thing is dimensioned nDim x nDim x nDim.
//
// Created by JMO, Thu Jul 17 17:05:58 PDT 2008
//----------------------------------------------------------------------------//
#ifndef __Spheral_GeomThirdRankTensor_hh__
#define __Spheral_GeomThirdRankTensor_hh__

#include <iostream>

#include "Geometry/RankNTensor.hh"
#include "Geometry/GeomThirdRankTensor_fwd.hh"

namespace Spheral {

template<int nDim>
class GeomThirdRankTensor : public RankNTensor<nDim, 3, GeomThirdRankTensor<nDim>> {
  using BaseType = RankNTensor<nDim, 3, GeomThirdRankTensor<nDim>>;

public:
  //--------------------------- Public Interface ---------------------------//
  using size_type = typename BaseType::size_type;

  SPHERAL_HOST_DEVICE static constexpr GeomThirdRankTensor zero() { return GeomThirdRankTensor<nDim>(); }

  // Constructors.
  SPHERAL_HOST_DEVICE GeomThirdRankTensor() = default;
  SPHERAL_HOST_DEVICE explicit GeomThirdRankTensor(const double val);
  SPHERAL_HOST_DEVICE GeomThirdRankTensor(const GeomThirdRankTensor& rhs) = default;

  // Assignment.
  SPHERAL_HOST_DEVICE GeomThirdRankTensor& operator=(const GeomThirdRankTensor& rhs) = default;
  SPHERAL_HOST_DEVICE GeomThirdRankTensor& operator=(const double rhs);

  // Access the elements by indicies.
  SPHERAL_HOST_DEVICE double operator()(const size_type i, const size_type j, const size_type k) const;
  SPHERAL_HOST_DEVICE double& operator()(const size_type i, const size_type j, const size_type k);

private:
  //--------------------------- Private Interface ---------------------------//
  using RankNTensor<nDim, 3, GeomThirdRankTensor>::mElements;
};

}

#include "GeomThirdRankTensorInline.hh"

#endif
