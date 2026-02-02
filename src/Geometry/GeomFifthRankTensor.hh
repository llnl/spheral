//---------------------------------Spheral++----------------------------------//
// GeomFifthRankTensor -- Geometric Tensor (rank 5) Class.
//
// This is a very simple, limited functionality rank 5 tensor.  Assumes
// the thing is dimensioned nDim x nDim x nDim x nDim.
//
// Created by JMO, Tue Oct 13 15:24:56 PDT 2015
//----------------------------------------------------------------------------//
#ifndef __Spheral_GeomFifthRankTensor_hh__
#define __Spheral_GeomFifthRankTensor_hh__

#include <iostream>

#include "Geometry/RankNTensor.hh"
#include "Geometry/GeomFifthRankTensor_fwd.hh"

namespace Spheral {

template<int nDim>
class GeomFifthRankTensor : public RankNTensor<nDim, 5, GeomFifthRankTensor<nDim>> {
  using BaseType = RankNTensor<nDim, 5, GeomFifthRankTensor<nDim>>;

public:
  //--------------------------- Public Interface ---------------------------//
  using size_type = typename BaseType::size_type;

  SPHERAL_HOST_DEVICE static constexpr GeomFifthRankTensor zero() { return GeomFifthRankTensor<nDim>(); }

  // Constructors.
  SPHERAL_HOST_DEVICE GeomFifthRankTensor() = default;
  SPHERAL_HOST_DEVICE explicit GeomFifthRankTensor(const double val);
  SPHERAL_HOST_DEVICE GeomFifthRankTensor(const GeomFifthRankTensor& rhs) = default;

  // Assignment.
  SPHERAL_HOST_DEVICE GeomFifthRankTensor& operator=(const GeomFifthRankTensor& rhs) = default;
  SPHERAL_HOST_DEVICE GeomFifthRankTensor& operator=(const double rhs);

  // Access the elements by indicies.
  SPHERAL_HOST_DEVICE double operator()(const size_type i, const size_type j, const size_type k, const size_type m, const size_type n) const;
  SPHERAL_HOST_DEVICE double& operator()(const size_type i, const size_type j, const size_type k, const size_type m, const size_type n);

private:
  //--------------------------- Private Interface ---------------------------//
  using RankNTensor<nDim, 5, GeomFifthRankTensor>::mElements;
};

}

#include "GeomFifthRankTensorInline.hh"

#endif
