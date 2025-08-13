//---------------------------------Spheral++----------------------------------//
// GeomFourthRankTensor -- Geometric Tensor (rank 4) Class.
//
// This is a very simple, limited functionality rank 4 tensor.  Assumes
// the thing is dimensioned nDim x nDim x nDim x nDim.
//
// Created by JMO, Sun Oct 11 22:10:34 PDT 2015
//----------------------------------------------------------------------------//
#ifndef __Spheral_GeomFourthRankTensor_hh__
#define __Spheral_GeomFourthRankTensor_hh__

#include <iostream>

#include "Geometry/RankNTensor.hh"
#include "Geometry/GeomFourthRankTensor_fwd.hh"

namespace Spheral {

template<int nDim>
class GeomFourthRankTensor : public RankNTensor<nDim, 4, GeomFourthRankTensor<nDim> > {
  using BaseType = RankNTensor<nDim, 4, GeomFourthRankTensor<nDim>>;

public:
  //--------------------------- Public Interface ---------------------------//
  using size_type = typename BaseType::size_type;
  static constexpr size_type numElements = BaseType::numElements;

  // Useful static member data.
  static const GeomFourthRankTensor zero;

  // Constructors.
  SPHERAL_HOST_DEVICE GeomFourthRankTensor() = default;
  SPHERAL_HOST_DEVICE explicit GeomFourthRankTensor(const double val);
  SPHERAL_HOST_DEVICE GeomFourthRankTensor(const GeomFourthRankTensor& rhs) = default;

  // Assignment.
  SPHERAL_HOST_DEVICE GeomFourthRankTensor& operator=(const GeomFourthRankTensor& rhs) = default;
  SPHERAL_HOST_DEVICE GeomFourthRankTensor& operator=(const double rhs);

  // Access the elements by indicies.
  SPHERAL_HOST_DEVICE double operator()(const size_type i, const size_type j, const size_type k, const size_type m) const;
  SPHERAL_HOST_DEVICE double& operator()(const size_type i, const size_type j, const size_type k, const size_type m);

private:
  //--------------------------- Private Interface ---------------------------//
  using RankNTensor<nDim, 4, GeomFourthRankTensor>::mElements;
};

template<int nDims> const GeomFourthRankTensor<nDims> GeomFourthRankTensor<nDims>::zero = GeomFourthRankTensor<nDims>(0.0);

}

#include "GeomFourthRankTensorInline.hh"

#endif
