//---------------------------------Spheral++----------------------------------//
// invertRankNTensor -- Overloaded method to compute the inverse of our
// specialized rank N tensors.
//
// Created by JMO, Mon Oct 12 14:16:42 PDT 2015
//----------------------------------------------------------------------------//

#include <iostream>
#include <algorithm>
#include <assert.h>
#include "Eigen/Dense"

#include "invertRankNTensor.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

// General form.
template<typename Tensor>
Tensor
invertRankNTensor(const Tensor& tensor) {
  assert(Tensor::nrank() % 2 == 0);

  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixType;
  const unsigned dim = (unsigned)std::pow(double(Tensor::nDimensions()), int(Tensor::nrank()/2));
  MatrixType A(dim, dim);

  // Note our convention for how rank N tensors work -- last index changing the fastest.
  // We use this knowledge of the 1D packing order here.
  std::copy(tensor.begin(), tensor.end(), A.data());

  // std::cout << "Input tensor: " << tensor << std::endl
  //           << "Eigen matrix: " << A << std::endl;

  // Check if this thing is invertable.
  VERIFY2(std::abs(A.determinant()) > 1.0e-20, "invertRankNTensor : input appears to be singular.");

  // Invert and load the return type.
  const MatrixType Ainv = A.inverse();
  Tensor result;
  std::copy(Ainv.data(), Ainv.data() + Ainv.size(), result.begin());
  return result;
}

}

