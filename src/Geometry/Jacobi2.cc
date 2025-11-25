// A very simple implementation of the iterative Jacobi rotation method for 
// diagonalizing a real symmetric matrix. 

#include "Jacobi2.hh"
#include "Dimension.hh"

#include <iostream>
using std::abs;

namespace Spheral {

template<typename Dimension>
inline
double sumOffDiagonalElements(const typename Dimension::SymTensor& D) {
  typedef typename Dimension::Tensor Tensor;
  double result = 0.0;
  for (typename Tensor::size_type row = 1; row != Dimension::nDim - 1; ++row) {
    for (typename Tensor::size_type column = row; column != Dimension::nDim; ++column) {
      result += abs(D(row, column));
    }
  }
  return result;
}

template<typename Dimension>
int jacobiDiagonalize(const typename Dimension::SymTensor& A,
                      typename Dimension::Tensor& eigenVectors,
                      typename Dimension::Vector& eigenValues,
                      const double convergenceThreshold,
                      const int maxSweeps) {

  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Tensor::size_type size_type;
  const int nDim = Dimension::nDim;

  // Begin by copying the input tensor to a local version, which we will
  // iteratively rotate until it is diagonal.
  SymTensor D = A;

  // Initialize the eigen vectors tensor as the identity tensor.
  eigenVectors.Identity();

  // Prepare to sweep over the matrix, elimnating each off diagonal element
  // in turn.
  int nSweeps = 0;
  int nRotations = 0;
  double error = sumOffDiagonalElements<Dimension>(D);
  const double epsilon = 1e-10;
  while (nSweeps < maxSweeps && error > convergenceThreshold) {
    const double elementThreshold = 0.2*error/(nDim*nDim);

    // Loop over each of the off diagonal elements.
    for (size_type row = 0; row != nDim - 1; ++row) {
      for (size_type column = row + 1; column != nDim; ++column) {

        // We carry out this rotation only if this element is above
        // a critical threshold.
        const double a = D(row, row);
        const double b = D(row, column);
        const double c = D(column, column);
        bool rotate = true;
        if (nSweeps < 4 && 
            abs(b) < elementThreshold) rotate = false;
        if (nSweeps > 3 &&
            (abs(b) < epsilon*abs(a) && abs(b) < epsilon*abs(c))) {
          rotate = false;
          D(row, column) = 0.0;
        }

        // Generate the rotation matrix appropriate to eliminate this
        // element.
        if (rotate) {
          const double theta = 0.5*atan2(2.0*b, c - a);
          const double xhat = cos(theta);
          const double yhat = sin(theta);
          Tensor R = Tensor::one();
          R(row, row) = xhat;
          R(row, column) = -yhat;
          R(column, row) = yhat;
          R(column, column) = xhat;

          // Apply this rotational transform to the D tensor.
          D.rotationalTransform(R);

          // Accumulate the local rotation on the eigen vectors.
          eigenVectors *= R.Transpose();
          ++nRotations;
        }
      }
    }

    // Calculate the new error in the off diagonal terms.
    error = sumOffDiagonalElements<Dimension>(D);
    ++nSweeps;
  }

//    // The eigen vectors tensor is currently in row form, so switch to
//    // column vectors.
//    eigenVectors = eigenVectors.Transpose();

  // The eigen vectors tensor is already set, but we need to set the Vector
  // of eigen values to the diagonal of D.
  eigenValues = D.diagonalElements();

  // Return the number of rotations that were required.
  return nRotations;
}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
template int jacobiDiagonalize< Dim<3> >(const Dim<3>::SymTensor& A,
                                         Dim<3>::Tensor& eigenVectors,
                                         Dim<3>::Vector& eigenValues,
                                         const double convergenceThreshold,
                                         const int maxSweeps);

}
