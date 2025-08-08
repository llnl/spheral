//---------------------------------Spheral++----------------------------------//
// Kernel -- The interpolation kernel for use in smoothed field estimates.
//
// Created by JMO, Thu Jul 29 19:43:35 PDT 1999
//----------------------------------------------------------------------------//
#ifndef __Spheral_Kernel_hh__
#define __Spheral_Kernel_hh__

#include "config.hh"

namespace Spheral {

template<typename Dimension, typename Descendant>
class Kernel {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  // Cast as the descendent type.
  SPHERAL_HOST_DEVICE Descendant& asDescendant() const;

  //======================================================================
  // Return the kernel weight
  SPHERAL_HOST_DEVICE double operator()(const double& etaij, const Scalar& Hdet) const;
  SPHERAL_HOST_DEVICE double operator()(const Vector& etaj, const Vector& etai, const Scalar& Hdet) const;

  //======================================================================
  // Return the gradient value for a given normalized distance or position.
  SPHERAL_HOST_DEVICE double grad(const double& etaij, const Scalar& Hdet) const;
  SPHERAL_HOST_DEVICE double grad(const Vector& etaj, const Vector& etai, const Scalar& Hdet) const;

  //======================================================================
  // Return the second derivative of the kernel for a given normalized distance
  //  or position.
  SPHERAL_HOST_DEVICE double grad2(const double& etaij, const Scalar& Hdet) const;
  SPHERAL_HOST_DEVICE double grad2(const Vector& etaj, const Vector& etai, const Scalar& Hdet) const;

  //======================================================================
  // Get the volume normalization constant.
  SPHERAL_HOST_DEVICE double volumeNormalization() const;

  // Get the extent of the kernel (the cutoff distance in eta over which the
  // kerel is non-zero.
  SPHERAL_HOST_DEVICE double kernelExtent() const;

  // We also require that all kernels provide their inflection point, i.e., the
  // point at which their gradient maxes out and starts rolling over.
  SPHERAL_HOST_DEVICE double inflectionPoint() const;

  // Call the descendent Kernel implementations to get the real values.
  // All Kernels are required to define the "kernelValue", "gradValue",
  // and "grad2Value" methods, with the same call signatures 
  // as these functions.
  SPHERAL_HOST_DEVICE double kernelValue(double etaij, const double Hdet) const;
  SPHERAL_HOST_DEVICE double gradValue(double etaij, const double Hdet) const;
  SPHERAL_HOST_DEVICE double grad2Value(double etaij, const double Hdet) const;

protected:
  //--------------------------- Protected Interface ---------------------------//
  // Descendant Kernel classes are allowed (in fact required!) to set the 
  // volume normalization.
  SPHERAL_HOST_DEVICE void setVolumeNormalization(double volumeNormalization);
  SPHERAL_HOST_DEVICE void setKernelExtent(double extent);
  SPHERAL_HOST_DEVICE void setInflectionPoint(double x);

// private:
//   //--------------------------- Private Interface ---------------------------//
  double mVolumeNormalization = 0.0;
  double mKernelExtent = 0.0;
  double mInflectionPoint = 0.0;

};

}

#include "KernelInline.hh"

#endif
