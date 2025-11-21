//---------------------------------Spheral++----------------------------------//
// QuadraticInterpolatorView
//
// View class for the QuadraticInterpolator class
// Created by LDO, Wed Oct 29 15:00:00 PST 2025
//----------------------------------------------------------------------------//
#ifndef __Spheral_QuadraticInterpolatorView__
#define __Spheral_QuadraticInterpolatorView__

#include "chai/ManagedArray.hpp"
#include "config.hh"

namespace Spheral {

class QuadraticInterpolatorView {
public:
  using ContainerType = typename chai::ManagedArray<double>;
  //--------------------------- Public Interface ---------------------------//
  // Constructors, destructors
  SPHERAL_HOST_DEVICE QuadraticInterpolatorView() = default;

  // Comparisons
  SPHERAL_HOST_DEVICE bool operator==(const QuadraticInterpolatorView& rhs) const;

  // Interpolate for the y value
  SPHERAL_HOST_DEVICE double operator()(const double x) const;
  SPHERAL_HOST_DEVICE double prime(const double x) const;    // First derivative
  SPHERAL_HOST_DEVICE double prime2(const double x) const;   // Second derivative
  // Index access
  SPHERAL_HOST_DEVICE double operator[](const size_t i) const;

  // Same as above, but use a pre-computed table position (from lowerBound)
  SPHERAL_HOST_DEVICE double operator()(const double x, const size_t i0) const;
  SPHERAL_HOST_DEVICE double prime(const double x, const size_t i0) const;    // First derivative
  SPHERAL_HOST_DEVICE double prime2(const double x, const size_t i0) const;   // Second derivative

  // Return the lower bound index in the table for the given x coordinate
  SPHERAL_HOST_DEVICE size_t lowerBound(const double x) const;

  // Allow read access the internal data representation
  SPHERAL_HOST_DEVICE size_t size() const  { return 3*(mN1 + 1u); }
  SPHERAL_HOST_DEVICE double xmin() const  { return mXmin; }
  SPHERAL_HOST_DEVICE double xmax() const  { return mXmax; }
  SPHERAL_HOST_DEVICE double xstep() const { return mXstep; }
  SPHERAL_HOST_DEVICE double* data() const { return mcoeffs.data(); }

  void move(chai::ExecutionSpace space)    { mcoeffs.move(space); }

protected:
  //--------------------------- Private Interface --------------------------//
  // Member data
  size_t mN1 = 0u;
  double mXmin = 0.;
  double mXmax = 0.;
  double mXstep = 0.;
  ContainerType mcoeffs;
};
}

#include "QuadraticInterpolatorViewInline.hh"

#endif
