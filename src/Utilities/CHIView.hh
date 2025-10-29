//---------------------------------Spheral++----------------------------------//
// CHIView
//
// View class of the CubicHermiteInterpolator class
//
// Created by LDO, Wed Oct 29 15:18:00 PST 2025
//----------------------------------------------------------------------------//
#ifndef __Spheral_CHIView__
#define __Spheral_CHIView__

#include "chai/ManagedArray.hpp"
#include "config.hh"

#include <cstddef>
#include <vector>

namespace Spheral {
class CHIView {
public:
  using ContainerType = typename chai::ManagedArray<double>;
  //--------------------------- Public Interface ---------------------------//
  // Constructors, destructors
  SPHERAL_HOST_DEVICE CHIView() = default;
  // Comparisons
  SPHERAL_HOST_DEVICE bool operator==(const CHIView& rhs) const;

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

  // Compute the Hermite basis functions
  SPHERAL_HOST_DEVICE double h00(const double t) const { return (2.0*t - 3.0)*t*t + 1.0; }
  SPHERAL_HOST_DEVICE double h10(const double t) const { return (t - 2.0)*t*t + t; }
  SPHERAL_HOST_DEVICE double h01(const double t) const { return (3.0 - 2.0*t)*t*t; }
  SPHERAL_HOST_DEVICE double h11(const double t) const { return (t - 1.0)*t*t; }

  // Allow read access the internal data representation
  SPHERAL_HOST_DEVICE size_t size() const              { return mN; }
  SPHERAL_HOST_DEVICE double xmin() const              { return mXmin; }
  SPHERAL_HOST_DEVICE double xmax() const              { return mXmax; }
  SPHERAL_HOST_DEVICE double xstep() const             { return mXstep; }
  SPHERAL_HOST_DEVICE double* data() const             { return mVals.data(); }

  void move(chai::ExecutionSpace space)                { mVals.move(space); }

protected:
  //--------------------------- Protected Interface --------------------------//
  // Member data
  size_t mN = 0u;
  double mXmin = 0.;
  double mXmax = 0.;
  double mXstep = 0.;
  ContainerType mVals;
};
}
#include "CHIViewInline.hh"

#endif
