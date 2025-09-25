//---------------------------------Spheral++----------------------------------//
// CubicHermiteInterpolator
//
// An (optionally monotonic) form of cubic Hermite interpolation.
//
// Created by JMO, Fri Apr  1 14:22:04 PDT 2022
//----------------------------------------------------------------------------//
#ifndef __Spheral_CubicHermiteInterpolator__
#define __Spheral_CubicHermiteInterpolator__

#include "chai/ManagedArray.hpp"
#include "config.hh"

#include <cstddef>
#include <vector>

namespace Spheral {
class CHIView : public chai::CHAICopyable {
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
  SPHERAL_HOST_DEVICE double h00(const double x) const;
  SPHERAL_HOST_DEVICE double h10(const double x) const;
  SPHERAL_HOST_DEVICE double h01(const double x) const;
  SPHERAL_HOST_DEVICE double h11(const double x) const;

  // Allow read access the internal data representation
  SPHERAL_HOST_DEVICE size_t size() const;                        // The number of tabulated values
  SPHERAL_HOST_DEVICE double xmin() const;                        // Minimum x coordinate for table
  SPHERAL_HOST_DEVICE double xmax() const;                        // Maximum x coordinate for table
  SPHERAL_HOST_DEVICE double xstep() const;                       // delta x between tabulated values
  SPHERAL_HOST void move(chai::ExecutionSpace space) { mVals.move(space); }
  SPHERAL_HOST_DEVICE double* data() const { return mVals.data(); }
  SPHERAL_HOST_DEVICE void touch(chai::ExecutionSpace space) {
    mVals.registerTouch(space);
  }

protected:
  //--------------------------- Protected Interface --------------------------//
  // Member data
  size_t mN = 0u;
  double mXmin = 0.;
  double mXmax = 0.;
  double mXstep = 0.;
  ContainerType mVals;
};

class CubicHermiteInterpolator : public CHIView {
public:
  //--------------------------- Public Interface ---------------------------//
  // Constructors, destructors
  template<typename Func>
  CubicHermiteInterpolator(const double xmin,
                           const double xmax,
                           const size_t n,
                           const Func& F);
  template<typename Func, typename GradFunc>
  CubicHermiteInterpolator(const double xmin,
                           const double xmax,
                           const size_t n,
                           const Func& F,
                           const GradFunc& Fgrad);
  CubicHermiteInterpolator(const double xmin,
                           const double xmax,
                           const std::vector<double>& values);
  CubicHermiteInterpolator(const CubicHermiteInterpolator& rhs);
  CubicHermiteInterpolator& operator=(const CubicHermiteInterpolator& rhs);
  CubicHermiteInterpolator() = default;
  ~CubicHermiteInterpolator();

  // (Re)initialize after construction, same options as construction
  template<typename Func>
  void initialize(const double xmin,
                  const double xmax,
                  const size_t n,
                  const Func& F);
  template<typename Func, typename GradFunc>
  void initialize(const double xmin,
                  const double xmax,
                  const size_t n,
                  const Func& F,
                  const GradFunc& Fgrad);
  void initialize(const double xmin,
                  const double xmax,
                  const std::vector<double>& yvals);

  // Force interpolation to be monotonic (may introduce structure between tabulated points)
  void makeMonotonic();

  CHIView view() { return static_cast<CHIView>(*this); }

  template<typename F> inline
  void setUserCallback(F&& extension) {
    mVals.setUserCallback(getNPLCallback(std::forward<F>(extension)));
  }

protected:
  template<typename F>
  auto getNPLCallback(F callback) {
    return [callback](
      const chai::PointerRecord * record,
      chai::Action action,
      chai::ExecutionSpace space) {
             callback(record, action, space);
           };
  }
private:
  //--------------------------- Private Interface --------------------------//
  // Initialize the gradient at the interpolation points based on the tabulated
  // interpolation values
  std::vector<double> mVec;
  void initializeGradientKnots();
  void initView();
};
}
#include "CubicHermiteInterpolatorInline.hh"

#endif
