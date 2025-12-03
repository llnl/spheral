//---------------------------------Spheral++----------------------------------//
// CubicHermiteInterpolator
//
// An (optionally monotonic) form of cubic Hermite interpolation.
//
// Created by JMO, Fri Apr  1 14:22:04 PDT 2022
//----------------------------------------------------------------------------//
#ifndef __Spheral_CubicHermiteInterpolator__
#define __Spheral_CubicHermiteInterpolator__

#include "CubicHermiteInterpolatorView.hh"
#include "chai/ManagedArray.hpp"
#include "chai/config.hpp"
#include "config.hh"

#include <cstddef>
#include <vector>

namespace Spheral {
class CubicHermiteInterpolator : public CubicHermiteInterpolatorView {
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

  CubicHermiteInterpolatorView view() { return static_cast<CubicHermiteInterpolatorView>(*this); }

#ifndef CHAI_DISABLE_RM
  template<typename F> inline
  void setUserCallback(F&& extension) {
    mVals.setUserCallback(getNPLCallback(std::forward<F>(extension)));
  }
#endif

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
