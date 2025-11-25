//---------------------------------Spheral++----------------------------------//
// QuadraticInterpolator
//
// Encapsulates the algorithm and data for parabolic interpolation in 1D
// Assumes the results is interpolated as y_interp = a + b*x + c*x^2
//
// Created by JMO, Fri Dec  4 14:28:08 PST 2020
//----------------------------------------------------------------------------//
#ifndef __Spheral_QuadraticInterpolator__
#define __Spheral_QuadraticInterpolator__

#include "QuadraticInterpolatorView.hh"
#include "chai/ManagedArray.hpp"
#include "config.hh"

#include <cstddef>
#include <vector>

namespace Spheral {

class QuadraticInterpolator : public QuadraticInterpolatorView {
public:
  template<typename Func>
  QuadraticInterpolator(double xmin, double xmax, size_t n, const Func& F);
  QuadraticInterpolator(double xmin, double xmax, const std::vector<double>& yvals);
  QuadraticInterpolator() = default;
  ~QuadraticInterpolator();
  QuadraticInterpolator(const QuadraticInterpolator& rhs);
  QuadraticInterpolator& operator=(const QuadraticInterpolator& rhs);

  // Initialize after construction, either with a function or tabulated values
  template<typename Func>
  void initialize(double xmin, double xmax, size_t n, const Func& f);
  void initialize(double xmin, double xmax, const std::vector<double>& yvals);

  QuadraticInterpolatorView view() { return static_cast<QuadraticInterpolatorView>(*this); }

  template<typename F> inline
  void setUserCallback(F&& extension) {
    mcoeffs.setUserCallback(getNPLCallback(std::forward<F>(extension)));
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
  std::vector<double> mVec;
  void initView();
};
}

#include "QuadraticInterpolatorInline.hh"

#endif
