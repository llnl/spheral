//---------------------------------Spheral++----------------------------------//
// PySolverFunction -- interface override of SolverFunction to simplify use
// with Python
// 
// Created by JMO, Wed Mar 12 13:50:41 PDT 2025
//----------------------------------------------------------------------------//
#ifndef __Spheral_PySolverFunction__
#define __Spheral_PySolverFunction__

#include "Solvers/SolverFunction.hh"

namespace Spheral {

class PySolverFunction: public SolverFunction {

public:
  //--------------------------- Public Interface ---------------------------//
  // Constructors, destructor
  PySolverFunction(size_t numUnknowns): SolverFunction(numUnknowns) {}
  virtual ~PySolverFunction() = default;

  // Functor method to be overriden with the function to compute residuals
  virtual void __call__(std::vector<double>& residuals,            // residual output array
                        const std::vector<double>& x) const = 0;   // input array of current unknowns

  // Override the required C++ method
  virtual void operator()(std::vector<double>& residuals,
                          const std::vector<double>& x) const override { this->__call__(residuals, x); }
};

}

#endif
