//---------------------------------Spheral++----------------------------------//
// SolverFunction -- functor abstract class defining interface for solvers
// calling the function to get residuals for minimization
// 
// Created by JMO, Mon Mar 10 15:54:37 PDT 2025
//----------------------------------------------------------------------------//
#ifndef __Spheral_SolverFunction__
#define __Spheral_SolverFunction__

#include <vector>

namespace Spheral {

class SolverFunction {

public:
  //--------------------------- Public Interface ---------------------------//
  // Constructors, destructor
  SolverFunction(size_t numUnknowns = 0u): mNumUnknowns(numUnknowns) {}
  virtual ~SolverFunction() = default;

  // Functor method to be overriden with the function to compute residuals
  virtual void operator()(std::vector<double>& residuals,            // residual output array
                          const std::vector<double>& x) const = 0;   // input array of current unknowns

  // Accessors
  size_t numUnknowns()               const { return mNumUnknowns; }
  void   numUnknowns(const size_t x)       { mNumUnknowns = x; }

private:
  //---------------------------  Private Interface ---------------------------//
  size_t mNumUnknowns;
};

}

#endif
