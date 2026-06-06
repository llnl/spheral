//---------------------------------Spheral++----------------------------------//
// KINSOL -- Wrapper around the Sundials KINSOL non-linear solver
// 
// Created by JMO, Mon Mar 10 15:54:37 PDT 2025
//----------------------------------------------------------------------------//
#ifndef __Spheral_KINSOL__
#define __Spheral_KINSOL__

#include "Solvers/SolverFunction.hh"

#include <sundials/sundials_types.h> /* defs. of sunrealtype, sunindextype      */

#include <functional>

namespace Spheral {

class KINSOL {

public:
  //--------------------------- Public Interface ---------------------------//
  // Constructors, destructor
  KINSOL();
  virtual ~KINSOL() = default;

  // Solve based on an initial guess stored in an input array.  On exit the solution
  // is returned in the input arrray.
  virtual size_t solve(SolverFunction& func,
                       std::vector<double>& x);

  // Accessors
  int globalstrategy()                    const { return mglobalstrategy; }
  sunrealtype fnormtol()                  const { return mfnormtol; }
  sunrealtype scsteptol()                 const { return mscsteptol; }
  long int numMaxIters()                  const { return mNumMaxIters; }

  void globalstrategy(const int x)              { mglobalstrategy = x; }
  void fnormtol(const sunrealtype x)            { mfnormtol = x; }
  void scsteptol(const sunrealtype x)           { mscsteptol = x; }
  void numMaxIters(const long int x)            { mNumMaxIters = x; }

private:
  //---------------------------  Private Interface ---------------------------//
  SUNContext mctx;
  int mglobalstrategy;
  sunrealtype mfnormtol, mscsteptol;
  long int mNumMaxIters;
};

}

#endif
