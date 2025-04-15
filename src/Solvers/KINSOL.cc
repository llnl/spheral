//---------------------------------Spheral++----------------------------------//
// KINSOL -- Wrapper around the Sundials KINSOL non-linear solver
// 
// Created by JMO, Mon Mar 10 15:54:37 PDT 2025
//----------------------------------------------------------------------------//

#include "Solvers/KINSOL.hh"
#include "Distributed/Communicator.hh"
#include "Distributed/allReduce.hh"

#include <kinsol/kinsol.h> /* access to KINSOL func., consts.      */
#include <sundials/sundials_nvector.h>
#include <nvector/nvector_parallel.h> /* access to parallel (MPI) N_Vector            */
#include <sundials/sundials_dense.h> /* use generic dense solver in precond. */
#include <sunlinsol/sunlinsol_spbcgs.h> /* access to SPBCGS SUNLinearSolver     */
#include <sunlinsol/sunlinsol_spfgmr.h> /* access to SPFGMR SUNLinearSolver     */
#include <sunlinsol/sunlinsol_spgmr.h> /* access to SPGMR SUNLinearSolver      */
#include <sunlinsol/sunlinsol_sptfqmr.h> /* access to SPTFQMR SUNLinearSolver    */

namespace Spheral {

namespace { // anonymous
//------------------------------------------------------------------------------
// Adapter static function so KINSOL can call a SolverFunction functor
// object as though its a static C function
//------------------------------------------------------------------------------
inline
int SpheralKINSOLfunc(N_Vector xSD, N_Vector fSD, void *user_data) {
  SolverFunction* fptr = reinterpret_cast<SolverFunction*>(user_data);
  const auto n = fptr->numUnknowns();
  CHECK(size_t(NV_LOCLENGTH_P(xSD)) == n and size_t(NV_LOCLENGTH_P(fSD)) == n);
  std::vector<double> x(n), res(n, 0.0);
  for (auto i = 0u; i < n; ++i) x[i] = NV_Ith_P(xSD, i);
  (*fptr)(res, x);
  for (auto i = 0u; i < n; ++i) NV_Ith_P(fSD, i) = res[i];
  return 0;
}

/*
 * Print final statistics contained in iopt
 */

static void PrintFinalStats(void* kmem)
{
  long int nni, nfe, nli, npe, nps, ncfl, nfeSG;
  int flag;
  CONTRACT_VAR(flag);

  flag = KINGetNumNonlinSolvIters(kmem, &nni);
  flag = KINGetNumFuncEvals(kmem, &nfe);

  flag = KINGetNumLinIters(kmem, &nli);
  flag = KINGetNumPrecEvals(kmem, &npe);
  flag = KINGetNumPrecSolves(kmem, &nps);
  flag = KINGetNumLinConvFails(kmem, &ncfl);
  flag = KINGetNumLinFuncEvals(kmem, &nfeSG);

  printf("Final Statistics.. \n");
  printf("NumNLSolves    = %5ld    NumLinIter = %5ld\n", nni, nli);
  printf("NumFuncEval    = %5ld    NumLinFevl = %5ld\n", nfe, nfeSG);
  printf("NumPrecSolv    = %5ld    NumPreEval = %5ld     numLinConvFail = %5ld\n", nps, npe, ncfl);
  printf("\n================================================================================\n\n");
}

}  // anonymous
  
//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
KINSOL::KINSOL():
  mctx(nullptr),
  mglobalstrategy(KIN_NONE),
  mfnormtol(1e-8),
  mscsteptol(1e-13),
  mNumMaxIters(200) {
  auto& comm = Communicator::communicator();
  SUNContext_Create(comm, &mctx);
}

//------------------------------------------------------------------------------
// Solve the given function based on an initial guess
//------------------------------------------------------------------------------
size_t
KINSOL::solve(SolverFunction& func,
              std::vector<double>& initialGuess) {

  // Check the size of our system
  const auto nloc = func.numUnknowns();
  const auto nglob = allReduce(nloc, SPHERAL_OP_SUM);
  REQUIRE2(initialGuess.size() == nloc, "KINSOL ERROR: must specify initial guess of (local, global) dimension (" << nloc << ", " << nglob << ")");

  // Grab KINSOL memory
  void* mkmem = KINCreate(mctx);
  VERIFY2(mkmem != nullptr, "KINSOL error in KINCreate");
  int flag = KINSetFuncNormTol(mkmem, mfnormtol);
  VERIFY2(flag == KIN_SUCCESS, "KINSOL error setting fnormtol: " << flag);
  flag = KINSetScaledStepTol(mkmem, mscsteptol);
  VERIFY2(flag == KIN_SUCCESS, "KINSOL error setting scsteptol: " << flag);
  flag = KINSetNumMaxIters(mkmem, mNumMaxIters);
  VERIFY2(flag == KIN_SUCCESS, "KINSOL error setting max iterations: " << flag);

  // Initial residuals
  std::vector<double> residuals0(nloc);
  func(residuals0, initialGuess);
  CHECK(residuals0.size() == nloc);

  // Prepare the vectors for unknowns (and scalings)
  auto& comm = Communicator::communicator();
  N_Vector mXvec = N_VMake_Parallel(comm, nloc, nglob, initialGuess.data(), mctx);  // Reusing memory in our input std::vector
  N_Vector mUscale = N_VNew_Parallel(comm, nloc, nglob, mctx);
  N_Vector mFscale = N_VNew_Parallel(comm, nloc, nglob, mctx);
  for (auto i = 0u; i < nloc; ++i) {
    NV_Ith_P(mUscale, i) = std::max(std::max(1e-3, mscsteptol), std::abs(initialGuess[i])); // std::min(10.0, std::max(1.0, safeInvVar(initialGuess[i], std::max(1e-3, 10.0*mfnormtol))));
    NV_Ith_P(mFscale, i) = std::max(std::max(1e-3, mscsteptol), std::abs(residuals0[i]));
  }

  // Initialize KINSOL, with x as a template for size of the problem.  We also have
  // to do something horrible, and pass our functional evaluation method
  // as a void* in order to reconstitute the object for querying on the other
  // side of the KINSOL C calls.
  flag = KINInit(mkmem, SpheralKINSOLfunc, mXvec);
  VERIFY2(flag == KIN_SUCCESS, "KINSOL error in KINInit: " << flag);
  flag = KINSetUserData(mkmem, &func);
  VERIFY2(flag == KIN_SUCCESS, "KINSOL error in KINSetUserData pointing to func functor: " << flag);

  // Attach a linear solver (GMRES used here)
  int maxl = 15;
  auto LS = SUNLinSol_SPGMR(mXvec, SUN_PREC_NONE, maxl, mctx);
  VERIFY2(LS != nullptr, "KINSOL error constructing GMRES linear solver");
  flag = KINSetLinearSolver(mkmem, LS, NULL);
  VERIFY2(flag == KIN_SUCCESS, "KINSOL error constructing GMRES linear solver: " << flag);

  // Set the maximum number of restarts
  int maxlrst = 2;
  flag = SUNLinSol_SPGMRSetMaxRestarts(LS, maxlrst);
  VERIFY2(flag == KIN_SUCCESS, "KINSOL error setting restart: " << flag);

  // Call KINSol to solve the problem
  flag = KINSol(mkmem,           /* KINSol memory block */
                mXvec,           /* initial guess on input; solution vector */
                mglobalstrategy, /* global strategy choice */
                mUscale,         /* scaling vector, for the variable x */
                mFscale);        /* scaling vector for function values fval */

  // Load the final solution into the initalGuess vector, whether it converged or not
  for (auto i = 0u; i < nloc; ++i) initialGuess[i] = NV_Ith_P(mXvec, i);

  // How many non-linear solves were required?
  long int nni;
  KINGetNumNonlinSolvIters(mkmem, &nni);
  PrintFinalStats(mkmem);

  // Clean up
  KINFree(&mkmem);
  SUNLinSolFree(LS);
  N_VDestroy(mXvec);
  N_VDestroy(mUscale);
  N_VDestroy(mFscale);

  // Return the number of iterations required to solve
  return size_t(nni);

  // Return the final status
  // cerr << "KINSOL << " << flag << " " << KIN_SUCCESS << endl;
  // return (flag >= KIN_SUCCESS);
}

}
