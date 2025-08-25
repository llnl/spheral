//---------------------------------Spheral++----------------------------------//
// HypreFunctions
//
// Convenience functions for Hypre
//----------------------------------------------------------------------------//
#ifndef __Spheral_HypreFunctions_hh__
#define __Spheral_HypreFunctions_hh__

#include <memory>
#include <string>

struct hypre_IJMatrix_struct;
struct hypre_IJVector_struct;
struct hypre_Solver_struct;

namespace Spheral {

class HypreOptions;

// These functions do not check for safety, so make sure objects are initialized
// Here is an example of how these might be called in order:
//   rhs = createVector(...);
//   lhs = createVector(...);
//   mat = createMatrix(...);
//   setMatrixValues(...);
//   assembleMatrix(...);
//   sol = createSolver(...);
//   prec = createPreconditioner(...);
//   setupSolver(...);
//   setVectorValues(...); // rhs
//   setVectorValues(...); // lhs
//   solveSystem(...);
//   getVectorValues(...); // lhs
struct HypreFunctions {
  typedef struct hypre_IJMatrix_struct MatrixType;
  typedef struct hypre_IJVector_struct VectorType;
  typedef struct hypre_Solver_struct SolverType;
  
  // Get Hypre vectors and matrices that self-destruct
  static std::shared_ptr<VectorType> createVector(const unsigned numLocVal,
                                                  const unsigned firstGlobInd);
  static std::shared_ptr<MatrixType> createMatrix(const unsigned numLocVal,
                                                  const unsigned firstGlobInd,
                                                  const unsigned* numValsPerRow);

  // Reinitialize matrix before setting new values
  static void initializeMatrix(std::shared_ptr<MatrixType> matrix);
  
  // Set values in the matrix
  static void setMatrixValues(const unsigned numRows,
                              const unsigned* numColsPerRow,
                              const unsigned* globRowInd,
                              const unsigned* globColInd,
                              const double* colVal,
                              std::shared_ptr<HypreOptions> opt,
                              std::shared_ptr<MatrixType> matrix);
  
  // Set or get the values of a Hypre vector
  static void setVectorValues(const unsigned numVals,
                               const unsigned firstGlobInd,
                               const double* val,
                               std::shared_ptr<VectorType> hypreVector);
  static void getVectorValues(const unsigned numVals,
                               const unsigned firstGlobInd,
                               double* val,
                               std::shared_ptr<VectorType> hypreVector);

  // Assemble the matrix
  static void assembleMatrix(std::shared_ptr<MatrixType> matrix);

  // Get solver and preconditioner that self-destructs
  static std::shared_ptr<SolverType> createSolver(std::shared_ptr<HypreOptions> opt);
  static std::shared_ptr<SolverType> createPreconditioner(std::shared_ptr<HypreOptions> opt,
                                                          std::shared_ptr<SolverType> solver);
  
  // Set remaining solver options and set up preconditioner
  static void setupSolver(std::shared_ptr<HypreOptions> opt,
                          std::shared_ptr<SolverType> solver,
                          std::shared_ptr<MatrixType> matrix,
                          std::shared_ptr<VectorType> lhs,
                          std::shared_ptr<VectorType> rhs);

  // Solve the linear system
  // Return iterations and residual norm
  static std::pair<int, double> solveSystem(std::shared_ptr<HypreOptions> opt,
                                            std::shared_ptr<SolverType> solver,
                                            std::shared_ptr<MatrixType> matrix,
                                            std::shared_ptr<VectorType> lhs,
                                            std::shared_ptr<VectorType> rhs,
                                            std::string description);

  // Multiply the linear system
  static void multiplySystem(std::shared_ptr<HypreOptions> opt,
                             std::shared_ptr<MatrixType> matrix,
                             std::shared_ptr<VectorType> lhs,
                             std::shared_ptr<VectorType> rhs,
                             std::string description);
}; // end struct HypreFunctions

} // end namespace Spheral

#endif
