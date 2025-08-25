//---------------------------------Spheral++----------------------------------//
// HypreLinearSolver
//
// Represents a solver for a linear system of equations
//----------------------------------------------------------------------------//
#ifndef __Spheral_HypreLinearSolver_hh__
#define __Spheral_HypreLinearSolver_hh__

#include <memory>
#include <vector>

#include "LinearSolver.hh"

struct hypre_IJMatrix_struct;
struct hypre_IJVector_struct;
struct hypre_Solver_struct;

namespace Spheral {

class HypreOptions;
template<typename DataType> class IncrementalStatistic;

class HypreLinearSolver : public LinearSolver {
public:
  typedef struct hypre_IJMatrix_struct MatrixType;
  typedef struct hypre_IJVector_struct VectorType;
  typedef struct hypre_Solver_struct SolverType;

  using LinearSolver::Component;
  
  // Constructor
  HypreLinearSolver(std::shared_ptr<HypreOptions> options);

  // Create matrices and vectors
  virtual void initialize(const unsigned numLocVal,
                          const unsigned firstGlobInd,
                          const unsigned* numValsPerRow) override;

  // Begin the matrix fill
  virtual void beginFill() override;
  
  // Set values in matrix
  virtual void setMatRow(const unsigned numVals,     // Number of values in this row
                         const unsigned globRowInd,  // Global index of this row
                         const unsigned* globColInd, // [numVals] Global indices for columns
                         const double* colVal) override;  // [numVals] Values
  virtual void setMatRows(const unsigned numRows,        // Number of rows
                          const unsigned* numColsPerRow, // [numRows] Number of columns for each of these rows
                          const unsigned* globRowInd,    // [numRows] Global indices for these rows
                          const unsigned* globColInd,    // [numColsPerRow[numRows]] Global column indies for each value
                          const double* colVal) override;     // [numColsPerRow[numRows]] Values

  // Assemble matrix after fill
  virtual void assemble() override;
  
  // Create solver/preconditioner
  virtual void finalize() override;

  // Set/get values
  virtual void set(const Component c,
                   const unsigned numVals,      // Number of vector values to set 
                   const unsigned firstGlobInd, // First global index for values
                   const double* val) override; // [numVals] Values
  virtual void get(const Component c,
                   const unsigned numVals,      // Number of vector values to set 
                   const unsigned firstGlobInd, // First global index for values
                   double* val) override;       // [numVals] Values
  
  // Solve the system of equations in place
  virtual void solve() override;

  // Multiply the matrix by a vector
  virtual void multiply() override;

  // Set the description
  virtual void setDescription(std::string desc) override;
  
  // Return statistics
  virtual std::vector<std::shared_ptr<IncrementalStatistic<double>>> statistics() const override;
  
private:
  // Which stages have we completed?
  bool mInitialized;
  bool mFilling;
  bool mAssembled;
  bool mFinalized;
  
  // Data
  std::shared_ptr<HypreOptions> mHypreOptions;
  std::shared_ptr<MatrixType> mHypreMatrix;
  std::shared_ptr<VectorType> mHypreLHS;
  std::shared_ptr<VectorType> mHypreRHS;
  std::shared_ptr<SolverType> mHypreSolver;
  std::shared_ptr<SolverType> mHyprePreconditioner;

  // Iteration statistics
  std::shared_ptr<IncrementalStatistic<double>> mIterationStatistics;
  std::shared_ptr<IncrementalStatistic<double>> mFinalResidualStatistics;
}; // end class HypreLinearSolver

} // end namespace Spheral

#include "HypreLinearSolverInline.hh"

#endif
