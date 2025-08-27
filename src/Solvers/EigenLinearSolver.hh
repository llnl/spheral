//---------------------------------Spheral++----------------------------------//
// EigenLinearSolver
//
// Solver linear system using Eigen
// Only works on one processor
//----------------------------------------------------------------------------//
#ifndef __Spheral_EigenLinearSolver_hh__
#define __Spheral_EigenLinearSolver_hh__

#include <memory>
#include <vector>

#include "Eigen/OrderingMethods"
#include "Eigen/Sparse"
#include "Eigen/SparseLU"
#include "Eigen/SparseQR"

#include "EigenOptions.hh"
#include "LinearSolver.hh"

namespace Spheral {

class EigenLinearSolver : public LinearSolver {
public:
  typedef Eigen::SparseMatrix<double> MatrixType;
  typedef Eigen::Triplet<double> DataType;
  typedef Eigen::VectorXd VectorType;
  typedef Eigen::SparseLU<MatrixType, Eigen::COLAMDOrdering<int>> LUSolverType;
  typedef Eigen::SparseQR<MatrixType, Eigen::COLAMDOrdering<int>> QRSolverType;

  using LinearSolver::Component;
  
  // Constructor
  EigenLinearSolver(std::shared_ptr<EigenOptions> options);

  // Create matrices and vectors
  virtual void initialize(const unsigned numLocVal,
                          const unsigned firstGlobInd,
                          const unsigned* numValsPerRow) override;

  // Begin the matrix fill
  virtual void beginFill() override {}
  
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
  
private:
  // Have we initialized things correctly?
  bool mGraphChangedSinceFill;
  bool mMatrixChangedSinceFactorization;
  
  // Input data
  std::shared_ptr<EigenOptions> mOptions;
  
  // Eigen data
  MatrixType mMatrix;
  LUSolverType mLUSolver;
  QRSolverType mQRSolver;
  VectorType mRhs;
  VectorType mLhs;
  
}; // end class EigenLinearSolver

} // end namespace Spheral

#endif
