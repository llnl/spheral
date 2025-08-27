//---------------------------------Spheral++----------------------------------//
// LinearSolver
//
// Represents a solver for a linear system of equations
// A given type only needs to fill in the pure virtual methods, and the rest
// is done through the setMap and setData functions. 
//----------------------------------------------------------------------------//
#ifndef __Spheral_LinearSolver_hh__
#define __Spheral_LinearSolver_hh__

#include <memory>
#include <string>
#include <vector>

namespace Spheral {

template<class DataType> class IncrementalStatistic;

/* To set up the class for a solve (each time the matrix changes)
     initialize(...);
     for (auto i = 0; i < numLocalRows; ++i) {
       setMatRow(...);
     }
     assemble(...);
     finalize(...);
   To solve or multiply (as many times as needed after initialization)
     setLHS(...); // Set guess
     setRHS(...); // Set source
     solve();     // Solve in place
     getLHS(...); // Get solution

*/
class LinearSolver {
public:
  // Specify which part of the system to apply set/get to
  enum Component {
    LHS,
    RHS,
  };
    
  // Constructor
  LinearSolver();

  // Create matrices and vectors
  virtual void initialize(const unsigned numLocVal,
                          const unsigned firstGlobInd,
                          const unsigned* numValsPerRow) = 0;

  // Begin matrix fill
  virtual void beginFill() = 0;

  // Set values in matrix
  virtual void setMatRow(const unsigned numVals,     // Number of values in this row
                         const unsigned globRowInd,  // Global index of this row
                         const unsigned* globColInd, // [numVals] Global indices for columns
                         const double* colVal) = 0;  // [numVals] Values
  virtual void setMatRows(const unsigned numRows,        // Number of rows
                          const unsigned* numColsPerRow, // [numRows] Number of columns for each of these rows
                          const unsigned* globRowInd,    // [numRows] Global indices for these rows
                          const unsigned* globColInd,    // [numColsPerRow[numRows]] Global column indies for each value
                          const double* colVal) = 0;     // [numColsPerRow[numRows]] Values

  // Assemble matrix after fill
  virtual void assemble() = 0;
  
  // Create solver/preconditioner
  virtual void finalize() = 0;

  // Set/get values in vectors
  virtual void set(const Component c,                // Which vector to apply this to
                   const unsigned numVals,      // Number of vector values to set 
                   const unsigned firstGlobInd, // First global index for values
                   const double* val) = 0;      // [numVals] Values
  virtual void get(const Component c,                // Which vector to apply this to
                   const unsigned numVals,      // Number of vector values to set 
                   const unsigned firstGlobInd, // First global index for values
                   double* val) = 0;            // [numVals] Values

  // Solve the system of equations in place
  virtual void solve() = 0;

  // Multiply the matrix by a vector
  virtual void multiply() = 0;

  // Set/get description
  virtual void setDescription(std::string desc);
  virtual std::string getDescription() const;

  // Return statistics if available, empty vector if not
  virtual std::vector<std::shared_ptr<IncrementalStatistic<double>>> statistics() const;

protected:
  std::string mDescription;
}; // end class LinearSolver

} // end namespace Spheral

#include "LinearSolverInline.hh"

#endif
