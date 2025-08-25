//---------------------------------Spheral++----------------------------------//
// EigenLinearSolver
//
// Represents a solver for a linear system of equations
//----------------------------------------------------------------------------//

#include "EigenLinearSolver.hh"

#include <map>
#include <numeric> // for std::accumulate
#include "Utilities/DBC.hh"
#include "Distributed/Process.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
EigenLinearSolver::
EigenLinearSolver(std::shared_ptr<EigenOptions> options) :
  mGraphChangedSinceFill(true),
  mMatrixChangedSinceFactorization(true),
  mOptions(options) {
  this->setDescription("EigenLinearSolver");
}

//------------------------------------------------------------------------------
// Initialize the matrices and vectors
//------------------------------------------------------------------------------
void
EigenLinearSolver::
initialize(const unsigned numLocVal,
           const unsigned firstGlobInd,
           const unsigned* numValsPerRow) {
  VERIFY(Process::getTotalNumberOfProcesses() == 1);

  mMatrix.setZero();
  mMatrix.resize(numLocVal, numLocVal);
  mRhs.resize(numLocVal);
  mLhs.resize(numLocVal);

  mGraphChangedSinceFill = true;
}

//------------------------------------------------------------------------------
// Set matrix values
//------------------------------------------------------------------------------
void
EigenLinearSolver::
setMatRow(const unsigned numVals,
          const unsigned globRowInd,
          const unsigned* globColInd,
          const double* colVal) {
  VERIFY(globRowInd < mMatrix.rows());
  for (auto j = 0u; j < numVals; ++j) {
    mMatrix.coeffRef(globRowInd, globColInd[j]) = colVal[j];
  }
  mMatrixChangedSinceFactorization = true;
}

void
EigenLinearSolver::
setMatRows(const unsigned numRows,
           const unsigned* numColsPerRow,
           const unsigned* globRowInd,
           const unsigned* globColInd,
           const double* colVal) {
  VERIFY(numRows > 0);
  auto index = 0u;
  for (auto i = 0u; i < numRows; ++i) {
    setMatRow(numColsPerRow[i], globRowInd[i], globColInd + index, colVal + index);
    index += numColsPerRow[i];
  }
  mMatrixChangedSinceFactorization = true;
}

//------------------------------------------------------------------------------
// Assemble matrix after fill
//------------------------------------------------------------------------------
void
EigenLinearSolver::
assemble() {
  mMatrix.makeCompressed();
  
  mGraphChangedSinceFill = false;
  mMatrixChangedSinceFactorization = true;
}

//------------------------------------------------------------------------------
// Create solver and preconditioner
//------------------------------------------------------------------------------
void
EigenLinearSolver::
finalize() {
  VERIFY(!mGraphChangedSinceFill);

  if (mOptions->qr) {
    mQRSolver.analyzePattern(mMatrix);
    mQRSolver.factorize(mMatrix);
  }
  else {
    mLUSolver.analyzePattern(mMatrix);
    mLUSolver.factorize(mMatrix);
  }

  mMatrixChangedSinceFactorization = false;
}

//------------------------------------------------------------------------------
// Set/get values in LHS and RHS
//------------------------------------------------------------------------------
void
EigenLinearSolver::
set(const Component c,
    const unsigned numVals,
    const unsigned firstGlobInd,
    const double* val) {
  switch (c) {
  case RHS:
    for (auto i = 0u; i < numVals; ++i) {
      mRhs(firstGlobInd + i) = val[i];
    }
    break;
  case LHS:
    for (auto i = 0u; i < numVals; ++i) {
      mLhs(firstGlobInd + i) = val[i];
    }
    break;
  }
}

void
EigenLinearSolver::
get(const Component c,
    const unsigned numVals,
    const unsigned firstGlobInd,
    double* val) {
  switch (c) {
  case RHS:
    for (auto i = 0u; i < numVals; ++i) {
      val[i] = mRhs(firstGlobInd + i);
    }
    break;
  case LHS:
    for (auto i = 0u; i < numVals; ++i) {
      val[i] = mLhs(firstGlobInd + i);
    }
    break;
  }
}

//------------------------------------------------------------------------------
// Solve the system of equations in place
//------------------------------------------------------------------------------
void
EigenLinearSolver::
solve() {
  VERIFY(!mMatrixChangedSinceFactorization);
  if (mOptions->qr) {
    mLhs = mQRSolver.solve(mRhs);
  }
  else {
    mLhs = mLUSolver.solve(mRhs);
  }
}

//------------------------------------------------------------------------------
// Multiply by the matrix
//------------------------------------------------------------------------------
void
EigenLinearSolver::
multiply() {
  mLhs = mMatrix * mRhs;
}

} // end namespace Spheral
