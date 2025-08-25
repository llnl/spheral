//---------------------------------Spheral++----------------------------------//
// HypreLinearSolver
//
// Represents a solver for a linear system of equations
//----------------------------------------------------------------------------//

#include "HypreLinearSolver.hh"

#include "Solvers/HypreFunctions.hh"
#include "Solvers/HypreOptions.hh"
#include "Solvers/IncrementalStatistic.hh"
#include "Utilities/DBC.hh"
#include "Utilities/Timer.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
HypreLinearSolver::
HypreLinearSolver(std::shared_ptr<HypreOptions> options) :
  mInitialized(false),
  mFilling(false),
  mAssembled(false),
  mFinalized(false),
  mHypreOptions(options),
  mIterationStatistics(
    std::make_shared<IncrementalStatistic<double>>(options->meanIterationsGuess,
                                                   "HypreLinearSolver iterations",
                                                   options->printIterations)), // print
  mFinalResidualStatistics(
    std::make_shared<IncrementalStatistic<double>>(options->toleranceL2,
                                                   "HypreLinearSolver residual",
                                                   false)) { // print
  this->setDescription("HypreLinearSolver");
}

//------------------------------------------------------------------------------
// Initialize the matrices and vectors
//------------------------------------------------------------------------------
void
HypreLinearSolver::
initialize(const unsigned numLocVal,
           const unsigned firstGlobInd,
           const unsigned* numValsPerRow) {
  TIME_SCOPE("HypreLinearSolver::initialize");
  mInitialized = false;
  mAssembled = false;
  mFinalized = false;

  mHypreLHS = HypreFunctions::createVector(numLocVal, firstGlobInd);
  mHypreRHS = HypreFunctions::createVector(numLocVal, firstGlobInd);
  mHypreMatrix = HypreFunctions::createMatrix(numLocVal, firstGlobInd, numValsPerRow);
  mHypreSolver = HypreFunctions::createSolver(mHypreOptions);
  mHyprePreconditioner = HypreFunctions::createPreconditioner(mHypreOptions, mHypreSolver);

  mInitialized = true;
}

void
HypreLinearSolver::
beginFill() {
  TIME_SCOPE("HypreLinearSolver::beginFill");
  VERIFY(mInitialized);
  HypreFunctions::initializeMatrix(mHypreMatrix);
  mFilling = true;
}

//------------------------------------------------------------------------------
// Set matrix values
//------------------------------------------------------------------------------
void
HypreLinearSolver::
setMatRow(const unsigned numVals,
          const unsigned globRowInd,
          const unsigned* globColInd,
          const double* colVal) {
  TIME_SCOPE("HypreLinearSolver::setMatRow");
  VERIFY(mFilling);
  HypreFunctions::setMatrixValues(1, &numVals, &globRowInd,
                                  globColInd, colVal,
                                  mHypreOptions, mHypreMatrix);
}

void
HypreLinearSolver::
setMatRows(const unsigned numRows,
           const unsigned* numColsPerRow,
           const unsigned* globRowInd,
           const unsigned* globColInd,
           const double* colVal) {
  TIME_SCOPE("HypreLinearSolver::setMatRows");
  VERIFY(mFilling);
  HypreFunctions::setMatrixValues(numRows, numColsPerRow, globRowInd,
                                  globColInd, colVal,
                                  mHypreOptions, mHypreMatrix);
}

//------------------------------------------------------------------------------
// Assemble matrix after fill
//------------------------------------------------------------------------------
void
HypreLinearSolver::
assemble() {
  TIME_SCOPE("HypreLinearSolver::assemble");
  VERIFY(mFilling);
  HypreFunctions::assembleMatrix(mHypreMatrix);
  mFilling = false;
  mAssembled = true;
}

//------------------------------------------------------------------------------
// Create solver and preconditioner
//------------------------------------------------------------------------------
void
HypreLinearSolver::
finalize() {
  TIME_SCOPE("HypreLinearSolver::finalize");
  VERIFY(mAssembled);
  HypreFunctions::setupSolver(mHypreOptions, mHypreSolver, mHypreMatrix, mHypreLHS, mHypreRHS);
  mFinalized = true;
}

//------------------------------------------------------------------------------
// Set/get values in LHS and RHS
//------------------------------------------------------------------------------
void
HypreLinearSolver::
set(const Component c,
    const unsigned numVals,
    const unsigned firstGlobInd,
    const double* val) {
  TIME_SCOPE("HypreLinearSolver::set");
  VERIFY(mInitialized);
  switch (c) {
  case RHS:
    HypreFunctions::setVectorValues(numVals, firstGlobInd, val, mHypreRHS);
    break;
  case LHS:
    HypreFunctions::setVectorValues(numVals, firstGlobInd, val, mHypreLHS);
    break;
  }
}

void
HypreLinearSolver::
get(const Component c,
    const unsigned numVals,
    const unsigned firstGlobInd,
    double* val) {
  TIME_SCOPE("HypreLinearSolver::get");
  VERIFY(mInitialized);
  switch (c) {
  case RHS:
    HypreFunctions::getVectorValues(numVals, firstGlobInd, val, mHypreRHS);
    break;
  case LHS:
    HypreFunctions::getVectorValues(numVals, firstGlobInd, val, mHypreLHS);
    break;
  }
}

//------------------------------------------------------------------------------
// Solve the system of equations in place
//------------------------------------------------------------------------------
void
HypreLinearSolver::
solve() {
  TIME_SCOPE("HypreLinearSolver::solve");
  VERIFY(mFinalized);
  auto [it, norm] = HypreFunctions::solveSystem(mHypreOptions, mHypreSolver, mHypreMatrix, mHypreLHS, mHypreRHS, this->getDescription());
  mIterationStatistics->add(it);
  mFinalResidualStatistics->add(norm);
}

//------------------------------------------------------------------------------
// Multiply by the matrix
//------------------------------------------------------------------------------
void
HypreLinearSolver::
multiply() {
  TIME_SCOPE("HypreLinearSolver::multiply");
  VERIFY(mAssembled);
  HypreFunctions::multiplySystem(mHypreOptions, mHypreMatrix, mHypreLHS, mHypreRHS, this->getDescription());
}

//------------------------------------------------------------------------------
// Set the description, including for the stats
//------------------------------------------------------------------------------
inline
void
HypreLinearSolver::
setDescription(std::string desc) {
  this->mDescription = desc;
  mIterationStatistics->setName(desc + " iter");
  mFinalResidualStatistics->setName(desc + " norm");
}


} // end namespace Spheral
