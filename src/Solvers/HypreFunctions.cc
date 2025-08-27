//---------------------------------Spheral++----------------------------------//
// HypreFunctions
//
// Convenience functions for Hypre
//----------------------------------------------------------------------------//

#include "HypreFunctions.hh"

#include <map>
#include <numeric> // for std::iota

#include "_hypre_parcsr_ls.h"
#include "HYPRE.h"
#include "HYPRE_IJ_mv.h"

#include "Distributed/Communicator.hh"
#include "Solvers/HypreOptions.hh"
#include "Solvers/IncrementalStatistic.hh"
#include "Utilities/DBC.hh"
#include "Utilities/Timer.hh"

static_assert(sizeof(int) == sizeof(HYPRE_Int),
              "if this fails, need to refactor to convert to HYPRE_Int");
static_assert(sizeof(int) == sizeof(HYPRE_BigInt),
              "if this fails, need to refactor to convert to HYPRE_BigInt");

namespace Spheral {

//------------------------------------------------------------------------------
// Create Hypre vector that should self-destruct when it goes out of scope.
//------------------------------------------------------------------------------
std::shared_ptr<HypreFunctions::VectorType>
HypreFunctions::
createVector(const unsigned numLocVal,
             const unsigned firstGlobInd) {
  TIME_FUNCTION;
  const auto lastGlobInd = static_cast<int>(firstGlobInd + numLocVal) - 1;

  int hypreStatus;
  VectorType* tempVector;

  hypreStatus = HYPRE_IJVectorCreate(Communicator::communicator(),
                                     firstGlobInd,
                                     lastGlobInd,
                                     &tempVector);
  VERIFY2(hypreStatus == 0,
          "HYPRE_IJVectorCreate, error code " + std::to_string(hypreStatus));
  
  hypreStatus = HYPRE_IJVectorSetObjectType(tempVector, HYPRE_PARCSR);
  VERIFY2(hypreStatus == 0,
          "HYPRE_IJVectorSetObjectType, error code " + std::to_string(hypreStatus));

  hypreStatus = HYPRE_IJVectorInitialize(tempVector);
  VERIFY2(hypreStatus == 0,
          "HYPRE_IJVectorInitialize, error code " + std::to_string(hypreStatus));
  
  return std::shared_ptr<VectorType>(tempVector, HYPRE_IJVectorDestroy);
}

//------------------------------------------------------------------------------
// Create Hypre matrix that should self-destruct when it goes out of scope.
//------------------------------------------------------------------------------
std::shared_ptr<HypreFunctions::MatrixType>
HypreFunctions::
createMatrix(const unsigned numLocVal,
             const unsigned firstGlobInd,
             const unsigned* numValsPerRow) {
  TIME_FUNCTION;
  const auto lastGlobInd = static_cast<int>(firstGlobInd + numLocVal) - 1;

  int hypreStatus;
  MatrixType* tempMatrix;
  hypreStatus = HYPRE_IJMatrixCreate(Communicator::communicator(),
                                     firstGlobInd,
                                     lastGlobInd,
                                     firstGlobInd,
                                     lastGlobInd,
                                     &tempMatrix);
  VERIFY2(hypreStatus == 0,
          "HYPRE_IJMatrixCreate, error code " + std::to_string(hypreStatus));

  hypreStatus = HYPRE_IJMatrixSetObjectType(tempMatrix, HYPRE_PARCSR);
  VERIFY2(hypreStatus == 0,
          "HYPRE_IJMatrixSetObjectType, error code " + std::to_string(hypreStatus));
  
  hypreStatus = HYPRE_IJMatrixSetRowSizes(tempMatrix, (const int*)numValsPerRow);
  VERIFY2(hypreStatus == 0,
          "HYPRE_IJMatrixSetRowSizes, error code " + std::to_string(hypreStatus));
  
  hypreStatus = HYPRE_IJMatrixInitialize(tempMatrix);
  VERIFY2(hypreStatus == 0,
          "HYPRE_IJMatrixInitialize, error code " + std::to_string(hypreStatus));
  
  return std::shared_ptr<MatrixType>(tempMatrix, HYPRE_IJMatrixDestroy);
}

//------------------------------------------------------------------------------
// Create Hypre matrix that should self-destruct when it goes out of scope.
//------------------------------------------------------------------------------
void
HypreFunctions::
initializeMatrix(std::shared_ptr<MatrixType> matrix) {
  TIME_FUNCTION;
  int hypreStatus;
  hypreStatus = HYPRE_IJMatrixInitialize(matrix.get());
  VERIFY2(hypreStatus == 0,
          "HYPRE_IJMatrixInitialize, error code " + std::to_string(hypreStatus));
}

//------------------------------------------------------------------------------
// Put values into the matrix
//------------------------------------------------------------------------------
void
HypreFunctions::
setMatrixValues(const unsigned numRows,
                const unsigned* numColsPerRow,
                const unsigned* globRowInd,
                const unsigned* globColInd,
                const double* colVal,
                std::shared_ptr<HypreOptions> opt,
                std::shared_ptr<MatrixType> matrix) {
  TIME_FUNCTION;
  if (numRows == 0) {
    return;
  }
  BEGIN_CONTRACT_SCOPE
  {
    auto index = 0u;
    for (auto i = 0u; i < numRows; ++i) {
      for (auto j = 0u; j < numColsPerRow[i]; ++j, ++index) {
        CHECK2(std::isfinite(colVal[index]),
               "NaN found in input to Hypre matrix at index " + std::to_string(globRowInd[i]) + ", " + std::to_string(globColInd[index]));
      }
    }
  }
  END_CONTRACT_SCOPE

  int hypreStatus;
  if (opt->addToValues) {
    hypreStatus = HYPRE_IJMatrixAddToValues(matrix.get(),
                                            numRows,
                                            (int*)numColsPerRow,
                                            (const int*)globRowInd,
                                            (const int*)globColInd,
                                            colVal);
  }
  else {
    hypreStatus = HYPRE_IJMatrixSetValues(matrix.get(),
                                          numRows,
                                          (int*)numColsPerRow,
                                          (const int*)globRowInd,
                                          (const int*)globColInd,
                                          colVal);
  }
  VERIFY2(hypreStatus == 0,
          "HYPRE_IJMatrixSetValues, error code " + std::to_string(hypreStatus));
}

//------------------------------------------------------------------------------
// Put values into hypre vector
//------------------------------------------------------------------------------
void
HypreFunctions::
setVectorValues(const unsigned numVals,
                const unsigned firstGlobInd,
                const double* val,
                std::shared_ptr<VectorType> vec) {
  TIME_FUNCTION;
  if (numVals == 0) {
    return;
  }
  BEGIN_CONTRACT_SCOPE
  {
    for (auto i = 0u; i < numVals; ++i) {
      CHECK2(std::isfinite(val[i]), "NaN found in input to Hypre vector at index " + std::to_string(firstGlobInd + i));
    }
  }
  END_CONTRACT_SCOPE
  
  std::vector<int> globalIndices(numVals);
  std::iota(globalIndices.begin(), globalIndices.end(), firstGlobInd);
  
  int hypreStatus;
  hypreStatus = HYPRE_IJVectorSetValues(vec.get(),
                                        numVals,
                                        &globalIndices[0],
                                        val);
  VERIFY2(hypreStatus == 0,
          "HYPRE_IJVectorSetValues, error code " + std::to_string(hypreStatus));
  
}

//------------------------------------------------------------------------------
// Get values from hypre vector
//------------------------------------------------------------------------------
void
HypreFunctions::
getVectorValues(const unsigned numVals,
                const unsigned firstGlobInd,
                double* val,
                std::shared_ptr<VectorType> vec) {
  TIME_FUNCTION;
  if (numVals == 0) {
    return;
  }
  
  std::vector<int> globalIndices(numVals);
  std::iota(globalIndices.begin(), globalIndices.end(), firstGlobInd);
  
  int hypreStatus;
  hypreStatus = HYPRE_IJVectorGetValues(vec.get(),
                                        numVals,
                                        &globalIndices[0],
                                        val);
  VERIFY2(hypreStatus == 0,
          "HYPRE_IJVectorGetValues, error code " + std::to_string(hypreStatus));

  BEGIN_CONTRACT_SCOPE
  {
    for (auto i = 0u; i < numVals; ++i) {
      CHECK2(std::isfinite(val[i]), "NaN found in Hypre output at index " + std::to_string(firstGlobInd + i));
    }
  }
  END_CONTRACT_SCOPE
}

void
HypreFunctions::
assembleMatrix(std::shared_ptr<MatrixType> matrix) {
  TIME_FUNCTION;
  int hypreStatus;
  hypreStatus = HYPRE_IJMatrixAssemble(matrix.get());
  VERIFY2(hypreStatus == 0,
          "HYPRE_IJMatrixAssemble, error code " + std::to_string(hypreStatus));
}

//------------------------------------------------------------------------------
// Create Hypre solver that should self-destruct when it goes out of scope.
// Right now, the constants and type (GMRES) are hardcoded. At the least, the
// kDim should be changed to a variable. 
//------------------------------------------------------------------------------
std::shared_ptr<HypreFunctions::SolverType>
HypreFunctions::
createSolver(std::shared_ptr<HypreOptions> opt) {
  TIME_FUNCTION;
  int hypreStatus;
  SolverType* tempSolver;
  hypreStatus = HYPRE_ParCSRGMRESCreate(Communicator::communicator(), &tempSolver);
  VERIFY2(hypreStatus == 0,
          "HYPRE_ParCSRGMRESCreate, error code " + std::to_string(hypreStatus));

  hypreStatus = HYPRE_ParCSRGMRESSetKDim(tempSolver, opt->kDim);
  VERIFY2(hypreStatus == 0,
          "HYPRE_PARCSRGMRESSetKDim, error code " + std::to_string(hypreStatus));

  hypreStatus = HYPRE_ParCSRGMRESSetMinIter(tempSolver, opt->minIters);
  VERIFY2(hypreStatus == 0,
          "HYPRE_ParCSRGMRESSetMinIter, error code " + std::to_string(hypreStatus));

  hypreStatus = HYPRE_ParCSRGMRESSetLogging(tempSolver, opt->logLevel);
  VERIFY2(hypreStatus == 0,
          "HYPRE_ParCSRGMRESSetLogging, error code " + std::to_string(hypreStatus));

  hypreStatus = HYPRE_ParCSRGMRESSetPrintLevel(tempSolver, opt->printLevel);
  VERIFY2(hypreStatus == 0,
          "HYPRE_ParCSRMGRESSetPrintLevel, error code " + std::to_string(hypreStatus));

  return std::shared_ptr<SolverType>(tempSolver, HYPRE_ParCSRGMRESDestroy);
}

//------------------------------------------------------------------------------
// Create Hypre preconditioner that should self-destruct when it goes out of
// scope.
// The constants and type (AMG) are hardcoded to start out with.
//------------------------------------------------------------------------------
std::shared_ptr<HypreFunctions::SolverType>
HypreFunctions::
createPreconditioner(std::shared_ptr<HypreOptions> opt,
                     std::shared_ptr<SolverType> solver) {
  TIME_FUNCTION;
  CHECK(solver);

  int hypreStatus;
  SolverType* tempPreconditioner;

  switch (opt->preconditionerType) {
  case HypreOptions::HyprePreconditionerType::NoPreconditioner:
    return std::shared_ptr<SolverType>();
  case HypreOptions::HyprePreconditionerType::AMGPreconditioner:
    hypreStatus = HYPRE_BoomerAMGCreate(&tempPreconditioner);
    VERIFY2(hypreStatus == 0,
            "HYPRE_BoomerAMGCreate, error code " + std::to_string(hypreStatus));
    
    hypreStatus = HYPRE_BoomerAMGSetCoarsenType(tempPreconditioner, opt->coarsenTypeAMG);
    VERIFY2(hypreStatus == 0,
            "HYPRE_BoomerAMGSetCoarsenType, error code " + std::to_string(hypreStatus));

    hypreStatus = HYPRE_BoomerAMGSetMeasureType(tempPreconditioner, opt->measure_type);
    VERIFY2(hypreStatus == 0,
            "HYPRE_BoomerAMGSetMeasureType, error code " + std::to_string(hypreStatus));

    hypreStatus = HYPRE_BoomerAMGSetTol(tempPreconditioner, opt->pcTol);
    VERIFY2(hypreStatus == 0,
            "HYPRE_BoomerAMGSetTol, error code " + std::to_string(hypreStatus));

    hypreStatus
      = HYPRE_BoomerAMGSetStrongThreshold(tempPreconditioner, opt->strongThresholdAMG);
    VERIFY2(hypreStatus == 0,
            "HYPRE_BoomerAMGSetStrongThreshold, error code " + std::to_string(hypreStatus));

    hypreStatus = HYPRE_BoomerAMGSetMaxRowSum(tempPreconditioner, opt->maxRowSumAMG);
    VERIFY2(hypreStatus == 0,
            "HYPRE_BoomerAMGSetMaxRowSum, error code " + std::to_string(hypreStatus));

    if (opt->interpTypeAMG >= 0) {
      hypreStatus = HYPRE_BoomerAMGSetInterpType(tempPreconditioner, opt->interpTypeAMG);
      VERIFY2(hypreStatus == 0,
              "HYPRE_BoomerAMGSetInterpType, error code " + std::to_string(hypreStatus));
    }

    if (opt->aggNumLevelsAMG >= 0) {
      hypreStatus
        = HYPRE_BoomerAMGSetAggNumLevels(tempPreconditioner, opt->aggNumLevelsAMG);
      VERIFY2(hypreStatus == 0,
              "HYPRE_BoomerAMGSetAggNumLevels, error code " + std::to_string(hypreStatus));
    }
    if (opt->aggInterpTypeAMG >= 0) {
      hypreStatus
        = HYPRE_BoomerAMGSetAggInterpType(tempPreconditioner, opt->aggInterpTypeAMG);
      VERIFY2(hypreStatus == 0,
              "HYPRE_BoomerAMGSetAggInterpType, error code " + std::to_string(hypreStatus));
    }
    if (opt->pMaxElmtsAMG >= 0){
      hypreStatus = HYPRE_BoomerAMGSetPMaxElmts(tempPreconditioner, opt->pMaxElmtsAMG);
      VERIFY2(hypreStatus == 0,
              "HYPRE_BoomerAMGSetPMaxElmts, error code " + std::to_string(hypreStatus));
    }

    hypreStatus = HYPRE_BoomerAMGSetTruncFactor(tempPreconditioner, opt->truncFactorAMG);
    VERIFY2(hypreStatus == 0,
            "HYPRE_BoomerAMGSetTruncFactor, error code " + std::to_string(hypreStatus));

    if (opt->printLevelAMG == -1) {
      opt->printLevelAMG = opt->printLevel;
    }
    hypreStatus = HYPRE_BoomerAMGSetPrintLevel(tempPreconditioner, opt->printLevelAMG);
    VERIFY2(hypreStatus == 0,
            "HYPRE_BoomerAMGSetPrintLevel, error code " + std::to_string(hypreStatus));

    hypreStatus = HYPRE_BoomerAMGSetLogging(tempPreconditioner, opt->logLevelAMG);
    VERIFY2(hypreStatus == 0,
            "HYPRE_BoomerAMGSetLogging, error code " + std::to_string(hypreStatus));

    hypreStatus = HYPRE_BoomerAMGSetMinIter(tempPreconditioner, opt->minItersAMG);
    VERIFY2(hypreStatus == 0,
            "HYPRE_BoomerAMGSetMinIter, error code " + std::to_string(hypreStatus));

    hypreStatus = HYPRE_BoomerAMGSetMaxIter(tempPreconditioner, opt->maxItersAMG);
    VERIFY2(hypreStatus == 0,
            "HYPRE_BoomerAMGSetMaxIter, error code " + std::to_string(hypreStatus));

    hypreStatus = HYPRE_BoomerAMGSetCycleType(tempPreconditioner, opt->cycleType);
    VERIFY2(hypreStatus == 0,
            "HYPRE_BoomerAMGSetCycleType, error code " + std::to_string(hypreStatus));

    hypreStatus = HYPRE_BoomerAMGSetRelaxWt(tempPreconditioner, opt->relaxWeightAMG);
    VERIFY2(hypreStatus == 0,
            "HYPRE_BoomerAMGSetRelaxWt, error code " + std::to_string(hypreStatus));

    if (opt->relaxTypeCoarseAMG == -1) {
      opt->relaxTypeCoarseAMG = opt->relaxTypeAMG;
    }

    hypreStatus
      = HYPRE_BoomerAMGSetCycleRelaxType(tempPreconditioner,
                                         opt->relaxTypeAMG, 1);
    VERIFY2(hypreStatus == 0,
            "HYPRE_BoomerAMGSetCycleRelaxType, error code " + std::to_string(hypreStatus));

    hypreStatus
      = HYPRE_BoomerAMGSetCycleRelaxType(tempPreconditioner,
                                         opt->relaxTypeAMG, 2);
    VERIFY2(hypreStatus == 0,
            "HYPRE_BoomerAMGSetCycleRelaxType, error code " + std::to_string(hypreStatus));

    hypreStatus
      = HYPRE_BoomerAMGSetCycleRelaxType(tempPreconditioner,
                                         opt->relaxTypeCoarseAMG, 3);
    VERIFY2(hypreStatus == 0,
            "HYPRE_BoomerAMGSetCycleRelaxType, error code " + std::to_string(hypreStatus));

    hypreStatus
      = HYPRE_BoomerAMGSetCycleNumSweeps(tempPreconditioner,
                                         opt->cycleNumSweepsAMG, 1);
    VERIFY2(hypreStatus == 0,
            "HYPRE_BoomerAMGSetCycleNumSweeps, error code " + std::to_string(hypreStatus));

    hypreStatus
      = HYPRE_BoomerAMGSetCycleNumSweeps(tempPreconditioner,
                                         opt->cycleNumSweepsAMG, 2);
    VERIFY2(hypreStatus == 0,
            "HYPRE_BoomerAMGSetCycleNumSweeps, error code " + std::to_string(hypreStatus));

    if (opt->relaxTypeCoarseAMG == 19
        || opt->relaxTypeCoarseAMG == 29
        || opt->relaxTypeCoarseAMG == 9) {
      opt->cycleNumSweepsCoarseAMG = 1;
    }
    else if (opt->cycleNumSweepsCoarseAMG == -1) {
      opt->cycleNumSweepsCoarseAMG = opt->cycleNumSweepsAMG;
    }
    hypreStatus = HYPRE_BoomerAMGSetCycleNumSweeps(tempPreconditioner,
                                                   opt->cycleNumSweepsCoarseAMG, 3);
    VERIFY2(hypreStatus == 0,
            "HYPRE_BoomerAMGSetCycleNumSweeps, error code " + std::to_string(hypreStatus));

    hypreStatus = HYPRE_BoomerAMGSetMaxLevels(tempPreconditioner, opt->maxLevelsAMG);
    VERIFY2(hypreStatus == 0,
            "HYPRE_BoomerAMGSetMaxLevels, error code " + std::to_string(hypreStatus));

    hypreStatus = HYPRE_ParCSRGMRESSetPrecond(solver.get(),
                                              HYPRE_BoomerAMGSolve,
                                              HYPRE_BoomerAMGSetup,
                                              tempPreconditioner);
    VERIFY2(hypreStatus == 0,
            "HYPRE_ParCSRGMRESSetPrecond, error code " + std::to_string(hypreStatus));
      
    // Set up preconditioner
    hypreStatus = HYPRE_BoomerAMGSetSetupType(tempPreconditioner, 1);
    VERIFY2(hypreStatus == 0,
            "HYPRE_BoomerAMGSetSetupType, error code " + std::to_string(hypreStatus));
  
    return std::shared_ptr<SolverType>(tempPreconditioner, HYPRE_BoomerAMGDestroy);
  case HypreOptions::HyprePreconditionerType::ILUPreconditioner:
    hypreStatus = HYPRE_EuclidCreate(Communicator::communicator(), &tempPreconditioner);
    VERIFY2(hypreStatus == 0,
            "HYPRE_BoomerAMGDestroy, error code " + std::to_string(hypreStatus));
    
    hypreStatus = HYPRE_EuclidSetLevel(tempPreconditioner, opt->factorLevelILU);
    VERIFY2(hypreStatus == 0,
            "HYPRE_EuclidSetLevel, error code " + std::to_string(hypreStatus));

    if (opt->printLevelILU == -1) {
      opt->printLevelILU = opt->printLevel;
    }
    hypreStatus = HYPRE_EuclidSetStats(tempPreconditioner, opt->printLevelILU);
    VERIFY2(hypreStatus == 0,
            "HYPRE_EuclidSetStats, error code " + std::to_string(hypreStatus));

    if (opt->useILUT) {
      hypreStatus = HYPRE_EuclidSetILUT(tempPreconditioner, opt->dropToleranceILU);
    }
    else {
      hypreStatus = HYPRE_EuclidSetSparseA(tempPreconditioner, opt->dropToleranceILU);
      VERIFY2(hypreStatus == 0,
              "HYPRE_EuclidSetILUT, error code " + std::to_string(hypreStatus));
    }

    hypreStatus = HYPRE_EuclidSetRowScale(tempPreconditioner, opt->rowScaleILU);
    VERIFY2(hypreStatus == 0,
            "HYPRE_EuclidSetRowScale, error code " + std::to_string(hypreStatus));
    
    hypreStatus = HYPRE_ParCSRGMRESSetPrecond(solver.get(),
                                              HYPRE_EuclidSolve,
                                              HYPRE_EuclidSetup,
                                              tempPreconditioner);
    VERIFY2(hypreStatus == 0,
            "HYPRE_ParCSRGMRESSetPrecond, error code " + std::to_string(hypreStatus));
    
    return std::shared_ptr<SolverType>(tempPreconditioner, HYPRE_EuclidDestroy);
  default:
    VERIFY2(false, "incorrect preconditioner type");
    return std::shared_ptr<SolverType>(); // So compiler doesn't complain
  }
}

void
HypreFunctions::
setupSolver(std::shared_ptr<HypreOptions> opt,
            std::shared_ptr<SolverType> solver,
            std::shared_ptr<MatrixType> matrix,
            std::shared_ptr<VectorType> lhs,
            std::shared_ptr<VectorType> rhs) {
  TIME_FUNCTION;
  // Get storage in CRS format
  HYPRE_ParCSRMatrix hypreParMatrix;
  HYPRE_ParVector hypreParLHS;
  HYPRE_ParVector hypreParRHS;

  int hypreStatus;
  hypreStatus = HYPRE_IJMatrixGetObject(matrix.get(),
                                        (void **)&hypreParMatrix);
  VERIFY2(hypreStatus == 0,
          "HYPRE_IJMatrixGetObject, error code " + std::to_string(hypreStatus));

  hypreStatus = HYPRE_IJVectorGetObject(lhs.get(),
                                        (void **)&hypreParLHS);
  VERIFY2(hypreStatus == 0,
          "HYPRE_IJVectorGetObject, error code " + std::to_string(hypreStatus));

  hypreStatus = HYPRE_IJVectorGetObject(rhs.get(),
                                        (void **)&hypreParRHS);
  VERIFY2(hypreStatus == 0,
          "HYPRE_IJVectorGetObject, error code " + std::to_string(hypreStatus));
  
  hypreStatus = HYPRE_ParCSRGMRESSetLogging(solver.get(), opt->logLevel);
  VERIFY2(hypreStatus == 0,
          "HYPRE_IJMatrixPrint, error code " + std::to_string(hypreStatus));

  hypreStatus = HYPRE_ParCSRGMRESSetPrintLevel(solver.get(), opt->printLevel);
  VERIFY2(hypreStatus == 0,
          "HYPRE_ParCSRGMRESSetPrintLevel, error code " + std::to_string(hypreStatus));

  hypreStatus
    = HYPRE_ParCSRGMRESSetMaxIter(solver.get(), opt->maxNumberOfIterations);
  VERIFY2(hypreStatus == 0,
          "HYPRE_ParCSRGMRESSetMaxIter, error code " + std::to_string(hypreStatus));

  hypreStatus = HYPRE_ParCSRGMRESSetTol(solver.get(), opt->toleranceL2);
  VERIFY2(hypreStatus == 0,
          "HYPRE_ParCSRGMRESSetTol, error code " + std::to_string(hypreStatus));

  if (opt->useRobustTolerance) {
    hypreStatus = HYPRE_GMRESSetSkipRealResidualCheck(solver.get(), 1);
    VERIFY2(hypreStatus == 0,
            "HYPRE_GMRESSetSkipRealResidualCheck, error code " + std::to_string(hypreStatus));
  }

  hypreStatus
    = HYPRE_ParCSRGMRESSetAbsoluteTol(solver.get(), opt->absoluteTolerance);
  VERIFY2(hypreStatus == 0,
          "HYPRE_ParCSRGMRESSetAbsoluteTol, error code " + std::to_string(hypreStatus));

  hypreStatus = HYPRE_ParCSRGMRESSetup(solver.get(),
                                       hypreParMatrix,
                                       hypreParRHS,
                                       hypreParLHS);
  VERIFY2(hypreStatus == 0,
          "ParCSRGMRESSetup, error code " + std::to_string(hypreStatus));
}

//------------------------------------------------------------------------------
// Solve the linear system
//------------------------------------------------------------------------------
static int solveCall = 0;
std::pair<int, double>
HypreFunctions::
solveSystem(std::shared_ptr<HypreOptions> opt,
            std::shared_ptr<SolverType> solver,
            std::shared_ptr<MatrixType> matrix,
            std::shared_ptr<VectorType> lhs,
            std::shared_ptr<VectorType> rhs,
            std::string description) {
  TIME_FUNCTION;
  // Get storage in CRS format
  HYPRE_ParCSRMatrix hypreParMatrix;
  HYPRE_ParVector hypreParLHS;
  HYPRE_ParVector hypreParRHS;

  int hypreStatus;
  hypreStatus = HYPRE_IJMatrixGetObject(matrix.get(),
                                        (void **)&hypreParMatrix);
  VERIFY2(hypreStatus == 0,
          "HYPRE_IJMatrixGetObject, error code " + std::to_string(hypreStatus));

  hypreStatus = HYPRE_IJVectorGetObject(lhs.get(),
                                        (void **)&hypreParLHS);
  VERIFY2(hypreStatus == 0,
          "HYPRE_IJVectorGetObject, error code " + std::to_string(hypreStatus));

  hypreStatus = HYPRE_IJVectorGetObject(rhs.get(),
                                        (void **)&hypreParRHS);
  VERIFY2(hypreStatus == 0,
          "HYPRE_IJVectorGetObject, error code " + std::to_string(hypreStatus));
  
  if (opt->saveLinearSystem) {
    const auto matName = "HypreMatrix_" + description + "_" + std::to_string(solveCall);
    const auto rhsName = "HypreRHS_" + description + "_" + std::to_string(solveCall);
    const auto initName = "HypreInit_" + description + "_" + std::to_string(solveCall);
    HYPRE_IJMatrixPrint(matrix.get(), matName.c_str());
    HYPRE_IJVectorPrint(rhs.get(), rhsName.c_str());
    HYPRE_IJVectorPrint(lhs.get(), initName.c_str());
  }
  
  int solveStatus = HYPRE_ParCSRGMRESSolve(solver.get(),
                                           hypreParMatrix,
                                           hypreParRHS,
                                           hypreParLHS);
  int numIterations = 0;
  if (solveStatus) {
    VERIFY2(!HYPRE_CheckError(solveStatus, HYPRE_ERROR_GENERIC),
            "INF or NaN in Matrix, RHS, or initial guess (" + description + ")");
    VERIFY2(!HYPRE_CheckError(solveStatus, HYPRE_ERROR_MEMORY),
            "HYPRE ParCSRGMRESSetupcannot allocate memory (" + description + ")");
    VERIFY2(HYPRE_CheckError(solveStatus, HYPRE_ERROR_CONV),
            "GMRES failed with (unknown) solver status (" + description + "):" << solveStatus);
    HYPRE_ClearError(HYPRE_ERROR_CONV);
    numIterations = opt->maxNumberOfIterations;
  }
  else {
    hypreStatus
      = HYPRE_ParCSRGMRESGetNumIterations(solver.get(), &numIterations);
    VERIFY2(hypreStatus == 0,
            "HYPRE_ParCSRGMRESGetNumIterations, error code " + std::to_string(hypreStatus));
  }

  double finalResNorm;
  hypreStatus = HYPRE_ParCSRGMRESGetFinalRelativeResidualNorm(solver.get(),
                                                              &finalResNorm);
  VERIFY2(hypreStatus == 0,
          "HYPRE_ParCSRGMRESGetFinalRelativeResidualNorm, error code " + std::to_string(hypreStatus));

  if (opt->saveLinearSystem) {
    const auto lhsName = "HypreLHS_" + description + "_" + std::to_string(solveCall);
    HYPRE_IJVectorPrint(lhs.get(), lhsName.c_str());
    ++solveCall;
  }

  if (numIterations >= opt->maxNumberOfIterations) {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(2);
    oss << "Hypre solver (" << description << ") did not converge, norm: " << finalResNorm;
    const auto message = oss.str();
    if (opt->quitIfDiverged) {
      VERIFY2(false, message);
    }
    else {
      if (opt->warnIfDiverged && Process::getRank() == 0) {
        std::cerr << message << std::endl;
      }
    }
  }

  return std::make_pair(numIterations, finalResNorm);
}

//------------------------------------------------------------------------------
// Multiply the linear system
//------------------------------------------------------------------------------
void
HypreFunctions::
multiplySystem(std::shared_ptr<HypreOptions> opt,
               std::shared_ptr<MatrixType> matrix,
               std::shared_ptr<VectorType> lhs,
               std::shared_ptr<VectorType> rhs,
               std::string description) {
  TIME_FUNCTION;
  // Get storage in CRS format
  HYPRE_ParCSRMatrix hypreParMatrix;
  HYPRE_ParVector hypreParLHS;
  HYPRE_ParVector hypreParRHS;

  int hypreStatus;
  hypreStatus = HYPRE_IJMatrixGetObject(matrix.get(),
                                        (void **)&hypreParMatrix);
  VERIFY2(hypreStatus == 0,
          "HYPRE_IJMatrixGetObject, error code " + std::to_string(hypreStatus));

  hypreStatus = HYPRE_IJVectorGetObject(lhs.get(),
                                        (void **)&hypreParLHS);
  VERIFY2(hypreStatus == 0,
          "HYPRE_IJVectorGetObject, error code " + std::to_string(hypreStatus));

  hypreStatus = HYPRE_IJVectorGetObject(rhs.get(),
                                        (void **)&hypreParRHS);
  VERIFY2(hypreStatus == 0,
          "HYPRE_IJVectorGetObject, error code " + std::to_string(hypreStatus));

  if (opt->saveLinearSystem) {
    const auto matName = "HypreMatrix_" + description + "_" + std::to_string(solveCall);
    const auto initName = "HypreLHS_" + description + "_" + std::to_string(solveCall);
    HYPRE_IJMatrixPrint(matrix.get(), matName.c_str());
    HYPRE_IJVectorPrint(lhs.get(), initName.c_str());
  }
  
  // Perform matrix multiplication, y = \alpha * A x + \beta y
  hypreStatus = HYPRE_ParCSRMatrixMatvec(1.0, // \alpha
                                         hypreParMatrix, // A
                                         hypreParLHS, // x
                                         0.0, // \beta
                                         hypreParRHS); // y
  
  if (opt->saveLinearSystem) {
    const auto rhsName = "HypreRHS_" + description + "_" + std::to_string(solveCall);
    HYPRE_IJVectorPrint(rhs.get(), rhsName.c_str());
    ++solveCall;
  }
  
  VERIFY2(hypreStatus == 0,
          "HYPRE_ParCSRMatrixMatvec, error code " + std::to_string(hypreStatus));
}

} // end namespace Spheral
