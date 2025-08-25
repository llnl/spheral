//---------------------------------Spheral++----------------------------------//
// HypreOptions
//
// Holds the options for Hypre solvers and preconditioners
//----------------------------------------------------------------------------//
#ifndef __Spheral_HypreOptions_hh__
#define __Spheral_HypreOptions_hh__

namespace Spheral {

struct HypreOptions {
  HypreOptions() { }

  // Quit program if solution did not converge
  bool quitIfDiverged = false;
  bool warnIfDiverged = true;

  // Sum all row values that share an index
  bool sumDuplicates = true;
  bool addToValues = false;
  
  // Shared options
  int logLevel = 0;
  int printLevel = 0;
  
  // Options used in GMRES solve
  bool saveLinearSystem = false;
  bool printIterations = false;
  int meanIterationsGuess = 5;
  int maxNumberOfIterations = 100;
  double toleranceL2 = 1.e-14;
  int useRobustTolerance = 1;
  double absoluteTolerance = 0.0;
  
  // Options used in GMRES initialization
  int kDim = 10;
  int minIters = 1;

  // General preconditioner options
  enum HyprePreconditionerType {
    NoPreconditioner,
    AMGPreconditioner,
    ILUPreconditioner
  };
  HyprePreconditionerType preconditionerType = HyprePreconditionerType::AMGPreconditioner;

  // Options used in BoomerAMG initialization
  int measure_type = 0;
  double pcTol = 0.;
  int minItersAMG = 1;
  int maxItersAMG = 1;
  int cycleType = 1;
  int coarsenTypeAMG = 6;
  double strongThresholdAMG = 0.25;
  double maxRowSumAMG = 0.5;
  int interpTypeAMG = -1;
  int aggNumLevelsAMG = -1;
  int aggInterpTypeAMG = -1;
  int pMaxElmtsAMG = -1;
  double truncFactorAMG = 0.0;
  int printLevelAMG = -1;
  int logLevelAMG = 0;
  double relaxWeightAMG = 1.0;
  int relaxTypeAMG = 8;
  int relaxTypeCoarseAMG = 19;
  int cycleNumSweepsAMG = 1;
  int cycleNumSweepsCoarseAMG = -1;
  int maxLevelsAMG = 25;

  // Options used in Euclid (ILU) initialization
  bool useILUT = false;
  int factorLevelILU = 1;
  int rowScaleILU = 0;
  int printLevelILU = -1;
  double dropToleranceILU = 0;
  
  // Postprocessing
  bool zeroOutNegativities = false;
}; // end struct HypreOptions

} // end namespace Spheral

#endif
