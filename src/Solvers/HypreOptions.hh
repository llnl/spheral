//---------------------------------Spheral++----------------------------------//
// HypreOptions
//
// Holds the options for Hypre solvers and preconditioners
//----------------------------------------------------------------------------//
#ifndef __Spheral_HypreOptions_hh__
#define __Spheral_HypreOptions_hh__

namespace Spheral {

enum class HyprePreset {
  Default,    // Don't change anything
  OldDefault, // Historical defaults
  Robust,     // Make more robust
  Fast        // Make faster
};

enum class HyprePreconditionerType {
  NoPreconditioner,   // Use GMRES without preconditioning
  AMGPreconditioner,  // BoomerAMG (default, best for diffusion problems)
  ILUPreconditioner   // Euclid ILU (for debugging or fallback)
};

struct HypreOptions {
  // Constructor
  HypreOptions(HyprePreset preset = HyprePreset::Default) {
    if (preset != HyprePreset::Default) {
      // Avoid recursion
      setPreset(preset);
    }
  }

  // ---------------------------------------------------------------------------
  // General options
  // ---------------------------------------------------------------------------
  
  // If true, abort execution when the linear solver fails to converge.
  bool quitIfDiverged = false;

  // If true, emit a warning when the solver fails to converge.
  bool warnIfDiverged = true;

  // If true, sum all row values that share an index
  bool addToValues = false;
  
  // 0 = no logging; set to 1+ for increasing verbosity.
  int logLevel = 0;

  // 0 = no printing; set to 1+ for increasing verbosity.
  int printLevel = 0;
  
  // ---------------------------------------------------------------------------
  // GMRES solve options
  // ---------------------------------------------------------------------------

  // Save the linear system (A, b, x) to disk for debugging.
  bool saveLinearSystem = false;

  // Print information after each iteration.
  bool printIterations = false;

  // Only used for keeping statistics.
  int meanIterationsGuess = 5;

  // Maximum number of GMRES iterations.
  int maxNumberOfIterations = 100;

  // Relative L2 convergence tolerance.
  double toleranceL2 = 1.e-8;

  // If enabled, use Hypre's robust convergence criterion.
  // Recommended for ill-conditioned systems.
  int useRobustTolerance = 1;

  // Absolute tolerance. Used with an OR along with the L2 tolerance above. 
  double absoluteTolerance = 1.e-12;
  
  // ---------------------------------------------------------------------------
  // GMRES initialization
  // ---------------------------------------------------------------------------

  // Krylov subspace dimension.
  int kDim = 10;

  // Minimum number of GMRES iterations before convergence is tested.
  int minIters = 1;

  // ---------------------------------------------------------------------------
  // Preconditioner selection
  // ---------------------------------------------------------------------------

  HyprePreconditionerType preconditionerType = HyprePreconditionerType::AMGPreconditioner;

  // ---------------------------------------------------------------------------
  // BoomerAMG options
  // ---------------------------------------------------------------------------

  // Strength-of-connection measure.
  // 0 = classical (based on matrix entries).
  // Appropriate for scalar diffusion operators.
  int measure_type = 0;

  // AMG convergence tolerance.
  // For preconditioning, this should be zero.
  double pcTol = 0.;

  // Minimum and maximum AMG iterations per application.
  // For a preconditioner, both should be 1 (single V-cycle).
  int minItersAMG = 1;
  int maxItersAMG = 1;

  // Multigrid cycle type.
  // 1 = V-cycle (default and recommended).
  // 2 = W-cycle (more expensive; rarely needed).
  int cycleType = 1;

  // Coarsening algorithm.
  // 6 = Falgout (better for structured problems)
  // 8 = PMIS (recommended for irregular graphs and large stencils).
  // 10 = HMIS (more conservative, may improve convergence)
  int coarsenTypeAMG = 8;

  // Strength-of-connection threshold.
  // Lower values retain more algebraic couplings.
  // Decrease to ~0.1 if not converging. Increase to ~0.7 for speed. 
  double strongThresholdAMG = 0.25;

  // Maximum allowed sum of weak connections relative to strong ones.
  // Used to prevent pathological strength graphs.
  // Values near 1.0 retain most couplings; smaller values prune aggressively.
  double maxRowSumAMG = 0.9;

  // Interpolation type.
  // 6 = extended classical interpolation (recommended for PMIS/HMIS).
  int interpTypeAMG = 6;

  // Number of aggressive coarsening levels.
  // -1 disables aggressive coarsening.
  int aggNumLevelsAMG = -1;

  // Interpolation type used during aggressive coarsening.
  // Only relevant if aggNumLevelsAMG >= 0.
  int aggInterpTypeAMG = -1;

  // Maximum number of nonzeros per row in interpolation operator.
  // Controls operator complexity.
  // Decrease to ~4 for more performance or increase to ~8 for more stability.
  int pMaxElmtsAMG = 6;

  // Truncation factor for interpolation coefficients.
  // 0.0 disables truncation by magnitude.
  // Increase to 0.1-0.2 to reduce interpolation complexity. 
  double truncFactorAMG = 0.0;

  // AMG-specific print verbosity.
  // Increase to > 0 for printing
  int printLevelAMG = -1;

  // AMG-specific logging level.
  // Increase to > 0 for logging
  int logLevelAMG = 0;

  // Relaxation weight for smoothers.
  // Usually 1.0 for GS-based methods.
  double relaxWeightAMG = 1.0;

  // Relaxation (smoother) type on fine and intermediate levels.
  // In testing, 8 and 16 work, but 88 and 89 do not. 
  // 8 = L1-scaled hybrid symmetric Gauss-Seidel (default)
  // 16 = Chebyshev (also works)
  int relaxTypeAMG = 8;

  // Relaxation type on the coarsest level.
  // -1 = Use value from relaxTypeAMG
  // 9 = Gaussian elimination
  // 19 = Gaussian elimination (old version)
  int relaxTypeCoarseAMG = 9;

  // Number of pre- and post-smoothing sweeps.
  // Typically 1 is sufficient.
  int cycleNumSweepsAMG = 1;

  // Number of coarse-level sweeps.
  // -1 defaults to cycleNumSweepsAMG unless overridden by solver type.
  int cycleNumSweepsCoarseAMG = -1;

  // Maximum number of AMG levels.
  // Acts as a safety cap.
  int maxLevelsAMG = 25;

  // Number of DOFs per node for systems of PDEs.
  // 1 = scalar (default). Set to 3 for coupled vector systems.
  int numFunctionsAMG = 1;

  // Nodal coarsening for systems of PDEs (requires numFunctionsAMG > 1).
  // 0 = unknown-based (coarsens each function independently).
  // 1 = Frobenius norm of nodal blocks (recommended for coupled systems).
  // 2 = sum of absolute values in each block.
  // 3 = largest element in each block.
  // 4 = row-sum norm.
  // 6 = sum of all values in each block.
  int nodalAMG = 1;

  // ---------------------------------------------------------------------------
  // Euclid (ILU) options
  // ---------------------------------------------------------------------------

  // Use ILUT instead of ILU(k).
  // Generally not useful except for debugging.
  bool useILUT = false;

  // Fill level for ILU(k).
  int factorLevelILU = 1;

  // Enable row scaling prior to ILU factorization.
  int rowScaleILU = 0;

  // Euclid print verbosity.
  // -1 inherits from global printLevel.
  int printLevelILU = -1;

  // Drop tolerance for ILUT or sparse ILU.
  double dropToleranceILU = 0;
  
  // ---------------------------------------------------------------------------
  // Postprocessing
  // ---------------------------------------------------------------------------

  // If true, clamp negative solution values to zero after solve.
  bool zeroOutNegativities = false;

  // ---------------------------------------------------------------------------
  // Set presets
  // ---------------------------------------------------------------------------
  
  void setPreset(HyprePreset preset) {
    // Reset to member defaults before applying preset
    *this = HypreOptions(HyprePreset::Default);

    switch (preset) {
    case HyprePreset::Default:
      break;

    case HyprePreset::OldDefault:
      // Historical defaults
      measure_type              = 0;
      pcTol                     = 0.0;
      minItersAMG               = 1;
      maxItersAMG               = 1;
      cycleType                 = 1;
      coarsenTypeAMG            = 6;
      strongThresholdAMG        = 0.25;
      maxRowSumAMG              = 0.5;
      interpTypeAMG             = -1;
      aggNumLevelsAMG           = -1;
      aggInterpTypeAMG          = -1;
      pMaxElmtsAMG              = -1;
      truncFactorAMG            = 0.0;
      printLevelAMG             = -1;
      logLevelAMG               = 0;
      relaxWeightAMG            = 1.0;
      relaxTypeAMG              = 8;
      relaxTypeCoarseAMG        = 19;
      cycleNumSweepsAMG         = 1;
      cycleNumSweepsCoarseAMG   = -1;
      maxLevelsAMG              = 25;
      break;

    case HyprePreset::Robust:
      // Improve convergence
      kDim                      = 20;    // Larger Krylov space
      cycleType                 = 2;     // W-cycle for better coarse convergence
      strongThresholdAMG        = 0.20;  // More strong connections
      maxRowSumAMG              = 1.0;   // Disable diagonal-dominance check
      pMaxElmtsAMG              = 8;     // More interpolation elements
      truncFactorAMG            = 0.0;   // No truncation
      break;

    case HyprePreset::Fast:
      // Improve speed where convergence isn't an issue
      toleranceL2               = 1.e-6; // Looser tolerance
      absoluteTolerance         = 1.e-10; 
      strongThresholdAMG        = 0.50;  // Aggressive coarsening
      maxRowSumAMG              = 0.5;   // Increase number of rows considered diagonally dominant
      aggNumLevelsAMG           = 1;     // One level of aggressive coarsening
      pMaxElmtsAMG              = 3;     // Fewer interpolation elements
      truncFactorAMG            = 0.15;  // Truncate weak entries
      break;
    }
  }
}; // end struct HypreOptions

} // end namespace Spheral

#endif
