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

  // ***************************************************************************
  // BoomerAMG options
  // ***************************************************************************
  struct AMG {
    // Strength-of-connection measure.
    // 0 = classical (based on matrix entries).
    int measureType = 0;

    // AMG convergence tolerance. For preconditioning, this should be zero.
    double tol = 0.;

    // Minimum and maximum AMG iterations per application.
    // For a preconditioner, both should be 1 (single V-cycle).
    int minIters = 1;
    int maxIters = 1;

    // Multigrid cycle type. 1 = V-cycle, 2 = W-cycle.
    int cycleType = 1;

    // Coarsening algorithm.
    // 6 = Falgout (structured), 8 = PMIS (irregular graphs), 10 = HMIS (conservative).
    int coarsenType = 8;

    // Strength-of-connection threshold.
    double strongThreshold = 0.25;

    // Maximum allowed sum of weak connections relative to strong ones.
    double maxRowSum = 0.9;

    // Interpolation type. 6 = extended classical (recommended for PMIS/HMIS).
    int interpType = 6;

    // Number of aggressive coarsening levels. -1 disables.
    int aggNumLevels = -1;

    // Interpolation type during aggressive coarsening. -1 = default.
    int aggInterpType = -1;

    // Maximum nonzeros per row in interpolation. -1 = default.
    int pMaxElmts = 6;

    // Truncation factor for interpolation coefficients. 0.0 = disabled.
    double truncFactor = 0.0;

    // AMG-specific print verbosity. -1 inherits from global printLevel.
    int printLevel = -1;

    // AMG-specific logging level.
    int logLevel = 0;

    // Relaxation weight for smoothers.
    double relaxWeight = 1.0;

    // Relaxation (smoother) type on fine and intermediate levels.
    // 8 = L1-scaled hybrid symmetric Gauss-Seidel (default), 16 = Chebyshev.
    int relaxType = 8;

    // Relaxation type on the coarsest level. -1 = use relaxType, 9 = Gaussian elimination.
    int relaxTypeCoarse = 9;

    // Number of pre- and post-smoothing sweeps.
    int cycleNumSweeps = 1;

    // Number of coarse-level sweeps. -1 = use cycleNumSweeps.
    int cycleNumSweepsCoarse = -1;

    // Maximum number of AMG levels.
    int maxLevels = 25;

    // Number of DOFs per node for systems of PDEs. 1 = scalar.
    int numFunctions = 1;

    // Nodal coarsening for systems (requires numFunctions > 1).
    // 0 = unknown-based, 1 = Frobenius norm, 2 = sum abs, 3 = largest, 4 = row-sum, 6 = sum all.
    int nodal = 1;
  };

  // ***************************************************************************
  // Euclid (ILU) options
  // ***************************************************************************
  struct ILU {
    // Use ILUT instead of ILU(k).
    bool useILUT = false;

    // Fill level for ILU(k).
    int factorLevel = 1;

    // Enable row scaling prior to ILU factorization.
    int rowScale = 0;

    // Euclid print verbosity. -1 inherits from global printLevel.
    int printLevel = -1;

    // Drop tolerance for ILUT or sparse ILU.
    double dropTolerance = 0;
  };

  // Sub-struct members
  AMG amg;
  ILU ilu;
  
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
      amg.measureType            = 0;
      amg.tol                    = 0.0;
      amg.minIters               = 1;
      amg.maxIters               = 1;
      amg.cycleType              = 1;
      amg.coarsenType            = 6;
      amg.strongThreshold        = 0.25;
      amg.maxRowSum              = 0.5;
      amg.interpType             = -1;
      amg.aggNumLevels           = -1;
      amg.aggInterpType          = -1;
      amg.pMaxElmts              = -1;
      amg.truncFactor            = 0.0;
      amg.printLevel             = -1;
      amg.logLevel               = 0;
      amg.relaxWeight            = 1.0;
      amg.relaxType              = 8;
      amg.relaxTypeCoarse        = 19;
      amg.cycleNumSweeps         = 1;
      amg.cycleNumSweepsCoarse   = -1;
      amg.maxLevels              = 25;
      break;

    case HyprePreset::Robust:
      // Improve convergence
      kDim                       = 20;    // Larger Krylov space
      amg.cycleType              = 2;     // W-cycle for better coarse convergence
      amg.strongThreshold        = 0.20;  // More strong connections
      amg.maxRowSum              = 1.0;   // Disable diagonal-dominance check
      amg.pMaxElmts              = 8;     // More interpolation elements
      amg.truncFactor            = 0.0;   // No truncation
      break;

    case HyprePreset::Fast:
      // Improve speed where convergence isn't an issue
      toleranceL2                = 1.e-6; // Looser tolerance
      absoluteTolerance          = 1.e-10; 
      amg.strongThreshold        = 0.50;  // Aggressive coarsening
      amg.maxRowSum              = 0.5;   // Increase number of rows considered diagonally dominant
      amg.aggNumLevels           = 1;     // One level of aggressive coarsening
      amg.pMaxElmts              = 3;     // Fewer interpolation elements
      amg.truncFactor            = 0.15;  // Truncate weak entries
      break;
    }
  }
}; // end struct HypreOptions

} // end namespace Spheral

#endif
