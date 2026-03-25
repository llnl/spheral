from PYB11Generator import *

HyprePreconditionerType = PYB11enum(("NoPreconditioner",
                                     "AMGPreconditioner",
                                     "ILUPreconditioner"),
                                    export_values = False,
                                    doc = "Preconditioners available for Hypre")

HyprePreset = PYB11enum(("Default",
                         "OldDefault",
                         "Robust",
                         "Fast"),
                        export_values = False,
                        doc = "Hypre presets for speed or robustness")

@PYB11holder("std::shared_ptr")
class HypreOptions:
    def pyinit(self,
               preset = ("HyprePreset", "HyprePreset::Default")):
        "Holds the options for Hypre preconditioners and solver"

    # General options
    quitIfDiverged = PYB11readwrite(doc="Abort execution when the linear solver fails to converge")
    warnIfDiverged = PYB11readwrite(doc="Emit a warning when the solver fails to converge")
    addToValues = PYB11readwrite(doc="Sum all row values that share an index")
    logLevel = PYB11readwrite(doc="Logging verbosity (0 = off)")
    printLevel = PYB11readwrite(doc="Print verbosity (0 = off)")

    # GMRES solve options
    saveLinearSystem = PYB11readwrite(doc="Save the linear system (A, b, x) to disk for debugging")
    printIterations = PYB11readwrite(doc="Print information after each iteration")
    meanIterationsGuess = PYB11readwrite(doc="Initial guess for mean iterations (statistics only)")
    maxNumberOfIterations = PYB11readwrite(doc="Maximum number of GMRES iterations")
    toleranceL2 = PYB11readwrite(doc="Relative L2 convergence tolerance")
    useRobustTolerance = PYB11readwrite(doc="Use Hypre robust convergence criterion for ill-conditioned systems")
    absoluteTolerance = PYB11readwrite(doc="Absolute tolerance (OR with L2 tolerance)")

    # GMRES initialization
    kDim = PYB11readwrite(doc="Krylov subspace dimension")
    minIters = PYB11readwrite(doc="Minimum GMRES iterations before convergence test")

    # Preconditioner selection
    preconditionerType = PYB11readwrite(doc="Preconditioner type (AMG, ILU, or none)")

    #---------------------------------------------------------------------------
    # AMG nested type — BoomerAMG options
    #---------------------------------------------------------------------------
    class AMG:
        def pyinit(self):
            "BoomerAMG preconditioner options"

        measureType = PYB11readwrite(doc="Strength-of-connection measure (0 = classical)")
        tol = PYB11readwrite(doc="AMG convergence tolerance (0 for preconditioning)")
        minIters = PYB11readwrite(doc="Minimum AMG iterations per application")
        maxIters = PYB11readwrite(doc="Maximum AMG iterations per application")
        cycleType = PYB11readwrite(doc="Multigrid cycle type (1 = V-cycle, 2 = W-cycle)")
        coarsenType = PYB11readwrite(doc="Coarsening algorithm (6 = Falgout, 8 = PMIS, 10 = HMIS)")
        strongThreshold = PYB11readwrite(doc="Strength-of-connection threshold")
        maxRowSum = PYB11readwrite(doc="Maximum allowed sum of weak connections relative to strong ones")
        interpType = PYB11readwrite(doc="Interpolation type (6 = extended classical)")
        aggNumLevels = PYB11readwrite(doc="Number of aggressive coarsening levels (-1 = disabled)")
        aggInterpType = PYB11readwrite(doc="Interpolation type during aggressive coarsening (-1 = default)")
        pMaxElmts = PYB11readwrite(doc="Maximum nonzeros per row in interpolation (-1 = default)")
        truncFactor = PYB11readwrite(doc="Truncation factor for interpolation coefficients (0.0 = disabled)")
        printLevel = PYB11readwrite(doc="AMG print verbosity (-1 inherits from global)")
        logLevel = PYB11readwrite(doc="AMG logging level")
        relaxWeight = PYB11readwrite(doc="Relaxation weight for smoothers")
        relaxType = PYB11readwrite(doc="Relaxation type on fine/intermediate levels (8 = L1 GS, 16 = Chebyshev)")
        relaxTypeCoarse = PYB11readwrite(doc="Relaxation type on coarsest level (-1 = use relaxType, 9 = Gaussian elimination)")
        cycleNumSweeps = PYB11readwrite(doc="Number of pre- and post-smoothing sweeps")
        cycleNumSweepsCoarse = PYB11readwrite(doc="Number of coarse-level sweeps (-1 = use cycleNumSweeps)")
        maxLevels = PYB11readwrite(doc="Maximum number of AMG levels")
        numFunctions = PYB11readwrite(doc="DOFs per node for systems of PDEs (1 = scalar)")
        nodal = PYB11readwrite(doc="Nodal coarsening for systems (requires numFunctions > 1)")

    #---------------------------------------------------------------------------
    # ILU nested type — Euclid ILU options
    #---------------------------------------------------------------------------
    class ILU:
        def pyinit(self):
            "Euclid ILU preconditioner options"

        useILUT = PYB11readwrite(doc="Use ILUT instead of ILU(k)")
        factorLevel = PYB11readwrite(doc="Fill level for ILU(k)")
        rowScale = PYB11readwrite(doc="Enable row scaling prior to ILU factorization")
        printLevel = PYB11readwrite(doc="Euclid print verbosity (-1 inherits from global)")
        dropTolerance = PYB11readwrite(doc="Drop tolerance for ILUT or sparse ILU")

    #---------------------------------------------------------------------------
    # Top-level sub-struct members
    #---------------------------------------------------------------------------
    amg = PYB11readwrite(doc="BoomerAMG preconditioner options")
    ilu = PYB11readwrite(doc="Euclid ILU preconditioner options")

    # Postprocessing
    zeroOutNegativities = PYB11readwrite(doc="Clamp negative solution values to zero after solve")

